# -*- coding: utf-8 -*-
"""

"""
import biosteam as bst
from .property_package import STRAP_chemicals_outline, create_property_package, create_property_package_MSW
from .process_settings import GWP as GWP_key, load_process_settings
from .systems import (
    create_single_layer_batch_separation_system, 
    create_multilayer_batch_separation_system,
    create_STRAPMSW_system
)
from .tea import create_baseline_tea
from .data import price_distributions_2023 as dist
from .data.lca_characterization_factors import indicators, set_CFs as set_GWPCF
CFs = indicators['GWP']
from chaospy import distributions as shape
from plastics import strap
from biosteam.utils import CABBI_colors, GG_colors, colors
from scipy.optimize import minimize
import flexsolve as flx
import thermosteam as tmo
import numpy as np
import os

__all__ = (
    'BaselineSTRAPProcess',
    'STRAPMSWProcess',
    'STRAPProcessPE',
    'MagnetRecovery',
    'define_solvent',
    'define_dissolution',
    'define_precipitation',
)

# %% Old STRAP-MSW compositional analysis
# import biosteam as bst

# class ChemicalData:
#     __slots__ = ('titer', 'conversion', 'theoretical_yield', 'chemical')
    
#     def __init__(self, ID, titer, conversion, theoretical_yield):
#         self.titer = titer
#         self.conversion = conversion
#         self.theoretical_yield = theoretical_yield
#         if ID is not None: self.chemical = bst.Chemical(ID, db='BioSTEAM')
        
#     def __getattr__(self, name):
#         return getattr(self.chemical, name)

# # Glucose, xylose, total
# glucan = bst.Chemical('glucan', db='BioSTEAM')
# xylan = bst.Chemical('xylan', db='BioSTEAM')
# glucose = bst.Chemical('glucose')
# xylose = bst.Chemical('xylose')
# glucan = 31.3 * (glucan.MW / glucose.MW)
# xylan = 5.2 * (xylan.MW / xylose.MW)
# glucose = ChemicalData('glucose', 34.4, 17.2 / 100, 54.9 / 100)
# xylose = ChemicalData('xylose', 6.8, 3.4 / 100, 65.8 / 100)
# sugar = ChemicalData(None, 41.2, 20.6 / 100, 50.6 / 100)
# total_material = sugar.titer / sugar.conversion # per L
# total_glucan = glucose.titer / glucose.theoretical_yield * glucan.MW / glucose.MW
# total_xylan = xylose.titer / xylose.theoretical_yield * xylan.MW / xylan.MW
# biogenic_content = 0.781 
# target_polymer_fraction = 0.766 # Of the plastic
# total_without_target_plastic = 1 - (1 - biogenic_content) * target_polymer_fraction
# biogenic_content_without_PEPP = biogenic_content / total_without_target_plastic
# total_biomaterial = total_material * biogenic_content_without_PEPP
# total_lignin_and_others = total_biomaterial - total_glucan - total_xylan

# %%

bst.System.strict_convergence = False
kg_per_ton = 907.18474
kg_per_MT = 1000
L_per_gal = 3.7854
ethanol_kg_per_gal = 2.98668849
ethanol_gal_per_kg = 1. / ethanol_kg_per_gal
ethanol_L_per_kg = ethanol_gal_per_kg * L_per_gal
ethanol_kg_per_L = 1. / ethanol_L_per_kg

# https://www.recyclingpyrolysisplant.com/FAQ/pyrolysis_plant/pyrolysis-oil-applications-68.html
pyrolysis_oil_density_kg_per_m3 = 820.5 # kg / m3
pyrolysis_oil_price_range_USD_per_L = np.array([0.455, 0.78]) # USD / L
pyrolysis_oil_price_range_USD_per_kg = 1000 * pyrolysis_oil_price_range_USD_per_L / pyrolysis_oil_density_kg_per_m3
# https://pubs.acs.org/doi/10.1021/acssuschemeng.9b04763

# https://resource-recycling.com/plastics/2024/05/15/recycled-plastic-prices-continue-to-climb-higher/
recycled_PP_price_USD_per_ton = 0.06 * 2.20462 * 907.185 # Recycled PP
recycled_HDPE_price_USD_per_ton = 0.3244 * 2.20462 * 907.185 # Recycled natural HDPE

# Bloomberg original source: https://www.statista.com/statistics/1171074/price-high-density-polyethylene-forecast-globally/
HDPE_price_range = (837., 1211.)

def update_chemicals_outline(plastic, solvent):
    if plastic not in STRAP_chemicals_outline:
        STRAP_chemicals_outline.extend([
            bst.ChemicalDraft(
                plastic, # Model generic plastic film as PET.
                formula='C10H8O4',
                search_db=False,
                phase='s',
                rho=1380, # kg / m3,
                Tm=523,
                Tb=623,
                Cp=1,
                default=True,
                LHV= 21285 * 192.16812,
            ),
            bst.ChemicalDraft(
                plastic + 'oligomer', # Model generic dissolved resin as hexene
                search_ID='1-Hexene',
                CAS=plastic + 'oligomer',
            ),
        ])
    if (solvent not in STRAP_chemicals_outline 
        and solvent not in solvent_mixture_names):
        STRAP_chemicals_outline.append(solvent)

solvent_mixtures = [
]

solvent_mixture_names = set(['DMSOWater'])

def define_solvent(
        name, chemicals, composition, wt=True
    ):
    solvent_mixtures.append(
        (name, chemicals, composition, wt)
    )
    solvent_mixture_names.add(name)

def default_plastic_solvent_pair(plastic, solvent):
    update_chemicals_outline(plastic, solvent)
    define_dissolution(plastic, solvent, override=False)
    define_precipitation(plastic, solvent, override=False)

def define_dissolution(
        plastic: str,
        solvent: str,
        capacity: float=0.05,
        solvent_content: float=0.5,
        T: float=130 + 273.15,
        tau: float=0.5,
        override: bool=True,
    ):
    name = f'{plastic}_{solvent}_dissolution'
    if not override and hasattr(strap.dissolution_steps, name): return 
    def f():
        return strap.dissolution_steps.DissolutionStep(
            plastic, plastic + 'oligomer', solvent, 
            tmo.Reaction(f'{plastic} -> {plastic}oligomer', plastic, X=1.0, basis='wt'), 
            capacity, solvent_content, T, tau,
        )
    setattr(strap.dissolution_steps, name, f)
    f.__name__ = name

def define_precipitation(
        plastic: str,
        solvent: str,
        solubility: float=0,
        precipitate_solvent_content: float=0.8,
        screw_press_solvent_content: float=0.4,
        T: float=308.15,
        tau: float=0.5,
        T_condensation: float=None, # Not actually used in new configuration
        override: bool=True,
    ):
    name = f'{plastic}_{solvent}_precipitation'
    if not override and hasattr(strap.precipitation_steps, name): return 
    def f():
        return strap.precipitation_steps.PrecipitationStep(
            solvent,
            plastic,
            plastic + 'oligomer',
            solubility,
            precipitate_solvent_content,
            screw_press_solvent_content,
            T,
            tau,
            T_condensation,
        )
    setattr(strap.precipitation_steps, name, f)
    f.__name__ = name
    
# TODO: 
#   Set up feed composition (biogenic material + plastic)
#   Set up Xylene as a solvent (dissolution/precipitation temperature, recovery, solvent/plastic ratio)
#   Simulate and test STRAP system
#   Set up biogenic material to ethanol process (no pretreatment, preliminary hydrolysis/ferm performance)

class STRAPMSWProcess(bst.ProcessModel):
    """
    Create a model for a solvent targeted precipitation and dissolution process.
    The dissolution and precipitation steps default to PE.
    
    Examples
    --------
    >>> from plastics import strap
    >>> pm = strap.STRAPMSWProcess(simulate=False, scenario='baseline')
    >>> pm.system.diagram(kind='cluster', number=True)
    >>> assumptions, results = pm.baseline()
    >>> results
    Ethanol  GWP [kg*CO2e/L]   0.00158
    -        IRR [%]            0.0285
             MSP [USD/kg]         1.41
    dtype: float64
    
    >>> assumptions
    MSW                  Tipping fee [USD/ton]                     58.5
                         Biogenic content [wt %]                  78.1
    Plastics             Solute content [wt %]                  0.0015
    MSW                  Processing capacity [ton/yr]          1.62e+05
    Resin                Price [USD/ton]                            200
    Solvent              Price [USD/kg]                            2.17
    Specification        IRR [%]                                   0.15
    Polymer              Mass fraction                            0.766
    Centrifuged plastic  Solvent content [%]                         50
    Precipitate          Solvent content [%]                        0.8
    Dissolution          Solvent capacity [wt %]                      2
    Ethanol              Price [USD/L]                            0.358
    RIN D3               Price [USD/RIN]                          0.534
    Natural gas          Price [USD/m3]                           0.167
    Electricity          Price [USD/kWh]                         0.0689
    Cellulase            Price [USD/kg]                           0.212
                         Cellulase loading [wt % cellulose]       0.02
    Saccharification     Solids loading [wt % solids]            0.667
    Cellulase            GWP [kg*CO2e/kg]                          8.07
    Saccharification     Glucose yield [%]                         54.9
                         Xylose yield [%]                          65.8
    Cofermenation        Xylose to ethanol yield [%]               71.1
                         Glucose to ethanol yield [%]              71.1
    Cofermentation       Ethanol productivity [g/L/h]             0.172
                         Ethanol titer [g/L]                       8.27
    dtype: float64
    
    >>> from plastics import strap
    >>> pm = strap.STRAPMSWProcess(simulate=False, scenario='potential')
    >>> assumptions, results = pm.baseline()
    >>> results
    Ethanol  GWP [kg*CO2e/L]   0.00101
    -        IRR [%]             0.332
             MSP [USD/kg]       -0.912
    dtype: float64
    
    >>> assumptions
    MSW                  Tipping fee [USD/ton]                     102
                         Biogenic content [wt %]                 78.1
    Plastics             Solute content [wt %]                 0.0015
    MSW                  Processing capacity [ton/yr]          7.7e+05
    Resin                Price [USD/ton]                           200
    Solvent              Price [USD/kg]                           2.17
    Specification        IRR [%]                                  0.15
    Polymer              Mass fraction                           0.766
    Centrifuged plastic  Solvent content [%]                        50
    Precipitate          Solvent content [%]                       0.8
    Dissolution          Solvent capacity [wt %]                     2
    Ethanol              Price [USD/L]                           0.358
    RIN D3               Price [USD/RIN]                         0.534
    Natural gas          Price [USD/m3]                          0.167
    Electricity          Price [USD/kWh]                        0.0689
    Cellulase            Price [USD/kg]                          0.212
                         Cellulase loading [wt % cellulose]      0.02
    Saccharification     Solids loading [wt % solids]           0.667
    Cellulase            GWP [kg*CO2e/kg]                         8.07
    Saccharification     Glucose yield [%]                          90
                         Xylose yield [%]                           90
    Cofermenation        Xylose to ethanol yield [%]                85
                         Glucose to ethanol yield [%]               95
    Cofermentation       Ethanol productivity [g/L/h]              1.5
                         Ethanol titer [g/L]                        54
    dtype: float64
    
    """
    cache = {}
    area_hatches = {
        'MRF': '',
        'STRAP': '',
        'EtOH': '',
        'WWT': '',
        'Facilities': '',
    }
    area_colors = {
        'STRAP': GG_colors.blue, 
        'EtOH': GG_colors.purple,
        'WWT': GG_colors.green,
        'Facilities': GG_colors.red,
        'MRF': CABBI_colors.brown,
    }
    # Tipping fee based on: https://erefdn.org/2022-msw-landfill-tipping-fees/
    # 2022 survay, mean +- std
    tipping_fees_by_region = {
        'nationwide': (58.47, 35.71),
        'pacific': (69.02, 50.79), 
        'northeast': (75.92, 22.98),
        'mountains/plains': (50.84, 29.74), 
        'midwest': (62.02, 39.01), 
        'southeast': (44.75, 14.59),
        'southcentral': (48.70, 26.54), 
    }
    tipping_fees_by_state = {
        # pacific
        'AK': (186.49, 126.63),
        'AZ': (43.88, 13.03),
        'CA': (60.56, 20.43),
        'HI': (107, 6),
        'ID': (41.64, 22.64),
        'NV': (36.11, 5.78),
        'OR': (64.82, 26.58),
        'WA': (99.66, 43.95),
        # northeast
        'CT': (np.nan, np.nan),
        'DE': (85, 0),
        'ME': (107.67, 18.61), 
        'MD': (70.91, 11.36),
        'MA': (102.50, 10.61),
        'NH': (80.75, 6.01),
        'NJ': (81.54, 6.24), 
        'NY': (61.40, 27.27),
        'PA': (89.72, 16.73),
        'RI': (115, 0),
        'VT': (np.nan, np.nan),
        'VA': (59.89, 18.40),
        'WV': (43.90, 3.50),
        # mountains/plains
        'CO': (45.30, 23.47), 
        'MT': (28.20, 5.31),
        'ND': (47.86, 2.01),
        'SD': (54.88, 8.95),
        'UT': (34.38, 4.65),
        'WY': (86.55, 39),
        # midwest
        'IL': (75.00, 43.13), 
        'IN': (41.60, 26.13),
        'IA': (57.81, 28.81),
        'KS': (54.66, 15.24),
        'MI': (52.38, 31.51),
        'MN': (83.66, 31.41),
        'MO': (125.47, 98.17),
        'NE': (49.97, 17.44),
        'OH': (56.63, 24.73),
        'WI': (54.50, 10.39), 
        # southeast
        'AL': (37.32, 13.04),
        'FL': (50.06, 15.26),
        'GA': (44.88, 14.14),
        'KY': (51.61, 27.12),
        'MS': (35.17, 15.32),
        'NC': (44.72, 7.34),
        'SC': (49.20, 11.95),
        'TN': (42.69, 11.30), 
        # southcentral
        'AR': (66.02, 40.16),
        'LA': (35.58, 10.72),
        'NM': (39.23, 6.44),
        'OK': (62.68, 27.41),
        'TX': (45.45, 25.26),
    }
    regions = {
        'pacific': ('AK', 'AZ', 'CA', 'HI', 'ID', 'NV', 'OR', 'WA'), 
        'northeast': ('CT', 'DE', 'ME', 'MD', 'MA', 'NH', 'NJ', 'NY', 'PA', 'RI', 'VT', 'VA', 'WV'),
        'mountains/plains': ('CO', 'MT', 'ND', 'SD', 'UT', 'WY'), 
        'midwest': ('IL', 'IN', 'IA', 'KS', 'MI', 'MN', 'MO', 'NE', 'OH', 'WI'), 
        'southeast': ('AL', 'FL', 'GA', 'KY', 'MS', 'NC', 'SC', 'TN'), 
        'southcentral': ('AR', 'LA', 'NM', 'OK', 'TX'),
    }
    # tons/year: https://erefdn.org/2022-msw-landfill-tipping-fees/
    landfill_sizes = dict(
        small=31434,
        medium=162287,
        large=770041,
    )
    # Baseline is based on a large MRF in MA https://www.recyclingtoday.com/article/largest-material-recovery-facilities-on-the-rebound/
    # Theoretical is 19% of material sent to a large landfill: https://erefdn.org/2022-msw-landfill-tipping-fees/
    # 19% of material sent to a MRF can be assumed to end up as MRF residues: https://www.sciencedirect.com/science/article/pii/S0956053X24006408
    MRF_residue_availabilities = dict(
        baseline=173250 * 0.19,
        theoretical=770041 * 0.19,
    )
    
    baseline_resin_price = np.mean(HDPE_price_range)
    resin_price_range = HDPE_price_range
    
    # Baseline fermentation performance
    # TODO: Update with values from compositional analysis
    experimental_glucose_yield = (44.7, 1e-6) # %
    experimental_xylose_yield = (38.1, 1e-6) # %
    experimental_sugars_to_ethanol_yield = (74.7, 1.3) # %
    # experimental_xylose_to_ethanol_yield = (74.7, 1.3) # %
    # experimental_glucose_to_ethanol_yield = (74.7, 1.3)
    experimental_cofermentation_ethanol_productivity = (0.523, 0.043)
    experimental_cofermentation_ethanol_titer = (11.5, 0.2)
    # experimental_glucose_yield = (54.9, 13.2) # %
    # experimental_xylose_yield = (65.8, 3.9) # %
    # experimental_xylose_to_ethanol_yield = (71.1, 2.1) # %
    # experimental_glucose_to_ethanol_yield = (71.1, 2.1)
    # experimental_cofermentation_ethanol_productivity = (0.1722, 0.0192)
    # experimental_cofermentation_ethanol_titer = (8.26514, 0.922)
    
    # NREL fermentation performance
    NREL_glucose_yield = 90
    NREL_xylose_yield = 90
    NREL_sugars_to_ethanol_yield = (85 * 3.91 + 95 * 6.97) / (3.91 + 6.97)
    # NREL_xylose_to_ethanol_yield = 85
    # NREL_glucose_to_ethanol_yield = 95
    NREL_cofermentation_ethanol_productivity = 1.5
    NREL_cofermentation_ethanol_titer = 54.0
    
    # Prospective fermentation performance
    # Mid point between baseline and NREL performance.
    # prospective_glucose_yield = 0.5 * (NREL_glucose_yield + experimental_glucose_yield[0])
    # prospective_xylose_yield = 0.5 * (NREL_xylose_yield + experimental_xylose_yield[0])
    # prospective_xylose_to_ethanol_yield = 0.5 * (NREL_xylose_to_ethanol_yield + experimental_xylose_to_ethanol_yield[0])
    # prospective_glucose_to_ethanol_yield = 0.5 * (NREL_glucose_to_ethanol_yield + experimental_glucose_to_ethanol_yield[0])
    # prospective_cofermentation_ethanol_productivity = 0.5 * (NREL_cofermentation_ethanol_productivity + experimental_cofermentation_ethanol_productivity[0])
    # prospective_cofermentation_ethanol_titer = 0.5 * (NREL_cofermentation_ethanol_titer + experimental_cofermentation_ethanol_titer[0])
    
    # # Default values based on finding % change in all parameters to get MSP equal to the market price.
    # prospective_glucose_yield = 69.7
    # prospective_xylose_yield = 77.0
    # prospective_xylose_to_ethanol_yield = 77.0
    # prospective_glucose_to_ethanol_yield = 81.2
    # prospective_cofermentation_ethanol_productivity = 0.731
    # prospective_cofermentation_ethanol_titer = 27.5
    
    # # Default values based on minimizing % change in parameters to get MSP equal to the market price.
    # prospective_glucose_yield = 54.9 + (72.88128026385473 - 54.9)
    # prospective_xylose_yield = 65.8 + (70.34379753877613 - 65.8)
    # prospective_xylose_to_ethanol_yield = 71.1 + (72.16367656868835 - 71.1)
    # prospective_glucose_to_ethanol_yield = 71.1 + (77.37170892159932 - 71.1)
    # prospective_cofermentation_ethanol_productivity = 0.1722 + (0.5225144196586275 - 0.1722)
    # prospective_cofermentation_ethanol_titer = 8.26514 + (34.03267574550374 - 8.26514)
    
    @classmethod
    def get_tipping_fee(cls, location):
        if location in cls.tipping_fees_by_region:
            return cls.tipping_fees_by_region[location]
        elif location in cls.tipping_fees_by_state:
            return cls.tipping_fees_by_state[location]
        else:
            raise ValueError(f"invalid location {location!r}; location must be a state or a region")
    
    class Scenario:
        solvent: str = '# Solvent used to separate the target plastic.'
        target_plastic: str = '# The polymer layer being dissolved.'
        target_plastic_percent: float = '# Fraction of target plastic in feedstock [%].'
        biogenic_material_percent: float = '# Fraction of biogenic material in feedstock [%].'
        precipitation_configuration: str = 'integrated heat transfer', "# Must be either 'solvent mixing' or 'integrated heat transfer'."
        fermentation_performance: str = 'experimental', "# Must be either 'full range', 'experimental', 'prospective', or 'NREL'."
        processing_capacity: str = 'baseline', "# Must be either 'full range', 'baseline', or 'theoretical'."
        tipping_fee_location: str = 'nationwide', "# Must be either a state (e.g., 'WI'), 'nationwide', or by region (e.g., 'northeast', 'midwest')"
        resin_product_price: str = 'MRF waste plastic' "# Must be either 'MRF waste plastic', 'none', or 'full range'"
        preprocessing: bool = False, '# Whether or not to include preprocessing'
        name: str = None, '# Shorthand name of scenario'
        hand_sorted_percent: float = 10.2, '# Fraction of hand sorted material in MSW'
        metals_percent: float = 2.4, '# Fraction of metals in MSW by wt'
        unders_percent: float = 6.42, '# Fraction of unders in MSW by wt'
        overs_percent: float = 2.98, '# Fraction of overs in MSW by wt'
        @property
        def bulk_plastic_percent(self):
            return 100 - self.target_plastic_percent - self.biogenic_material_percent
    
    @property
    def name(self):
        scenario = self.scenario
        return (
            "STRAPMSW_"
            f"{scenario.fermentation_performance}_fermentation_performance_"
            f"{scenario.processing_capacity}_processing_capacity_"
            f"{scenario.tipping_fee_location}_tipping_fee"
        )
    
    # Ideas for figures
    # 1. Monte Carlo simulations baseline vs. prospective scenarios.
    # 2. USA map, IRR by state tipping fees, prospective scenarios.
    # 3. Monte Carlo across process capacity for prospective scenario.
    # 4. SI: Monte Carlo and Spearman's rho full landscape of scenarios.
    
    @classmethod
    def run_monte_carlo(cls, scenario=None, N=2000):
        return strap.run_monte_carlo(N, ProcessModel=cls, scenario=scenario)
    
    @classmethod
    def as_scenario(cls, scenario):
        name = scenario
        scenario = scenario.replace('-', '_').replace(' ', '_')
        preprocessing = 'verification' in scenario
        if scenario == 'all':
            fermentation_performance = processing_capacity = resin_product_price = 'full range'
            tipping_fee = 'nationwide'
        elif scenario == 'verification_baseline':
            fermentation_performance = 'experimental'
            processing_capacity = 'full range'
            tipping_fee = 'nationwide'
            resin_product_price = 'HDPE'
        elif scenario == 'verification_NREL':
            fermentation_performance = 'NREL'
            processing_capacity = 'full range'
            tipping_fee = 'nationwide'
            resin_product_price = 'HDPE'
        elif scenario in {'experimental', 'NREL'}:
            fermentation_performance = scenario
            processing_capacity = 'baseline'
            tipping_fee = '*MA'
            resin_product_price = 'HDPE'
        elif scenario == 'baseline':
            fermentation_performance = 'experimental'
            processing_capacity = 'baseline'
            tipping_fee = '*MA'
            resin_product_price = 'HDPE'
        elif scenario == 'potential':
            fermentation_performance = 'NREL'
            processing_capacity = 'baseline'
            tipping_fee = '*MA'
            resin_product_price = 'HDPE'
        else:
            raise ValueError(f'invalid scenario {scenario!r}')
        
        # Composition from experimental data:
        # PEPP = 2.54 / 15.13
        # other plastic = 0.62 / 10.08 * (1 - PEPP)
        # biogenic = 1 - (PET + PEPP)
        return cls.Scenario(
            'Xylene', 'PEPP', 16.79, 78.09, 'integrated heat transfer',
            fermentation_performance, processing_capacity, tipping_fee,
            resin_product_price, preprocessing, name
        )
    
    @classmethod
    def default_scenario(cls):
        return cls.as_scenario('all')
    
    def create_thermo(self):
        # Assume Michigan 100g experiment
        chemicals = create_property_package_MSW()
        scenario = self.scenario
        chemicals.set_alias('Extract', 'Extractives')
        chemicals.define_group(
            'Plastic', 
            [scenario.target_plastic, 'BulkPlastic'], 
            [scenario.target_plastic_percent, 
             scenario.bulk_plastic_percent],
            wt=True
        )
        # TODO: Update biogenic material composition with full data.
        biogenic_material_composition = [
            ('Ash', 11.28),
            ('Extractives', 2.95 + 1.3),
            ('Lignin', 34.64),
            ('SolubleLignin', 2.41),
            ('Glucan', 37.7 + 0.89 + 2.7), # Includes 0.89 Galactan and 2.7 Mannan
            ('Xylan', 6.83),
            ('Acetate', 0.28),
        ]
        chemicals.define_group(
            'BiogenicMaterial', 
            *zip(*biogenic_material_composition),
            wt=True
        )
        chemicals.define_group(
            'Solutes', 
            ['Solubles', 'Protein', 'Acetate'], 
            [1, 0, 0], 
            wt=True
        )
        chemicals.define_group(
            'MSW_Residue',
            ['Ash', 'Extractives', 'Lignin', 
             'SolubleLignin', 'Glucan', 'Xylan', 
             'Acetate', 'Solubles', 'Protein', 'Acetate',
             scenario.target_plastic, 'BulkPlastic'] 
        )
        chemicals.define_group(
            'NonResidue',
            ['HandSorted', 'Metal', 'Overs', 'Unders'] 
        )
        # TODO: Why is Xylan in the same amount as Glucan according data?
        # polymer + water -> monomer
        # xylose_titer = 6.8 # g / L final
        # glucose_titer = 34.4 # g / L final
        # xylose_yield = 3.4 # % wt
        # glucose_yield = 17.2 # % wt
        # xylose = xylose_titer / (xylose_yield/100) # g / L initial
        # glucose = glucose_titer / (glucose_yield/100)
        # print(glucose / xylose) -> 1.0
        
        return chemicals
        
    def create_system(self):
        load_process_settings()
        scenario = self.scenario
        self.dissolution_step = dissolution_step = getattr(strap.dissolution_steps, f'{scenario.target_plastic }_{scenario.solvent}_dissolution')()
        self.precipitation_step = precipitation_step = getattr(strap.precipitation_steps, f'{scenario.target_plastic }_{scenario.solvent}_precipitation')()
        system = create_STRAPMSW_system(
            preprocessing=scenario.preprocessing,
            dissolution_step=dissolution_step,
            precipitation_step=precipitation_step,
            precipitation_configuration=scenario.precipitation_configuration,
        )
        feedstock = system.ins[0]
        feedstock.register_alias('feedstock', override=True, safe=False)
        system.flowsheet.M401.ID = 400
        if scenario.preprocessing:
            self.areas = [
                'MRF',
                'STRAP',
                'EtOH',
                'WWT',
                'Facilities',
            ]
        else:
            self.areas = [
                'STRAP',
                'EtOH',
                'WWT',
                'Facilities',
            ]
        self.unit_groups = unit_groups = bst.UnitGroup.group_by_area(system.units)
        for i, j in zip(unit_groups, self.areas): i.name = j
        # for HXN_group in unit_groups:
        #     if HXN_group.name == 'HXN':
        #         HXN_group.filter_savings = False # Allow negative values in heat utilities
        #         HXN = HXN_group.units[0]
        #         try:
        #             assert isinstance(HXN, bst.HeatExchangerNetwork)
        #         except: breakpoint()
        # self.HXN = HXN
        # HXN.acceptable_energy_balance_error = 0.02
        # HXN.replace_unit_heat_utilities = False
        # HXN.force_ideal_thermo = True
        # HXN.cache_network = True
        # HXN.avoid_recycle = True
        self.tea = create_baseline_tea(system, years=30)
        CO2_bulk_plastic = self.chemicals.BulkPlastic.atoms['C'] * 44.01
        self.direct_nonbiogenic_emissions = lambda: self.lignin.imol['BulkPlastic'] * CO2_bulk_plastic + self.natural_gas.F_mol * 44.01
        system.define_process_impact(
            key=GWP_key,
            name='Direct non-biogenic emissions',
            basis='kg',
            inventory=self.direct_nonbiogenic_emissions,
            CF=1.,
        )
        system.set_tolerance(
            rmol=1e-6, mol=1e-6, maxiter=200, subsystems=True
        )
        
        @system.add_specification()
        def adjust_MSW_composition():
            if not scenario.preprocessing: return
            feedstock = self.feedstock
            total_flow = feedstock.F_mass
            non_residues_fraction = (
                self.hand_sorted_content 
                + self.metals_content
                + self.overs_content
                + self.unders_content
            )
            residues_fraction = 100 - non_residues_fraction
            feedstock.imass['NonResidue'] = 0 
            feedstock.F_mass = residues_fraction * total_flow / 100
            feedstock.imass['HandSorted', 'Metal', 'Overs', 'Unders'] = total_flow * np.array([
                self.hand_sorted_content,
                self.metals_content,
                self.overs_content,
                self.unders_content,
            ]) / 100
        
        minimum_fraction_burned = 0
        maximum_fraction_burned = 0.99
        @system.add_bounded_numerical_specification(
            x0=minimum_fraction_burned, x1=maximum_fraction_burned, 
            xtol=1e-6, ytol=100, maxiter=300,
        )
        def adjust_fraction_to_boiler(fraction_burned):
            # Returns energy consumption at given fraction processed (not sent to boiler).
            splitter = self.incineration_splitter
            splitter.split[:] = 1 - fraction_burned
            system.simulate(design_and_cost=True)
            excess = self.BT._excess_electricity_without_fuel
            if fraction_burned == minimum_fraction_burned and excess > 0:
                splitter.neglect_natural_gas_streams = False # No need to neglect
                return 0 # No need to burn bagasse
            elif fraction_burned == maximum_fraction_burned and excess < 0: 
                splitter.neglect_natural_gas_streams = False # Cannot be neglected
                return 0 # Cannot satisfy energy demand even at 99% sent to boiler (or minimum fraction processed)
            else:
                splitter.neglect_natural_gas_streams = True
                return excess
        
        self.adjust_fraction_to_boiler = adjust_fraction_to_boiler
        
        @system.add_specification()
        def assume_negligible_natural_gas_streams():
            splitter = self.incineration_splitter
            if splitter.neglect_natural_gas_streams: self.natural_gas.empty()
        
        return system
    
    def create_model(self):
        system = self.system
        self.dissolution_vessel.minimum_solvent_to_feed = 2.0
        bst.settings.register_credit('Ethanol RIN D3', dist.mean_RIN_D3_price * ethanol_L_per_kg)
        ethanol_storage = self.ethanol.source
        ethanol_storage.define_credit('Ethanol RIN D3', self.ethanol)
        bst.PowerUtility.set_CF(
            GWP_key, # [kg*CO2*eq / kWhr] From GREET; NG-Fired Simple-Cycle Gas Turbine CHP Plant 
            0.36
        )
        self.natural_gas.set_CF(
            GWP_key,
            0.33, # Natural gas from shell conventional recovery, GREET; includes non-biogenic emissions
        )
        self.solvent.set_CF(
            GWP_key,
            0.8199, # GREET; Mixed xylenes production from catalytic reforming of naphtha
        )
        
        # Set non-negligible characterization factors
        stream_names = ('FGD_lime', 'cellulase', 'urea', 'caustic', 
                        'HCl', 'NaOH', 'denaturant', 'lime', 'H3PO4')
        
        get_stream = lambda ID: getattr(self, ID) if hasattr(self, ID) else bst.MockStream(ID) # Retrieve mock-stream in case stream does not exist in configuration
        (FGD_lime, cellulase, urea, caustic, 
         HCl, NaOH, denaturant, lime, H3PO4) = [get_stream(i) for i in stream_names]
        
        def get_renamed_stream(search, ID):
            if hasattr(self, search):
                stream = getattr(self, search)
                stream.ID = ID
            else:
                stream = bst.MockStream(ID)
            return stream
        
        NaOCl = get_renamed_stream('naocl_R602', 'NaOCl')
        citric_acid = get_renamed_stream('citric_R602', 'citric_acid')
        bisulfite = get_renamed_stream('bisulfite_R602', 'bisulfite')
        polishing_filter_air = get_renamed_stream('air_R603', 'polishing_filter_air')
        polishing_filter_vent = get_renamed_stream('vent_R603', 'polishing_filter_vent')
        set_GWPCF(NaOCl, 'NaOCl')
        set_GWPCF(citric_acid, 'citric acid')
        set_GWPCF(bisulfite, 'bisulfite')
        set_GWPCF(H3PO4, 'H3PO4')
        set_GWPCF(lime, 'lime', dilution=0.046) # Diluted with water
        set_GWPCF(denaturant, 'gasoline')
        set_GWPCF(FGD_lime, 'lime', dilution=0.451)
        set_GWPCF(cellulase, 'cellulase', dilution=0.05) 
        set_GWPCF(urea, 'urea')
        set_GWPCF(caustic, 'NaOH', 0.5)
        set_GWPCF(HCl, 'HCl')
        set_GWPCF(NaOH, 'NaOH')
        scenario = self.scenario
        dissolution_step = self.dissolution_step
        precipitation_step = self.precipitation_step
        model = bst.Model(system)
        self.key_fermentation_parameters = []
        self.key_contextual_parameters = []
        self.nonkey_parameters =[]
        self.tea_parameters = []
        self.lca_parameters = []
        self.general_parameters = []
        
        def fermentation(p):
            self.key_fermentation_parameters.append(p)
            self.nonkey_parameters.remove(p)
            return p
        
        def contextual(p):
            self.key_contextual_parameters.append(p)
            self.nonkey_parameters.remove(p)
            return p
        
        def metric(f=None, *args, **kwargs):
            if f is None: return lambda f: metric(f, *args, **kwargs) 
            return model.metric(f, *args, **kwargs)
            
        def parameter(f=None, *args, group=None, **kwargs):
            if f is None:
                return lambda f: parameter(f, *args, group=group, **kwargs)
            else:
                p = model.parameter(f, *args, **kwargs)
                if group is None: group = 'general'
                general = self.general_parameters
                tea = self.tea_parameters
                lca = self.lca_parameters
                if group == 'general': 
                    general.append(p)
                    tea.append(p)
                    lca.append(p)
                elif group == 'tea': 
                    tea.append(p)
                elif group == 'lca': 
                    lca.append(p)
                elif group == 'none': pass
                else: raise ValueError('invalind parameter category')
                self.nonkey_parameters.append(p)
                return p
        
        def uniform(lb, ub, *args, **kwargs):
            return parameter(*args, distribution=shape.Uniform(lb, ub), bounds=(lb, ub), **kwargs)
        
        def default(baseline, *args, **kwargs):
            lb = 0.75*baseline
            ub = 1.25*baseline
            if baseline < 0: lb, ub = ub, lb
            return parameter(*args, distribution=shape.Uniform(lb, ub), bounds=(lb, ub),
                             baseline=baseline, **kwargs)
        
        def default_gwp(baseline, *args, **kwargs):
            lb = 0.90*baseline
            ub = 1.10*baseline
            return parameter(*args, distribution=shape.Uniform(lb, ub), bounds=(lb, ub), 
                             baseline=baseline, **kwargs)
        
        def triangular(lb, mid, ub, *args, **kwargs):
            return parameter(*args, distribution=shape.Triangle(lb, mid, ub), bounds=(lb, ub), **kwargs)
        
        def gaussian(mean, std, bounds=None, *args, **kwargs):
            normal = shape.Normal(mean, std)
            lb = mean - 2 * std
            ub = mean + 2 * std
            distribution = shape.Trunc(normal, lower=lb, upper=ub)
            return parameter(
                *args, distribution=distribution, 
                bounds=(lb, ub) if bounds is None else bounds, baseline=mean, **kwargs
            )
        
        def constant(constant, *args, **kwargs):
            lb = constant
            ub = constant + (constant * 1e-12 if constant > 0 else + 1e-12)
            return parameter(*args, 
                distribution=shape.Uniform(lb, ub), 
                bounds=(lb, ub), 
                baseline=constant, 
                group='none', 
                **kwargs
            )
        
        feedstock_flow = lambda: system.operating_hours * self.feedstock.F_mass / kg_per_MT # MT / y
        ethanol_flow = lambda: system.operating_hours * self.ethanol.F_mass * ethanol_L_per_kg # L / y
        resin_flow = lambda: system.operating_hours * self.resin.F_mass # kg / y
        electricity = lambda: system.operating_hours * system.power_utility.rate
        
        @metric(units='kg/MT')
        def resin_production():
            return resin_flow() / feedstock_flow()
        
        @metric(units='L/MT')
        def ethanol_production():
            return ethanol_flow() / feedstock_flow()
        
        @metric(units='kWh/MT')
        def electricity_production():
            return -electricity() / feedstock_flow()
           
        @metric(name='Carbon intensity', units='kg*CO2e/GGE')
        def GWP_energy():
            return system.get_property_allocated_impact(
                key=GWP_key, name='energy', basis='GGE',
                products=[self.resin, self.ethanol]
            ) # kg-CO2e / kJ
        
        @metric(units='kg*CO2e/kWh', name='Carbon intensity', element='electricity')
        def GWP_electricity():
            GWP = GWP_energy()
            return bst.units_of_measure.convert(GWP, '1/GGE', '1/kWh')
        
        @metric(units='kg*CO2e/L', name='Carbon intensity', element='ethanol')
        def GWP_ethanol():
            GWP = GWP_energy()
            return GWP / 1.5 / L_per_gal # per liter of ethanol
        
        @metric(units='kg*CO2e/kg', name='Carbon intensity', element='polymer resin')
        def GWP_polymer_resin():
            GWP = GWP_energy()
            mass_flow = self.resin.F_mass # kg / hr
            energy_flow = self.resin.get_property('LHV', 'GGE/hr')
            GGE_per_kg = energy_flow / mass_flow
            return GWP * GGE_per_kg
        
        self.GWP_resin = GWP_polymer_resin
        
        @metric(units='%')
        def IRR():
            try:
                return 100 * self.tea.solve_IRR(bounds=(0, 1))
            except:
                return 0
        
        @metric(units='10^6 USD')
        def TCI(): 
            return self.tea.TCI / 1e6
        
        @metric(units='USD/kg')
        def MSP(): # Price of Pyrolysis oil- https://www.sciencedirect.com/science/article/pii/S0360544218303220?casa_token=3dIJ06zwXScAAAAA:0w1izrLvVhSg1b1ZoY3h1F7PbyFBnjeTj-wqXIDtmr_EZ3BSzcV7sNcQkuYG8X4wHLlsipHAow#tbl2fnb
            return self.tea.solve_price(self.resin)
        
        location = scenario.tipping_fee_location
        if location[0] == '*':
            tipping_fee = constant(
                self.get_tipping_fee(location[1:])[0], 
                element='MSW',
                units='USD/ton',
            )
        else:
            tipping_fee = gaussian(
                *self.get_tipping_fee(location), 
                element='MSW',
                units='USD/ton',
                group='tea',
                hook=lambda x: max(x, 0)
            )
            
        baseline = scenario.biogenic_material_percent
        @uniform(
            baseline - 5, baseline + 5,
            baseline=baseline,
            element='MSW',
            units='wt %',
        )
        def set_biogenic_content(biogenic_content):
            biogenic_content /= 100
            total_flow = self.feedstock.F_mass - self.feedstock.imass['NonResidue']
            self.feedstock.imass['BiogenicMaterial'] = total_flow * biogenic_content
            self.feedstock.imass['Plastic'] = total_flow * (1. - biogenic_content)
        
        if scenario.preprocessing:
            @contextual
            @tipping_fee
            def set_MSW_tipping_fee(tipping_fee):
                self.MSW_tipping_fee = tipping_fee
                price = -tipping_fee / 907.185
                self.feedstock.price = price + self.manual_sorting_cost / 1000
                self.overs.price = self.unders.price = price
            
            self.manual_sorting_cost = 33
            
            # From jupyter notebook for preprocessing
            @contextual
            @parameter(baseline=33, bounds=(33, 44), element='Preprocess', units='USD/MT')
            def set_MSW_manual_sorting_cost(manual_sorting_cost):
                self.manual_sorting_cost = manual_sorting_cost
                self.feedstock.price = -self.MSW_tipping_fee / 907.185 + manual_sorting_cost / 1000
            
            baseline = scenario.hand_sorted_percent
            @parameter(
                bounds=(baseline * 0.5, baseline * 1.5),
                baseline=baseline,
                distribution='triangular',
                element='MSW',
                units='wt %',
            )
            def set_hand_sorted_content(hand_sorted):
                self.hand_sorted_content = hand_sorted
            
            baseline = scenario.metals_percent
            @parameter(
                bounds=(baseline * 0.5, baseline * 1.5),
                baseline=baseline,
                distribution='triangular',
                element='MSW',
                units='wt %',
            )
            def set_metals_content(metals):
                self.metals_content = metals
            
            baseline = scenario.unders_percent
            @parameter(
                bounds=(baseline * 0.1, baseline * 1.1),
                baseline=baseline,
                distribution='triangular',
                element='MSW',
                units='wt %',
            )
            def set_unders_content(unders):
                self.unders_content = unders
            
            baseline = scenario.overs_percent
            @parameter(
                bounds=(baseline * 0.05, baseline * 1.05),
                baseline=baseline,
                distribution='triangular',
                element='MSW',
                units='wt %',
            )
            def set_overs_content(overs):
                self.overs_content = overs
            
            # TODO: check assumption with Srikar
            @contextual
            @parameter(baseline=0.0435, bounds=(0, 0.087), element='Preprocess', units='USD/ton')
            def set_hand_sorted_material_price(hand_sorted_material_price):
                self.sorted_material.price = hand_sorted_material_price / 907.185
                
            # TODO: the ‘eddy current- non ferrous metals’ at $0.06-0.15/lb
            # the magnet separator ‘ferrous metals’ at  $0.0394/lb.
            # source: https://www.scrapmetalbuyers.com/current-prices)
            @contextual
            @parameter(bounds=(0.132, 0.33), element='Preprocess', units='USD/kg')
            def set_nonferrous_metal_price(nonferrous_metal_price):
                self.nonferrous_metal.price = nonferrous_metal_price
            
            @contextual
            @parameter(baseline=0.087, bounds=(0.087 * 0.5, 0.087 * 1.15), element='Preprocess', units='USD/kg')
            def set_ferrous_metal_price(ferrous_metal_price):
                self.ferrous_metal.price = ferrous_metal_price
            
        else:
            @contextual
            @tipping_fee
            def set_MSW_tipping_fee(tipping_fee):
                self.feedstock.price = -tipping_fee / 907.185
        
        @uniform(
            0.001, 0.01,
            baseline=0.005,
            element='MSW',
            units='wt % plastic',
        )
        def set_solute_content(solute_content):
            solute_content = solute_content / 100
            plastics = self.feedstock.imass['Plastic']
            self.feedstock.imass['Solutes'] = 0
            self.feedstock.imass['Plastic'] = plastics * (1 - solute_content)
            self.feedstock.imass['Solutes'] = plastics * solute_content
        
        # @uniform(
        #     1, 3,
        #     baseline=1,
        #     element='Feedstock',
        #     units='hr',
        # )
        # def set_cycle_time(cycle_time):
        #     self.adsorption_column.t_cycle = cycle_time
            
        # @uniform(
        #     3.6, 14.4,
        #     baseline=14.4,
        #     element='Feedstock',
        #     units='m/hr', 
        # )
        # def set_superficial_velocity(superficial_velocity):
        #     self.adsorption_column.superficial_velocity = superficial_velocity
        
        # @uniform(
        #     1100, 1400,
        #     baseline=1270,
        #     element='Feedstock',
        #     units='L/g',
        # )
        # def set_K_Langmuir(K_Langmuir):
        #     self.adsorption_column.KL = K_Langmuir
            
        # @uniform(
        #     2, 18,
        #     baseline=16.6,
        #     element='Feedstock',
        #     units='1/hr',
        # )
        # def set_k_MT(k_MT):
        #     self.adsorption_column.k = k_MT
        
        baseline, theoretical = self.MRF_residue_availabilities.values()
        match scenario.processing_capacity:
            case 'full range':
                processing_capacity = uniform(
                    baseline, theoretical,
                    element='MSW',
                    units='ton/yr',
                    baseline=baseline,
                )
            case 'baseline':
                processing_capacity = constant(
                    baseline, 
                    element='MSW',
                    units='ton/yr',
                )
            case 'theoretical':
                processing_capacity = constant(
                    theoretical, 
                    element='MSW',
                    units='ton/yr',
                )
        
        @contextual
        @processing_capacity
        def set_processing_capacity(processing_capacity):
            feedstock = self.feedstock
            self.system.rescale(
                feedstock, 
                processing_capacity * 907.185 / self.tea.operating_hours / feedstock.F_mass
            )
            
        match scenario.resin_product_price:
            case 'full range':
                resin_product_price = parameter(
                    bounds=self.resin_price_range, 
                    element='resin',
                    units='USD/ton',
                )
            case 'HDPE':
                resin_product_price = constant(
                    self.baseline_resin_price, 
                    element='resin',
                    units='USD/ton',
                )
            case 'none':
                resin_product_price = constant(
                    0, 
                    element='resin',
                    units='USD/ton',
                )
            case 'MRF recycled HDPE':
                resin_product_price = constant(
                    recycled_HDPE_price_USD_per_ton, 
                    element='resin',
                    units='USD/ton',
                )
            case 'MRF recycled PP':
                resin_product_price = constant(
                    recycled_PP_price_USD_per_ton, 
                    element='resin',
                    units='USD/ton',
                )
        
        @resin_product_price
        def set_resin_product_price(price):
            self.resin.price = price / 907.185
            
        baseline = 2.17
        @default(
            baseline=baseline,
            element='Solvent',
            units='USD/kg',
            group='tea',
        )
        def set_solvent_price(price):
            self.solvent.price = price
        
        # @parameter(
        #     baseline=0.10,
        #     units='%',
        #     distribution=shape.Uniform(0.1, 0.2), 
        #     element='specification',
        #     group='tea',
        # )
        # def set_IRR(IRR):
        #     self.tea.IRR = IRR
        self.tea.IRR = 0.1 # Consistent with other waste reducing processes (NREL).
        
        baseline = scenario.target_plastic_percent / (100 - scenario.biogenic_material_percent)
        @parameter(
            baseline=baseline,
            element='feedstock',
            distribution=shape.Uniform(baseline - 0.1, baseline + 0.1),
            units='wt % plastic',
        )
        def set_extracted_polymer_fraction(extracted_polymer_fraction):
            s = self.feedstock
            F_mass = s.imass['Plastic']
            plastic = self.dissolution_step.plastic
            s.imass['BulkPlastic'] = F_mass * (1 - extracted_polymer_fraction)
            s.imass[plastic] = extracted_polymer_fraction * F_mass
        
        def solvent_content(baseline, *args, **kwargs):
            bounds = (baseline - 10, baseline + 10)
            return parameter(*args, bounds=bounds, baseline=baseline,
                             distribution=shape.Uniform(*bounds), 
                             units='%', **kwargs)
        
        @uniform(
            100 * 0.4, 100 * 0.6,
            element='centrifuged plastic',
        )
        def set_centrifuged_plastic_solvent_content(solvent_content):
            dissolution_step.solvent_content = solvent_content / 100
        
        @solvent_content(
            100 * precipitation_step.precipitate_solvent_content,
            element='precipitate',
        )
        def set_precipitate_solvent_content(solvent_content):
            precipitation_step.precipitate_solvent_content = solvent_content / 100
        
        @parameter(
            element='Extracted polymer to solvent', 
            units='wt % solvent',
            distribution='triangular', bounds=(2, 10),
            baseline=5,
        )
        def set_extracted_polymer_to_solvent_ratio(ratio):
            dissolution_step.capacity = ratio / 100
        
        # USDA ERS historical price data without EPA RIN prices
        @parameter(distribution=dist.ethanol_no_RIN_price_distribution, element='Ethanol', 
                   baseline=dist.mean_ethanol_no_RIN_price, units='USD/L',
                   group='tea')
        def set_ethanol_price(price): # Triangular distribution fitted over the past 10 years Sep 2009 to Nov 2020
            self.ethanol.price = price * ethanol_L_per_kg
        
        @parameter(distribution=dist.RIN_D3_price_distribution, element='RIN D3', 
                   baseline=dist.mean_RIN_D3_price, units='USD/RIN',
                   group='tea')
        def set_RIN_D3_price(price): # Triangular distribution fitted over the past 10 years Sep 2009 to Nov 2020
            bst.settings.register_credit('Ethanol RIN D3', price * ethanol_L_per_kg)
        
        # set_GWPCF(NaOCl, 'NaOCl')
        # set_GWPCF(citric_acid, 'citric acid')
        # set_GWPCF(bisulfite, 'bisulfite')
        # set_GWPCF(lime, 'lime', dilution=0.046) # Diluted with water
        # set_GWPCF(denaturant, 'gasoline')
        # set_GWPCF(FGD_lime, 'lime', dilution=0.451)
        # set_GWPCF(cellulase, 'cellulase', dilution=0.05) 
        # set_GWPCF(urea, 'urea')
        # set_GWPCF(caustic, 'NaOH', 0.5)
        # set_GWPCF(HCl, 'HCl')
        # set_GWPCF(NaOH, 'NaOH')
        
        # DO NOT DELETE:
        # natural_gas.phase = 'g'
        # natural_gas.set_property('T', 60, 'degF')
        # natural_gas.set_property('P', 14.73, 'psi')
        # original_value = natural_gas.imol['CH4']
        # natural_gas.imass['CH4'] = 1 
        # V_ng = natural_gas.get_total_flow('m3/hr')
        # natural_gas.imol['CH4'] = original_value
        V_ng = 1.473318463076884 # Natural gas volume at 60 F and 14.73 psi [m3 / kg]
        
        # https://www.eia.gov/energyexplained/natural-gas/prices.php
        @parameter(distribution=dist.natural_gas_price_distribution, element='Natural gas', units='USD/m3',
                    baseline=4.73 * 35.3146667/1e3, group='tea')
        def set_natural_gas_price(price): 
            self.BT.natural_gas_price = price * V_ng
        
        @parameter(distribution=dist.electricity_price_distribution, units='USD/kWh',
                   element='electricity', baseline=dist.mean_electricity_price,
                   group='tea')
        def set_electricity_price(price): 
            bst.PowerUtility.price = price
        
        @default(0.212 * 150 / 20, units='USD/kg', element='cellulase', group='tea')
        def set_cellulase_price(price):
            self.cellulase.price = price
        
        glucan_fraction_in_biogenic = 0.282 # TODO: Update once real fraction is known
        biogenic_content = scenario.biogenic_material_percent
        glucan_fraction = glucan_fraction_in_biogenic * biogenic_content / (biogenic_content + scenario.target_plastic_percent)
        self.cellulase_mixer.enzyme_loading = 0.02668 # Simply assume 4x the amount use in NREL (experimental conditions are too unoptimized)
        # self.cellulase_mixer.enzyme_loading = 0.05 * 150 / 20 / glucan_fraction #  5 % v/w-MSW at 150g/L to % w of Glucan at 20 g/L
        self.cellulase_mixer.enzyme_concentration = 0.150 # in experiments, 150  g cellulase / 1000g cellulose mixture at industrial scale
        MSW = self.cellulase_mixer.ins[0]
        self.cellulase_mixer.loading_basis = lambda: MSW.imass['Glucan']
        experimental_cellulase_loading = 100 * self.cellulase_mixer.enzyme_loading
        NREL_cellulase_loading = 100 * 0.00667
        @default_gwp(CFs['cellulase'], name='GWP', 
                     element='Cellulase', units='kg*CO2e/kg', group='lca')
        def set_cellulase_GWP(value):
            self.cellulase.characterization_factors['GWP'] = value * 0.150
            
        match scenario.fermentation_performance:
            case 'full range': # Lower bound experimental data, upper bound NREL
                def f(name):
                    mean, std = getattr(self, f'experimental_{name}')
                    ub = getattr(self, f'NREL_{name}')
                    lb = mean - 2 * std
                    if ub < lb: lb, ub = ub, lb
                    return dict(lb=lb, ub=ub, baseline=mean)
                
                cellulase_loading = uniform(NREL_cellulase_loading, experimental_cellulase_loading, units='wt % cellulose', element='cellulase')
                glucose_yield = uniform(**f('glucose_yield'), units='%', element='Saccharification')
                xylose_yield = uniform(**f('xylose_yield'), units='%', element='Saccharification')
                sugars_to_ethanol_yield = uniform(**f('sugars_to_ethanol_yield'), units='%', element='Cofermenation')
                cofermentation_ethanol_productivity = uniform(**f('cofermentation_ethanol_productivity'), units='g/L/h', element='Cofermentation')
                cofermentation_ethanol_titer = uniform(**f('cofermentation_ethanol_titer'), units='g/L', element='Cofermentation')
                solids_loading = uniform(0.09524, 0.2, units='wt % solids', element='Saccharification')
            case 'experimental': # Based on experimental data
                cellulase_loading = constant(experimental_cellulase_loading, units='wt % cellulose', element='cellulase')
                glucose_yield = gaussian(*self.experimental_glucose_yield, units='%', element='Saccharification', bounds=[54.7 - 2 * 13.2, 95])
                xylose_yield = gaussian(*self.experimental_xylose_yield, units='%', element='Saccharification', bounds=[65.8 - 2 * 3.9, 95])
                # TODO: Ask experimental collaborators numbers for xylose and glucose to ethanol (not just overall)
                sugars_to_ethanol_yield = gaussian(*self.experimental_sugars_to_ethanol_yield, units='%', element='Cofermenation', bounds=[71.1 - 2 * 2.1, 85])
                # TODO: Ask experimental collaborators
                cofermentation_ethanol_productivity = gaussian(*self.experimental_cofermentation_ethanol_productivity, units='g/L/h', element='Cofermentation', bounds=[0.1722 - 2 * 0.0192, 1.5])
                cofermentation_ethanol_titer = gaussian(*self.experimental_cofermentation_ethanol_titer, units='g/L', element='Cofermentation', bounds=[8.26514 - 2 * 0.922, 60])
                solids_loading = constant(0.09524, units='wt % solids', element='Saccharification')
            case 'NREL': # Based on NREL assumptions
                cellulase_loading = constant(NREL_cellulase_loading, units='wt % cellulose', element='cellulase')
                glucose_yield = constant(self.NREL_glucose_yield, units='%', element='Saccharification')
                xylose_yield = constant(self.NREL_xylose_yield, units='%', element='Saccharification')
                sugars_to_ethanol_yield = constant(self.NREL_sugars_to_ethanol_yield, units='%', element='Cofermenation')
                cofermentation_ethanol_productivity = constant(self.NREL_cofermentation_ethanol_productivity, units='g/L/h', element='Cofermentation')
                cofermentation_ethanol_titer = constant(self.NREL_cofermentation_ethanol_titer, units='g/L', element='Cofermentation')
                solids_loading = constant(0.2, units='wt % solids', element='Saccharification')
            case 'prospective': # Based on NREL assumptions
                f = lambda name: getattr(self, f'prospective_{name}')
                cellulase_loading = constant(f('cellulase_loading'), units='wt % cellulose', element='cellulase')
                glucose_yield = constant(f('glucose_yield'), units='%', element='Saccharification')
                xylose_yield = constant(f('xylose_yield'), units='%', element='Saccharification')
                sugars_to_ethanol_yield = constant(f('sugars_to_ethanol_yield'), units='%', element='Cofermenation')
                cofermentation_ethanol_productivity = constant(f('cofermentation_ethanol_productivity'), units='g/L/h', element='Cofermentation')
                cofermentation_ethanol_titer = constant(f('cofermentation_ethanol_titer'), units='g/L', element='Cofermentation')
                solids_loading = constant(0.2, units='wt % solids', element='Saccharification')
        
        @cellulase_loading
        def set_cellulase_loading(cellulase_loading):
            self.cellulase_mixer.enzyme_loading = cellulase_loading / 100
        
        @solids_loading
        def set_solids_loading(solids_loading):
            self.cellulase_mixer.solids_loading = self.cellulase_mixer.insoluble_solids_loading = solids_loading
            
        @fermentation
        @glucose_yield
        def set_glucose_yield(glucose_yield):
            self.hydrolysis.reactions[0].X = glucose_yield / 100
        
        @fermentation
        @xylose_yield
        def set_xylose_yield(xylose_yield):
            self.hydrolysis.reactions[2].X = xylose_yield / 100
         
        # TODO: Ask experimental collaborators numbers for xylose and glucose to ethanol (not just overall)
        @fermentation
        @sugars_to_ethanol_yield
        def set_sugars_to_ethanol_yield(sugars_to_ethanol_yield):
            seed_train = self.seed_train
            seed_splitter = self.seed_splitter
            fermentor = self.fermentor
            sugars_to_ethanol_yield *= 0.01
            split = seed_splitter.split.mean()
            
            # Xylose
            X1 = split * seed_train.reactions.X[1]
            X2 = split * seed_train.reactions.X[3]
            X3 = (sugars_to_ethanol_yield - X1) / (1 - X1 - X2)
            X_excess = X3 * 1.0526 - 1
            if X_excess > 0.: X3 = 0.95 - 1e-16 # Maximum
            fermentor.cofermentation.X[1] = X3
            fermentor.cofermentation.X[3] = X3 * 0.0526 # 95% towards ethanol, the other 5% goes towards cell mass
        
            # Glucose
            X1 = split * seed_train.reactions.X[0]
            X2 = split * seed_train.reactions.X[2]
            X3 = (sugars_to_ethanol_yield - X1) / (1 - X1 - X2)
            X_excess = X3 * 1.0526 - 1
            if X_excess > 0.: X3 = 0.95 - 1e-16 # Maximum
            fermentor.cofermentation.X[0] = X3
            fermentor.cofermentation.X[2] = X3 * 0.0526 # 95% towards ethanol, the other 5% goes towards cell mass    
        
        @fermentation
        @cofermentation_ethanol_productivity
        def set_cofermentation_ethanol_productivity(ethanol_productivity):
            self.fermentor.productivity = ethanol_productivity
        
        @fermentation
        @cofermentation_ethanol_titer
        def set_cofermentation_ethanol_titer(ethanol_titer):
            self.fermentor.titer = ethanol_titer
        
        return model

    @classmethod
    def get_distribution_tables(cls):
        pm_all = cls(scenario='all', simulate=False)
        tables = {
            'Full problem space': strap.get_distributions(pm_all.parameters, save=False)
        }
        # Key contextual table
        pm_worst = cls(scenario='worst', simulate=False)
        tables['Contextual scenario subspaces'] = table = strap.get_distributions(pm_worst.key_contextual_parameters, save=False)
        table.rename(columns={'Baseline': 'Worst'}, inplace=True)
        table.drop(columns={'Distribution'}, inplace=True)
        for scenario in ('baseline', 'best'):
            pm = cls(scenario=scenario, simulate=False)
            baseline = strap.get_distributions(pm.key_contextual_parameters, save=False)['Baseline']
            table[scenario.capitalize()] = baseline
        # Fermentation performance table
        pm_experimental = cls(scenario='experimental', simulate=False)
        tables['Fermentation scenario subspaces'] = table = strap.get_distributions(pm_experimental.key_fermentation_parameters, save=False)
        table.rename(columns={'Distribution': 'Experimental'}, inplace=True)
        table.drop(columns={'Baseline'}, inplace=True)
        for scenario in ('prospective', 'NREL'):
            pm = cls(scenario=scenario, simulate=False)
            baseline = strap.get_distributions(pm.key_fermentation_parameters, save=False)['Baseline']
            table[scenario if scenario.isupper() else scenario.capitalize()] = baseline
        results_folder = strap.simulation.results_folder
        for name, table in tables.items():
            file = os.path.join(results_folder, name.replace(' ', '_') + '.xlsx')
            table.to_excel(file)
            table.index.name = '#'
        return tables
            
    def prospective_fermentation_performance(self, overall=True, update=True):
        from warnings import filterwarnings
        filterwarnings('ignore')
        names = (
            'glucose_yield',
            'xylose_yield',
            'xylose_to_ethanol_yield',
            'glucose_to_ethanol_yield',
            'cofermentation_ethanol_productivity',
            'cofermentation_ethanol_titer',
        )
        parameters = [getattr(self, f'set_{i}') for i in names]
        setters = [i.setter for i in parameters]
        experimental_values = np.array([getattr(self, f'experimental_{i}')[0] for i in names])
        NREL_values = np.array([getattr(self, f'NREL_{i}') for i in names])
        
        if overall:
            def f(x):
                for f, baseline, NREL in zip(setters, experimental_values, NREL_values):
                    value = baseline + x * (NREL - baseline)
                    f(value) 
                self.system.simulate()
                return 100 * self.tea.IRR - self.IRR()
            
            prospective_improvements = flx.IQ_interpolation(f, 0, 1, xtol=1e-3)
        else:
            def objective(xs):
                return (xs * xs).sum()
            
            experimental_stds = [getattr(self, f'experimental_{i}')[1] for i in names]
            lbs = [mean - 2 * std for mean, std in zip(experimental_values, experimental_stds)]
            def constraint(xs):
                for f, x, baseline, NREL, lb in zip(setters, xs, experimental_values, NREL_values, lbs):
                    value = baseline + x * (NREL - baseline)
                    if value < lb: value = lb
                    f(value) 
                self.system.simulate()
                return 100 * self.tea.IRR - self.IRR()
            
            N = len(setters)
            result = minimize(
                objective, np.zeros(N), 
                method='COBYQA', bounds=N * [(0, 1)], 
                constraints=[{'type': 'eq', 'fun': constraint}],
                tol=1e-2
            )
            prospective_improvements = result.x
        
        prospectives = experimental_values + prospective_improvements * (NREL_values - experimental_values)
        prospectives = {name: prospective for name, prospective in zip(names, prospectives)}
        if update:
            for name, value in prospectives.items(): setattr(self, f'prospective_{name}', value)
        if overall:
            return prospectives
        else:
            return result, prospectives

    def preprocessing_production_cost(self):
        units = [*self.unit_groups[0].units]
        installed_equipment_cost = sum([i.installed_cost for i in units])
        depreciation = installed_equipment_cost / 10 / (330 * 24)
        utility_cost = sum([i.utility_cost for i in units if i.utility_cost])
        manual_sorting_cost = self.MRF_residues.F_mass * self.manual_sorting_cost / 1000
        MSW_residue = self.MRF_residues.F_mass
        production_cost = (manual_sorting_cost + utility_cost + depreciation) / MSW_residue
        columns = ('Value', 'Units')
        index = []
        values = []
        def add(a, b, value, units):
            index.append((a, b))
            values.append((value, units))
        add('MRF residues', '-', self.MRF_residues.F_mass * self.tea.operating_hours / 1e6, '10^3 MT/y') 
        add('preprocessing cost', 'manual sorting cost', 1000 * manual_sorting_cost / MSW_residue, 'USD/MT')
        add('preprocessing cost', 'depreciation', 1000 * depreciation / MSW_residue, 'USD/MT')
        add('preprocessing cost', 'utility cost', 1000 * utility_cost / MSW_residue, 'USD/MT')
        add('preprocessing cost', 'total', production_cost * 1000, 'USD/kg')
        import pandas as pd
        return pd.DataFrame(values, index=pd.MultiIndex.from_tuples(index), columns=columns)

    def get_production_cost_contribution(self):
        import biosteam as bst
        feedstock = self.feedstock
        name = self.scenario.name
        self.system.rescale(
            feedstock, 
            150 * 1e6 / self.tea.operating_hours / feedstock.F_mass
        ) # 150 10^3 MT/yr
        self.system.simulate()
        MSP = self.MSP()
        production = self.PEPP.F_mass * self.tea.operating_hours
        opex_contributions = {}
        VOC_table = self.tea.VOC_table(products=[self.PEPP])
        VOC_table = VOC_table['Cost [MM$/yr]'] * 1e6 / production
        materials = VOC_table['Raw materials']
        key_raw_materials = ('MRF residues', 'Cellulase')
        values = []
        def account(x, production=False):
            if production: x = -x
            values.append(x)
            return x
        for i in key_raw_materials: opex_contributions[i] = account(materials[i])
        # opex_contributions['FOC'] = self.tea.FOC / production
        # breakpoint()
        opex_contributions['Electricity production'] = account(VOC_table['Co-products & credits', 'Electricity production'], production=True)
        products = VOC_table['Co-products & credits']
        opex_contributions['Ethanol & RIN D3'] = account(products['Ethanol'] + products['Ethanol RIN D3'], production=True)
        VOC = VOC_table.loc['Raw materials'].sum() - VOC_table.loc['Other utilities & fees'].sum() - VOC_table.loc['Co-products & credits'].sum()
        opex_contributions['Other'] = VOC - sum(values)
        OPEX_contribution = sum(opex_contributions.values())
        CAPEX_contribution = MSP - OPEX_contribution
        for group in self.unit_groups: 
            try:
                group.autofill_metrics(
                    shorthand=False, 
                    installed_cost=True,
                    cooling_duty=False,
                    heating_duty=False,
                    electricity_consumption=False,
                    electricity_production=False,
                    material_cost=False
                )
            except:
                pass
        CAPEX_by_area = bst.UnitGroup.df_from_groups(
            self.unit_groups, fraction=True,
            scale_fractions_to_positive_values=True,
        )
        CAPEX_contributions = CAPEX_by_area.values * CAPEX_contribution
        capex_contributions = {}
        for i, j in zip(self.areas, CAPEX_contributions): 
            if 'EtOH' in i: i = 'EtOH prod.'
            if i in capex_contributions:
                capex_contributions[i] += float(j) / 100 
            else:
                capex_contributions[i] = float(j) / 100
            
        import pandas as pd
        import matplotlib.pyplot as plt
        from matplotlib.pyplot import colormaps
        import thermosteam as tmo
        bst.set_figure_size(aspect_ratio=1.2, width='half')
        tmo.utils.set_font(12)
        total = pd.DataFrame({'Total': [MSP]}, index=['Total'])
        CAPEX = pd.DataFrame({'CAPEX': [*capex_contributions.values()]}, index=[*capex_contributions.keys()])
        OPEX = pd.DataFrame({'OPEX': [*opex_contributions.values()]}, index=[*opex_contributions.keys()])
        df = pd.concat([total, CAPEX, OPEX], axis=1, sort=False)
        colors = {}
        colors['Total'] = 'black'
        area_names = [*self.area_colors.keys()]
        area_names[area_names.index('EtOH')] = 'EtOH prod.'
        other_names = [i for i in df.index if i not in area_names and i != 'Total']
        for x in [area_names, other_names]:
            n = 0
            for i in x: 
                if i not in colors:
                    colors[i] = colormaps['Set1'](n)
                    n += 1
        ax = df.T.plot.bar(stacked=True, rot=0)
        plt.ylabel(r'MSP [USD$\cdot$kg$^{-1}$]')
        plt.ylim(-2, 2)
        plt.axhline(y=self.baseline_resin_price, color='gray', linestyle='--')
        plt.axhline(y=0, color='k', linestyle='--')
        legend = ax.get_legend()
        legend.set_bbox_to_anchor((1, 1))
        plt.subplots_adjust(left=0.12, right=0.70)
        for i in ('svg', 'png'):
            plt.savefig(f'MSP_contributions_{name}.{i}', dpi=900, transparent=True)
        return opex_contributions, capex_contributions
        
        
    def plot_map(self, indicator):
        import plotly.express as px
        tipping_fees = self.tipping_fees_by_state
        if hasattr(self, indicator): indicator = getattr(self, indicator)
        locations = []
        values = []
        for location, (mean, std) in tipping_fees.items():
            if np.isnan(mean): continue
            self.set_MSW_tipping_fee(mean)
            locations.append(location)
            values.append(indicator())
        fig = px.choropleth(
            locations=locations,
            locationmode='USA-states',
            scope='usa',
            color=values,
            color_continuous_scale='viridis',
            title='STRAP-MSW profitability by state'
        )
        fig.add_scattergeo(
            locations = locations,
            locationmode='USA-states', 
            text = locations, # [f"{i}\n{int(j)}" for i, j in zip(locations, values)],
            mode = 'text',
        )
        fig.update_layout(coloraxis_colorbar_title_text = indicator.name_with_units)
        fig.show()


class BaselineSTRAPProcess(bst.ProcessModel):
    """
    Create a model for a solvent targeted precipitation and dissolution process.
    The dissolution and precipitation steps default to PE.
    
    Examples
    --------
    >>> from plastics.strap import BaselineSTRAPProcess
    >>> pm = BaselineSTRAPProcess(simulate=False)
    >>> pm.system.diagram(kind='cluster', number=True)
    >>> pm.system.simulate()
    >>> assumptions, results = pm.baseline()
    >>> assumptions
    Natural gas          Price [USD/m3]                0.167
    Feedstock            Processing capacity [MT/yr]   5e+03
                         Price [USD/kg]                 0.01
    -                    IRR [%]                        0.15
    Solvent              Price [USD/kg]                 2.17
    Polymer              Mass fraction                   0.5
    Centrifuged plastic  Solvent content [%]              50
    Plastic              Feedstock distance [km]         500
    Solvent              Solvent loss [%]                0.1
    Dissolution          Temperature [K]                 368
    Precipitation        Temperature [K]                 308
    Dissolution          Solvent capacity [wt %]           3
    dtype: float64
    
    >>> results
    -  GWP [kg*CO2e/kg]   1.54
       FFC [MJ/kg]        1.87
       MSP [USD/kg]       2.95
    dtype: float64
    
    """
    @bst.scenario
    class Scenario:
        solvent: str|tuple[str, ...] = '# Solvent used to separate the target plastic'
        target_plastic: str|tuple[str, ...] = '# The polymer layer being dissolved'
        target_plastic_percent: float|tuple[float, ...] = '# Fraction in feedstock [%]'
        processing_capacity: float = 5000, '# Feedstock flow rate [MT-plastic/yr]'
        sell_leftover_plastic: bool = False, '# Whether the MSP will include all products'
        burn_leftover_plastic: bool = True, '# Produce heat and power from leftover plastic'
        facilities: bool = True, '# On-site cooling tower, heat and power generation'
        precipitation_temperature_format: str = 'constant', "# Use 'drop' for % temperature drop to solvent melting point. Use 'constant' to set in Kelvin."
        precipitation_configuration: str = 'integrated heat transfer', "# Must be either 'solvent mixing' or 'integrated heat transfer'."
        turbogenerator: bool = True, '# On-site electricity generation'
        
        @property
        def multistep(self):
            return not isinstance(self.target_plastic, str)
        
        @property    
        def N_steps(self):
            return len(self.target_plastic) if self.multistep else 1

    @classmethod
    def get_scenarios(cls):
        return tuple(cls._scenarios.values())
    
    @classmethod
    def get_scenario(cls, scenario):
        return cls._scenarios[scenario] 
    
    @classmethod
    def set_scenario(cls, scenario):
        cls._scenarios[f'{scenario.target_plastic}/{scenario.solvent}'] = scenario
    
    _scenarios = [
        Scenario('THF', 'PC', 65, 250, False, False, False),
        Scenario('Toluene', 'PE', 50, 5000, False, True, True),
        Scenario(('Toluene', 'DMSOWater'), ('PE', 'EVOH'), (50, 3.22), 5000, False, True, True),
        Scenario('Xylene', 'PE', 90, 5e3, False, False, False),
        Scenario('Xylene', 'HDPE', 14, 325, False, False, False)
    ]
    _scenarios = {
        f'{i.target_plastic}/{i.solvent}': i for i in _scenarios
    }
    
    @property
    def name(self):
        scenario = self.scenario
        if scenario.multistep:
            name = f"{'_'.join(scenario.solvent)}_{'_'.join(scenario.target_plastic)}"
            if scenario.sell_leftover_plastic:
                name += "_sell_leftover_plastic"
            elif scenario.burn_leftover_plastic: 
                name += "_burn_leftover_plastic"
            if scenario.facilities: 
                name += "_with_facilities"
        else:
            name = f"{scenario.solvent}_{scenario.target_plastic}"
            if scenario.sell_leftover_plastic:
                name += "_sell_leftover_plastic"
            elif scenario.burn_leftover_plastic: 
                name += "_burn_leftover_plastic"
            if scenario.facilities: 
                name += "_with_facilities"
        return name
    
    @classmethod
    def default_scenario(cls):
        return cls._scenarios['PE/Toluene']
        
    @classmethod
    def as_scenario(cls, scenario):
        if isinstance(scenario, (str, tuple)):
            return cls._scenarios[scenario]
        else:
            raise TypeError('invalid scenario type')
    
    def create_thermo(self):
        scenario = self.scenario
        if scenario.multistep:
            for i, j in zip(scenario.target_plastic, scenario.solvent):
                default_plastic_solvent_pair(i, j)
        else:
            default_plastic_solvent_pair(scenario.target_plastic, scenario.solvent)
        chemicals = create_property_package()
        for i in solvent_mixtures: chemicals.define_group(*i)
        return chemicals
    
    def create_system(self):
        scenario = self.scenario
        chemicals = self.chemicals
        load_process_settings()
        if scenario.multistep:
            dissolution_step = tuple([
                getattr(strap.dissolution_steps, f'{i}_{j}_dissolution')()
                for i, j in zip(scenario.target_plastic, scenario.solvent)
            ])
            precipitation_step = tuple([
                getattr(strap.precipitation_steps, f'{i}_{j}_precipitation')()
                for i, j in zip(scenario.target_plastic, scenario.solvent)
            ])
        else:
            dissolution_step = getattr(strap.dissolution_steps, f'{scenario.target_plastic }_{scenario.solvent}_dissolution')()
            precipitation_step = getattr(strap.precipitation_steps, f'{scenario.target_plastic }_{scenario.solvent}_precipitation')()
        facilities = scenario.facilities
        target_plastic_percent = scenario.target_plastic_percent
        if scenario.multistep:
            bulk_plastic_percent = 100 - sum(target_plastic_percent) # PET
            chemicals.define_group('Plastic', [*scenario.target_plastic, 'BulkPlastic'], [*target_plastic_percent, bulk_plastic_percent], wt=True)
        else:
            bulk_plastic_percent = 100 - target_plastic_percent # PET
            chemicals.define_group('Plastic', [scenario.target_plastic, 'BulkPlastic'], [target_plastic_percent, bulk_plastic_percent], wt=True)
        chemicals.define_group('Solutes', ['Minerals', 'Solubles'], [0.8, 0.2], wt=True)
        if scenario.multistep:
            system = create_multilayer_batch_separation_system(
                dissolution_steps=dissolution_step, 
                precipitation_steps=precipitation_step,
                shred=True,
                facilities=scenario.facilities,
                core_facilities=True,
                turbogenerator=scenario.turbogenerator,
                precipitation_configuration=scenario.precipitation_configuration,
            )
            for i, step in zip(system.subsystems, dissolution_step):
                i.ID = f"{step.plastic}/{step.solvent}"
        else:
            system = create_single_layer_batch_separation_system(
                dissolution_step=dissolution_step,
                precipitation_step=precipitation_step,
                facilities=facilities,
                relative_molar_tolerance=1e-6, 
                molar_tolerance=1e-2,
                method='fixed-point',
                burn_leftover_plastic=scenario.burn_leftover_plastic,
                core_facilities=True,
                precipitation_configuration=scenario.precipitation_configuration,
                turbogenerator=scenario.turbogenerator,
            )
            for i, step in zip(system.subsystems, [dissolution_step]):
                i.ID = f"{step.plastic}/{step.solvent}"
        if scenario.multistep:
            self.dissolution_steps = dissolution_step
            self.precipitation_steps = precipitation_step
        else:
            self.dissolution_step = dissolution_step
            self.precipitation_step = precipitation_step
        self.tea = create_baseline_tea(system)
        self.direct_nonbiogenic_emissions = lambda: self.emissions.imass['CO2'] * system.operating_hours
        system.define_process_impact(
            key=GWP_key,
            name='Direct non-biogenic emissions',
            basis='kg',
            inventory=self.direct_nonbiogenic_emissions,
            CF=1.,
        )
        system.define_process_impact(
            key='FFC',
            name='Direct non-biogenic emissions',
            basis='kg',
            inventory=lambda: 0,
            CF=1.,
        )
        
        #add water consumption
        system.define_process_impact(
            key='WU',
            name='Direct non-biogenic emissions',
            basis='kg',
            inventory=lambda: 0,
            CF=1.,
        )
        
        #add human toxicity, cancer
        system.define_process_impact(
            key='HTC',
            name='Human Toxicity, cancer',
            basis='kg',
            inventory=lambda: 0,
            CF=1.,
        )
        
        #add human toxicity, non-cancer
        system.define_process_impact(
            key='HTNC',
            name='Human Toxicity, non-cancer',
            basis='kg',
            inventory=lambda: 0,
            CF=1.,
        )
        
        #acidification
        system.define_process_impact(
            key='ACD',
            name='Acidification',
            basis='kg',
            inventory=lambda: 0,
            CF=1.,
        )
        
        #freshwater ecotoxicity
        system.define_process_impact(
            key='ETOX',
            name='Ecotoxicity',
            basis='kg',
            inventory=lambda: 0,
            CF=1.,
        )
        
        #ozone depletion
        system.define_process_impact(
            key='OZD',
            name='Ozone depletion',
            basis='kg',
            inventory=lambda: 0,
            CF=1.,
        )
        
        #Photochemical Ozone Creation Potential
        system.define_process_impact(
            key='POCP',
            name='Photochemical Ozone Creation Potential',
            basis='kg',
            inventory=lambda: 0,
            CF=1.,
        )
        
        system.set_tolerance(mol=1e-6, rmol=1e-9, T=1e-6, rT=1e-9, subsystems=True, maxiter=200)      
        return system
        
    def create_model(self):
        scenario = self.scenario
        if not scenario.turbogenerator: self.BT = self.B
        facilities = scenario.facilities
        processing_capacity = scenario.processing_capacity
        target_plastic_percent = scenario.target_plastic_percent
        system = self.system
        if scenario.sell_leftover_plastic:
            products = system.outs
            
            #TO DO: if sell_leftover_pastic=False, set leftover_plastic.price = 0
            
        else:
            _, *products = system.outs
        
        self.products = products
        model = bst.Model(system)
        parameter = model.parameter
        metric = model.metric
        self.tea_parameters = []
        def tea_param(f):
            self.tea_parameters.append(f)
            return f
        
        self.lca_parameters = []
        def lca_param(f):
            self.lca_parameters.append(f)
            return f
        
        self.general_parameters = []
        def gen_param(f):
            self.tea_parameters.append(f)
            self.lca_parameters.append(f)
            self.general_parameters.append(f)
            return f
        
        #global warming potential
        @metric(units='kg*CO2e/kg')
        def GWP():
            if scenario.burn_leftover_plastic:
                GWP_material = system.get_total_feeds_impact(GWP_key)
                GWP_electricity_production = CFs['Electricity'] * (self.system.get_electricity_production() - self.system.get_electricity_consumption())
                GWP_emissions = system.get_process_impact(GWP_key) # kg CO2 eq. / y
                GWP_total = GWP_material + GWP_emissions - GWP_electricity_production # kg CO2 eq. / y
                return GWP_total / (self.PE_resin.F_mass * self.tea.operating_hours)
            else:
                GWP = system.get_property_allocated_impact(
                    key=GWP_key, name='mass', basis='kg',
                    products=products
                ) # kg-CO2e / kg
                if GWP < 0: breakpoint()
            return GWP
        
        #fossil fuel consumption
        @metric(units='MJ/kg')
        def FFC():
            if scenario.burn_leftover_plastic:
                FFC_material = system.get_total_feeds_impact('FFC')
                FFC_electricity_production = CFs['Electricity'] * (self.system.get_electricity_production() - self.system.get_electricity_consumption())
                #FFC_emissions = system.get_process_impact('FFC') # kg CO2 eq. / y
                FFC_total = FFC_material  - FFC_electricity_production # kg CO2 eq. / y
                return FFC_total / (self.PE_resin.F_mass * self.tea.operating_hours)
            else:
                FFC = system.get_property_allocated_impact(
                    key='FFC', name='mass', basis='kg',
                    products=[i for i in products if 'resin' in i.ID or 'product' in i.ID]
                ) # MJ / kg
                if FFC < 0: breakpoint()
            return FFC
        
        #water usage
        @metric(units='m3/kg')
        def WU():
            try:
                if scenario.burn_leftover_plastic:
                    WU_material = system.get_total_feeds_impact('WU')
                    WU_electricity_production = CFs['Electricity'] * (self.system.get_electricity_production() - self.system.get_electricity_consumption())
                    #FFC_emissions = system.get_process_impact('FFC') # kg CO2 eq. / y
                    WU_total = WU_material  - WU_electricity_production # kg CO2 eq. / y
                    return WU_total / (self.PE_resin.F_mass * self.tea.operating_hours)
                else:
                    WU = system.get_property_allocated_impact(
                        key='WU', name='mass', basis='kg',
                        products=[i for i in products if 'resin' in i.ID or 'product' in i.ID]
                    ) # MJ / kg
                    if WU < 0: breakpoint()
            except:
                WU=None
            return WU
        
        #human toxicity - CANCER
        @metric(units='CTUh/kg')
        def HTC():
            if scenario.burn_leftover_plastic:
                HTC_material = system.get_total_feeds_impact('HTC')
                HTC_electricity_production = CFs['Electricity'] * (self.system.get_electricity_production() - self.system.get_electricity_consumption())
                #FFC_emissions = system.get_process_impact('FFC') # kg CO2 eq. / y
                HTC_total = HTC_material  - HTC_electricity_production # kg CO2 eq. / y
                return HTC_total / (self.PE_resin.F_mass * self.tea.operating_hours)
            else:
                HTC = system.get_property_allocated_impact(
                    key='HTC', name='mass', basis='kg',
                    products=[i for i in products if 'resin' in i.ID or 'product' in i.ID]
                ) # MJ / kg
                if HTC < 0: breakpoint()
            return HTC
        
        #human toxicity - non CANCER
        @metric(units='CTUh/kg')
        def HTNC():
            if scenario.burn_leftover_plastic:
                HTNC_material = system.get_total_feeds_impact('HTNC')
                HTNC_electricity_production = CFs['Electricity'] * (self.system.get_electricity_production() - self.system.get_electricity_consumption())
                #FFC_emissions = system.get_process_impact('FFC') # kg CO2 eq. / y
                HTNC_total = HTNC_material  - HTNC_electricity_production # kg CO2 eq. / y
                return HTNC_total / (self.PE_resin.F_mass * self.tea.operating_hours)
            else:
                HTNC = system.get_property_allocated_impact(
                    key='HTNC', name='mass', basis='kg',
                    products=[i for i in products if 'resin' in i.ID or 'product' in i.ID]
                ) # MJ / kg
                if HTNC < 0: breakpoint()
            return HTNC
        
        #acidification
        @metric(units='MOL H+ eq/kg')
        def ACD():
            if scenario.burn_leftover_plastic:
                ACD_material = system.get_total_feeds_impact('ACD')
                ACD_electricity_production = CFs['Electricity'] * (self.system.get_electricity_production() - self.system.get_electricity_consumption())
                #FFC_emissions = system.get_process_impact('FFC') # kg CO2 eq. / y
                ACD_total = ACD_material  - ACD_electricity_production # kg CO2 eq. / y
                return ACD_total / (self.PE_resin.F_mass * self.tea.operating_hours)
            else:
                ACD = system.get_property_allocated_impact(
                    key='ACD', name='mass', basis='kg',
                    products=[i for i in products if 'resin' in i.ID or 'product' in i.ID]
                ) # MJ / kg
                if ACD < 0: breakpoint()
            return ACD
        
        #ecotoxicity
        @metric(units='CTU eq/kg')
        def ETOX():
            if scenario.burn_leftover_plastic:
                ETOX_material = system.get_total_feeds_impact('ETOX')
                ETOX_electricity_production = CFs['Electricity'] * (self.system.get_electricity_production() - self.system.get_electricity_consumption())
                #FFC_emissions = system.get_process_impact('FFC') # kg CO2 eq. / y
                ETOX_total = ETOX_material  - ETOX_electricity_production # kg CO2 eq. / y
                return ETOX_total / (self.PE_resin.F_mass * self.tea.operating_hours)
            else:
                ETOX = system.get_property_allocated_impact(
                    key='ETOX', name='mass', basis='kg',
                    products=[i for i in products if 'resin' in i.ID or 'product' in i.ID]
                ) # MJ / kg
                if ETOX < 0: breakpoint()
            return ETOX
        
        #ozone depletion
        @metric(units='kg CFC11 eq/kg')
        def OZD():
            if scenario.burn_leftover_plastic:
                OZD_material = system.get_total_feeds_impact('OZD')
                OZD_electricity_production = CFs['Electricity'] * (self.system.get_electricity_production() - self.system.get_electricity_consumption())
                #FFC_emissions = system.get_process_impact('FFC') # kg CO2 eq. / y
                OZD_total = OZD_material  - OZD_electricity_production # kg CO2 eq. / y
                return OZD_total / (self.PE_resin.F_mass * self.tea.operating_hours)
            else:
                OZD = system.get_property_allocated_impact(
                    key='OZD', name='mass', basis='kg',
                    products=[i for i in products if 'resin' in i.ID or 'product' in i.ID]
                ) # MJ / kg
                if OZD < 0: breakpoint()
            return OZD
        
        #photochemical ozone creation potential
        @metric(units='kg CFC11 eq/kg')
        def POCP():
            if scenario.burn_leftover_plastic:
                POCP_material = system.get_total_feeds_impact('POCP')
                POCP_electricity_production = CFs['Electricity'] * (self.system.get_electricity_production() - self.system.get_electricity_consumption())
                #FFC_emissions = system.get_process_impact('FFC') # kg CO2 eq. / y
                POCP_total = POCP_material  - POCP_electricity_production # kg CO2 eq. / y
                return POCP_total / (self.PE_resin.F_mass * self.tea.operating_hours)
            else:
                POCP = system.get_property_allocated_impact(
                    key='POCP', name='mass', basis='kg',
                    products=[i for i in products if 'resin' in i.ID or 'product' in i.ID]
                ) # MJ / kg
                if POCP < 0: breakpoint()
            return POCP
        
        @metric(units='USD/kg')
        def MSP():
            return self.tea.solve_price(products)
        
        V_ng = 1.473318463076884 # Natural gas volume at 60 F and 14.73 psi [m3 / kg]
        
        # https://www.eia.gov/energyexplained/natural-gas/prices.php
        # @parameter(analysis='Sobol', group='tea')
        if facilities:
            @tea_param
            @parameter(distribution=dist.natural_gas_price_distribution, element='Natural gas', units='USD/m3',
                       baseline=4.73 * 35.3146667/1e3)
            def set_natural_gas_price(price): 
                self.BT.natural_gas_price = price * V_ng
    
        
        # Processing capacity is entirely arbitrary for now
        @parameter(
            bounds=(processing_capacity * 0.5, processing_capacity * 2),
            element='Feedstock',
            units='MT/yr',
            baseline=processing_capacity,
        )
        def set_processing_capacity(processing_capacity):
            self.feedstock.F_mass = processing_capacity * 1000 / self.tea.operating_hours
        
        # Feedstock price will be equal to transportation cost.
        # TODO: Base estimate to transportation cost on availability of 
        # post-industrial plastic per area.
        
        @tea_param
        @parameter(
            baseline=0.01,
            element='Feedstock',
            units='USD/kg',
            distribution=shape.Uniform(0, 0.02)
        )
        def set_feedstock_price(price):
            self.feedstock.price = price
        
        @gen_param
        @parameter(
            baseline = 500, #km
            element = 'feedstock',
            units='km',
            distribution=shape.Uniform(20, 2000)
        )
        def set_feedstock_distance(distance):
            GWP = 'GWP'
            FFC = 'FFC'
            WU = 'WU'
            HTC = 'HTC'
            HTNC = 'HTNC'
            ETOX = 'ETOX'
            ACD = 'ACD'
            OZD = 'OZD'
            POCP = 'POCP'
            #from Articulated lorry transport, Total weight 12-14 t, mix Euro 0-5, consumption mix, to consumer, diesel driven, Euro 0 - 5 mix, cargo, 12 - 14t gross weight / 9.3t payload capacity - RNA
            self.plastic.set_CF(GWP, distance * 0.083 * 1/1000, ) #distance (km) * .08 kg co2 / (t*km) * t/1000kg EF 2.0
            self.plastic.set_CF(FFC, distance * 1.08828 * 1/1000, ) #distance (km) * .08 kg co2 / (t*km) * t/1000kg EF 2.0
            self.plastic.set_CF(WU, distance * 0.00799 * 1/1000, ) #distance (km) * .08 kg co2 / (t*km) * t/1000kg EF 2.0
            self.plastic.set_CF(HTC, distance * 1.25436e-9 * 1/1000, )
            self.plastic.set_CF(HTNC, distance * 2.884e-9 * 1/1000, )
            self.plastic.set_CF(ETOX, distance * 0.02214 * 1/1000, )
            self.plastic.set_CF(ACD, distance * 0.00065 * 1/1000, )
            self.plastic.set_CF(OZD, distance * 2.70507e-13 * 1/1000, )
            self.plastic.set_CF(POCP, distance * 0.0006 * 1/1000, )
        
        @tea_param
        @parameter(
            element='Cashflow analysis',
            baseline=0.10,
            units='%',
            distribution=shape.Uniform(0.1, 0.2)
        )
        def set_IRR(IRR):
            self.tea.IRR = IRR
        
        if scenario.multistep:
            baseline = 100 * self.dissolution_steps[0].solvent_content
            @gen_param
            @parameter(
                distribution=shape.Uniform(baseline - 10, baseline + 10),
                element='centrifuged plastic', units='%',
            )
            def set_centrifuged_plastic_solvent_content(solvent_content):
                for i in self.dissolution_steps:
                    i.solvent_content = solvent_content / 100
            
            @gen_param
            @parameter(
                baseline=0.005,
                element='Solvent',
                units='%',
                distribution=shape.Uniform(0.001 * 100, 0.01 * 100)
            )
            def set_solvent_loss(solvent_loss):
                for i in self.dissolution_steps:
                    getattr(self, i.solvent + '_loss').split[:] = solvent_loss / 100
            
            first_step, *other_steps = self.dissolution_steps
            self.target_plastics_ratio = target_plastics_ratio = {i.plastic: 0 for i in other_steps}    
            
            def create_parameters(dissolution_step, precipitation_step):
                
                def step_parameter(*args, element=None, **kwargs):
                    return lambda f: _step_parameter(f, *args, element=element, **kwargs)
                    
                def _step_parameter(f, *args, element, **kwargs):
                    element = f'{dissolution_step.plastic} step-{element}'
                    f.__name__ = f.__name__.replace('set_', f'set_{dissolution_step.plastic}_')
                    return parameter(f, *args, element=element, **kwargs)
                
                baseline = 2.17
                @tea_param
                @step_parameter(
                    baseline=baseline,
                    element='solvent',
                    units='USD/kg',
                    distribution=shape.Uniform(0.5 * baseline, 1.5 * baseline),
                )
                def set_solvent_price(price):
                    getattr(self, dissolution_step.solvent).price = price
                
                if dissolution_step is not first_step:
                    @gen_param
                    @step_parameter(
                        baseline=target_plastic_percent[self.dissolution_steps.index(dissolution_step)] / target_plastic_percent[0],
                        element='polymer',
                        distribution=shape.Uniform(0.1, 1),
                    )
                    def set_polymer_ratio(ratio):
                        self.target_plastics_ratio[dissolution_step.plastic] = ratio
                
                @gen_param
                @step_parameter(
                    element='dissolution', units='wt %',
                    distribution=shape.Uniform(1, 5),
                    baseline=dissolution_step.capacity * 100,
                )
                def set_dissolution_capacity(solvent_capacity):
                    dissolution_step.capacity = solvent_capacity / 100
                
                chemicals = bst.settings.chemicals
                solvent = chemicals[dissolution_step.solvent]
                if isinstance(solvent, list):
                    T = dissolution_step.T
                    Tmax = min([i.Tb for i in solvent]) - 5
                    Tmin = max(
                        max([i.Tm for i in solvent]) + 5, 265
                    )
                    if T > Tmax: T = Tmax - 1
                else:
                    T = dissolution_step.T
                    Tmax = solvent.Tb - 5
                    Tmin = max(solvent.Tm + 5, 265)
                    if T > Tmax: T = Tmax - 1
                
                @step_parameter(
                    element='dissolution', units='K', 
                    distribution=shape.Triangle(Tmin, T, Tmax)
                )
                def set_dissolution_temperature(temperature):
                    dissolution_step.T = temperature
                
                if scenario.precipitation_temperature_format == 'drop':
                    @gen_param
                    @step_parameter(
                        element='precipitation', units='%',
                        distribution=shape.Uniform(50, 100),
                    )
                    def set_precipitation_temperature_drop(temperature_drop):
                        T = dissolution_step.T
                        precipitation_step.T = (
                            T - temperature_drop / 100 * (T - Tmin)
                        )
                elif scenario.precipitation_temperature_format == 'constant':
                    T_precipitation = precipitation_step.T
                    @gen_param
                    @step_parameter(
                        element='precipitation', units='K',
                        baseline=T_precipitation,
                        distribution=shape.Uniform(T_precipitation - 5, T_precipitation + 5),
                    )
                    def set_precipitation_temperature(temperature):
                        precipitation_step.T = temperature
                
            
            for i in range(scenario.N_steps):
                create_parameters(
                    dissolution_step=self.dissolution_steps[i],
                    precipitation_step=self.precipitation_steps[i]
                )
                
            @gen_param
            @parameter(
                baseline=sum(target_plastic_percent) / 100.,
                element='target polymer',
                distribution=shape.Uniform(0.3, 0.9)
            )
            def set_polymer_mass_fraction(mass_fraction):
                self.target_polymer_mass_fraction = mass_fraction
        else:
            @gen_param
            @parameter(
                baseline=target_plastic_percent / 100.,
                element='target polymer',
                distribution=shape.Uniform(0.3, 0.9)
            )
            def set_polymer_mass_fraction(mass_fraction):
                s = self.feedstock
                F_mass = s.F_mass
                plastic = self.dissolution_step.plastic
                s.imass[plastic] = 0
                other_composition = s.mass / s.F_mass
                s.mass = other_composition * F_mass * (1 - mass_fraction)
                s.imass[plastic] = mass_fraction * F_mass
            
            def solvent_content(baseline, *args, **kwargs):
                bounds = (baseline - 10, baseline + 10)
                return parameter(*args, bounds=bounds, baseline=baseline,
                                 distribution=shape.Uniform(*bounds), 
                                 units='%', **kwargs)
            
            @gen_param
            @solvent_content(
                100 * self.dissolution_step.solvent_content,
                element='centrifuged plastic'
            )
            def set_centrifuged_plastic_solvent_content(solvent_content):
                self.dissolution_step.solvent_content = solvent_content / 100
            
            baseline = 2.17
            @tea_param
            @parameter(
                baseline=baseline,
                element='Solvent',
                units='USD/kg',
                distribution=shape.Uniform(0.5 * baseline, 1.5 * baseline)
            )
            def set_solvent_price(price):
                self.solvent.price = price
            
            @gen_param
            @parameter(
                baseline=0.001 * 100,
                element='solvent',
                units='%',
                distribution=shape.Uniform(0.0001 * 100, 0.002 * 100)
            )
            def set_solvent_loss(solvent_loss):
                self.solvent_loss.split[:] = solvent_loss / 100.
            
            @gen_param
            @parameter(
                element='dissolution', units='wt %',
                distribution=shape.Uniform(1, 5),
                baseline=self.dissolution_step.capacity * 100,
            )
            def set_dissolution_capacity(solvent_capacity):
                self.dissolution_step.capacity = solvent_capacity / 100
            
            chemicals = bst.settings.chemicals
            solvent = chemicals[self.dissolution_step.solvent]
            T = self.dissolution_step.T
            Tmax = solvent.Tb - 5
            Tmin = max(solvent.Tm + 5, 265)
            if T > Tmax: T = Tmax - 1
            @parameter(
                element='dissolution', units='K',
                distribution=shape.Triangle(Tmin, T, Tmax),
                baseline=T,
            )
            def set_dissolution_temperature(temperature):
                self.dissolution_step.T = temperature
            
            if scenario.precipitation_temperature_format == 'drop':
                @gen_param
                @parameter(
                    element='precipitation', units='%',
                    distribution=shape.Uniform(50, 100),
                )
                def set_precipitation_temperature_drop(temperature_drop):
                    T = self.dissolution_step.T
                    self.precipitation_step.T = (
                        T - temperature_drop / 100 * (T - Tmin)
                    )
            elif scenario.precipitation_temperature_format == 'constant':
                T_precipitation = self.precipitation_step.T
                @gen_param
                @parameter(
                    element='precipitation', units='K',
                    baseline=T_precipitation,
                    distribution=shape.Uniform(T_precipitation - 5, T_precipitation + 5),
                )
                def set_precipitation_temperature(temperature):
                    self.precipitation_step.T = temperature
        
        # @gen_param
        # @solvent_content(
        #     precipitation_step.centrifuge_solvent_content,
        #     element='centrifuged precipitate'
        # )
        # def set_centrifuged_precipitate_solvent_content(solvent_content):
        #     precipitation_step.centrifuge_solvent_content = solvent_content
        
        # @gen_param
        # @solvent_content(
        #     precipitation_step.screw_press_solvent_content,
        #     element='screw pressed precipitate'
        # )
        # def set_screw_press_solvent_content(solvent_content):
        #     precipitation_step.screw_press_solvent_content = solvent_content
        
        # chemicals = bst.settings.chemicals
        # solvent = chemicals[dissolution_step.solvent]
        # solvent.Psat.method = 'BOILING_CRITICAL'
        # # This line resets the extrapolation coefficients
        # solvent.Psat.extrapolation_coeffs.clear()
        
        # @gen_param
        # @parameter(
        #     baseline=solvent.Tb,
        #     element='Solvent', units='K',
        #     distribution=shape.Uniform(solvent.Tb - 25, solvent.Tb + 25),
        # )
        # def set_boiling_point(normal_boiling_point):
        #     solvent.Tb = normal_boiling_point
        #     # This line resets the extrapolation coefficients
        #     solvent.Psat.extrapolation_coeffs.clear()
        if scenario.multistep:
            @system.add_specification(simulate=True)
            def adjust_composition():
                s = self.feedstock
                F_mass = s.F_mass
                IDs = [i.plastic for i in self.dissolution_steps]
                s.imass[IDs] = 0
                other_composition = s.mass / s.F_mass
                s.mass = other_composition * F_mass * (1 - self.target_polymer_mass_fraction)
                composition = np.array([1, *self.target_plastics_ratio.values()], dtype=float)
                composition /= composition.sum()
                s.imass[IDs] = self.target_polymer_mass_fraction * F_mass * composition
            
        self.load_model(model)
        for i in ('emissions', 'natural_gas', 'makeup_water', 'cooling_tower_makeup_water'):
            if not hasattr(self, i): setattr(self, i, bst.Stream(i))
        if facilities:
            self.natural_gas.set_CF(
                GWP_key,
                0.33, # Natural gas from shell conventional recovery, GREET; includes non-biogenic emissions
            )
            self.natural_gas.set_CF(
                'FFC',
                51, # [MJ / kg NG] From Open-LCA Environmental Footprint 2.0
            )
            
            #TODO: ADD NG WATER USE
        # TODO: Adjust solvent             CF accordingly
        if scenario.multistep:
            for i in self.dissolution_steps:
                solvent = getattr(self, i.solvent)
                solvent.set_CF(
                    GWP_key,
                    0.8199, # GREET; Mixed xylenes production from catalytic reforming of naphtha
                )
                solvent.set_CF(
                    'FFC',
                    54, # GREET; Mixed xylenes production from catalytic reforming of naphtha
                )
        else:
            self.solvent.set_CF(
                GWP_key,
                0.8199, # GREET; Mixed xylenes production from catalytic reforming of naphtha
            )
            self.solvent.set_CF(
                'FFC',
                54, # GREET; Mixed xylenes production from catalytic reforming of naphtha
            )
        return model

class STRAPProcessPE(bst.ProcessModel):
    """
    Create a model for a solvent targeted precipitation and dissolution process.
    The dissolution and precipitation steps default to PE.
    
    """
    cache = {}
    def __new__(
            cls,
            simulate=True,
            dissolution_step=None,
            precipitation_step=None,
            solvent=None,
            burn_leftover_plastic=False,
            
        ):
        
        chemicals = create_property_package()
        bst.settings.set_thermo(chemicals)
        
        if solvent is None:
             solvent = 'Xylene'
             
        load_process_settings()
        
        
        if dissolution_step is None:
            dissolution_step = strap.dissolution_steps.PE_Xylene_dissolution()
        if precipitation_step is None:
            precipitation_step = strap.precipitation_steps.PE_Xylene_precipitation()
        
        
        dissolution_step.solvent = precipitation_step.solvent = solvent
        dissolution_step.T = 378 #K or 105 C
        precipitation_step.T = 323 #K
        key = (dissolution_step, precipitation_step)
        if key in cls.cache: return cls.cache[key]
        self = super().__new__(cls)
        self.flowsheet = bst.Flowsheet()
        bst.main_flowsheet.set_flowsheet(self.flowsheet)
        PE_fraction = 0.98
        chemicals.define_group('Plastic', ['PE', 'PET'], [PE_fraction, 1-PE_fraction], wt=True)
        chemicals.define_group('Solutes', ['Minerals', 'Solubles'], [0.8, 0.2], wt=True)
        system = create_single_layer_batch_separation_system(
                dissolution_step=dissolution_step,
                precipitation_step=precipitation_step,
                facilities=True,
                burn_leftover_plastic=burn_leftover_plastic,
                molar_tolerance=1e-9,
                relative_molar_tolerance=1e-6
        )
        products = list(system.outs)
        self.products = products
        self.name = 'STRAP-B'
        self.dissolution_step = dissolution_step
        self.precipitation_step = precipitation_step
        self.tea = create_baseline_tea(system)
        
        system.define_process_impact(
            key=GWP_key,
            name='Direct non-biogenic emissions',
            basis='kg',
            inventory=lambda: self.emissions.imass['CO2'] * system.operating_hours,
            CF=1.,
        )
        
        
        model = bst.Model(system)
        parameter = model.parameter
        metric = model.metric
        self.tea_parameters = []
        def tea_param(f):
            self.tea_parameters.append(f)
            return f
        
        self.lca_parameters = []
        def lca_param(f):
            self.lca_parameters.append(f)
            return f
        
        self.general_parameters = []
        def gen_param(f):
            self.tea_parameters.append(f)
            self.lca_parameters.append(f)
            self.general_parameters.append(f)
            return f
        
        @metric(units='kg*CO2e/kg', element='product')
        def GWP():
            GWP = system.get_property_allocated_impact(
                key=GWP_key, name='mass', basis='kg',
                products=[i for i in system.outs if 'resin' in i.ID or 'product' in i.ID]
            ) # kg-CO2e / kg
            if GWP < 0: breakpoint()
            return GWP
        
        @metric(units='USD/kg')
        def MSP():
            return self.tea.solve_price(self.PE_resin)
        
        # Feedstock price will be equal to transportation cost.
        # TODO: Base estimate to transportation cost on availability of 
        # post-industrial plastic per area.
        
        @tea_param
        @parameter(
            baseline=0.035,
            element='Feedstock',
            units='USD/kg',
            distribution=shape.Uniform(0.01, 0.05)
        )
        def set_feedstock_price(price):
            self.feedstock.price = price
        
        # @tea_param
        # @parameter(
        #     baseline=0,
        #     element='Coproduct',
        #     units='USD/kg',
        #     distribution=shape.Uniform(0, 1.2)
        # )
        # def set_coproduct_price(price):
        #     self.coproduct.price = price
        
        baseline = 2.17
        @tea_param
        @parameter(
            baseline=baseline,
            element='Solvent',
            units='USD/kg',
            distribution=shape.Uniform(0.5 * baseline, 1.5 * baseline)
        )
        def set_solvent_price(price):
            self.solvent.price = price
        
        @tea_param
        @parameter(
            baseline=0.15,
            units='%',
            distribution=shape.Uniform(0.1, 0.2)
        )
        def set_IRR(IRR):
            self.tea.IRR = IRR
        
        @parameter(
            bounds=(300, 700),
            element='MSW',
            units='kg/hr',
            baseline=500,
        )
        def set_processing_capacity(processing_capacity):
            self.plastic.F_mass = processing_capacity
        
        
        @gen_param
        @parameter(
            baseline=PE_fraction,
            element='polymer',
            distribution=shape.Uniform(0.1, 0.9)
        )
        def set_polymer_mass_fraction(mass_fraction):
            s = self.feedstock
            F_mass = s.F_mass
            plastic = self.dissolution_step.plastic
            s.imass[plastic] = 0
            other_composition = s.mass / s.F_mass
            s.mass = other_composition * F_mass * (1 - mass_fraction)
            s.imass[plastic] = mass_fraction * F_mass
        
        def solvent_content(baseline, *args, **kwargs):
            bounds = (baseline - 10, baseline + 10)
            return parameter(*args, bounds=bounds, baseline=baseline,
                             distribution=shape.Uniform(*bounds), 
                             units='%', **kwargs)
        
        @gen_param
        @solvent_content(
            100 * dissolution_step.solvent_content,
            element='centrifuged plastic'
        )
        def set_centrifuged_plastic_solvent_content(solvent_content):
            dissolution_step.solvent_content = solvent_content / 100
        
        # @gen_param
        # @solvent_content(
        #     precipitation_step.centrifuge_solvent_content,
        #     element='centrifuged precipitate'
        # )
        # def set_centrifuged_precipitate_solvent_content(solvent_content):
        #     precipitation_step.centrifuge_solvent_content = solvent_content
        
        @gen_param
        @solvent_content(
            100 * precipitation_step.screw_press_solvent_content,
            element='screw pressed precipitate'
        )
        def set_screw_press_solvent_content(solvent_content):
            precipitation_step.screw_press_solvent_content = solvent_content / 100
        
        chemicals = bst.settings.chemicals
        solvent = chemicals[dissolution_step.solvent]
        T = dissolution_step.T
        solvent.Psat.method = 'BOILING_CRITICAL'
        # This line resets the extrapolation coefficients
        solvent.Psat.extrapolation_coeffs.clear()
        
        @gen_param
        @parameter(
            baseline=solvent.Tb,
            element='Solvent', units='K',
            distribution=shape.Uniform(solvent.Tb - 25, solvent.Tb + 25),
        )
        def set_boiling_point(normal_boiling_point):
            solvent.Tb = normal_boiling_point
            # This line resets the extrapolation coefficients
            solvent.Psat.extrapolation_coeffs.clear()
        
        Tmax = solvent.Tb - 5
        Tmin = max(solvent.Tm + 15, 265)
        
        @gen_param
        @parameter(
            element='Dissolution', units='K',
            distribution=shape.Triangle(Tmin, T, Tmax)
        )
        def set_dissolution_temperature(temperature):
            dissolution_step.T = temperature
        
        @gen_param
        @parameter(
            element='Precipitation', units='%',
            distribution=shape.Uniform(0, 100),
        )
        def set_precipitation_temperature_drop(temperature_drop):
            T = dissolution_step.T
            precipitation_step.T = (
                T - temperature_drop / 100 * (T - Tmin)
            )
        
        @gen_param
        @parameter(
            element='Dissolution', units='wt %',
            distribution=shape.Uniform(1, 10),
        )
        def set_dissolution_capacity(solvent_capacity):
            dissolution_step.capacity = solvent_capacity / 100
        
        self.load_system(system)
        self.BT.ins[0] = self.solvent_loss.outs[0]
        self.system.update_configuration()
        self.load_model(model)
        for i in model.parameters:
            if i.baseline is not None: i.setter(i.baseline)
        if simulate: system.simulate()
        self.natural_gas.set_CF(
            GWP_key,
            0.33, # Natural gas from shell conventional recovery, GREET; includes non-biogenic emissions
        )
        cls.cache[key] = self
        return self
                        
    
  
    
  
# Neodymium Magnet Recovery
class MagnetHandSorting(bst.Unit):
    N_shifts = 3
    N_workers = 3
    
    def _init(self, N_workers,N_shifts):
        self.N_workers = N_workers
        self.N_shifts = N_shifts
        self._cost()
        self.total_salary = 5e4*N_shifts*N_workers*1.6
        
    def _cost(self):
        worker_salary = 5e4
        self.total_salary = worker_salary*self.N_shifts*self.N_workers*1.6
    
    def results(self):
        import pandas as pd
        results = {
            "Parameters": ['Workers', 'Shifts (8hrs/shift)','Scaling(Vacation etc.)','Total Labor Cost'],
            "Value":[self.N_workers, self.N_shifts, 1.6, self.total_salary]
            }
        return pd.DataFrame(results)



class MagnetRecovery(bst.ProcessModel):
    """
    create a model for using STRAP to recover Neodymium magnets
    
    E.g. 
    >>> from plastics.strap import MagnetRecovery
    >>> pm = MagnetRecovery(simulate = False)
    
    
    
    .... another things 
    """
    
    
    
    @bst.scenario
    class Scenario:
        #copy paste from above
       solvent: str|tuple[str, ...] = '# Solvent used to separate the target plastic'
       target_plastic: str|tuple[str, ...] = '# The polymer layer being dissolved'
       target_plastic_percent: float|tuple[float, ...] = '# Fraction in feedstock [%]'
       processing_capacity: float = 5000, '# Feedstock flow rate [MT-plastic/yr]'
       sell_leftover_plastic: bool = False, '# Whether the MSP will include all products'
       burn_leftover_plastic: bool = True, '# Produce heat and power from leftover plastic'
       facilities: bool = True, '# On-site cooling tower, heat and power generation'
       precipitation_temperature_format: str = 'constant', "# Use 'drop' for % temperature drop to solvent melting point. Use 'constant' to set in Kelvin."
       precipitation_configuration: str = 'integrated heat transfer', "# Must be either 'solvent mixing' or 'integrated heat transfer'."
       turbogenerator: bool = True, '# On-site electricity generation'
       
       @property
       def N_steps(self):
           return 1
    
    
    @classmethod
    def get_scenario(cls, scenario):
        return cls._scenarios[scenario]
    
    
    @classmethod
    def set_scenario(cls, scenario):
        cls._scenarios[f'{scenario.target_plastic}/{scenario.solvent}'] = scenario
        
    _scenario =[
        Scenario('Xylene', 'HDPE', 14, 325, False, False, False)
            ]
    _scenario = {
            f'{i.target_plastic}/{i.solvent}': i for i in _scenario
            }
        
    @property
    def name(self):
        scenario = self.scenario
        
        name = f"{scenario.solvent}_{scenario.target_plastic}"
        if scenario.sell_leftover_plastic:
            name += "_sell_leftover_plastic"
        elif scenario.burn_leftover_plastic:
            name += "_burn_leftover_plastic"
        if scenario.facilities:
            name += "_with_facilities"
            
        return name
    
    @classmethod
    def default_scenario(cls):
        return cls._scenario['HDPE/Xylene'] 
    
    @classmethod
    def as_scenario(cls,scenario):
        if isinstance(scenario,(str,tuple)):
            return cls._scenarios[scenario]
        
        else:
            raise TypeError('invalid scenario')
    
    def create_thermo(self):
        scenario = self.scenario
        
        default_plastic_solvent_pair(scenario.target_plastic, scenario.solvent)
        chemicals = create_property_package()
        for i in solvent_mixtures: chemicals.define_groups(*i)
        return chemicals
    
    def create_system(self):
        scenario = self.scenario
        chemicals = self.chemicals
        load_process_settings()
        
        #changing CEPCI
        bst.settings.CEPCI = 836.9
        
        dissolution_step = getattr(strap.dissolution_steps, f'{scenario.target_plastic }_{scenario.solvent}_dissolution')()
        precipitation_step = getattr(strap.precipitation_steps,f'{scenario.target_plastic}_{scenario.solvent}_precipitation')()
        
        facilities = scenario.facilities
        target_plastic_percent = scenario.target_plastic_percent
        
        bulk_plastic_percent = 100 - target_plastic_percent
        chemicals.define_group('Plastic', [scenario.target_plastic, 'NdFeB'], [target_plastic_percent, bulk_plastic_percent], wt=True)
        
        chemicals.define_group('Solutes', ['Minerals', 'Solubles'], [0.8, 0.2], wt=True)
        
        system = create_single_layer_batch_separation_system(
            dissolution_step=dissolution_step,
            precipitation_step=precipitation_step,
            facilities=facilities,
            relative_molar_tolerance=1e-6, 
            molar_tolerance=1e-2,
            method='fixed-point',
            burn_leftover_plastic=scenario.burn_leftover_plastic,
            core_facilities=True,
            precipitation_configuration=scenario.precipitation_configuration,
            turbogenerator=scenario.turbogenerator,
        )
        for i, step in zip(system.subsystems, [dissolution_step]):
            i.ID = f"{step.plastic}/{step.solvent}"
        
        self.dissolution_step= dissolution_step
        self.precipitation_step = precipitation_step
        
        #making object of Handsorting class
        self.HS = MagnetHandSorting(ins = system.ins[0],N_workers=3,N_shifts=3)  #change workers and shifts
        
        u = system.flowsheet.unit
           
        u.T1.ins[0] = self.HS.outs[0]
        u.U3.ins[0] = u.T1.outs[0]
        
        
        
        self.HS.outs[0].ID = 'Impellers'         #changing the sorted feed of interest to impellers
        #process.M2.outs[0].ID = 'NdFeB Magnets'
        u.U9.outs[0].ID = 'HDPE resins'

        #creating a splitter for magnets
        self.S_mag = bst.Splitter(split=1)
        self.S_mag.isplit['NdFeB'] =  0
        
        self.S_mag.ins[0] = u.T4.outs[0]
        self.S_mag.outs[0] = u.P3.ins[0]
        
        u.P3.outs[0] = u.U6.ins[0]
        
        #making it neat and clean, removing empty stream from M3 mixer
        u.M3.ins[:] = [u.M3.ins[i] for i in (0,2,3)]
        
        system.outs[:] = [s for s in system.outs if 'M2' not in s.ID]
        
        system.units.remove(u.U1)
        system.units.remove(u.U4)
        system.units.remove(u.U5)
        system.units.remove(u.M2)
        system.units.append(self.HS)
        system.units.append(self.S_mag)
        
        system.update_configuration()
        
        
        self.tea = create_baseline_tea(system)
        
        self.direct_nonbiogenic_emissions = lambda: self.emissions.imass['CO2'] * system.operating_hours
        system.define_process_impact(
            key=GWP_key,
            name='Direct non-biogenic emissions',
            basis='kg',
            inventory=self.direct_nonbiogenic_emissions,
            CF=1.,
        )
        system.define_process_impact(
            key='FFC',
            name='Direct non-biogenic emissions',
            basis='kg',
            inventory=lambda: 0,
            CF=1.,
        )
        
        #add water consumption
        system.define_process_impact(
            key='WU',
            name='Direct non-biogenic emissions',
            basis='kg',
            inventory=lambda: 0,
            CF=1.,
        )
        
        #add human toxicity, cancer
        system.define_process_impact(
            key='HTC',
            name='Human Toxicity, cancer',
            basis='kg',
            inventory=lambda: 0,
            CF=1.,
        )
        
        #add human toxicity, non-cancer
        system.define_process_impact(
            key='HTNC',
            name='Human Toxicity, non-cancer',
            basis='kg',
            inventory=lambda: 0,
            CF=1.,
        )
        
        #acidification
        system.define_process_impact(
            key='ACD',
            name='Acidification',
            basis='kg',
            inventory=lambda: 0,
            CF=1.,
        )
        
        #freshwater ecotoxicity
        system.define_process_impact(
            key='ETOX',
            name='Ecotoxicity',
            basis='kg',
            inventory=lambda: 0,
            CF=1.,
        )
        
        #ozone depletion
        system.define_process_impact(
            key='OZD',
            name='Ozone depletion',
            basis='kg',
            inventory=lambda: 0,
            CF=1.,
        )
        
        #Photochemical Ozone Creation Potential
        system.define_process_impact(
            key='POCP',
            name='Photochemical Ozone Creation Potential',
            basis='kg',
            inventory=lambda: 0,
            CF=1.,
        )
        
        system.set_tolerance(mol=1e-6, rmol=1e-9, T=1e-6, rT=1e-9, subsystems=True, maxiter=200)      
        self.HDPE_resins = u.U9.outs[0]
        self.HDPE_resins.ID = 'HDPE resins'
        
        
        self.NdFeB_Magnets = self.S_mag.outs[1]
        self.NdFeB_Magnets.ID = 'neodymium magnets'
        
        #self.__dict__['products'] = [self.HDPE_resins,self.NdFeB_Magnets]
        
        return system
    
    
    def create_model(self):
        
        
        scenario = self.scenario
        if not scenario.turbogenerator: self.BT = self.B
        facilities = scenario.facilities
        processing_capacity = scenario.processing_capacity
        target_plastic_percent = scenario.target_plastic_percent
        system = self.system
        if scenario.sell_leftover_plastic:
            products = system.outs
            
            #TO DO: if sell_leftover_pastic=False, set leftover_plastic.price = 0
            
        else:
            _, *products = system.outs
        
        self.products = products = [self.HDPE_resins, self.NdFeB_Magnets]             # changed from products
        model = bst.Model(system)
        parameter = model.parameter
        metric = model.metric
        self.tea_parameters = []
        def tea_param(f):
            self.tea_parameters.append(f)
            return f
        
        self.lca_parameters = []
        def lca_param(f):
            self.lca_parameters.append(f)
            return f
        
        self.general_parameters = []
        def gen_param(f):
            self.tea_parameters.append(f)
            self.lca_parameters.append(f)
            self.general_parameters.append(f)
            return f
        
        '''#updating labor cost
        @self.system.add_specification
        def update_labor():
            strap_labor = 580000 # 2 x operators x 4.8 x $50K, 1 engineer x $100k
            self.tea.labor_cost = strap_labor + self.HS.total_salary or 0
            self.tea.operating_days = 328.5
        '''
        
        
        #global warming potential
        @metric(units='kg*CO2e/kg')
        def GWP():
            if scenario.burn_leftover_plastic:
                GWP_material = system.get_total_feeds_impact(GWP_key)
                GWP_electricity_production = CFs['Electricity'] * (self.system.get_electricity_production() - self.system.get_electricity_consumption())
                GWP_emissions = system.get_process_impact(GWP_key) # kg CO2 eq. / y
                GWP_total = GWP_material + GWP_emissions - GWP_electricity_production # kg CO2 eq. / y
                return GWP_total / (self.HDPE_resin.F_mass * self.tea.operating_hours)
            else:
                GWP = system.get_property_allocated_impact(
                    key=GWP_key, name='mass', basis='kg',
                    products=products
                ) # kg-CO2e / kg
                if GWP < 0: breakpoint()
            return GWP
        
        #fossil fuel consumption
        @metric(units='MJ/kg')
        def FFC():
            if scenario.burn_leftover_plastic:
                FFC_material = system.get_total_feeds_impact('FFC')
                FFC_electricity_production = CFs['Electricity'] * (self.system.get_electricity_production() - self.system.get_electricity_consumption())
                #FFC_emissions = system.get_process_impact('FFC') # kg CO2 eq. / y
                FFC_total = FFC_material  - FFC_electricity_production # kg CO2 eq. / y
                return FFC_total / (self.PE_resin.F_mass * self.tea.operating_hours)
            else:
                FFC = system.get_property_allocated_impact(
                    key='FFC', name='mass', basis='kg',
                    products=[i for i in products if 'resin' in i.ID or 'product' in i.ID]
                ) # MJ / kg
                if FFC < 0: breakpoint()
            return FFC
        
        #water usage
        @metric(units='m3/kg')
        def WU():
            try:
                if scenario.burn_leftover_plastic:
                    WU_material = system.get_total_feeds_impact('WU')
                    WU_electricity_production = CFs['Electricity'] * (self.system.get_electricity_production() - self.system.get_electricity_consumption())
                    #FFC_emissions = system.get_process_impact('FFC') # kg CO2 eq. / y
                    WU_total = WU_material  - WU_electricity_production # kg CO2 eq. / y
                    return WU_total / (self.PE_resin.F_mass * self.tea.operating_hours)
                else:
                    WU = system.get_property_allocated_impact(
                        key='WU', name='mass', basis='kg',
                        products=[i for i in products if 'resin' in i.ID or 'product' in i.ID]
                    ) # MJ / kg
                    if WU < 0: breakpoint()
            except:
                WU=None
            return WU
        
        #human toxicity - CANCER
        @metric(units='CTUh/kg')
        def HTC():
            if scenario.burn_leftover_plastic:
                HTC_material = system.get_total_feeds_impact('HTC')
                HTC_electricity_production = CFs['Electricity'] * (self.system.get_electricity_production() - self.system.get_electricity_consumption())
                #FFC_emissions = system.get_process_impact('FFC') # kg CO2 eq. / y
                HTC_total = HTC_material  - HTC_electricity_production # kg CO2 eq. / y
                return HTC_total / (self.PE_resin.F_mass * self.tea.operating_hours)
            else:
                HTC = system.get_property_allocated_impact(
                    key='HTC', name='mass', basis='kg',
                    products=[i for i in products if 'resin' in i.ID or 'product' in i.ID]
                ) # MJ / kg
                if HTC < 0: breakpoint()
            return HTC
        
        #human toxicity - non CANCER
        @metric(units='CTUh/kg')
        def HTNC():
            if scenario.burn_leftover_plastic:
                HTNC_material = system.get_total_feeds_impact('HTNC')
                HTNC_electricity_production = CFs['Electricity'] * (self.system.get_electricity_production() - self.system.get_electricity_consumption())
                #FFC_emissions = system.get_process_impact('FFC') # kg CO2 eq. / y
                HTNC_total = HTNC_material  - HTNC_electricity_production # kg CO2 eq. / y
                return HTNC_total / (self.PE_resin.F_mass * self.tea.operating_hours)
            else:
                HTNC = system.get_property_allocated_impact(
                    key='HTNC', name='mass', basis='kg',
                    products=[i for i in products if 'resin' in i.ID or 'product' in i.ID]
                ) # MJ / kg
                if HTNC < 0: breakpoint()
            return HTNC
        
        #acidification
        @metric(units='MOL H+ eq/kg')
        def ACD():
            if scenario.burn_leftover_plastic:
                ACD_material = system.get_total_feeds_impact('ACD')
                ACD_electricity_production = CFs['Electricity'] * (self.system.get_electricity_production() - self.system.get_electricity_consumption())
                #FFC_emissions = system.get_process_impact('FFC') # kg CO2 eq. / y
                ACD_total = ACD_material  - ACD_electricity_production # kg CO2 eq. / y
                return ACD_total / (self.PE_resin.F_mass * self.tea.operating_hours)
            else:
                ACD = system.get_property_allocated_impact(
                    key='ACD', name='mass', basis='kg',
                    products=[i for i in products if 'resin' in i.ID or 'product' in i.ID]
                ) # MJ / kg
                if ACD < 0: breakpoint()
            return ACD
        
        #ecotoxicity
        @metric(units='CTU eq/kg')
        def ETOX():
            if scenario.burn_leftover_plastic:
                ETOX_material = system.get_total_feeds_impact('ETOX')
                ETOX_electricity_production = CFs['Electricity'] * (self.system.get_electricity_production() - self.system.get_electricity_consumption())
                #FFC_emissions = system.get_process_impact('FFC') # kg CO2 eq. / y
                ETOX_total = ETOX_material  - ETOX_electricity_production # kg CO2 eq. / y
                return ETOX_total / (self.PE_resin.F_mass * self.tea.operating_hours)
            else:
                ETOX = system.get_property_allocated_impact(
                    key='ETOX', name='mass', basis='kg',
                    products=[i for i in products if 'resin' in i.ID or 'product' in i.ID]
                ) # MJ / kg
                if ETOX < 0: breakpoint()
            return ETOX
        
        #ozone depletion
        @metric(units='kg CFC11 eq/kg')
        def OZD():
            if scenario.burn_leftover_plastic:
                OZD_material = system.get_total_feeds_impact('OZD')
                OZD_electricity_production = CFs['Electricity'] * (self.system.get_electricity_production() - self.system.get_electricity_consumption())
                #FFC_emissions = system.get_process_impact('FFC') # kg CO2 eq. / y
                OZD_total = OZD_material  - OZD_electricity_production # kg CO2 eq. / y
                return OZD_total / (self.PE_resin.F_mass * self.tea.operating_hours)
            else:
                OZD = system.get_property_allocated_impact(
                    key='OZD', name='mass', basis='kg',
                    products=[i for i in products if 'resin' in i.ID or 'product' in i.ID]
                ) # MJ / kg
                if OZD < 0: breakpoint()
            return OZD
        
        #photochemical ozone creation potential
        @metric(units='kg CFC11 eq/kg')
        def POCP():
            if scenario.burn_leftover_plastic:
                POCP_material = system.get_total_feeds_impact('POCP')
                POCP_electricity_production = CFs['Electricity'] * (self.system.get_electricity_production() - self.system.get_electricity_consumption())
                #FFC_emissions = system.get_process_impact('FFC') # kg CO2 eq. / y
                POCP_total = POCP_material  - POCP_electricity_production # kg CO2 eq. / y
                return POCP_total / (self.PE_resin.F_mass * self.tea.operating_hours)
            else:
                POCP = system.get_property_allocated_impact(
                    key='POCP', name='mass', basis='kg',
                    products=[i for i in products if 'resin' in i.ID or 'product' in i.ID]
                ) # MJ / kg
                if POCP < 0: breakpoint()
            return POCP
        
        @metric(units='USD/kg')
        def MSP():
            return self.tea.solve_price(products)
        
        V_ng = 1.473318463076884 # Natural gas volume at 60 F and 14.73 psi [m3 / kg]
        
        # https://www.eia.gov/energyexplained/natural-gas/prices.php
        # @parameter(analysis='Sobol', group='tea')
        if facilities:
            @tea_param
            @parameter(distribution=dist.natural_gas_price_distribution, element='Natural gas', units='USD/m3',
                       baseline=4.73 * 35.3146667/1e3)
            def set_natural_gas_price(price): 
                self.BT.natural_gas_price = price * V_ng
    
        
        # Processing capacity is entirely arbitrary for now
        @parameter(
            bounds=(processing_capacity * 0.5, processing_capacity * 2),
            element='Feedstock',
            units='MT/yr',
            baseline=processing_capacity,
        )
        def set_processing_capacity(processing_capacity):
            self.feedstock.F_mass = processing_capacity * 1000 / self.tea.operating_hours
        
        # Feedstock price will be equal to transportation cost.
        # TODO: Base estimate to transportation cost on availability of 
        # post-industrial plastic per area.
        
        @tea_param
        @parameter(
            baseline=0.01,
            element='Feedstock',
            units='USD/kg',
            distribution=shape.Uniform(0, 0.02)
        )
        def set_feedstock_price(price):
            self.feedstock.price = price
        
        @gen_param
        @parameter(
            baseline = 500, #km
            element = 'feedstock',
            units='km',
            distribution=shape.Uniform(20, 2000)
        )
        def set_feedstock_distance(distance):
            GWP = 'GWP'
            FFC = 'FFC'
            WU = 'WU'
            HTC = 'HTC'
            HTNC = 'HTNC'
            ETOX = 'ETOX'
            ACD = 'ACD'
            OZD = 'OZD'
            POCP = 'POCP'
            #from Articulated lorry transport, Total weight 12-14 t, mix Euro 0-5, consumption mix, to consumer, diesel driven, Euro 0 - 5 mix, cargo, 12 - 14t gross weight / 9.3t payload capacity - RNA
            self.plastic.set_CF(GWP, distance * 0.083 * 1/1000, ) #distance (km) * .08 kg co2 / (t*km) * t/1000kg EF 2.0
            self.plastic.set_CF(FFC, distance * 1.08828 * 1/1000, ) #distance (km) * .08 kg co2 / (t*km) * t/1000kg EF 2.0
            self.plastic.set_CF(WU, distance * 0.00799 * 1/1000, ) #distance (km) * .08 kg co2 / (t*km) * t/1000kg EF 2.0
            self.plastic.set_CF(HTC, distance * 1.25436e-9 * 1/1000, )
            self.plastic.set_CF(HTNC, distance * 2.884e-9 * 1/1000, )
            self.plastic.set_CF(ETOX, distance * 0.02214 * 1/1000, )
            self.plastic.set_CF(ACD, distance * 0.00065 * 1/1000, )
            self.plastic.set_CF(OZD, distance * 2.70507e-13 * 1/1000, )
            self.plastic.set_CF(POCP, distance * 0.0006 * 1/1000, )
        
        @tea_param
        @parameter(
            element='Cashflow analysis',
            baseline=0.10,
            units='%',
            distribution=shape.Uniform(0.1, 0.2)
        )
        def set_IRR(IRR):
            self.tea.IRR = IRR
        
        @gen_param
        @parameter(
            baseline=target_plastic_percent / 100.,
            element='target polymer',
            distribution=shape.Uniform(0.3, 0.9)
        )
        def set_polymer_mass_fraction(mass_fraction):
            s = self.feedstock
            F_mass = s.F_mass
            plastic = self.dissolution_step.plastic
            s.imass[plastic] = 0
            other_composition = s.mass / s.F_mass
            s.mass = other_composition * F_mass * (1 - mass_fraction)
            s.imass[plastic] = mass_fraction * F_mass
        
        def solvent_content(baseline, *args, **kwargs):
            bounds = (baseline - 10, baseline + 10)
            return parameter(*args, bounds=bounds, baseline=baseline,
                             distribution=shape.Uniform(*bounds), 
                             units='%', **kwargs)
        
        @gen_param
        @solvent_content(
            100 * self.dissolution_step.solvent_content,
            element='centrifuged plastic'
        )
        def set_centrifuged_plastic_solvent_content(solvent_content):
            self.dissolution_step.solvent_content = solvent_content / 100
        
        baseline = 2.17
        @tea_param
        @parameter(
            baseline=baseline,
            element='Solvent',
            units='USD/kg',
            distribution=shape.Uniform(0.5 * baseline, 1.5 * baseline)
        )
        def set_solvent_price(price):
            self.solvent.price = price
        
        @gen_param
        @parameter(
            baseline=0.001 * 100,
            element='solvent',
            units='%',
            distribution=shape.Uniform(0.0001 * 100, 0.002 * 100)
        )
        def set_solvent_loss(solvent_loss):
            self.solvent_loss.split[:] = solvent_loss / 100.
        
        @gen_param
        @parameter(
            element='dissolution', units='wt %',
            distribution=shape.Uniform(1, 5),
            baseline=self.dissolution_step.capacity * 100,
        )
        def set_dissolution_capacity(solvent_capacity):
            self.dissolution_step.capacity = solvent_capacity / 100
        
        chemicals = bst.settings.chemicals
        solvent = chemicals[self.dissolution_step.solvent]
        T = self.dissolution_step.T
        Tmax = solvent.Tb - 5
        Tmin = max(solvent.Tm + 5, 265)
        if T > Tmax: T = Tmax - 1
        @parameter(
            element='dissolution', units='K',
            distribution=shape.Triangle(Tmin, T, Tmax),
            baseline=T,
        )
        def set_dissolution_temperature(temperature):
            self.dissolution_step.T = temperature
        
        if scenario.precipitation_temperature_format == 'drop':
            @gen_param
            @parameter(
                element='precipitation', units='%',
                distribution=shape.Uniform(50, 100),
            )
            def set_precipitation_temperature_drop(temperature_drop):
                T = self.dissolution_step.T
                self.precipitation_step.T = (
                    T - temperature_drop / 100 * (T - Tmin)
                )
        elif scenario.precipitation_temperature_format == 'constant':
            T_precipitation = self.precipitation_step.T
            @gen_param
            @parameter(
                element='precipitation', units='K',
                baseline=T_precipitation,
                distribution=shape.Uniform(T_precipitation - 5, T_precipitation + 5),
            )
            def set_precipitation_temperature(temperature):
                self.precipitation_step.T = temperature
            
            self.load_model(model)
            for i in ('emissions', 'natural_gas', 'makeup_water', 'cooling_tower_makeup_water'):
                if not hasattr(self, i): setattr(self, i, bst.Stream(i))
            if facilities:
                self.natural_gas.set_CF(
                    GWP_key,
                    0.33, # Natural gas from shell conventional recovery, GREET; includes non-biogenic emissions
                )
                self.natural_gas.set_CF(
                    'FFC',
                    51, # [MJ / kg NG] From Open-LCA Environmental Footprint 2.0
                )
                
            #CF for xylenes
            self.solvent.set_CF(GWP_key,0.8199)
            self.solvent.set_CF('FFC',54)
            
            #self.products[:] = [self.HDPE_resins, self.NdFeB_Magnets]
            
        return model



    