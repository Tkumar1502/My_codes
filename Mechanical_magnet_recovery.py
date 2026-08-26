# -*- coding: utf-8 -*-
"""
Created on Wed Jun  3 14:19:56 2026

@author: tlnu
"""

import biosteam as bst
from plastics.strap.property_package import STRAP_chemicals_outline, create_property_package, create_property_package_MSW
from biosteam import units
from warnings import filterwarnings
import numpy as np
from biosteam import Unit
import pandas as pd
from biosteam.units.decorators import cost
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines


filterwarnings('ignore')

#from biosteam import settings
custom_chemicals = create_property_package()


Films = bst.Chemical('Films',search_ID = 'Polyethylene',phase = 's')
FittingsFilters = bst.Chemical('FittingsFilters', search_ID = 'Polypropylene', phase = 's')
#BrownSupport = bst.Chemical('BrownSupport', search_ID = 'Polypropylene', phase = 's')
BrownSupport = FittingsFilters.copy('BrownSupport')
#SiliconeTubings = bst.Chemical('SiliconeTubings', search_ID = 'Polypropylene', phase = 's')
SiliconeTubings = FittingsFilters.copy('SiliconeTubings')
paa = bst.Chemical('PeraceticAcid')
H2O2 = bst.Chemical('H2O2')

NaHSO3 = bst.Chemical('NaHSO3')
NaOH = bst.Chemical('NaOH')
Urea = bst.Chemical('Urea')
H3PO4 = bst.Chemical('H3PO4')
Na2SO4 = bst.Chemical('Na2SO4')
CH3COONa = bst.Chemical('CH3COONa')
#BiogenicResidue = bst.Chemical('BiogenicResidue', search_ID = 'Yeast', phase = 's')
BiogenicResidue = bst.Chemical('yeast').copy('BiogenicResidue')
#BiogenicResdue = 'l'

bst.settings.set_thermo(
    [custom_chemicals.NdFeB,
     custom_chemicals.HDPE,
     custom_chemicals.PeraceticAcid,
     custom_chemicals.HydrogenPeroxide,
     custom_chemicals.AceticAcid,
    'Water',
    Films,
    FittingsFilters,
    BrownSupport,
    SiliconeTubings,
    NaHSO3, NaOH, Urea, BiogenicResidue,H3PO4, Na2SO4, CH3COONa
    ])

HDPE =custom_chemicals.HDPE
NdFeB =custom_chemicals.NdFeB

#bst.settings.CEPCI = set a value
#bst.settings.heating_agents = set a value


chemicals = bst.settings.chemicals
chemicals.define_group(
    name = 'Impeller',
    IDs = ['HDPE', 'NdFeB'],
    composition = [0.14, 0.86],
    wt = True
    )

Impeller_composition = {
    HDPE.ID:0.14,
    NdFeB.ID: 0.86
    }




#%%setting the plant conditions
capacity = 1000*1976     #kg per year
operating_days = 309
N_shifts = 2 
operating_hours = operating_days*N_shifts*8

feed_per_hour = capacity/operating_hours

weight_per_bag = 4
water_per_bag = 50

fresh_water_utility = bst.Stream(ID = 'Fresh_water')
fresh_paa_utility = bst.Stream(ID = 'Peracetic Acid (15%)')

nahso3_utility = bst.Stream(ID = 'NaHSO3 (40%)')
naoh_utility = bst.Stream(ID = 'NaOH (50%)')


feed_composition = {'NdFeB': 0.1548,
        'HDPE': 0.0252,
        'Films': 0.10,
        'FittingsFilters': 0.35,
        'BrownSupport': 0.07,
        'SiliconeTubings': 0.28,
        'BiogenicResidue': 0.02
        }

feed = bst.Stream(ID = 'Bioreactor Bags', units = 'kg/hr',
                  **{chem: feed_per_hour*frac for chem, frac in feed_composition.items()})
bst.settings.CEPCI = 836.9      #updated to 2025

#%%defining units


# =============================================================================
# 1. DISINFECTION UNIT (QUENCH & PH CONTROL INCLUDED)
# =============================================================================

class Disinfection_Unit(bst.Unit):
    _N_ins = 5  # [bags_in, water_in, paa_in, nahso3_in, naoh_in]
    _N_outs = 3 # [clean_bags_out, evap_loss, wastewater]
    
    _units = {'Volume':'m3'}
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, soak_time=0.5, 
                 cycles_before_dump=4, N_workers=1, N_shifts=2):
        super().__init__(ID, ins, outs, thermo)
        self.soak_time = soak_time
        self.cycles_before_dump = cycles_before_dump
        self.N_workers = N_workers
        self.N_shifts = N_shifts
    
    def _run(self):
        bags_in, fresh_water, paa_in, nahso3_in, naoh_in = self.ins
        clean_bags_out, evap_loss, wastewater = self.outs
        
        # 1. Mass & Water Dynamics
        plastic_flow = bags_in.F_mass 
        vessel_water_pool = (plastic_flow / 4.0) * 70.0 * self.soak_time 
        
        hourly_evap_water = 5.0 
        hourly_dragout_water = plastic_flow * 0.1 
        hourly_purge_water = vessel_water_pool / (self.cycles_before_dump * self.soak_time)
        total_water_makeup = hourly_purge_water + hourly_evap_water + hourly_dragout_water
        
        # 2. Commercial PAA Dosing (15/22/16/47 formulation)
        comm_paa = total_water_makeup * 0.01
        paa_mass = comm_paa * 0.15
        h2o2_mass = comm_paa * 0.22
        aa_in_mass = comm_paa * 0.16
        paa_water = comm_paa * 0.47
        
        # 3. Quenching Stoichiometry (NaHSO3)
        active_nahso3 = 1.368 * paa_mass + 3.059 * h2o2_mass
        nahso3_water = (active_nahso3 / 0.40) - active_nahso3  # 40% commercial solution
        
        nahso4_gen = 1.579 * paa_mass + 3.530 * h2o2_mass
        total_aa = aa_in_mass + (0.790 * paa_mass)
        water_from_h2o2_quench = 0.530 * h2o2_mass
        
        # 4. Neutralization Stoichiometry (NaOH) -> Na2SO4 and CH3COONa
        naoh_for_nahso4 = 0.333 * nahso4_gen
        naoh_for_aa = 0.666 * total_aa
        active_naoh = naoh_for_nahso4 + naoh_for_aa
        naoh_water = (active_naoh / 0.50) - active_naoh        # 50% commercial solution

        na2so4_gen = 1.183 * nahso4_gen
        ch3coona_gen = 1.366 * total_aa
        water_from_neutralization = (0.150 * nahso4_gen) + (0.300 * total_aa)

        # --- Set Input Streams ---
        # Deduct water content entering from PAA, NaHSO3, and NaOH utilities
        fresh_water.empty()
        fresh_water.imass['Water'] = max(0.0, total_water_makeup - paa_water - nahso3_water - naoh_water)
        
        paa_in.empty()
        paa_in.imass['PeraceticAcid', 'HydrogenPeroxide', 'AceticAcid', 'Water'] = (
            paa_mass, h2o2_mass, aa_in_mass, paa_water
        )

        nahso3_in.empty()
        nahso3_in.imass['NaHSO3', 'Water'] = active_nahso3, nahso3_water

        naoh_in.empty()
        naoh_in.imass['NaOH', 'Water'] = active_naoh, naoh_water

        # --- Set Output Streams ---
        clean_bags_out.copy_like(bags_in)
        clean_bags_out.imass['BiogenicResidue'] = 0.0  # Completely clean off bags
        clean_bags_out.imass['Water'] = hourly_dragout_water
        
        evap_loss.empty()
        evap_loss.imass['Water'] = hourly_evap_water

        wastewater.empty()
        wastewater.imass['BiogenicResidue'] = bags_in.imass['BiogenicResidue']
        wastewater.imass['Na2SO4'] = na2so4_gen
        wastewater.imass['CH3COONa'] = ch3coona_gen
        
        # Closed Water Balance
        wastewater.imass['Water'] = (
            fresh_water.imass['Water'] +
            paa_water + 
            nahso3_water + 
            naoh_water + 
            water_from_h2o2_quench + 
            water_from_neutralization -
            hourly_evap_water -
            hourly_dragout_water
        )
    
    def _design(self):
        batch_bags = (self.ins[0].F_mass / 4.0) * self.soak_time
        vessel_volume = ((batch_bags * 70.0) / 1000.0) * 1.2 
        
        self.design_results['Vessels Required'] = 1
        self.design_results['Volume'] = vessel_volume
        self.power_utility.rate = 0.25 * vessel_volume

    @property
    def total_salary(self):
        return 70000.0 * self.N_shifts * self.N_workers * 1.6
    
    def _cost(self):
        V = self.design_results['Volume']
        base_cost = 31000.0 * (V / 4.2)**0.6
        self.purchase_costs['Disinfection System'] = base_cost + 15000.0 # Includes dosing skids
        self.installed_costs['Disinfection System'] = self.purchase_costs['Disinfection System'] * 2.0


# =============================================================================
# 2. WASTEWATER TREATMENT UNITS (MBBR vs SURCHARGE)
# =============================================================================

class Municipal_Sewer(bst.Unit):
    _N_ins = 1
    _N_outs = 1
    
    def _run(self):
        # Pass the water through the unit
        wastewater_in = self.ins[0]
        treated_water = self.outs[0]
        treated_water.copy_like(wastewater_in)
        
    def _cost(self):
        wastewater_in = self.ins[0]
        
        # 1. Calculate BOD in kg/hr
        aa_bod = wastewater_in.imass['AceticAcid'] * 0.78
        bio_bod = wastewater_in.imass['BiogenicResidue'] * 0.65
        bod_kg_hr = aa_bod + bio_bod
        
        # 2. Convert kg/hr to lbs/hr (1 kg = 2.20462 lbs)
        bod_lbs_hr = bod_kg_hr * 2.20462
        
        # 3. Calculate hourly penalty at $0.65/lb
        hourly_penalty = bod_lbs_hr * 0.65
        
        # 4. Explicitly tell BioSTEAM's TEA to track this cost
        self.add_OPEX = {'BOD_Surcharge': hourly_penalty}


class MBBR_Tank(bst.Unit):
    """If MBBR=True: Biological tank requiring Urea and Phosphoric Acid."""
    _N_ins = 3 # [Wastewater, Urea_in, H3PO4_in]
    _N_outs = 2 # [Clean_Water, Sludge]
    _units = {'Volume':'m3'}
    
    def _run(self):
        wastewater, urea_in, h3po4_in = self.ins
        clean_water, sludge = self.outs
        
        # 1. BOD Load
        aa_bod = wastewater.imass['AceticAcid'] * 0.78
        bio_bod = wastewater.imass['BiogenicResidue'] * 0.65
        self.bod_kg_hr = aa_bod + bio_bod
        
        # 2. Nutrient Dosing (BOD:N:P = 100:5:1)
        # 32.5% Urea solution
        pure_urea = (self.bod_kg_hr * 0.05) / 0.46
        urea_water = (pure_urea / 0.325) - pure_urea
        
        # 75% Phosphoric Acid solution
        pure_h3po4 = (self.bod_kg_hr * 0.01) / 0.316
        h3po4_water = (pure_h3po4 / 0.75) - pure_h3po4
        
        urea_in.empty()
        urea_in.imass['Urea', 'Water'] = pure_urea, urea_water
        
        h3po4_in.empty()
        h3po4_in.imass['H3PO4', 'Water'] = pure_h3po4, h3po4_water
        
        # 3. Treatment Outputs
        clean_water.empty()
        sludge.empty()
        
        clean_water.imass['Water'] = wastewater.imass['Water'] + urea_water + h3po4_water
        sludge.imass['BiogenicResidue'] = self.bod_kg_hr * 0.25

    def _design(self):
        # VOLR = 5.0 kg BOD / m3 / day = 0.208 kg BOD / m3 / hr
        v_working = self.bod_kg_hr / 0.208
        self.design_results['Volume'] = v_working * 1.20 # 20% freeboard
        self.power_utility.rate = 2 # kW blower

    def _cost(self):
        V = self.design_results['Volume']
        base_cost = 125000 * (V / 100.0)**0.6 
        self.purchase_costs['MBBR System'] = base_cost
        self.installed_costs['MBBR System'] = base_cost * 1.5


# =============================================================================
# 3. FACTORY ROUTING FUNCTION
# =============================================================================

def build_wastewater_treatment(wastewater_stream, MBBR=False):
    if MBBR:
        urea_feed = bst.Stream('Urea_Feed')
        h3po4_feed = bst.Stream('H3PO4_Feed')
        return MBBR_Tank('MBBR_System', 
                         ins=(wastewater_stream, urea_feed, h3po4_feed), 
                         outs=('Treated_Water', 'Bio_Sludge'))
    else:
        return Municipal_Sewer('Sewer_Discharge', 
                               ins=wastewater_stream, 
                               outs='Wastewater (to sewer)')


class HandSorting_Unit (bst.Unit):
    _N_ins=1
    _N_outs=2
    
    
    def __init__(self, ID= '', ins =None, outs=(),thermo=None,N_workers=2,N_shifts=2):
        super().__init__(ID,ins,outs,thermo)
        self.N_workers = N_workers
        self.N_shifts = N_shifts
        self.F_BM = 2
    def _run(self):
        feed = self.ins[0]
        other_components = self.outs[1]
        impellers = self.outs[0]
        
        #mass balance logic
        other_components.empty()
        impellers.empty()
        
        impellers.imass['NdFeB'] = feed.imass['NdFeB']
        impellers.imass['HDPE']= feed.imass['HDPE']
        
        other_components.imass['Films'] = feed.imass['Films']
        other_components.imass['FittingsFilters'] = feed.imass['FittingsFilters']
        other_components.imass['BrownSupport'] = feed.imass['BrownSupport']
        other_components.imass['SiliconeTubings'] = feed.imass['SiliconeTubings']
        other_components.imass['Water'] = feed.imass['Water']
     
        
    def _design(self):
        self.power_utility.rate = 0.25
        
    
    @property
    def total_salary(self):
        worker_salary = 5e4
        return worker_salary*self.N_shifts*self.N_workers*1.6        
        

    def _cost(self):
        self.power_utility.rate = 0.25                    
        self.baseline_purchase_costs['Conveyor'] = 15000
        self.purchase_costs['Conveyor'] = 15000
        self.installed_costs['Conveyor'] = 15000 * self.F_BM
    
    def results(self):
        import pandas as pd
        elec_power = self.power_utility.rate  # kW
        elec_cost = self.utility_cost  # USD/hr

        # Cost summaries
        conveyor_cost = self.purchase_costs.get("Conveyor", 0.0)
        total_purchase = sum(self.purchase_costs.values())
        total_installed = sum(self.installed_costs.values())

        # Construct MultiIndex key-value pairs
        data = [
            ("Electricity", "Power", "kW", elec_power),
            ("Electricity", "Cost", "USD/hr", elec_cost),
            ("Purchase cost", "Conveyor", "USD", conveyor_cost),
            ("Total purchase cost", "", "USD", total_purchase),
            ("Installed equipment cost", "", "USD", total_installed),
            ("Utility cost", "", "USD/hr", elec_cost),          #calculated multiplying utility rate(power)xelectricity rate
        ]

        # Convert to Pandas DataFrame with a 2-level MultiIndex on rows
        index = pd.MultiIndex.from_tuples(
            [(item[0], item[1]) for item in data], names=[self.line, ""]
        )

        df = pd.DataFrame(
            {"Units": [item[2] for item in data], self.ID: [item[3] for item in data]},
            index=index,
        )
        
        return df
       

class Lathe_Machine (bst.Unit):
    _N_ins = 1
    _N_outs = 2
    
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, power_kw=3, 
                 purchase_cost=15000, BM =2.0,N_workers=1,N_shifts=2) :
        super().__init__(ID, ins, outs, thermo)
        self.power_kw = power_kw          # kW rating of the industrial lathe
        self.custom_purchase_cost = purchase_cost # Capital cost ($) for TEA
        self.F_BM['Industrial Lathe'] = BM
        self.N_shifts = N_shifts
        self.N_workers= N_workers
        
    def _run(self):
        feed = self.ins[0]
        machined, shavings = self.outs
        
        # Mass balance logic
        machined.copy_like(feed)
        shavings.empty()
        
        shavings_mass = feed.imass['HDPE'] * 0.05
        machined.imass['HDPE'] -= shavings_mass
        shavings.imass['HDPE'] = shavings_mass
        
    def _design(self):
        self.power_utility.rate = self.power_kw 
        self.design_results['Power rating'] = self.power_kw
    
    @property
    def total_salary(self):
        worker_salary = 5e4
        return worker_salary*self.N_shifts*self.N_workers*1.6

    def _cost(self):
        self.baseline_purchase_costs['Industrial Lathe'] = self.custom_purchase_cost
        
        

class Arbor_Press(bst.Unit):
    _N_ins = 1
    _N_outs = 2
    def __init__(self, ID='', ins=None, outs=(), thermo=None, power_kw=0, purchase_cost=7000,
                 BM = 1.5,N_workers=1,N_shifts=2):
        super().__init__(ID, ins, outs, thermo)
        self.power_kw = power_kw
        self.custom_purchase_cost = purchase_cost
        self.F_BM['Arbor Press']= BM
        self.N_workers=N_workers
        self.N_shifts = N_shifts
        
    def _run(self):
        feed=self.ins[0]
        magnets = self.outs[0]
        waste_hdpe = self.outs[1]
        
        magnets.empty()
        waste_hdpe.empty()
        
        magnets.imass['NdFeB'] = feed.imass['NdFeB']
        waste_hdpe.imass['HDPE'] = feed.imass['HDPE']
        
    def _design(self):
        self.power_utility.rate = self.power_kw
        
    @property
    def total_salary(self):
        worker_salary = 5e4
        return worker_salary*self.N_shifts*self.N_workers*1.6
    
    def _cost(self):
        self.baseline_purchase_costs['Arbor Press'] = self.custom_purchase_cost
        self.purchase_costs['Arbor Press'] = self.custom_purchase_cost
        

#%%TEA parameters
class MechanicalRecyclingTEA(bst.TEA):
    """Custom TEA engine built for localized mechanical recycling."""
    def __init__(self, system, IRR, duration, depreciation,operating_days, income_tax, 
                 lang_factor, construction_schedule, labor_cost, fringe_benefits, 
                 property_tax, WC_over_FCI,property_insurance, supplies,maintenance,
                 administration):
        
        super().__init__(
            system, IRR, duration, depreciation, income_tax,
                         operating_days, lang_factor, construction_schedule,
                         startup_months=0, startup_FOCfrac=0, startup_VOCfrac=0,
                         startup_salesfrac=0, finance_interest=0, finance_years=0,
                         finance_fraction=0, WC_over_FCI=WC_over_FCI
        )
        self.labor_cost = labor_cost
        self.fringe_benefits = fringe_benefits
        self.property_tax = property_tax
        self.property_insurance = property_insurance
        self.supplies= supplies
        self.maintenance = maintenance
        self.administration = administration
        self.custom_osbl = None
        
        self.osbl_itemized = {
        "1. Site Acquisition (Pre-Built Facility)": 2250000.0,
        "2. Electrical Infrastructure Upgrade": 60000.0,
        "3. Fire Suppression & Life Safety Retrofit": 30000.0,
        "4. Medical, Health & OSHA Safety": 12000.0,
        "5. Compressed Air & Utility Piping Network": 35000.0,
        "6. Civil & Concrete Floor Reinforcement": 28000.0,
        "7. QC Analytical Station / Lab Space": 35000.0,
        "8. Maintenance Shop & Tool Crib": 25000.0,
    }
    
    # --- Direct Costs ---
    @property
    def ISBL(self):
        return self.installed_equipment_cost

    @property
    def OSBL(self):
        return sum(self.osbl_itemized.values())
    
    def OSBL_table(self,formatted=True):
        """Returns the itemized OSBL Cost Breakdown table (Cost Category & Cost)."""
        categories = list(self.osbl_itemized.keys()) + [
            "TOTAL DIRECT ASSET & OSBL COST"
        ]
        costs = list(self.osbl_itemized.values()) + [self.OSBL]
    
        if formatted:
          formatted_costs = [f"${c:,.0f}" for c in costs]
          df = pd.DataFrame(
              {"Cost Category": categories, "Cost ($)": formatted_costs}
          )
        else:
          df = pd.DataFrame({"Cost Category": categories, "Cost ($)": costs})
    
        return df.set_index("Cost Category")

    # --- BioSTEAM Engine Overrides ---
    def _DPI(self, installed_equipment_cost):
        """Direct Permanent Investment"""
        return self.ISBL + self.OSBL

    def _TDC(self, DPI):
        """Total Direct Cost"""
        return DPI

    def _FCI(self, TDC):
        """Fixed Capital Investment (Direct + Indirects)"""
        indirect_costs = 0.90 * TDC
        return TDC + indirect_costs

    def _FOC(self, FCI): 
        """Fixed Operating Costs"""
        return (FCI * (self.property_tax + self.property_insurance
                     + self.maintenance + self.administration)
                + self.labor_cost * (1 + self.fringe_benefits + self.supplies))

    # --- Convenience Properties ---
    @property
    def engineering_cost(self):
        return 0.50 * self.TDC

    @property
    def contingency(self):
        return 0.40 * self.TDC

    @property
    def TIC(self):
        return self.engineering_cost + self.contingency

    @property
    def WC(self):
        return self.WC_over_FCI * self.FCI
        

    @property
    def TCI(self):
        return self.FCI + self.WC

    # --- CAPEX Table Method ---
    def CAPEX_table(self):
        """Displays your simplified CAPEX table whenever process.tea.CAPEX_table() is called."""
        isbl = self.ISBL/1e6
        osbl = self.OSBL/1e6
        tdc = self.TDC/1e6
        eng = self.engineering_cost/1e6
        cont = self.contingency/1e6
        tic = self.TIC/1e6
        fci = self.FCI/1e6
        wc = self.WC/1e6
        tci = self.TCI/1e6
        
        index = [
            ('Direct costs', 'ISBL installed equipment cost'),
            ('Direct costs', 'OSBL installed equipment cost'),
            ('Total direct cost (TDC)', ''),
            ('Indirect costs', 'Engineering & field overhead'),
            ('Indirect costs', 'Contingency'),
            ('Total indirect cost (TIC)', ''),
            ('Fixed capital investment (FCI)', ''),
            ('Working capital (WC)', ''),
            ('Total capital investment (TCI)', '')
        ]

        notes = [
            '', 
            '', 
            '', 
            '50.0% of TDC', 
            '40.0% of TDC', 
            '', 
            'TDC + TIC', 
            '25% of FCI', 
            'FCI + WC'
        ]

        costs = [isbl, osbl, tdc, eng, cont, tic, fci, wc, tci]

        df = pd.DataFrame(
            {'Notes': notes, 'Cost [MM$]': [round(c, 3) for c in costs]},
            index=pd.MultiIndex.from_tuples(index)
        )
        return df


#%% assigning units
U1 = Disinfection_Unit('U1', ins=(feed, fresh_water_utility, fresh_paa_utility,
                                  nahso3_utility, naoh_utility), 
                       outs=('clean_rinsed_bags', 'accidental_spills_loss', 'disinfectant_solution'))
U2 = HandSorting_Unit('U2', ins = U1.outs[0], N_shifts =2, N_workers = 4)
U3 = Lathe_Machine('U3', ins = U2.outs[0],power_kw=3)

S1 = bst.Splitter('S1', ins=U3.outs[0], outs = ('to Arbor press 1', 'to Arbor press B'), split=0.5)
U4 = Arbor_Press('U4', ins=S1.outs[0])
U5 = Arbor_Press('U5',ins=S1.outs[1])
M1 = bst.Mixer('M1', ins = (U4.outs[0], U5.outs[0]), outs='Total NdFeB Magnets')
M2 = bst.Mixer('M2', ins = (U4.outs[1], U5.outs[1]), outs = 'Total Waste HDPE casings')
S2 = bst.Splitter('S2',ins=U2.outs[1],split=1)
S2.isplit['Water']=0

WM_Mixer = bst.Mixer('WM_Mixer', ins=(U1.outs[2], S2.outs[1]), 
                     outs =  'Combined_Wastewater')

#%%CHECK MBBR
MBBR = True
print("MBBR: ", MBBR)
WWT_System = build_wastewater_treatment(WM_Mixer.outs[0], MBBR)


U1.outs[0].ID = 'sterilized bags'
U2.outs[0].ID = 'impellers'
U2.outs[1].ID = 'Other components'
U3.outs[0].ID = 'trimmed Impellers'
U4.outs[0].ID = U5.outs[0].ID = 'NdFeB Magnets'
U4.outs[1].ID = U5.outs[1].ID = 'waste hdpe casings'
U3.outs[1].ID = 'waste hdpe shavings'
S2.outs[0].ID = 'other components'
S2.outs[1].ID = 'wastewater'

process = bst.System('mechanical_magnet_recovery', path = (U1, U2, U3,S1,U4,U5,M1,M2,S2,
                                                           WM_Mixer, WWT_System))

#%% wastewater BOD fee
def wastewater_bod_fee(stream):
    flow_rate_kg_hr = stream.F_mass
    if flow_rate_kg_hr == 0:
        return 0.0

    mass_fractions = {
        'Water': stream.imass['Water'] / flow_rate_kg_hr,
        'BiogenicResidue': stream.imass['BiogenicResidue'] / flow_rate_kg_hr,
        'Na2SO4': stream.imass['Na2SO4'] / flow_rate_kg_hr,
        'CH3COONa': stream.imass['CH3COONa'] / flow_rate_kg_hr
    }
    
    volumetric_rate = 10.13  # LASAN $/HCF
    bod_rate = 1.02         # LASAN $/lb BOD
    
    fBOD_CH3COONa = 0.85
    fBOD_BiogenicResidue = 0.65
    
    bod_factors = {
        'CH3COONa': (64.0 / 82.0) * fBOD_CH3COONa,
        'BiogenicResidue': (160.0 / 113.0) * fBOD_BiogenicResidue,
    }
    
    total_bod_kg_hr = sum(
        flow_rate_kg_hr * frac * bod_factors.get(comp, 0.0)
        for comp, frac in mass_fractions.items()
    )
    
    bod_lbs_hr = total_bod_kg_hr * 2.20462
    bod_surcharge = bod_lbs_hr * bod_rate
    
    return round (bod_surcharge/flow_rate_kg_hr, 4)

#%% inputs and outputs price
process.simulate()
process.diagram()

#products
NdFeB_magnets = M1.outs[0]
HDPE_casing = M2.outs[0] 
HDPE_shavings = U3.outs[1]
other_components = S2.outs[0]

#setting up prices
NdFeB_magnets.price = 100
HDPE_casing.price= -0.072
HDPE_shavings.price = -0.072
other_components.price = -0.072
          
#feed stock handling and transportation
feed.price = 0.25

#raw materials cost
fresh_water_utility.price = 0.0015
fresh_paa_utility.price = 8.8
 
#price the new utilities
nahso3_utility.price = 0.65
naoh_utility.price = 0.75

if isinstance(WWT_System, MBBR_Tank):
    WWT_System.ins[1].price = 0.60      #Urea feed
    WWT_System.ins[2].price = 0.75      #Phosphoric acid feed
    WWT_System.outs[0].price = -0.0036  #wastewater without bod
    WWT_System.outs[1].price = -0.072   #landfilling
else:
    WWT_System.outs[0].price = -(0.0036 + wastewater_bod_fee(WWT_System.outs[0])) #without mbbr

#adjusted total days for 2 shifts and 250 days, to make 4000 operating_hours
adjusted_days = operating_hours/24

# handsorting=2 workers, 1 worker lathe machine, 2 workers on 2 arbor press
total_labor_cost = U1.total_salary +  U2.total_salary + U4.total_salary + U5.total_salary + U3.total_salary

tea = MechanicalRecyclingTEA(
    system=process,
    IRR=0.15,
    duration=(2026, 2046),
    depreciation='MACRS7',
    income_tax=0.21,
    operating_days=adjusted_days,
    lang_factor=None,
    construction_schedule=(1.0,),
    WC_over_FCI=0.25,
    labor_cost=total_labor_cost,
    fringe_benefits=0.4,
    property_tax=0.001,
    property_insurance=0.005,
    supplies=0.05,
    maintenance=0.03,
    administration=0.005
)

#%%defining functions for TEA
def my_profit(processing_capacity, ndfeb_price, ndfeb_conc, feed_price, peraceticacid_price, freshwater_price,
              othercomponent_fee, wastewater_fee, hdpe_fee, worker_salary, nahso3_price, naoh_price,
              urea_price=0.60, h3po4_price=0.75):
    
    # 1. Clear intermediate stream prices to prevent BioSTEAM confusion
    WM_Mixer.outs[0].price = 0.0
    
    # 2. Apply Standard Prices
    capacity = processing_capacity * 1000
    NdFeB_magnets.price = ndfeb_price
    feed.price = feed_price
    fresh_paa_utility.price = peraceticacid_price
    fresh_water_utility.price = freshwater_price
    naoh_utility.price = naoh_price
    nahso3_utility.price = nahso3_price
    
    other_components.price = -othercomponent_fee
    HDPE_casing.price = -hdpe_fee
    HDPE_shavings.price = -hdpe_fee
    
    # 3. DYNAMIC WWT ROUTING FIX (Including Urea & H3PO4)
    if isinstance(WWT_System, MBBR_Tank):
        WWT_System.ins[1].price = urea_price
        WWT_System.ins[2].price = h3po4_price
        WWT_System.outs[0].price = -0.0036  # Base clean water discharge
        WWT_System.outs[1].price = -othercomponent_fee # Sludge disposal
    else:
        # If no MBBR, apply the heavy wastewater surcharge
        WWT_System.outs[0].price = -wastewater_fee
        
    # 4. Update Capacity and Labor
    N_shifts = 2
    tea.labor_cost = worker_salary * N_shifts * 6 * 1.6
    
    work_hours = capacity / 400
    feed_per_hour = capacity / work_hours
    tea.operating_hours = work_hours
    process.operating_hours = work_hours
    
    # 5. Mass Balance Update
    feed.imass['NdFeB'] = feed_per_hour * (0.18 * ndfeb_conc)
    feed.imass['HDPE'] = feed_per_hour * (0.18 * (1.0 - ndfeb_conc))
    feed.imass['Films'] = feed_per_hour * 0.10
    feed.imass['FittingsFilters'] = feed_per_hour * 0.36
    feed.imass['BrownSupport'] = feed_per_hour * 0.07
    feed.imass['SiliconeTubings'] = feed_per_hour * 0.29
    
    # Run the simulation once
    process.simulate()
    
    return tea.net_earnings


import scipy.optimize as opt
def scale_feedstock_flow():
    """Helper function to scale the feed stream for a target annual MT processing scale.

    Maintains fixed flow basis (400 kg/hr throughput). Defaults to baseline
    (1,976 MT/yr).
    """
    scale_mt_yr = 1976.0  # Baseline processing scale in MT/yr
    capacity_kg = scale_mt_yr * 1000.0
    
    feed_rate = 400
    # Calculate operating hours keeping feed_per_hour standard
    work_hours = (
        capacity_kg / feed_rate
    )  # Total hourly feed rate at baseline (1976 MT / 4000 hrs = 494 kg/hr)
    feed_per_hour = capacity_kg / work_hours if work_hours > 0 else 0

    # Update operating hours in TEA and System
    tea.operating_hours = work_hours
    process.operating_hours = work_hours

    # Scale component mass flow rates (kg/hr)
    for chem, frac in feed_composition.items():
        feed.imass[chem] = feed_per_hour * frac


# ==========================================
# 2. Capacity Curve Sweep Function
# ==========================================

def get_capacity_curve(min_scale_mt=50, max_scale_mt=2000):
    """
    Calculates capacity vs. net earnings for Mechanical Recovery.
    Preserves ALL low-scale simulation points + break-even + 250 MT steps.
    Includes a state reset guard to guarantee consistent baseline outputs.
    """
    # -------------------------------------------------------------------------
    # 0. STATE RESET GUARD: Restore original stream prices and baseline feed
    # -------------------------------------------------------------------------
    NdFeB_magnets.price = 100.0
    feed.price = 0.25
    fresh_paa_utility.price = 8.8
    fresh_water_utility.price = 0.0015
    nahso3_utility.price = 0.65
    naoh_utility.price = 0.75
    
    # Negative priced products / disposal fees
    HDPE_casing.price = -0.072
    HDPE_shavings.price = -0.072
    other_components.price = -0.072
    WM_Mixer.outs[0].price = -0.003
    
    # Restore baseline worker salary
    baseline_worker_salary = 50000.0
    N_shifts = 2
    total_workers = 6
    tea.labor_cost = baseline_worker_salary * N_shifts * total_workers * 1.6
    
    # Restore baseline feed composition 
    baseline_feed_composition = {
        'NdFeB': 0.1548,
        'HDPE': 0.0252,
        'Films': 0.10,
        'FittingsFilters': 0.35,
        'BrownSupport': 0.07,
        'SiliconeTubings': 0.28,
        'BiogenicResidue': 0.02
    }
    
    # 1. Reset baseline (1,976 MT/yr)
    scale_mt_yr = 1976.0
    capacity_kg = scale_mt_yr * 1000.0
    baseline_feed_rate = 494.0 # kg/hr equivalent for baseline operations
    work_hours = capacity_kg / baseline_feed_rate
    
    tea.operating_hours = work_hours
    process.operating_hours = work_hours
    
    for chem, frac in baseline_feed_composition.items():
        feed.imass[chem] = baseline_feed_rate * frac
        
    process.simulate()

    baseline_sales, baseline_voc, baseline_foc = tea.sales, tea.VOC, tea.FOC
    tax_rate = getattr(tea, "income_tax", 0.21)

    # Average annual depreciation
    cashflow_df = tea.get_cashflow_table()
    dep_col = [c for c in cashflow_df.columns if 'Depreciation' in c][0]
    dep_in_usd = cashflow_df[dep_col] * (1e6 if 'MM$' in dep_col else 1.0)
    baseline_dep = dep_in_usd[dep_in_usd > 0].mean()

    # 2. Profit simulator
    def simulate_profit(scale_mt):
        if scale_mt < 1976:
            ratio = scale_mt / 1976.0
            taxable_income = (baseline_sales * ratio) - (baseline_voc * ratio + baseline_foc + baseline_dep)
            return taxable_income * (1.0 - tax_rate)
        else:
            sim_capacity_kg = scale_mt * 1000.0
            sim_work_hours = sim_capacity_kg / baseline_feed_rate
            tea.operating_hours = process.operating_hours = sim_work_hours
            for chem, frac in baseline_feed_composition.items():
                feed.imass[chem] = baseline_feed_rate * frac
            process.simulate()
            return tea.net_earnings

    # 3. Find break-even point
    try:
        be_point = opt.root_scalar(simulate_profit, bracket=[min_scale_mt, 1976], method="bisect", xtol=0.1).root
    except ValueError:
        be_point = min_scale_mt

    # All points preserved
    low_scales = np.linspace(min_scale_mt, 250, 6)
    high_scales = np.arange(500, max_scale_mt + 1, 250)
    all_scales = np.unique(np.sort(np.concatenate(([be_point], low_scales, high_scales))))

    scale_list = [round(s, 1) for s in all_scales]
    profit_list = [round(simulate_profit(s) / 1e6, 2) for s in all_scales]

    return scale_list, profit_list, be_point


# %% best case vs worst case
def best_case():
    """
    Evaluates the Mechanical plant under the most favorable (optimistic) operating conditions
    derived from the upper/lower bounds in the sensitivity analysis tornado plot.
    """
    best_params = {
        'ndfeb_conc': 0.9,          # Max magnet concentration in impeller (0.9)
        'processing_capacity': 2965,               # Max processing capacity (2,965 tons/yr)
        'feed_price': 0.10,        # Min feedstock purchase price ($0.10/kg)
        'ndfeb_price': 150.0,           # Max magnet selling price ($150.00/kg)
        'worker_salary': 40.0,          # Min worker salary ($40k/yr)
        'freshwater_price': 0.0008,     # Min freshwater price ($0.0008/kg)
        'peraceticacid_price': 4.40,              # Min peracetic acid price ($4.40/kg)
        'wastewater_fee': 0.001,        # Min wastewater fee ($0.001/L)
        'othercomponent_fee': 0.01,     # Min disposal fee ($0.01/kg)
        'hdpe_fee': 0.01,                # Min HDPE disposal fee ($0.01/kg)
        'nahso3_price': 0.325,
        'naoh_price': 0.375
    }
    
    net_earnings = my_profit(**best_params)
    irr_val = tea.solve_IRR() * 100
    
    print(f"=== MECHANICAL PLANT BEST-CASE SCENARIO ===")
    print(f"Net Earnings: ${net_earnings / 1e6:.2f} MM/yr")
    print(f"Internal Rate of Return (IRR): {irr_val:.2f}%\n")
    
    return net_earnings, irr_val


def worst_case():
    """
    Evaluates the Mechanical plant under the most unfavorable (pessimistic) operating conditions
    derived from the upper/lower bounds in the sensitivity analysis tornado plot.
    """
    worst_params = {
        'ndfeb_conc': 0.8,          # Min magnet concentration in impeller (0.8)
        'processing_capacity': 990,                # Min processing capacity (990 tons/yr)
        'feed_price': 0.50,        # Max feedstock purchase price ($0.50/kg)
        'ndfeb_price': 50.0,            # Min magnet selling price ($50.00/kg)
        'worker_salary': 70.0,          # Max worker salary ($70k/yr)
        'freshwater_price': 0.002,      # Max freshwater price ($0.002/kg)
        'peraceticacid_price': 13.20,             # Max peracetic acid price ($13.20/kg)
        'wastewater_fee': 0.005,        # Max wastewater fee ($0.005/L)
        'othercomponent_fee': 0.10,     # Max disposal fee ($0.10/kg)
        'hdpe_fee': 0.10,                # Max HDPE disposal fee ($0.10/kg)
        'nahso3_price': 0.975,
        'naoh_price': 1.125
    }
    
    net_earnings = my_profit(**worst_params)
    irr_val = tea.solve_IRR() * 100
    
    print(f"=== MECHANICAL PLANT WORST-CASE SCENARIO ===")
    print(worst_params)
    print(f"Net Earnings: ${net_earnings / 1e6:.2f} MM/yr")
    print(f"Internal Rate of Return (IRR): {irr_val:.2f}%\n")
    
    return net_earnings, irr_val


def simulate_and_get_profit(use_mbbr):
    global process, tea, WWT_System, WM_Mixer, U1, U2, U3, S1, U4, U5, M1, M2, S2, feed
    global fresh_water_utility, fresh_paa_utility, nahso3_utility, naoh_utility
    global NdFeB_magnets, HDPE_casing, HDPE_shavings, other_components

    bst.main_flowsheet.clear()
    
    feed = bst.Stream(ID='Bioreactor Bags', units='kg/hr',
                       **{chem: feed_per_hour * frac for chem, frac in feed_composition.items()})
    fresh_water_utility = bst.Stream(ID='Fresh_water')
    fresh_paa_utility = bst.Stream(ID='Peracetic Acid (15%)')
    nahso3_utility = bst.Stream(ID='NaHSO3_utility')
    naoh_utility = bst.Stream(ID='NaOH_utility')

    U1 = Disinfection_Unit('U1', ins=(feed, fresh_water_utility, fresh_paa_utility,
                                      nahso3_utility, naoh_utility), 
                           outs=('clean_rinsed_bags', 'accidental_spill_loss', 'disinfectant_solution'))
    U2 = HandSorting_Unit('U2', ins=U1.outs[0], N_shifts=2, N_workers=4)
    U3 = Lathe_Machine('U3', ins=U2.outs[0], power_kw=3)

    S1 = bst.Splitter('S1', ins=U3.outs[0], outs=('to Arbor press 1', 'to Arbor press B'), split=0.5)
    U4 = Arbor_Press('U4', S1.outs[0])
    U5 = Arbor_Press('U5', ins=S1.outs[1])
    M1 = bst.Mixer('M1', ins=(U4.outs[0], U5.outs[0]), outs='Total NdFeB Magnets')
    M2 = bst.Mixer('M2', ins=(U4.outs[1], U5.outs[1]), outs='Total Waste HDPE casings')
    S2 = bst.Splitter('S2', ins=U2.outs[1], split=1)
    S2.isplit['Water'] = 0

    WM_Mixer = bst.Mixer('WM_Mixer', ins=(U1.outs[2], S2.outs[1]), outs='Combined_Wastewater')
    WWT_System = build_wastewater_treatment(WM_Mixer.outs[0], MBBR=use_mbbr)

    process = bst.System('mechanical_magnet_recovery', 
                         path=(U1, U2, U3, S1, U4, U5, M1, M2, S2, WM_Mixer, WWT_System))

    process.simulate()
    
    NdFeB_magnets = M1.outs[0]
    HDPE_casing = M2.outs[0] 
    HDPE_shavings = U3.outs[1]
    other_components = S2.outs[0]
    
    U3.outs[1].ID = 'Waste_HDPE_shavings'
    S2.outs[0].ID = 'Other_components'
    
    NdFeB_magnets.price = 100               
    HDPE_casing.price = -0.072            
    HDPE_shavings.price = -0.072            
    other_components.price = -0.072            

    feed.price = 0.25
    fresh_water_utility.price = 0.0015
    fresh_paa_utility.price = 8.8
    nahso3_utility.price = 0.65
    naoh_utility.price = 0.75

    if use_mbbr:
        WWT_System.ins[1].price = 0.60
        WWT_System.ins[2].price = 0.75
        WWT_System.outs[0].price = -0.0036
        WWT_System.outs[1].price = -0.072
    else:
        WWT_System.outs[0].price = -(0.0036 + wastewater_bod_fee(WWT_System.outs[0])) 

    total_labor_cost = U1.total_salary + U2.total_salary + U4.total_salary + U5.total_salary + U3.total_salary
    adjusted_days = operating_hours / 24

    tea = MechanicalRecyclingTEA(
        system=process,
        IRR=0.15,
        duration=(2026, 2046),
        depreciation='MACRS7',
        income_tax=0.21, 
        operating_days=adjusted_days,
        lang_factor=None,
        construction_schedule=(1.0,),
        WC_over_FCI=0.25,
        labor_cost=total_labor_cost,
        fringe_benefits=0.4,
        property_tax=0.001,
        property_insurance=0.005,
        supplies=0.05,
        maintenance=0.03,
        administration=0.005
    )

    return tea.net_earnings

def compare_wwt_profits(target_mbbr=MBBR):
    print("Simulating process WITH on-site MBBR...")
    profit_mbbr = simulate_and_get_profit(use_mbbr=True)
    
    print("Simulating process WITHOUT MBBR (Sewer Surcharge)...")
    profit_sewer = simulate_and_get_profit(use_mbbr=False)
    
    # Finalize state according to user target setting
    simulate_and_get_profit(use_mbbr=target_mbbr)

    print("-" * 45)
    print(f"Profit WITH MBBR:    ${profit_mbbr / 1e6:,.2f} / year")
    print(f"Profit WITHOUT MBBR: ${profit_sewer / 1e6:,.2f} / year")
    
    return {"MBBR_Profit": profit_mbbr, "Sewer_Profit": profit_sewer}

results = compare_wwt_profits()


#%%defining plotting functions
def tornado_plot():
    """Runs economic sensitivity analysis and generates tornado plots for both
    Net Earnings and IRR simultaneously to save computation time.
    """
    is_mbbr_active = isinstance(WWT_System, MBBR_Tank)
    
    # Fetch the correct baseline wastewater fee based on WWT configuration
    if is_mbbr_active:
        base_ww_fee = 0.0036
    else:
        base_ww_fee = 0.003 + wastewater_bod_fee(WWT_System.outs[0])

    # Core parameters
    base_params = {
        "processing_capacity": 1976,
        "ndfeb_price": 100,
        "ndfeb_conc": 0.86,
        "feed_price": 0.25,
        "peraceticacid_price": 8.8,
        "nahso3_price": 0.65,
        "naoh_price": 0.75,
        "freshwater_price": 0.0015,
        "hdpe_fee": 0.072,
        "wastewater_fee": base_ww_fee,
        "othercomponent_fee": 0.072,
        "worker_salary": 5e4,
    }

    ranges = {
        "processing_capacity": (990, 2965),
        "ndfeb_price": (50, 150),
        "ndfeb_conc": (0.80, 0.90),
        "feed_price": (0.1, 0.5),
        "peraceticacid_price": (4.4, 13.2),
        "nahso3_price": (0.325, 0.975),
        "naoh_price": (0.375, 1.125),
        "freshwater_price": (0.0008, 0.002),
        "hdpe_fee": (0.01, 0.1),
        "wastewater_fee": (0.5 * base_ww_fee, round(1.5 * base_ww_fee, 4)),
        "othercomponent_fee": (0.01, 0.1),
        "worker_salary": (4e4, 7e4),
    }

    label_map = {
        "processing_capacity": "Capacity",
        "ndfeb_price": "NdFeB Price ($/kg)",
        "ndfeb_conc": "NdFeB Concentration in Impeller",
        "feed_price": "Feedstock Price ($/kg)",
        "peraceticacid_price": "Peracetic Acid Price ($/kg)",
        "freshwater_price": "Freshwater Price ($/kg)",
        "nahso3_price": "NaHSO3 Price ($/kg)",
        "naoh_price": "NaOH Price ($/kg)",
        "hdpe_fee": "HDPE Disposal Fee ($/kg)",
        "wastewater_fee": "Wastewater Fee ($/L)",
        "othercomponent_fee": "Other Components Fee ($/kg)",
        "worker_salary": "Worker Salary (in k's)",
    }

    # Conditionally inject Urea and H3PO4 if MBBR is true
    if is_mbbr_active:
        base_params["urea_price"] = 0.60
        base_params["h3po4_price"] = 0.75
        
        ranges["urea_price"] = (0.30, 0.90)      # +/- 50%
        ranges["h3po4_price"] = (0.375, 1.125)   # +/- 50%
        
        label_map["urea_price"] = "Urea Price ($/kg)"
        label_map["h3po4_price"] = "H3PO4 Price ($/kg)"

    results = []
    results_irr = []

    # Get Baseline (Simulates once, grabs both metrics)
    base_profit_m = my_profit(**base_params) / 1e6
    base_irr_m = tea.solve_IRR() * 100

    for param, (low, high) in ranges.items():
        # Test Low Bound
        p_low = base_params.copy()
        p_low[param] = low
        low_profit_m = my_profit(**p_low) / 1e6
        low_irr = tea.solve_IRR() * 100

        # Test High Bound
        p_high = base_params.copy()
        p_high[param] = high
        high_profit_m = my_profit(**p_high) / 1e6
        high_irr = tea.solve_IRR() * 100

        display_low = low / 1000.0 if param == "worker_salary" else low
        display_high = high / 1000.0 if param == "worker_salary" else high

        # Append Profit Results
        results.append((
            label_map[param], low_profit_m, high_profit_m,
            abs(high_profit_m - low_profit_m), display_low, display_high,
        ))
        
        # Append IRR Results
        results_irr.append((
            label_map[param], low_irr, high_irr,
            abs(high_irr - low_irr), display_low, display_high,
        ))

    # Render Chart Function
    def _render_chart(data, base_val, title, x_label):
        data.sort(key=lambda x: x[3], reverse=True)
        names = [d[0] for d in data]
        lows = [d[1] for d in data]
        highs = [d[2] for d in data]
        p_lows = [d[4] for d in data]
        p_highs = [d[5] for d in data]

        plt.figure(figsize=(11, 7))
        
        # Handle inverted metrics (e.g., lower fee = higher profit)
        for i in range(len(names)):
            if lows[i] < base_val:
                plt.barh(names[i], base_val - lows[i], left=lows[i], color="red", edgecolor="black", alpha=0.8)
                plt.barh(names[i], highs[i] - base_val, left=base_val, color="green", edgecolor="black", alpha=0.8)
            else:
                plt.barh(names[i], lows[i] - base_val, left=base_val, color="green", edgecolor="black", alpha=0.8)
                plt.barh(names[i], base_val - highs[i], left=highs[i], color="red", edgecolor="black", alpha=0.8)

        plt.axvline(base_val, color="blue", linestyle="--", linewidth=2.0, zorder=3)

        x_min, x_max = plt.xlim()
        x_gap = 0.03 * (x_max - x_min)

        for i, (lo, hi) in enumerate(zip(p_lows, p_highs)):
            if lows[i] < base_val:
                plt.text(lows[i] - x_gap, i, f"{lo}", ha="right", va="center", fontsize=11, fontweight="bold", color="grey")
                plt.text(highs[i] + x_gap, i, f"{hi}", ha="left", va="center", fontsize=11, fontweight="bold", color="grey")
            else:
                plt.text(highs[i] - x_gap, i, f"{hi}", ha="right", va="center", fontsize=11, fontweight="bold", color="grey")
                plt.text(lows[i] + x_gap, i, f"{lo}", ha="left", va="center", fontsize=11, fontweight="bold", color="grey")

        all_points = lows + highs + [base_val]
        padding = (max(all_points) - min(all_points)) * 0.25
        plt.xlim(min(all_points) - padding, max(all_points) + padding)

        ax = plt.gca()
        ax.set_title(title, fontsize=15, fontweight="bold", pad=15)
        plt.xlabel(x_label, fontsize=13, fontweight="bold")
        plt.xticks(fontsize=11, fontweight="bold")
        plt.yticks(fontsize=11, fontweight="bold")
        
        for spine in ax.spines.values():
            spine.set_linewidth(2)

        red_patch = mpatches.Patch(color="red", label="Worse than Baseline")
        green_patch = mpatches.Patch(color="green", label="Better than Baseline")
        baseline_legend = mlines.Line2D([], [], color="blue", linestyle="--", label=f"Baseline: {base_val:.2f}")

        leg = plt.legend(handles=[red_patch, green_patch, baseline_legend], prop={"weight": "bold", "size": 11}, loc="best")
        leg.get_frame().set_edgecolor("black")
        leg.get_frame().set_linewidth(2)

        plt.tight_layout()
        plt.show()

    # Generate both plots
    _render_chart(results, base_profit_m, "Mechanical Method Sensitivity: Net Earnings", "Profit ($ Millions)")
    _render_chart(results_irr, base_irr_m, "Mechanical Method Sensitivity: IRR", "IRR (%)")

def plot_profit_vs_scale():
    """Plots Net Earnings vs. Processing Capacity.
    Automatically runs get_capacity_curve() internally with zero arguments.
    """
    
    # 1. Fetch capacity curve data internally
    scale_mt_list, profit_list, be_point = get_capacity_curve()
    
    print(f"==================================================")
    print(f"Exact Breakeven Capacity: {be_point:.1f} Metric Tonnes/yr")
    print(f"==================================================")

    impeller_frac = 1.0  # Default 100% processing scale basis

    # 2. Scale capacity list by impeller_frac
    scaled_scales = [s * impeller_frac for s in scale_mt_list]

    # 3. Deduplicate points that round to the same integer label
    dedup_data = {}
    for scale, profit in zip(scaled_scales, profit_list):
        label_key = int(round(scale))
        if label_key not in dedup_data:
            dedup_data[label_key] = profit

    x_labels = [str(k) for k in dedup_data.keys()]
    clean_profits = list(dedup_data.values())

    # 4. Setup figure
    plt.close("all")
    fig, ax = plt.subplots(figsize=(12, 6), dpi=300)

    colors = [
        "#8b0000" if p < 0 else plt.cm.Greens(0.3 + 0.7 * (p / max(clean_profits)))
        for p in clean_profits
    ]
    bars = ax.bar(x_labels, clean_profits, color=colors, width=0.5)
    ax.axhline(0, color="black", linestyle="--", linewidth=1)

    max_val = (
        max(abs(np.array(clean_profits))) if len(clean_profits) > 0 else 1.0
    )

    for bar, profit in zip(bars, clean_profits):
        yval = bar.get_height()
        if profit >= 0:
            va_align = "bottom"
            offset = 0.02 * max_val
            text_color = "#005a00"
        else:
            va_align = "top"
            offset = -0.02 * max_val
            text_color = "#8b0000"

        ax.text(
            bar.get_x() + bar.get_width() / 2.0,
            yval + offset,
            f"${profit:.2f}M",
            ha="center",
            va=va_align,
            fontsize=8,
            fontweight="bold",
            color=text_color,
        )

    title = "Mechanical Recovery Profitability"
    x_label = "Bioreactor Bags Feedstock [Metric Tonnes / yr]"

    ax.set_title(title, fontsize=14, fontweight="bold", pad=15)
    ax.set_xlabel(x_label, fontsize=11, fontweight="bold", labelpad=10)
    ax.set_ylabel(
        "Net Earnings [MM$ / yr]", fontsize=11, fontweight="bold", labelpad=10
    )

    ax.set_ylim(bottom=-5, top=25)
    plt.xticks(rotation=45, ha="right", fontweight="bold")
    plt.yticks(fontweight="bold")

    for spine in ax.spines.values():
        spine.set_linewidth(1.5)
        spine.set_color("black")

    plt.tight_layout()
    plt.show()

def plot_profit_line_mechanical():
    """
    Plots Net Earnings vs. Capacity for Mechanical Recovery.
    Uses multi-tiered vertical offsets for dense points (<= 250 MT) to guarantee
    clear spacing between adjacent labels.
    """
    scale_list, profit_list, be_point = get_capacity_curve()
    
    print(f"==================================================")
    print(f"Exact Breakeven Capacity: {be_point:.1f} Metric Tonnes/yr")
    print(f"==================================================")

    plt.figure(figsize=(12, 6.5), dpi=300)
    ax = plt.gca()

    # Plot data line & break-even reference line
    plt.plot(
        scale_list, 
        profit_list, 
        marker='o', 
        color='#005b96', 
        linewidth=2.5, 
        markersize=6, 
        label='Net Earnings'
    )
    plt.axhline(0, color='black', linewidth=1, linestyle='--')

    # No background grid
    ax.grid(False)

    y_range = max(profit_list) - min(profit_list)
    if y_range == 0: y_range = 1.0

    # MULTI-TIERED OFFSETS (Up-Low, Down-Deep, Up-High, Down-Shallow, ...)
    # Alternates both direction AND height to prevent horizontal label collisions
    dense_offsets = [0.09, -0.22, 0.17, -0.11, 0.09, -0.22, 0.17]
    low_scale_idx = 0

    # Annotations
    for x, y in zip(scale_list, profit_list):
        text_color = '#005a00' if y >= 0 else '#8b0000'
        
        if x <= 250:
            # Multi-tiered vertical callouts
            offset = dense_offsets[low_scale_idx % len(dense_offsets)] * y_range
            low_scale_idx += 1
            
            plt.annotate(
                f"${y:.2f}M",
                xy=(x, y),
                xytext=(x, y + offset),
                ha='center',
                va='center',
                fontsize=8,
                fontweight='bold',
                color=text_color,
                arrowprops=dict(
                    arrowstyle='-',
                    color='gray',
                    lw=0.8,
                    alpha=0.7
                )
            )
        else:
            # Direct top positioning for wide 250 MT steps
            plt.text(
                x, 
                y + 0.04 * y_range, 
                f"${y:.2f}M", 
                ha='center', 
                va='bottom', 
                fontsize=8.5, 
                fontweight='bold', 
                color=text_color
            )

    # Dynamic Y-Limits: Extra bottom padding so deep downward callout labels have ample space
    plt.ylim(bottom=min(profit_list) - 0.32 * y_range, top=max(profit_list) + 0.14 * y_range)

    # Ticks & Formatting
    x_ticks = [0, 250, 500, 750, 1000, 1250, 1500, 1750, 2000]
    plt.xticks(x_ticks, [f"{x:,}" for x in x_ticks], fontsize=10, fontweight='bold')
    plt.yticks(fontsize=10, fontweight='bold')

    plt.xlabel('Feedstock Capacity [Metric Tonnes / yr]', fontsize=11, fontweight='bold', labelpad=10)
    plt.ylabel('Net Earnings [MM$ / yr]', fontsize=11, fontweight='bold', labelpad=10)
    plt.title('Mechanical Recovery Plant Profitability vs. Scale', fontsize=13, fontweight='bold', pad=15)

    # Frame styling
    for spine in ax.spines.values():
        spine.set_linewidth(1.5)
        spine.set_color('black')

    plt.tight_layout()
    plt.show()

#%%LCA
#%% LCA Characterizatio n###
GWP = 'GWP'
FFC = 'FFC'
WU = 'WU'
HTC = 'HTC'
HTNC = 'HTNC'
ETOX = 'ETOX'
ACD = 'ACD'
OZD = 'OZD'
POCP = 'POCP'

#%%CFs 

# CF for peracetic acid
process.feeds[2].ID = "peracetic_acid (15%)"
process.feeds[2].set_CF(GWP, 0.0)

#CFs for freshwater
process.feeds[1].ID = 'fresh_water'
process.feeds[1].set_CF(GWP,0.0)

#CFs for the bioreactor bags
process.feeds[0].ID = 'bioreactor bags'
process.feeds[0].set_CF(GWP,0.0)

#CF for the power
bst.PowerUtility.set_CF('GWP',0.42)

#CF for NaHSO3
process.feeds[3].ID = "NaHSO3 (40%)"
process.feeds[3].set_CF(GWP, 0.0)

#CF for NaOH
process.feeds[4].ID = "NaOH (50%)"
process.feeds[4].set_CF(GWP, 0.0)

inventory_table = bst.report.lca_inventory_table(systems =[process],keys='GWP', items=process.products)