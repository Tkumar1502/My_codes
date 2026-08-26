import biosteam as bst
import importlib
import types
import pandas as pd
import numpy as np
import scipy.optimize as opt
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
from warnings import filterwarnings

bst.settings.ID_magic = False
from plastics import strap

importlib.reload(strap.property_package)
importlib.reload(strap.process_model)
filterwarnings('ignore')

original_u7_design = strap.units.Precipitator._design

# =============================================================================
# 1. GLOBAL CONFIGURATION & MASTER SWITCHES
# =============================================================================
#%% CHECK MBBR
# BRAND NEW VARIABLE: Overcomes IDE caching issues
ENABLE_MBBR_WWT = True     # True = On-site MBBR WWT; False = Municipal sewer discharge
print("MBBR Enabled:", ENABLE_MBBR_WWT)

capacity = 1976           # MT/yr feedstock capacity
bst.settings.CEPCI = 836.9
elec_price = bst.settings.electricity_price

N_operators = 1
N_plant_manager = 1
plant_manager_salary = 110000
operator_salary = 80000

feed_composition = {
    'NdFeB': 0.1548,
    'HDPE': 0.0252,
    'Films': 0.10,
    'FittingsFilters': 0.35,
    'BrownSupport': 0.07,
    'SiliconeTubings': 0.28,
    'Solutes': 0.001,
    'BiogenicResidue': 0.02
}

INDIVIDUAL_UNIT_POWER_KW = {
    'U10': 1.06, 'U3': 1.0, 'T2': 0.0, 'T3': 0.0, 'T4': 0.0,
    'P1': 0.01, 'P2': 0.3, 'P3': 0.10, 'U6': 1.0, 'U7': 0.0,
    'U8': 4.15, 'F1': 5.5, 'U9': 18.5, 'U2': 0.10, 'Vac_S': 0.75,
    'H1': 0.15, 'H2': 0.10, 'H3': 0.00, 'CT': 0.1
}

# =============================================================================
# 2. CUSTOM TEA CLASS
# =============================================================================
class SimplifiedTEA(bst.TEA):
    """Custom TEA class that overrides default BioSTEAM CAPEX calculations."""
    def __init__(self, system, IRR, duration, depreciation, operating_days, income_tax, 
                 lang_factor=None, construction_schedule=(1.0,), labor_cost=0.0, 
                 fringe_benefits=0.4, property_tax=0.001, WC_over_FCI=0.1, 
                 property_insurance=0.005, supplies=0.05, maintenance=0.03, 
                 administration=0.005, **kwargs):
        
        if not hasattr(system, 'installed_equipment_cost'):
            if hasattr(system, 'sys'): system = system.sys
            elif hasattr(system, 'system'): system = system.system
        
        super().__init__(system, IRR, duration, depreciation, income_tax, operating_days,
                         lang_factor, construction_schedule, startup_months=0, startup_FOCfrac=0,
                         startup_VOCfrac=0, startup_salesfrac=0, finance_interest=0, finance_years=0,
                         finance_fraction=0, WC_over_FCI=WC_over_FCI)
        
        self.labor_cost = labor_cost
        self.fringe_benefits = fringe_benefits
        self.property_tax = property_tax
        self.property_insurance = property_insurance
        self.supplies = supplies
        self.maintenance = maintenance
        self.administration = administration
        
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

    @property
    def ISBL(self): return self.installed_equipment_cost
    @property
    def OSBL(self): return sum(self.osbl_itemized.values())
    def _DPI(self, installed_equipment_cost): return self.ISBL + self.OSBL
    def _TDC(self, DPI): return DPI
    def _FCI(self, TDC): return TDC + 0.90 * TDC
    def _FOC(self, FCI): 
        return (FCI * (self.property_tax + self.property_insurance + self.maintenance + self.administration)
                + self.labor_cost * (1 + self.fringe_benefits + self.supplies))

    @property
    def engineering_cost(self): return 0.50 * self.TDC
    @property
    def contingency(self): return 0.40 * self.TDC
    @property
    def TIC(self): return self.engineering_cost + self.contingency
    @property
    def WC(self): return self.WC_over_FCI * self.FCI
    @property
    def TCI(self): return self.FCI + self.WC

    def OSBL_table(self, formatted=True):
        categories = list(self.osbl_itemized.keys()) + ["TOTAL DIRECT ASSET & OSBL COST"]
        costs = list(self.osbl_itemized.values()) + [self.OSBL]
        if formatted:
            return pd.DataFrame({"Cost Category": categories, "Cost ($)": [f"${c:,.0f}" for c in costs]}).set_index("Cost Category")
        return pd.DataFrame({"Cost Category": categories, "Cost ($)": costs}).set_index("Cost Category")

    def CAPEX_table(self):
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
        notes = ['', '', '', '50.0% of TDC', '40.0% of TDC', '', 'TDC + TIC', f'{self.WC_over_FCI * 100:.1f}% of FCI', 'FCI + WC']
        costs = [self.ISBL/1e6, self.OSBL/1e6, self.TDC/1e6, self.engineering_cost/1e6, self.contingency/1e6, self.TIC/1e6, self.FCI/1e6, self.WC/1e6, self.TCI/1e6]
        return pd.DataFrame({'Notes': notes, 'Cost [MM$]': [round(c, 3) for c in costs]}, index=pd.MultiIndex.from_tuples(index))


# =============================================================================
# 3. UNIFIED PROCESS BUILDER
# =============================================================================
def build_process(mbbr_flag=ENABLE_MBBR_WWT, proc_capacity=capacity):
    global process, tea, feedstock, NdFeB_Magnets, HDPE_resins, wastewater
    global othercomponents, freshwater, peracetic_acid, nahso3, naoh, solvent

    bst.main_flowsheet.clear()
    if hasattr(strap.MagnetRecovery, 'cache'):
        strap.MagnetRecovery.cache.clear()

    def safe_u7_design(self):
        total_vol = sum(i.F_vol for i in self.ins if not i.isempty())
        if total_vol < 1e-6:
            self.design_results.update({'Vessel volume': 0.0, 'Number of reactors': 0, 'Weight': 0.0, 'Diameter': 0.0, 'Length': 0.0})
        else:
            try: original_u7_design(self)
            except ZeroDivisionError:
                self.design_results.update({'Vessel volume': 0.0, 'Number of reactors': 0, 'Weight': 0.0, 'Diameter': 0.0, 'Length': 0.0})

    strap.units.Precipitator._design = safe_u7_design

    process = strap.MagnetRecovery(
        processing_capacity=proc_capacity,
        sell_leftover_plastic=True,
        simulate=False,
        mbbr=mbbr_flag
    )

    operating_hours = 309 * 16
    process.tea.operating_hours = operating_hours
    process.tea.operating_days = 309
    feed_per_hour = proc_capacity * 1000 / operating_hours

    feedstock = process.feedstock
    feedstock.empty()
    for chem, frac in feed_composition.items():
        feedstock.imass[chem] = feed_per_hour * frac

    strap_labor = N_plant_manager * plant_manager_salary + N_operators * operator_salary * 1.6 * 3
    total_labor_cost = strap_labor + process.U11.total_salary + process.U10.total_salary

    process.plastic.ID = 'Bioreactor bags'
    process.Vac_S.outs[0].ID = 'NdFeB_Magnets'
    process.Vac_S.outs[1].ID = 'Xylenes vapors'
    process.M2.disconnect()
    
    solvent = process.solvent
    solvent.ID = 'Xylenes'

    active_units = [unit for unit in process.system.units if unit.ID != 'M2']
    process.system = bst.System(ID='sys', path=active_units, facilities=process.system.facilities)

    process.s22.ID = 'Impurities'
    NdFeB_Magnets = process.Vac_S.outs[0]
    HDPE_resins = process.U9.outs[0]
    process.products[:] = [HDPE_resins, NdFeB_Magnets]

    HDPE_resins.price = 1.20
    NdFeB_Magnets.price = 100

    othercomponents = process.S3.outs[0]
    othercomponents.price = -0.072
    process.S3.outs[1].price = 0.0

    freshwater = process.U10.ins[1]
    peracetic_acid = process.U10.ins[2]
    nahso3 = process.U10.ins[3]
    naoh = process.U10.ins[4]

    feedstock.price = 0.25
    freshwater.price = 0.0015
    peracetic_acid.price = 8.8
    nahso3.price = 0.65
    naoh.price = 0.75
    process.set_solvent_price(0.85)

    def c_u3(self):
        self.baseline_purchase_costs.clear()
        self.baseline_purchase_costs['Vertical pressure vessel'] = 32352.9
        self.power_utility.rate = 0.0203 / elec_price
        self.BM = 3.4

    def c_t2(self):
        self.baseline_purchase_costs.clear()
        if hasattr(self, 'baseline_item_bms'): self.baseline_item_bms.clear()
        self.baseline_purchase_costs['Tank'] = 15000.0
        self.power_utility.rate = 0.0

    def c_t3(self):
        self.baseline_purchase_costs.clear()
        self.baseline_purchase_costs['Tank system'] = 15000.0
        self.F_BM['Tank system'] = 2.3
        self.power_utility.rate = 0.0

    def c_t4(self):
        self.baseline_purchase_costs.clear()
        if hasattr(self, 'baseline_item_bms'): self.baseline_item_bms.clear()
        self.baseline_purchase_costs['Tank'] = 20000.0
        self.power_utility.rate = 0.0

    def c_p1(self):
        self.baseline_purchase_costs.clear()
        self.baseline_purchase_costs['Pump'] = 8000.0
        self.power_utility.rate = 2.9e-6 / elec_price

    def c_p2(self):
        self.baseline_purchase_costs.clear()
        self.baseline_purchase_costs['Pump system'] = 8000.0
        self.F_BM['Pump system'] = 3.3
        self.power_utility.rate = 0.0034 / elec_price

    def c_p3(self):
        self.baseline_purchase_costs.clear()
        self.baseline_purchase_costs['Pump'] = 0.0
        self.power_utility.rate = 0.0270 / elec_price

    def c_u6(self):
        self.baseline_purchase_costs.clear()
        self.baseline_purchase_costs['combined U6 and P3'] = 250000.0
        self.F_BM['combined U6 and P3'] = 3.0
        self.power_utility.rate = 0.0021 / elec_price

    def c_u7(self):
        self.baseline_purchase_costs.clear()
        if hasattr(self, 'baseline_item_bms'): self.baseline_item_bms.clear()
        self.baseline_purchase_costs['Precipitator'] = 175000.0
        self.F_BM['Precipitator'] = 3.46
        self.heat_utilities.clear()
        self.power_utility.rate = 0.0133 / elec_price

    def c_u8(self):
        self.baseline_purchase_costs.clear()
        self.baseline_purchase_costs['Precipitator'] = 0.0
        self.power_utility.rate = 0.2888 / elec_price

    def c_u9(self):
        self.baseline_purchase_costs.clear()
        if hasattr(self, 'baseline_item_bms'): self.baseline_item_bms.clear()
        self.baseline_purchase_costs['Screw degasser'] = 250000.0
        self.F_BM['Screw degasser'] = 2.60
        self.power_utility.rate = 1.3061 / elec_price

    def c_h1(self):
        self.baseline_purchase_costs.clear()
        self.baseline_purchase_costs['Double pipe'] = 8000.0 * 1.776
        self.heat_utilities.clear()
        self.power_utility.rate = 0.0146 / elec_price

    def c_h2(self):
        self.baseline_purchase_costs.clear()
        if hasattr(self, 'baseline_item_bms'): self.baseline_item_bms.clear()
        self.baseline_purchase_costs['Chiller'] = 50000.0
        self.F_BM['Chiller'] = 2.5
        self.heat_utilities.clear()
        self.power_utility.rate = 4e-6 / elec_price

    def c_h3(self):
        self.baseline_purchase_costs.clear()
        self.baseline_purchase_costs['Double pipe'] = 125000.0 / 1.768
        self.power_utility.rate = 0.0

    def c_zero(self):
        self.baseline_purchase_costs.clear()
        self.power_utility.rate = 0.0

    process.U3._cost = types.MethodType(c_u3, process.U3)
    process.T2._cost = types.MethodType(c_t2, process.T2)
    process.T3._cost = types.MethodType(c_t3, process.T3)
    process.T4._cost = types.MethodType(c_t4, process.T4)
    process.P1._cost = types.MethodType(c_p1, process.P1)
    process.P2._cost = types.MethodType(c_p2, process.P2)
    process.P3._cost = types.MethodType(c_p3, process.P3)
    process.U6._cost = types.MethodType(c_u6, process.U6)
    process.U7._cost = types.MethodType(c_u7, process.U7)
    process.U8._cost = types.MethodType(c_u8, process.U8)
    process.U9._cost = types.MethodType(c_u9, process.U9)
    process.H1._cost = types.MethodType(c_h1, process.H1)
    process.H2._cost = types.MethodType(c_h2, process.H2)
    process.H3._cost = types.MethodType(c_h3, process.H3)
    process.CWP._cost = types.MethodType(c_zero, process.CWP)
    process.CT._cost = types.MethodType(c_zero, process.CT)

    for unit_id, target_kw in INDIVIDUAL_UNIT_POWER_KW.items():
        if hasattr(process, unit_id):
            u = getattr(process, unit_id)
            u.power_utility.consumption = target_kw

    tea = SimplifiedTEA(
        system=process.system,
        IRR=0.15,
        duration=(2025, 2045),
        depreciation='MACRS7',
        income_tax=0.21,
        operating_days=operating_hours / 24,
        labor_cost=total_labor_cost,
        WC_over_FCI=0.10
    )
    process.tea = tea

    process.system.simulate()
    
    volumetric_rate = 0.0036 

    if mbbr_flag:
        wastewater = process.MBBR_System.outs[0]
        wastewater.price = -volumetric_rate
        if hasattr(process, 'Sewer_Discharge'):
            process.Sewer_Discharge.outs[0].price = 0.0
    else:
        wastewater = process.Sewer_Discharge.outs[0]
        bod_fee = wastewater_bod_fee(wastewater)
        wastewater.price = -(volumetric_rate + bod_fee)

    return process


# =============================================================================
# 4. WWT COMPARISON & SENSITIVITY FUNCTIONS
# =============================================================================
def run_wwt_comparison(target_mbbr=ENABLE_MBBR_WWT):
    """Runs a side-by-side economic comparison between MBBR and Sewer Discharge."""
    print("Simulating process WITH on-site MBBR...")
    p_mbbr = build_process(mbbr_flag=True, proc_capacity=capacity)
    profit_mbbr = p_mbbr.tea.net_earnings

    print("Simulating process WITHOUT MBBR (Sewer Surcharge)...")
    p_sewer = build_process(mbbr_flag=False, proc_capacity=capacity)
    profit_sewer = p_sewer.tea.net_earnings

    build_process(mbbr_flag=target_mbbr, proc_capacity=capacity)

    print("=" * 55)
    print("       WASTEWATER TREATMENT COMPARISON RESULTS       ")
    print("=" * 55)
    print(f"Profit WITH MBBR:    ${profit_mbbr / 1e6:,.2f} MM/yr")
    print(f"Profit WITHOUT MBBR: ${profit_sewer / 1e6:,.2f} MM/yr")
    print("-" * 55)

    return {"MBBR_Profit": profit_mbbr, "Sewer_Profit": profit_sewer}


def profit(u3_ndfeb_ratio, capacity, solvent_loss, plastic_conc, solvent_price, feedstock_price, 
           NdFeB_price, HDPE_price, freshwater_price, paa_price, wastewater_fee, othercomponent_fee,
           naoh_price=0.75, nahso3_price=0.65, urea_price=0.50, h3po4_price=1.20, **kwargs):
    """Evaluates process net earnings given varying sensitivity parameters."""
    process.set_processing_capacity(capacity)
    process.set_solvent_loss(solvent_loss)
    process.set_dissolution_capacity(plastic_conc)
    process.set_solvent_price(solvent_price)
    
    feedstock.price = feedstock_price
    NdFeB_Magnets.price = NdFeB_price
    HDPE_resins.price = HDPE_price
    freshwater.price = freshwater_price
    peracetic_acid.price = paa_price
    naoh.price = naoh_price
    nahso3.price = nahso3_price
    wastewater.price = -wastewater_fee
    process.S3.outs[1].price = -wastewater_fee
    othercomponents.price = -othercomponent_fee
    
    if ENABLE_MBBR_WWT:
        if hasattr(process, 'urea'): process.urea.price = urea_price
        if hasattr(process, 'h3po4'): process.h3po4.price = h3po4_price

    work_hours = capacity * 1000 / (4 * 100)
    process.tea.operating_hours = work_hours
    feed_per_hour = capacity * 1000 / work_hours

    feedstock.imass['NdFeB'] = feed_per_hour * (0.18 * u3_ndfeb_ratio)
    feedstock.imass['HDPE'] = feed_per_hour * (0.18 * (1.0 - u3_ndfeb_ratio))
    feedstock.imass['Films'] = feed_per_hour * 0.10
    feedstock.imass['FittingsFilters'] = feed_per_hour * 0.36
    feedstock.imass['BrownSupport'] = feed_per_hour * 0.07
    feedstock.imass['SiliconeTubings'] = feed_per_hour * 0.29

    process.system.simulate()
    return process.tea.net_earnings


def tornado_plot():
    base_params = {
        'u3_ndfeb_ratio': 0.86, 'capacity': 1976, 'solvent_loss': 0.1, 'plastic_conc': 5.0, 
        'solvent_price': 0.85, 'feedstock_price': 0.25, 'NdFeB_price': 100, 'HDPE_price': 1.20,
        'freshwater_price': 0.002, 'wastewater_fee': 0.0036 + wastewater_bod_fee(), 
        'paa_price': 8.8, 'othercomponent_fee': 0.072,
        'naoh_price': 0.75, 'nahso3_price': 0.65
    }
    
    ranges = {
        'u3_ndfeb_ratio': (0.80, 0.90), 'capacity': (990, 2965), 'solvent_loss': (0.01, 1.0),
        'plastic_conc': (2.0, 7.5), 'solvent_price': (0.1, 1.6), 'feedstock_price': (0.1, 0.50),
        'NdFeB_price': (50.0, 150.0), 'HDPE_price': (1.00, 1.50), 'freshwater_price': (0.001, 0.01), 
        'wastewater_fee': (0.5*(0.0036 + wastewater_bod_fee()), 1.5*(0.0036 + wastewater_bod_fee())), 
        'paa_price': (4.4, 13.2), 'othercomponent_fee': (0.01, 0.1),
        'naoh_price': (0.375, 1.125), 'nahso3_price': (0.325, 0.975)
    }
    
    label_map = {
        'u3_ndfeb_ratio': 'NdFeB Concentration in Impeller', 'capacity': 'Processing Capacity (tons)',
        'solvent_loss': 'Solvent Loss Rate', 'plastic_conc': 'Feed-to-Solvent (wt%)',
        'solvent_price': 'Solvent Price ($/kg)', 'feedstock_price': 'Feedstock Price ($/kg)',
        'NdFeB_price': 'Neodymium Magnet Price ($/kg)', 'HDPE_price': 'HDPE Price ($/kg)',
        'freshwater_price': 'Freshwater Price ($/kg)', 'wastewater_fee': 'Wastewater Fee ($/L)',
        'paa_price': 'Peracetic Acid Price ($/kg)', 'othercomponent_fee': 'Othercomponent Fee ($/kg)',
        'naoh_price': 'Sodium Hydroxide (NaOH) Price ($/kg)',
        'nahso3_price': 'Sodium Bisulfite (NaHSO3) Price ($/kg)'
    }

    if ENABLE_MBBR_WWT:
        base_params.update({'urea_price': 0.50, 'h3po4_price': 1.20})
        ranges.update({
            'urea_price': (0.25, 0.75),
            'h3po4_price': (0.60, 1.80)
        })
        label_map.update({
            'urea_price': 'Urea Price ($/kg)',
            'h3po4_price': 'Phosphoric Acid (H3PO4) Price ($/kg)'
        })
    
    results, results_irr = [], []
    base_profit_m = profit(**base_params) / 1e6
    base_irr_m = my_irr(**base_params)
  
    for param, (low, high) in ranges.items():
        p_low = base_params.copy(); p_low[param] = low
        p_high = base_params.copy(); p_high[param] = high

        low_profit_m, high_profit_m = profit(**p_low) / 1e6, profit(**p_high) / 1e6
        low_irr, high_irr = my_irr(**p_low), my_irr(**p_high)
        
        results.append((label_map[param], low_profit_m, high_profit_m, abs(high_profit_m - low_profit_m), low, high))
        results_irr.append((label_map[param], low_irr, high_irr, abs(high_irr - low_irr), low, high))
        
    plot_it(results, base_profit_m, "Sensitivity: Net Earnings", "Profit ($ Millions)", "skyblue")
    plot_it(results_irr, base_irr_m, "Sensitivity: IRR", "IRR (%)", "teal")


def my_irr(**kwargs):
    """Evaluates process IRR given varying sensitivity parameters."""
    profit(**kwargs)
    return process.tea.solve_IRR() * 100


def wastewater_bod_fee(stream=None):
    """Calculates municipal BOD surcharge rate ($/kg wastewater)."""
    if stream is None:
        if not hasattr(process, 'Sewer_Discharge'):
            return 0.0
        stream = process.Sewer_Discharge.outs[0]

    flow_rate_kg_hr = stream.F_mass
    if flow_rate_kg_hr < 1e-6:
        return 0.0

    bod_rate = 1.02  # LASAN $/lb BOD
    bod_factors = {
        'CH3COONa': (64.0 / 82.0) * 0.85,
        'BiogenicResidue': (160.0 / 113.0) * 0.65,
    }

    total_bod_kg_hr = sum(
        stream.imass[comp] * factor 
        for comp, factor in bod_factors.items() 
        if comp in stream.chemicals
    )
    
    bod_lbs_hr = total_bod_kg_hr * 2.20462
    return round((bod_lbs_hr * bod_rate) / flow_rate_kg_hr, 4)


def best_case_scenario():
    best_params = {
        'u3_ndfeb_ratio': 0.9, 'capacity': 2965, 'solvent_loss': 0.01, 'plastic_conc': 7.5,
        'solvent_price': 0.10, 'feedstock_price': 0.10, 'NdFeB_price': 150.0, 'HDPE_price': 1.50,
        'freshwater_price': 0.001, 'paa_price': 4.40, 'wastewater_fee': 0.001, 'othercomponent_fee': 0.01,
        'naoh_price': 0.375, 'nahso3_price': 0.325
    }
    if ENABLE_MBBR_WWT:
        best_params.update({'urea_price': 0.25, 'h3po4_price': 0.60})

    net_earnings = profit(**best_params)
    irr_val = my_irr(**best_params)
    
    print(f"=== STRAP BEST-CASE SCENARIO (MBBR={ENABLE_MBBR_WWT}) ===")
    print("Parameters Used:")
    for param, val in best_params.items():
        print(f"  * {param}: {val}")
    print("-" * 50)
    print(f"Net Earnings: ${net_earnings / 1e6:.2f} MM/yr | IRR: {irr_val:.2f}%\n")
    return net_earnings, irr_val


def worst_case_scenario():
    worst_params = {
        'u3_ndfeb_ratio': 0.8, 'capacity': 990, 'solvent_loss': 1.0, 'plastic_conc': 2.0,
        'solvent_price': 1.60, 'feedstock_price': 0.50, 'NdFeB_price': 50.0, 'HDPE_price': 1.00,
        'freshwater_price': 0.01, 'paa_price': 13.20, 'wastewater_fee': 0.0273, 'othercomponent_fee': 0.10,
        'naoh_price': 1.125, 'nahso3_price': 0.975
    }
    if ENABLE_MBBR_WWT:
        worst_params.update({'urea_price': 0.75, 'h3po4_price': 1.80})

    net_earnings = profit(**worst_params)
    irr_val = my_irr(**worst_params)
    
    print(f"=== STRAP WORST-CASE SCENARIO (MBBR={ENABLE_MBBR_WWT}) ===")
    print("Parameters Used:")
    for param, val in worst_params.items():
        print(f"  * {param}: {val}")
    print("-" * 50)
    print(f"Net Earnings: ${net_earnings / 1e6:.2f} MM/yr | IRR: {irr_val:.2f}%\n")
    return net_earnings, irr_val

# =============================================================================
# 5. CAPACITY & SENSITIVITY PLOTTING
# =============================================================================
def get_capacity_curve(min_scale_mt=50, max_scale_mt=2000):
    process.set_processing_capacity(1976)
    process.system.simulate()
    
    baseline_sales = process.tea.sales
    baseline_voc = process.tea.VOC
    baseline_foc = process.tea.FOC
    tax_rate = getattr(process.tea, 'income_tax', 0.21)

    cashflow_df = process.tea.get_cashflow_table()
    dep_in_usd = cashflow_df['Depreciation [MM$]'] * 1e6
    baseline_depreciation = dep_in_usd[dep_in_usd > 0].mean()

    def simulate_profit(scale_mt):
        if scale_mt < 1976:
            scale_ratio = scale_mt / 1976.0
            sales = baseline_sales * scale_ratio
            voc = baseline_voc * scale_ratio
            taxable_income = sales - (voc + baseline_foc + baseline_depreciation)
            return taxable_income * (1.0 - tax_rate)
        else:
            process.set_processing_capacity(scale_mt)
            process.system.simulate()
            return process.tea.net_earnings

    try:
        sol = opt.root_scalar(simulate_profit, bracket=[min_scale_mt, 1976], method='bisect', xtol=0.1)
        exact_be = sol.root
    except ValueError:
        exact_be = min_scale_mt

    low_scales = np.linspace(min_scale_mt, 250, 6)
    high_scales = np.arange(500, max_scale_mt + 1, 250)
    all_scales = np.unique(np.sort(np.concatenate(([exact_be], low_scales, high_scales))))

    scale_mt_list = [round(s, 1) for s in all_scales]
    profit_mm_list = [round(simulate_profit(s) / 1e6, 2) for s in all_scales]

    return scale_mt_list, profit_mm_list


def plot_it(data, base_val, title, x_label, color):
    data.sort(key=lambda x: x[3], reverse=True)
    names = [d[0] for d in data]
    lows, highs = [d[1] for d in data], [d[2] for d in data]
    p_lows, p_highs = [d[4] for d in data], [d[5] for d in data]
    
    plt.figure(figsize=(10, 5))
    widths = [h - l for h, l in zip(highs, lows)]
    
    bars = plt.barh(names, widths, left=lows, color=color, edgecolor='black', alpha=0.8)
    plt.barh(names, [base_val - l for l in lows], left=lows, color='red', edgecolor='black', alpha=0.8)
    plt.barh(names, [h - base_val for h in highs], left=[base_val] * len(names), color='green', edgecolor='black', alpha=0.8)
    plt.axvline(base_val, color='blue', linestyle='--', linewidth=2.0, zorder=3)
    
    x_min, x_max = plt.xlim()
    x_gap = 0.03 * (x_max - x_min)
    for bar, lo, hi in zip(bars, p_lows, p_highs):
        y = bar.get_y() + bar.get_height() / 2
        plt.text(bar.get_x() - x_gap, y, f'{lo}', ha='right', va='center', fontsize=12, fontweight='bold', color='grey')
        plt.text(bar.get_x() + bar.get_width() + x_gap, y, f'{hi}', ha='left', va='center', fontsize=12, fontweight='bold', color='grey')

    all_points = lows + highs + [base_val]
    padding = (max(all_points) - min(all_points)) * 0.25
    plt.xlim(min(all_points) - padding, max(all_points) + padding)
    ax = plt.gca()
    ax.set_title(f'Cytiva STRAP {title}', fontsize=16, fontweight='bold', pad=15)
    plt.xlabel(x_label, fontsize=14, fontweight='bold')
    plt.xticks(fontsize=12, fontweight='bold')
    plt.yticks(fontsize=12, fontweight='bold')
    plt.grid(False)
    
    for spine in ax.spines.values(): spine.set_linewidth(2)
    red_patch = mpatches.Patch(color='red', label='Losing money')
    green_patch = mpatches.Patch(color='green', label='Gaining money')
    baseline_legend = mlines.Line2D([], [], color='black', linestyle='--', label=f'Baseline: {base_val:.2f}')
    
    leg = plt.legend(handles=[red_patch, green_patch, baseline_legend], loc='best')
    leg.get_frame().set_edgecolor('black')
    leg.get_frame().set_linewidth(2)
    plt.tight_layout()
    plt.show()


def plot_profit_vs_scale(scale_list=None, profit_list=None, impeller_frac=None):
    if scale_list is None or profit_list is None:
        build_process(mbbr_flag = ENABLE_MBBR_WWT, proc_capacity=capacity)
        scale_list, profit_list = get_capacity_curve()

    plot_scales = [round(s * impeller_frac, 1) for s in scale_list] if impeller_frac else scale_list
    x_title = f'Impellers [Metric Tonnes / yr]' if impeller_frac else 'Bioreactor Bags Feedstock [Metric Tonnes / yr]'

    plt.figure(figsize=(12, 6), dpi=300)
    ax = plt.gca()
    x_labels = [f"{int(round(s)):,}" for s in plot_scales]
    
    colors = [plt.cm.Greens(0.4 + 0.5 * (p / max(profit_list))) if p >= 0 else plt.cm.Reds(0.5 + 0.4 * (abs(p) / abs(min(profit_list)))) for p in profit_list]
    bars = plt.bar(x_labels, profit_list, color=colors, edgecolor='none', width=0.6)

    plt.axhline(0, color='black', linewidth=1, linestyle='--')
    plt.xlabel(x_title, fontsize=11, fontweight='bold', labelpad=10)
    plt.ylabel('Net Earnings [MM$ / yr]', fontsize=11, fontweight='bold', labelpad=10)
    plt.title('Cytiva STRAP Plant Profitability', fontsize=13, fontweight='bold', pad=15)
    plt.xticks(rotation=45, ha='right', fontsize=9, fontweight='bold')
    plt.yticks(fontsize=10, fontweight='bold')

    max_val = max(abs(np.array(profit_list))) if len(profit_list) > 0 else 1.0
    for bar, p in zip(bars, profit_list):
        yval = bar.get_height()
        offset = 0.02 * max_val if yval >= 0 else -0.04 * max_val
        plt.text(bar.get_x() + bar.get_width() / 2.0, yval + offset, f"${p:.2f}M", ha='center', va='bottom' if yval >= 0 else 'top', fontsize=8, fontweight='bold', color='#005a00' if p >= 0 else '#8b0000')

    for spine in ax.spines.values(): spine.set_linewidth(1.5)
    plt.ylim(bottom=min(-5, min(profit_list) * 1.15), top=max(25, max(profit_list) * 1.15))
    plt.tight_layout()
    plt.show()


def plot_profit_line():
    build_process(mbbr_flag=ENABLE_MBBR_WWT, proc_capacity=capacity)
    scale_list, profit_list = get_capacity_curve()

    plt.figure(figsize=(12, 6.5), dpi=300)
    ax = plt.gca()
    plt.plot(scale_list, profit_list, marker='o', color='#005b96', linewidth=2.5, markersize=6, label='Net Earnings')
    plt.axhline(0, color='black', linewidth=1, linestyle='--')
    ax.grid(False)

    y_range = max(profit_list) - min(profit_list) or 1.0
    dense_offsets = [0.09, -0.22, 0.18, -0.11, 0.09, -0.22, 0.18]
    low_scale_idx = 0

    for x, y in zip(scale_list, profit_list):
        text_color = '#005a00' if y >= 0 else '#8b0000'
        if x <= 250:
            offset = dense_offsets[low_scale_idx % len(dense_offsets)] * y_range
            low_scale_idx += 1
            plt.annotate(f"${y:.2f}M", xy=(x, y), xytext=(x, y + offset), ha='center', va='center', fontsize=8, fontweight='bold', color=text_color, arrowprops=dict(arrowstyle='-', color='gray', lw=0.8, alpha=0.7))
        else:
            plt.text(x, y + 0.04 * y_range, f"${y:.2f}M", ha='center', va='bottom', fontsize=8.5, fontweight='bold', color=text_color)

    x_ticks = [0, 250, 500, 750, 1000, 1250, 1500, 1750, 2000]
    plt.xticks(x_ticks, [f"{x:,}" for x in x_ticks], fontsize=10, fontweight='bold')
    plt.yticks(fontsize=10, fontweight='bold')
    plt.xlabel('Bioreactor Bags Feedstock [Metric Tonnes / yr]', fontsize=11, fontweight='bold', labelpad=10)
    plt.ylabel('Net Earnings [MM$ / yr]', fontsize=11, fontweight='bold', labelpad=10)
    plt.title('Cytiva STRAP Plant Profitability vs. Scale', fontsize=13, fontweight='bold', pad=15)

    for spine in ax.spines.values(): spine.set_linewidth(1.5)
    plt.ylim(bottom=min(profit_list) - 0.32 * y_range, top=max(profit_list) + 0.14 * y_range)
    plt.tight_layout()
    plt.show()


# =============================================================================
# 6. REPORTING & ITEMIZED COST SUMMARY
# =============================================================================
def itemized_cost():
    data = []
    for u in process.system.units:
        util_cost = getattr(u, "utility_cost", None)
        if util_cost is None:
            try:
                power_cost = u.power_utility.cost if hasattr(u, "power_utility") else 0
                heat_cost = sum([hu.cost for hu in u.heat_utilities]) if hasattr(u, "heat_utilities") else 0
                util_cost = (power_cost or 0) + (heat_cost or 0)
            except Exception:
                util_cost = 0

        utility_cost_yr = util_cost * process.tea.operating_hours
        p_cost = getattr(u, "purchase_cost", 0) or 0
        i_cost = getattr(u, "installed_cost", 0) or 0

        data.append({
            "Unit": u.ID,
            "Unit operation": getattr(u, "line", u.__class__.__name__),
            "Purchase cost (10^3 USD)": p_cost / 1e3,
            "Utility cost (10^3 USD/yr)": utility_cost_yr / 1e3,
            "Installed cost (10^3 USD)": i_cost / 1e3,
        })
    
    df_equipment = pd.DataFrame(data).sort_values(by="Unit").reset_index(drop=True)
    df_equipment.to_excel("Equipment_Summary_Table.xlsx", index=False)
    return df_equipment


# =============================================================================
# 7. INITIALIZE BASELINE MODEL
# =============================================================================
build_process(mbbr_flag=ENABLE_MBBR_WWT, proc_capacity=capacity)