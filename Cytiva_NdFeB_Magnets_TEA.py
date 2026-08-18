# -*- coding: utf-8 -*-
"""
Created on Tue Apr 21 15:27:00 2026

@author: tlnu
"""

import biosteam as bst
import importlib

bst.settings.ID_magic = False
from plastics import strap
importlib.reload(strap.property_package)
importlib.reload(strap.process_model)
from chaospy.distributions import Uniform
import pandas as pd
from warnings import filterwarnings
import numpy as np
import matplotlib.pyplot as plt
import os
from datetime import datetime
from matplotlib.ticker import PercentFormatter, FuncFormatter
from matplotlib.colors import TwoSlopeNorm, Normalize
import matplotlib.patches as mpatches
import matplotlib.lines as mlines


filterwarnings('ignore')

bst.main_flowsheet.clear()

if hasattr(strap.MagnetRecovery, 'cache'):
    strap.MagnetRecovery.cache.clear()

capacity = 1976

bst.settings.CEPCI =836.9
process = strap.MagnetRecovery(
    processing_capacity = capacity,      #tons
    sell_leftover_plastic = True,
    simulate=False)

#for STRAP plant
N_operators = 1
N_plant_manager = 1

plant_manager_salary = 110000
operator_salary = 80000


feed_composition = {'NdFeB': 0.1548,
        'HDPE': 0.0252,
        'Films': 0.10,
        'FittingsFilters': 0.36,
        'BrownSupport': 0.07,
        'SiliconeTubings': 0.29,
        'Solutes': 0.001
        }


#update the processing conditions if model starts behaving funky
process.tea.operating_days = 309
operating_hours = 309*16 #16hours/day for 2 shifts

#adjusting the total number of days for 2 shifts/day
#adjusted_days = operating_hours/24
process.tea.operating_hours = operating_hours
feed_per_hour = capacity*1000/process.tea.operating_hours

process.feedstock.empty()
for chem, frac in feed_composition.items():
    process.feedstock.imass[chem]=feed_per_hour*frac


strap_labor = N_plant_manager*plant_manager_salary + N_operators*operator_salary*1.6*3

labor_cost = strap_labor + process.U11.total_salary + process.U10.total_salary

process.plastic.ID =  'Bioreactor bags'
process.Vac_S.outs[0].ID = 'NdFeB_Magnets'
process.Vac_S.outs[1].ID = 'Xylenes vapors'

process.M2.disconnect()

process.solvent.ID = 'Xylenes'

active_units = [unit for unit in process.system.units if unit.ID != 'M2']
process.system  = bst.System(ID = 'sys', path =active_units, facilities = process.system.facilities)

#process.S2.outs[1].ID='NdFeB_Magnets'
process.s22.ID='Impurities'

process.NdFeB_Magnets = process.Vac_S.outs[0]
process.HDPE_resins = process.U9.outs[0]

products = [process.NdFeB_Magnets, process.HDPE_resins]
process.products[:] = [process.HDPE_resins, process.NdFeB_Magnets]

#products
process.HDPE_resins.price = 1.20
process.NdFeB_Magnets.price = 100

#wastes generated
wastewater = process.U10.outs[2]
othercomponents = process.S3.outs[0]
 
wastewater.price = -0.003
othercomponents.price = -0.072
process.S3.outs[1].price = -0.003

#raw_materials
bioreactor = process.feedstock
freshwater = process.U10.ins[1]
peracetic_acid = process.U10.ins[2]
process.set_solvent_price(0.85)


#setting up raw material price
bioreactor.price = 0.25
freshwater.price = 0.0015
peracetic_acid.price = 8.8







#other tea parameters
process.tea.IRR = 0.15





#process.set_processing_capacity(325)
process.system.simulate()





#adding a vacuum storage unit

process.system.diagram()

#TEA parameters
#process.set_IRR(0.15)



#%%
import types
import plastics.strap.units as strap_units

# 1. Grab the original design method from the class definition
original_u7_design = strap_units.Precipitator._design

# 2. Define the safe function accepting 'self'
def safe_u7_design(self):
    # Check if this specific unit's inlet is empty or has zero volumetric flow
    if self.ins[0].isempty() or self.ins[0].F_vol <1e-6:
        self.design_results['Vessel volume'] = 0.0
        self.design_results['Number of reactors'] = 0
    else:
        try:
            # Pass 'self' back to the original class method normally
            original_u7_design(self)
        except ZeroDivisionError:
            self.design_results['Vessel volume'] = 0.0
            self.design_results['Number of reactors'] = 0.0

# 3. Bind the function to the specific process.U7 instance as a method
process.U7._design = types.MethodType(safe_u7_design, process.U7)
process.system.simulate()

elec_price = bst.settings.electricity_price

# %% Updating power
# ==============================================================================
# 4. INDIVIDUAL UNIT POWER ASSIGNMENTS (No Correlations)
# ==============================================================================

# Dictionary of exact, standalone power ratings (kW) per unit operation
INDIVIDUAL_UNIT_POWER_KW = {
    'U10': 1.06,   # Pretreatment Granulator / Shredder
    'U3':  1.0,   # Dissolution Tank 1 Agitator
    'T2':  0.0,   # Dissolution Holding Tank 2 Agitator
    'T3':  0.0,   # Precipitation Tank 3 Agitator
    'T4':  0.0,   # Holding Tank 4 Agitator
    'P1':  0.01,   # Feed Pump 1
    'P2':  0.3,   # Pump 2
    'P3':  0.10,   # Solvent Pump 3
    'U6':  1.0,   # Centrifuge 1
    'U7':  0.0,   # Precipitator tank (in PE-film model file)
    'U8':  4.15,   # Dewatering Centrifuge / Separation
    'F1':  5.5,     # vacuum consistent with PE-film model file
    'U9':  18.5,  # Screw Degasser / Primary Driver
    'U2': 0.10,     #adsorption column
    'Vac_S': 0.75, # Vacuum Pump Package
    'H1':  0.15,   # Heat Exchanger 1
    'H2':  0.10,   # Chiller Unit
    'H3':  0.00,   # Heat Exchanger 3
    'CT': 0.1       #cooling tower
}

def apply_strap_power_corrections():
    """
    Sets each unit operation's electricity consumption (kW) separately
    using fixed, explicit ratings without any scaling correlations.
    """
    for unit_id, target_kw in INDIVIDUAL_UNIT_POWER_KW.items():
        if hasattr(process, unit_id):
            unit = getattr(process, unit_id)
            unit.power_utility.consumption = target_kw


def print_strap_power_summary():
    """Prints an itemized list of electrical power draw (kW) for every unit."""
    print("\n--- Itemized STRAP Unit Electrical Consumption (kW) ---")
    total_kw = 0.0
    for u in process.system.units:
        pwr = u.power_utility.consumption
        total_kw += pwr
        print(f"{u.ID:>10}: {pwr:.2f} kW")
    print("-" * 45)
    print(f"Total Plant Electrical Draw: {total_kw:.2f} kW\n")



apply_strap_power_corrections()
#%% updating units with price based on smallest size price as in MTU
# ------------------------------------------------------------------------------
# DISSOLUTION & TANKS
# ------------------------------------------------------------------------------

def custom_u3_cost():
    process.U3.baseline_purchase_costs.clear()
    process.U3.baseline_purchase_costs['Vertical pressure vessel'] = 32352.9
    process.U3.baseline_purchase_costs['Platform and ladders'] = 0
    process.U3.baseline_purchase_costs['Agitator - Agitator'] = 0
    
    # Target Utility Cost: $0.0203/hr (MTU dissolution vessel)
    process.U3.power_utility.rate = 0.0203 / elec_price
    process.U3.BM = 3.4
process.U3._cost = custom_u3_cost

def custom_t2_cost():
    process.T2.baseline_purchase_costs.clear()
    if hasattr(process.T2, 'baseline_item_bms'):
        process.T2.baseline_item_bms.clear()
    process.T2.baseline_purchase_costs['Tank'] = 15000.0
    process.T2.power_utility.rate = 0.0
process.T2._cost = custom_t2_cost

def custom_t3_cost():
    process.T3.baseline_purchase_costs.clear()
    process.T3.baseline_purchase_costs['Tank system'] = 15000.0
    process.T3.F_BM['Tank system'] = 2.3
    process.T3.power_utility.rate = 0.0
process.T3._cost = custom_t3_cost

def custom_t4_cost():
    process.T4.baseline_purchase_costs.clear()
    if hasattr(process.T4, 'baseline_item_bms'):
        process.T4.baseline_item_bms.clear()
    process.T4.baseline_purchase_costs['Tank'] = 20000.0
    process.T4.power_utility.rate = 0.0
process.T4._cost = custom_t4_cost


# ------------------------------------------------------------------------------
# PUMPS
# ------------------------------------------------------------------------------

def custom_p1_cost():
    process.P1.baseline_purchase_costs.clear()
    process.P1.baseline_purchase_costs['Pump'] = 8000.0
    # Target Utility Cost: $0.0000029/hr (intermittent feed pump)
    process.P1.power_utility.rate = 2.9e-6 / elec_price
process.P1._cost = custom_p1_cost

def custom_p2_cost():
    process.P2.baseline_purchase_costs.clear()
    process.P2.baseline_purchase_costs['Pump system'] = 8000.0
    process.P2.F_BM['Pump system'] = 3.3 
    # Target Utility Cost: $0.0034/hr (MTU pump_2 load)
    process.P2.power_utility.rate = 0.0034 / elec_price
process.P2._cost = custom_p2_cost

def custom_p3_cost():
    process.P3.baseline_purchase_costs.clear()
    process.P3.baseline_purchase_costs['Pump'] = 0.0
    # Target Utility Cost: $0.0270/hr (MTU solvent supply pump)
    process.P3.power_utility.rate = 0.0270 / elec_price
process.P3._cost = custom_p3_cost


# ------------------------------------------------------------------------------
# CENTRIFUGES & SEPARATION UNITS
# ------------------------------------------------------------------------------

def custom_u6_cost():
    process.U6.baseline_purchase_costs.clear()
    process.U6.baseline_purchase_costs['combined U6 and P3'] = 250000.0
    process.U6.F_BM['combined U6 and P3'] = 3.0
    # Target Utility Cost: $0.0021/hr (MTU centrifuge_1 load)
    process.U6.power_utility.rate = 0.0021 / elec_price
process.U6._cost = custom_u6_cost

def custom_u7_cost():
    process.U7.baseline_purchase_costs.clear()
    if hasattr(process.U7, 'baseline_item_bms'):
        process.U7.baseline_item_bms.clear()
    process.U7.baseline_purchase_costs['Precipitator'] = 175000.0
    process.U7.F_BM['Precipitator'] = 3.46
    # Target Utility Cost: $0.0133/hr (MTU precipitation tank)
    process.U7.heat_utilities.clear()
    process.U7.power_utility.rate = 0.0133 / elec_price
process.U7._cost = custom_u7_cost

def custom_u8_cost():
    process.U8.baseline_purchase_costs.clear()
    process.U8.baseline_purchase_costs['Precipitator'] = 0.0
    # Target Utility Cost: $0.2888/hr (MTU dewatering centrifuge)
    process.U8.power_utility.rate = 0.2888 / elec_price
process.U8._cost = custom_u8_cost

def custom_U9_cost():
    process.U9.baseline_purchase_costs.clear()
    if hasattr(process.U9, 'baseline_item_bms'):
        process.U9.baseline_item_bms.clear()
    process.U9.baseline_purchase_costs['Screw degasser'] = 250000.0
    process.U9.F_BM['Screw degasser'] = 2.60
    # Target Utility Cost: $1.3061/hr (MTU degasser load)
    process.U9.power_utility.rate = 1.3061 / elec_price
process.U9._cost = custom_U9_cost


# ------------------------------------------------------------------------------
# HEAT EXCHANGERS & CHILLERS
# ------------------------------------------------------------------------------

def custom_h1_cost():
    process.H1.baseline_purchase_costs.clear()
    process.H1.baseline_purchase_costs['Double pipe'] = 8000.0 * 1.776
    process.H1.heat_utilities.clear()                     # Clear default auto-utilities
    process.H1.power_utility.rate = 0.0146 / elec_price   # Target: $0.0146/hr
process.H1._cost = custom_h1_cost

def custom_h2_cost():
    process.H2.baseline_purchase_costs.clear()
    if hasattr(process.H2, 'baseline_item_bms'):
        process.H2.baseline_item_bms.clear()
    process.H2.baseline_purchase_costs['Chiller'] = 50000.0
    process.H2.F_BM['Chiller'] = 2.5
    process.H2.heat_utilities.clear()                     # Clear default auto-utilities
    process.H2.power_utility.rate = 4e-6 / elec_price     # Target: $0.000004/hr
process.H2._cost = custom_h2_cost

def custom_h3_cost():
    process.H3.baseline_purchase_costs.clear()
    process.H3.baseline_purchase_costs['Double pipe'] = 125000.0 / 1.768
    process.H3.power_utility.rate = 0.0
process.H3._cost = custom_h3_cost



# 2. Define zero-cost methods for CWP and CT
def custom_CWP_cost():
  process.CWP.baseline_purchase_costs.clear()
  process.CWP.power_utility.rate = 0.0


def custom_CT_cost():
  process.CT.baseline_purchase_costs.clear()
  process.CT.power_utility.rate = 0.0


# 3. Attach the cost functions
process.CWP._cost = custom_CWP_cost
process.CT._cost = custom_CT_cost


# Run simulation
process.system.simulate()
# Record total vendor baseline equipment purchase cost at 1,976 MT as an absolute floor
baseline_total_purchase_cost = sum([u.purchase_cost for u in process.system.units])




#%% functions for TEA
def profit(u3_ndfeb_ratio, capacity,solvent_loss,plastic_conc,solvent_price,feedstock_price, NdFeB_price, 
           HDPE_price, freshwater_price, paa_price,wastewater_fee,othercomponent_fee):
    
    
    process.set_processing_capacity(capacity)
    
    process.set_solvent_loss(solvent_loss)
    process.set_dissolution_capacity(plastic_conc)
    
    process.set_solvent_price(solvent_price)
    process.set_feedstock_price(feedstock_price)
    process.NdFeB_Magnets.price = NdFeB_price
    process.HDPE_resins.price = HDPE_price
    
    freshwater.price = freshwater_price
    peracetic_acid.price = paa_price
    
    wastewater.price = -wastewater_fee
    process.S3.outs[1].price = -wastewater_fee
    
    othercomponents.price = -othercomponent_fee
    

    process.tea.operating_hours = capacity*1000/(4*100)
    
    
    feed_per_hour = capacity * 1000 / process.tea.operating_hours
    
    
    process.feedstock.imass['NdFeB'] = feed_per_hour * (0.18 * u3_ndfeb_ratio)
    process.feedstock.imass['HDPE'] = feed_per_hour * (0.18 * (1.0 - u3_ndfeb_ratio))
    
    # Remaining baseline feedstock matrix components
    process.feedstock.imass['Films'] = feed_per_hour * 0.10
    process.feedstock.imass['FittingsFilters'] = feed_per_hour * 0.36
    process.feedstock.imass['BrownSupport'] = feed_per_hour * 0.07
    process.feedstock.imass['SiliconeTubings'] = feed_per_hour * 0.29
    
    process.system.simulate()
    return process.tea.net_earnings 



def my_irr(u3_ndfeb_ratio, capacity, solvent_loss, plastic_conc,solvent_price, feedstock_price, NdFeB_price, 
           HDPE_price, freshwater_price, paa_price,wastewater_fee,othercomponent_fee):
    
    process.set_processing_capacity(capacity)
    
    process.set_solvent_loss(solvent_loss)
    process.set_dissolution_capacity(plastic_conc)
    process.set_solvent_price(solvent_price)
    process.set_feedstock_price(feedstock_price)
    process.NdFeB_Magnets.price = NdFeB_price
    process.HDPE_resins.price = HDPE_price
    
    freshwater.price = freshwater_price
    peracetic_acid.price = paa_price
    
    wastewater.price = -wastewater_fee
    othercomponents.price = -othercomponent_fee
    
    feed_per_hour = capacity * 1000 / process.tea.operating_hours
    
    # Update composition map directly using your baseline dictionary values
    process.feedstock.imass['NdFeB'] = feed_per_hour * (0.18 * u3_ndfeb_ratio)
    process.feedstock.imass['HDPE'] = feed_per_hour * (0.18 * (1.0 - u3_ndfeb_ratio))
    
    # Remaining baseline feedstock matrix components
    process.feedstock.imass['Films'] = feed_per_hour * 0.10
    process.feedstock.imass['FittingsFilters'] = feed_per_hour * 0.36
    process.feedstock.imass['BrownSupport'] = feed_per_hour * 0.07
    process.feedstock.imass['SiliconeTubings'] = feed_per_hour * 0.29
    
    process.system.simulate()
    return process.tea.solve_IRR()*100

# %% Best case vs worst case
def best_case_scenario():
    """
    Evaluates the STRAP process under the most favorable (optimistic) operating conditions
    derived from the upper/lower bounds in the sensitivity analysis tornado plot.
    """
    best_params = {
        'u3_ndfeb_ratio': 0.9,          # Max magnet concentration in impeller (0.9)
        'capacity': 2965,               # Max processing capacity (2,965 tons/yr)
        'solvent_loss': 0.01,           # Min solvent loss rate (1%)
        'plastic_conc': 7.5,            # Max feed-to-solvent concentration (7.5 wt%)
        'solvent_price': 0.10,          # Min solvent price ($0.10/kg)
        'feedstock_price': 0.10,        # Min feedstock purchase price ($0.10/kg)
        'NdFeB_price': 150.0,           # Max magnet selling price ($150.00/kg)
        'HDPE_price': 1.50,             # Max HDPE selling price ($1.50/kg)
        'freshwater_price': 0.001,      # Min freshwater price ($0.001/kg)
        'paa_price': 4.40,              # Min peracetic acid price ($4.40/kg)
        'wastewater_fee': 0.001,        # Min wastewater fee ($0.001/L)
        'othercomponent_fee': 0.01      # Min disposal fee ($0.01/kg)
    }
    
    net_earnings = profit(**best_params)
    irr_val = my_irr(**best_params)
    
    print(f"=== STRAP BEST-CASE SCENARIO ===")
    print(f"Net Earnings: ${net_earnings / 1e6:.2f} MM/yr")
    print(f"Internal Rate of Return (IRR): {irr_val:.2f}%\n")
    
    return net_earnings, irr_val


def worst_case_scenario():
    """
    Evaluates the STRAP process under the most unfavorable (pessimistic) operating conditions
    derived from the upper/lower bounds in the sensitivity analysis tornado plot.
    """
    worst_params = {
        'u3_ndfeb_ratio': 0.8,          # Min magnet concentration in impeller (0.8)
        'capacity': 990,                # Min processing capacity (990 tons/yr)
        'solvent_loss': 1.0,            # Max solvent loss rate (100%)
        'plastic_conc': 2.0,            # Min feed-to-solvent concentration (2.0 wt%)
        'solvent_price': 1.60,          # Max solvent price ($1.60/kg)
        'feedstock_price': 0.50,        # Max feedstock purchase price ($0.50/kg)
        'NdFeB_price': 50.0,            # Min magnet selling price ($50.00/kg)
        'HDPE_price': 1.00,             # Min HDPE selling price ($1.00/kg)
        'freshwater_price': 0.01,       # Max freshwater price ($0.010/kg)
        'paa_price': 13.20,             # Max peracetic acid price ($13.20/kg)
        'wastewater_fee': 0.005,        # Max wastewater fee ($0.005/L)
        'othercomponent_fee': 0.10      # Max disposal fee ($0.10/kg)
    }
    
    net_earnings = profit(**worst_params)
    irr_val = my_irr(**worst_params)
    
    print(f"=== STRAP WORST-CASE SCENARIO ===")
    print(f"Net Earnings: ${net_earnings / 1e6:.2f} MM/yr")
    print(f"Internal Rate of Return (IRR): {irr_val:.2f}%\n")
    
    return net_earnings, irr_val
    




#%% simplified TEA definition

class SimplifiedTEA(bst.TEA):
    """
    Custom TEA class that overrides default BioSTEAM CAPEX calculations
    to use simplified percentages and matches MechanicalRecyclingTEA FOC logic.
    """
    def __init__(self, system, IRR, duration, depreciation, operating_days, income_tax, 
                 lang_factor=None, construction_schedule=(1.0,), labor_cost=0.0, 
                 fringe_benefits=0.4, property_tax=0.001, WC_over_FCI=0.1, 
                 property_insurance=0.005, supplies=0.05, maintenance=0.03, 
                 administration=0.005,  **kwargs):
        
        # --- Unwrap system if process container object (e.g., MagnetRecovery) is passed ---
        if not hasattr(system, 'installed_equipment_cost'):
            if hasattr(system, 'sys'):
                system = system.sys
            elif hasattr(system, 'system'):
                system = system.system
        
        # Call base bst.TEA __init__
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
            f'{self.WC_over_FCI * 100:.1f}% of FCI', 
            'FCI + WC'
        ]

        costs = [isbl, osbl, tdc, eng, cont, tic, fci, wc, tci]

        df = pd.DataFrame(
            {'Notes': notes, 'Cost [MM$]': [round(c, 3) for c in costs]},
            index=pd.MultiIndex.from_tuples(index)
        )
        return df

adjusted_days = process.tea.operating_hours / 24

process.tea = SimplifiedTEA(
    system=process.system,
    IRR=0.15,
    duration=(2025, 2045),
    depreciation='MACRS7',
    income_tax=0.21,
    operating_days=adjusted_days,
    labor_cost=labor_cost,
    WC_over_FCI=0.10
)



#%% Scale vs Profit
import scipy.optimize as opt

def get_capacity_curve(min_scale_mt=50, max_scale_mt=2000):
    """
    Calculates capacity vs. net earnings using dense points for low scales (< 250 MT)
    and clean 250 MT steps for higher scales (500 to 2000 MT).
    """
    # 1. Capture exact 1,976 MT baseline values
    process.set_processing_capacity(1976)
    process.system.simulate()
    
    baseline_sales = process.tea.sales
    baseline_voc = process.tea.VOC
    baseline_foc = process.tea.FOC
    tax_rate = getattr(process.tea, 'income_tax', 0.21)

    # 2. Extract baseline depreciation from cash flow table
    cashflow_df = process.tea.get_cashflow_table()
    dep_in_usd = cashflow_df['Depreciation [MM$]'] * 1e6
    baseline_depreciation = dep_in_usd[dep_in_usd > 0].mean()

    # 3. Income statement simulator
    def simulate_profit(scale_mt):
        if scale_mt < 1976:
            scale_ratio = scale_mt / 1976.0
            sales = baseline_sales * scale_ratio
            voc = baseline_voc * scale_ratio
            foc = baseline_foc
            dep = baseline_depreciation
            
            taxable_income = sales - (voc + foc + dep)
            earnings = taxable_income * (1.0 - tax_rate)
        else:
            process.set_processing_capacity(scale_mt)
            process.system.simulate()
            earnings = process.tea.net_earnings

        return earnings

    # 4. Find exact break-even point
    try:
        sol = opt.root_scalar(simulate_profit, bracket=[min_scale_mt, 1976], method='bisect', xtol=0.1)
        exact_be = sol.root
    except ValueError:
        exact_be = min_scale_mt

    # 5. Generate clean, evenly-spaced scale points
    low_scales = np.linspace(min_scale_mt, 250, 6)
    # Clean 250 MT step increments after 250 MT/yr
    high_scales = np.arange(500, max_scale_mt + 1, 250)  # [500, 750, 1000, 1250, 1500, 1750, 2000]
    
    all_scales = np.unique(np.sort(np.concatenate(([exact_be], low_scales, high_scales))))

    scale_mt_list = []
    profit_mm_list = []

    for scale_mt in all_scales:
        earnings = simulate_profit(scale_mt)
        scale_mt_list.append(round(scale_mt, 1))
        profit_mm_list.append(round(earnings / 1e6, 2))

    return scale_mt_list, profit_mm_list


#%% Tornado / Sensitivity Plot
def tornado_plot():
    """Generates sensitivity analysis tornado plots for Net Earnings and IRR."""
    base_params = {
        'u3_ndfeb_ratio': 0.86,
        'capacity': 1976,
        'solvent_loss': 0.1, 
        'plastic_conc': 5.0, 
        'solvent_price': 0.85, 
        'feedstock_price': 0.25,
        'NdFeB_price': 100,
        'HDPE_price': 1.20,
        'freshwater_price': 0.002,
        'wastewater_fee': 0.003,
        'paa_price': 8.8,
        'othercomponent_fee': 0.072
    }
    
    ranges = {
        'u3_ndfeb_ratio': (0.80, 0.90),
        'capacity': (990, 2965),
        'solvent_loss': (0.01, 1.0),
        'plastic_conc': (2.0, 7.5),
        'solvent_price': (0.1, 1.6),
        'feedstock_price': (0.1, 0.50),
        'NdFeB_price': (50.0, 150.0),
        'HDPE_price': (1.00, 1.50),
        'freshwater_price': (0.001, 0.01), 
        'wastewater_fee': (0.001, 0.005),
        'paa_price': (4.4, 13.2),
        'othercomponent_fee': (0.01, 0.1)
    }
    
    label_map = {
        'u3_ndfeb_ratio': 'NdFeB Concentration in Impeller',
        'capacity': 'Processing Capacity (tons)',
        'solvent_loss': 'Solvent Loss Rate',
        'plastic_conc': 'Feed-to-Solvent (wt%)',
        'solvent_price': 'Solvent Price ($/kg)',
        'feedstock_price': 'Feedstock Price ($/kg)',
        'NdFeB_price': 'Neodymium Magnet Price ($/kg)',
        'HDPE_price': 'HDPE Price ($/kg)',
        'freshwater_price': 'Freshwater Price ($/kg)',
        'wastewater_fee': 'Wastewater Fee ($/L)',
        'paa_price': 'Peracetic Acid Price ($/kg)',
        'othercomponent_fee': 'Othercomponent Fee ($/kg)'
    }
    
    results = []
    results_irr = []  
    
    base_profit_m = profit(**base_params) / 1e6
    base_irr_m = my_irr(**base_params)
  
    for param, (low, high) in ranges.items():
        p_low = base_params.copy(); p_low[param] = low
        p_high = base_params.copy(); p_high[param] = high

        low_profit_m = profit(**p_low) / 1e6
        high_profit_m = profit(**p_high) / 1e6
        
        low_irr = my_irr(**p_low)
        high_irr = my_irr(**p_high)
        
        results.append(
            (label_map[param], low_profit_m, high_profit_m, abs(high_profit_m - low_profit_m), low, high)
        )
        results_irr.append(
            (label_map[param], low_irr, high_irr, abs(high_irr - low_irr), low, high)
        )
        
    plot_it(results, base_profit_m, "Sensitivity: Net Earnings", "Profit ($ Millions)", "skyblue")
    plot_it(results_irr, base_irr_m, "Sensitivity: IRR", "IRR (%)", "teal")


#%% Plot Helper function for Tornado Plot
def plot_it(data, base_val, title, x_label, color):
    data.sort(key=lambda x: x[3], reverse=True)
    names = [d[0] for d in data]
    lows = [d[1] for d in data]
    highs = [d[2] for d in data]
    p_lows = [d[4] for d in data]
    p_highs = [d[5] for d in data]
    
    plt.figure(figsize=(10, 5))
    widths = [h - l for h, l in zip(highs, lows)]
    
    bars = plt.barh(
        names,
        widths,
        left=lows,
        color=color,
        edgecolor='black',
        alpha=0.8
    )

    # Low value to baseline = red
    plt.barh(
        names,
        [base_val - l for l in lows],
        left=lows,
        color='red',
        edgecolor='black',
        alpha=0.8
    )

    # Baseline to high value = green
    plt.barh(
        names,
        [h - base_val for h in highs],
        left=[base_val] * len(names),
        color='green',
        edgecolor='black',
        alpha=0.8
    )       
    
    plt.axvline(base_val, color='blue', linestyle='--', linewidth=2.0, zorder=3)
    
    x_min, x_max = plt.xlim()
    x_gap = 0.03 * (x_max - x_min)
    
    for bar, lo, hi in zip(bars, p_lows, p_highs):
        y = bar.get_y() + bar.get_height() / 2
    
        # Low value (left of bar)
        plt.text(
            bar.get_x() - x_gap,
            y,
            f'{lo}',
            ha='right',
            va='center',
            fontsize=12,
            fontweight='bold',
            color='grey',
            alpha=0.85
        )
    
        # High value (right of bar)
        plt.text(
            bar.get_x() + bar.get_width() + x_gap,
            y,
            f'{hi}',
            ha='left',
            va='center',
            fontsize=12,
            fontweight='bold',
            color='grey',
            alpha=0.85
        )

    all_points = lows + highs + [base_val]
    padding = (max(all_points) - min(all_points)) * 0.25
    plt.xlim(min(all_points) - padding, max(all_points) + padding)

    ax = plt.gca()
    ax.set_title('Cytiva STRAP Sensitivity Analysis', fontsize=16, fontweight='bold', pad=15)
    
    plt.xlabel(x_label, fontsize=14, fontweight='bold')
    plt.xticks(fontsize=12, fontweight='bold')
    plt.yticks(fontsize=12, fontweight='bold')
    plt.grid(False)
    
    for spine in ax.spines.values():
        spine.set_linewidth(2)

    red_patch = mpatches.Patch(color='red', label='Losing money')
    green_patch = mpatches.Patch(color='green', label='Gaining money')
    baseline_legend = mlines.Line2D([], [], color='black', linestyle='--', label=f'Baseline: {base_val:.2f}')
    
    leg = plt.legend(
        handles=[red_patch, green_patch, baseline_legend], 
        prop={'weight': 'bold', 'size': 12},
        loc='best'
    )
    leg.get_frame().set_edgecolor('black')
    leg.get_frame().set_linewidth(2)
    
    plt.tight_layout()
    plt.show()   


#%% Plot Profit vs Scale Bar Chart
def plot_profit_vs_scale(scale_list=None, profit_list=None, impeller_frac=None):
    """
    Plots a clean bar chart of Profit vs. Scale.
    If scale_list and profit_list are not supplied, it automatically runs get_capacity_curve().
    """
    if scale_list is None or profit_list is None:
        # 1. Reset baseline parameters that tornado_plot() alters
        # (Add/adjust stream/parameter names to match your script's variables)
        try:
            # Re-set key stream baseline prices
            if 'NdFeB' in globals(): process.NdFeB.price = 100.0  # $/kg
            if 'feedstock' in globals(): process.feedstock.price = 0.25 # $/kg
            if 'water' in globals(): process.water.price = 0.0015 # $/kg
            if 'paa' in globals(): process.paa.price = 8.80 # $/kg
            if 'HDPE' in globals(): process.HDPE.price = 1.20 # $/kg (or disposal fee for Mech)
            
            # Re-set baseline capacity (1,976 MT/yr)
            # 1,976,000 kg / 8000 operating hours = 247 kg/hr
            if 'feedstock' in globals(): process.feedstock.F_mass = 1976000.0 / 8000.0 

            # Re-simulate system to flush out mutated unit operations/TEA states
            if 'sys' in globals(): 
                process.sys.simulate()
            elif 'process' in globals() and hasattr(process, 'sys'):
                process.sys.simulate()

        except Exception as e:
            print(f"Warning: Could not automatically reset system parameters: {e}")

        # 2. Get clean capacity curve from baseline state
        scale_list, profit_list = get_capacity_curve()




    # Convert scale if plotting for impellers
    if impeller_frac is not None:
        plot_scales = [round(s * impeller_frac, 1) for s in scale_list]
        x_title = f'Impellers [Metric Tonnes / yr] ({impeller_frac*100:g}% of Bioreactor Bags)'
    else:
        plot_scales = scale_list
        x_title = 'Bioreactor Bags Feedstock [Metric Tonnes / yr]'

    plot_title = 'Cytiva STRAP Plant Profitability'

    plt.figure(figsize=(12, 6), dpi=300)
    ax = plt.gca()

    x_labels = [f"{int(round(s)):,}" for s in plot_scales]
    
    # Calculate relative color intensity based on profit height
    max_pos = max([p for p in profit_list if p > 0] or [1.0])
    max_neg = abs(min([p for p in profit_list if p < 0] or [-1.0]))
    
    colors = []
    for p in profit_list:
        if p >= 0:
            intensity = 0.4 + 0.5 * (p / max_pos)  # Scale between 0.4 and 0.9 in Greens map
            colors.append(plt.cm.Greens(intensity))
        else:
            intensity = 0.5 + 0.4 * (abs(p) / max_neg) # Scale between 0.5 and 0.9 in Reds map
            colors.append(plt.cm.Reds(intensity))

    # Draw bars directly with calculated color gradient shades and NO outlines
    bars = plt.bar(x_labels, profit_list, color=colors, edgecolor='none', width=0.6)

    # Zero line
    plt.axhline(0, color='black', linewidth=1, linestyle='--')
    
    plt.xlabel(x_title, fontsize=11, fontweight='bold', labelpad=10)
    plt.ylabel('Net Earnings [MM$ / yr]', fontsize=11, fontweight='bold', labelpad=10)
    plt.title(plot_title, fontsize=13, fontweight='bold', pad=15)
    
    plt.xticks(rotation=45, ha='right', fontsize=9, fontweight='bold')
    plt.yticks(fontsize=10, fontweight='bold')

    # Annotate bar values with dark green/red text
    max_val = max(abs(np.array(profit_list))) if len(profit_list) > 0 else 1.0
    for bar, profit in zip(bars, profit_list):
        yval = bar.get_height()
        va_align = 'bottom' if yval >= 0 else 'top'
        offset = 0.02 * max_val if yval >= 0 else -0.04 * max_val
        
        label = f"${profit:.2f}M"
        text_color = '#005a00' if profit >= 0 else '#8b0000'
        
        plt.text(
            bar.get_x() + bar.get_width() / 2.0, 
            yval + offset, 
            label, 
            ha='center', 
            va=va_align, 
            fontsize=8,
            fontweight='bold',
            color=text_color
        )

    # Style plot borders/spines
    for spine in ax.spines.values():
        spine.set_linewidth(1.5)
        spine.set_color('black')
        
    plt.ylim(bottom=-5, top=25)
    plt.tight_layout()
    plt.show()

#%%plot profit vs scale line
def plot_profit_line():
    """
    Plots Net Earnings vs. Scale for Cytiva STRAP Plant with multi-tiered callouts
    for dense low-scale points (<= 250 MT) to prevent label overlap.
    """
    # 1. Reset baseline parameters
    try:
        if 'NdFeB' in globals(): process.NdFeB.price = 100.0
        if 'feedstock' in globals(): process.feedstock.price = 0.25
        if 'water' in globals(): process.water.price = 0.0015
        if 'paa' in globals(): process.paa.price = 8.80
        if 'HDPE' in globals(): process.HDPE.price = 1.20
        if 'feedstock' in globals(): process.feedstock.F_mass = 1976000.0 / 8000.0 

        if 'sys' in globals(): 
            process.sys.simulate()
        elif 'process' in globals() and hasattr(process, 'sys'):
            process.sys.simulate()
    except Exception as e:
        print(f"Warning: Could not automatically reset system parameters: {e}")

    # 2. Get capacity curve data
    scale_list, profit_list = get_capacity_curve()

    # 3. Create Line Plot
    plt.figure(figsize=(12, 6.5), dpi=300)
    ax = plt.gca()

    # Plot data line
    plt.plot(
        scale_list, 
        profit_list, 
        marker='o', 
        color='#005b96', 
        linewidth=2.5, 
        markersize=6, 
        label='Net Earnings'
    )

    # Zero break-even reference line
    plt.axhline(0, color='black', linewidth=1, linestyle='--')

    # Remove background gridlines
    ax.grid(False)

    y_range = max(profit_list) - min(profit_list)
    if y_range == 0: y_range = 1.0
    
    # MULTI-TIERED OFFSETS: Alternates direction AND height (low-up, deep-down, high-up, shallow-down)
    # to guarantee horizontal and vertical separation between adjacent labels
    dense_offsets = [0.09, -0.22, 0.18, -0.11, 0.09, -0.22, 0.18]
    low_scale_idx = 0

    # 4. Annotations
    for x, y in zip(scale_list, profit_list):
        text_color = '#005a00' if y >= 0 else '#8b0000'
        
        if x <= 250:
            # Multi-tiered callouts with pointer lines for low capacities
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
            # Direct top annotations for neat round capacities (> 250 MT)
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

    # Set explicit round x-ticks
    x_ticks = [0, 250, 500, 750, 1000, 1250, 1500, 1750, 2000]
    plt.xticks(x_ticks, [f"{x:,}" for x in x_ticks], fontsize=10, fontweight='bold')
    plt.yticks(fontsize=10, fontweight='bold')

    # Labels and Title
    plt.xlabel('Bioreactor Bags Feedstock [Metric Tonnes / yr]', fontsize=11, fontweight='bold', labelpad=10)
    plt.ylabel('Net Earnings [MM$ / yr]', fontsize=11, fontweight='bold', labelpad=10)
    plt.title('Cytiva STRAP Plant Profitability vs. Scale', fontsize=13, fontweight='bold', pad=15)

    # Border spines
    for spine in ax.spines.values():
        spine.set_linewidth(1.5)
        spine.set_color('black')

    # Expanded bottom margin (-0.32 * y_range) to prevent deep callouts from crossing bottom frame
    plt.ylim(bottom=min(profit_list) - 0.32 * y_range, top=max(profit_list) + 0.14 * y_range)
    plt.tight_layout()
    plt.show()



#%% saving model results and assumptions in excel
tea = process.tea

assumptions_data = {
    'Model Assumptions': [
        'IRR',
        'duration',
        'depreciation',
        'income tax',
        'operating days',
        'lang factor',
        'construction schedule (FCI expenses / yr)',
        'startup months',
        'startup FOC fraction',
        'startup sales fraction',
        'startup VOC fraction',
        'WC over FCI',
        'finance interest',
        'finance years',
        'finance fraction',
        'labor cost',
        'labor burden',
        'property insurance',
        'maintenance',
        'steam power depreciation'
    ],
    'Value': [
        f"{getattr(tea, 'IRR', 0) * 100:.1f}%",
        f"{getattr(tea, 'duration', 0)} years",
        getattr(tea, 'depreciation', 'N/A'),
        f"{getattr(tea, 'income_tax', 0) * 100:.1f}%",
        getattr(tea, 'operating_days', 0),
        getattr(tea, 'lang_factor', None),
        getattr(tea, 'construction_schedule', 'N/A'),
        getattr(tea, 'startup_months', 0),
        getattr(tea, 'startup_FOC_fraction', 0),
        getattr(tea, 'startup_sales_fraction', 0),
        getattr(tea, 'startup_VOC_fraction', 0),
        getattr(tea, 'WC_over_FCI', 0),
        f"{getattr(tea, 'finance_interest', 0) * 100:.1f}%",
        getattr(tea, 'finance_years', 0),
        getattr(tea, 'finance_fraction', 0),
        f"${getattr(tea, 'labor_cost', 0):,.2f}",
        f"{getattr(tea, 'labor_burden', 0) * 100:.1f}%",
        f"{getattr(tea, 'property_insurance', 0) * 100:.2f}%",
        f"{getattr(tea, 'maintenance', 0) * 100:.1f}%",
        getattr(tea, 'steam_power_depreciation', 'N/A')
    ]
}

# Convert to DataFrame and export
df_assumptions = pd.DataFrame(assumptions_data)
#df_assumptions.to_excel('Model_Assumptions.xlsx', index=False)


#%%Itemized costs

def itemized_cost():
    data = []
    
    for u in process.system.units:
        # Safely extract utility_cost or calculate it from power and heat utilities
        util_cost = getattr(u, "utility_cost", None)
    
        if util_cost is None:
            # Fallback: calculate utility cost per hour from power and heat utilities if utility_cost is None
            try:
                power_cost = (
                    u.power_utility.cost if hasattr(u, "power_utility") else 0
                )
                heat_cost = (
                    sum([hu.cost for hu in u.heat_utilities])
                    if hasattr(u, "heat_utilities")
                    else 0
                )
                util_cost = (power_cost or 0) + (heat_cost or 0)
            except Exception:
                util_cost = 0
    
        # Calculate yearly utility cost
        utility_cost_yr = util_cost * process.tea.operating_hours
    
        # Safely get purchase and installed costs
        p_cost = getattr(u, "purchase_cost", 0) or 0
        i_cost = getattr(u, "installed_cost", 0) or 0
    
        data.append(
            {
                "Unit": u.ID,
                "Unit operation": getattr(u, "line", u.__class__.__name__),
                "Purchase cost (10^3 USD)": p_cost / 1e3,
                "Utility cost (10^3 USD/yr)": utility_cost_yr / 1e3,
                "Installed cost (10^3 USD)": i_cost / 1e3,
            }
        )
    
    # Convert to DataFrame
    df_equipment = pd.DataFrame(data)
    
    # Sort by Unit ID
    df_equipment = df_equipment.sort_values(by="Unit").reset_index(drop=True)
    
    # Save to Excel
    df_equipment.to_excel("Equipment_Summary_Table.xlsx", index=False)
    
    # Display
    df_equipment