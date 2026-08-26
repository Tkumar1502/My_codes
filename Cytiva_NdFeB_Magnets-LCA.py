# -*- coding: utf-8 -*-
"""
Created on Tue Apr 21 15:31:58 2026

@author: tlnu
"""

import biosteam as bst
import importlib

bst.settings.ID_magic = False
from plastics import strap
from chaospy.distributions import Uniform
import pandas as pd
from warnings import filterwarnings
import numpy as np
import matplotlib.pyplot as plt
import os
from datetime import datetime
from matplotlib.ticker import PercentFormatter, FuncFormatter
from matplotlib.colors import TwoSlopeNorm, Normalize



filterwarnings('ignore')

bst.main_flowsheet.clear()

if hasattr(strap.MagnetRecovery, 'cache'):
    strap.MagnetRecovery.cache.clear()

capacity = 1976
bst.settings.CEPCI = 836.9
process = strap.MagnetRecovery(
    processing_capacity = capacity,      #tons
    sell_leftover_plastic = True,
    simulate=False,
    mbbr = False)


feed_composition = {'NdFeB': 0.1548,
        'HDPE': 0.0252,
        'Films': 0.10,
        'FittingsFilters': 0.36,
        'BrownSupport': 0.07,
        'SiliconeTubings': 0.29,
        'Solutes': 0.001
        }

process.tea.operating_hours = 309*2*8       #days*shifts*hours
feed_per_hour = capacity*1000/process.tea.operating_hours

process.feedstock.empty()
for chem, frac in feed_composition.items():
    process.feedstock.imass[chem]=feed_per_hour*frac
    
    
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
  

process.dissolution_step.T = 383.15

process.tea.labor_cost = 580000 + process.HS.total_salary

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

process.HDPE_resins.price = 1.20
process.NdFeB_Magnets.price = 100

#add processing parameters again if flow rates (input) starts funky

#process.set_processing_capacity(325)    #somehow changing the order changes the input flow rate, add your processing capacity after declaring tea hours/days

process.system.update_configuration()



process.system.simulate()

process.system.diagram()



# %% UPDATED POWER RATES
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
    'CWP': 0.0,
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



#%%     LCA

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


#%%natural gas CF
'''
process.natural_gas.set_CF(GWP, 0.61146) # kg CO2 eq / kg NG
process.natural_gas.set_CF(FFC,50.937, )# 51MJ/kg from EF 2.0
process.natural_gas.set_CF(WU,0.01633, )
process.natural_gas.set_CF(HTC,2.708e-10, )
process.natural_gas.set_CF(HTNC,3.5399e-9 )
process.natural_gas.set_CF(ETOX, 0.02579)
process.natural_gas.set_CF(ACD,0.00130)
process.natural_gas.set_CF(OZD, 2.8956e-11)
process.natural_gas.set_CF(POCP,.00124)                        
'''
#%%add CFs for xylene
# EF 2.0; Xylene production, production mix, at plant, technology mix, 100% active substance
# solvent vapor trap incineration: EF 2.0 Waste incineration of inert material, production mix, at consumer, waste-to-energy plant with dry flue gas treatment, including transport and pre-treatment, inert material waste
process.solvent.set_CF(FFC,52.33 + 2.260, )#EF 2.0; Xylene production, production mix, at plant, technology mix, 100% active substance
process.solvent.set_CF(GWP,0.9383 + 0.1563, )
process.solvent.set_CF(WU,0.16893 + 0.12196, )
process.solvent.set_CF(HTC,2.6045e-8 + 3.56126e-10, )
process.solvent.set_CF(HTNC,8.74695e-8 + 8.0464e-9, )
process.solvent.set_CF(ETOX,0.8608 + 0.00900, )
process.solvent.set_CF(ACD,0.00463 + 0.00048, )
process.solvent.set_CF(OZD,2.78856e-9 + 2.8967e-11, )
process.solvent.set_CF(POCP,0.00335 + 0.00044, )

#%%CFs for water
process.makeup_water.set_CF(FFC,5.09/1000, )#MJ/kg
process.makeup_water.set_CF(GWP,0.33/1000, ) #kg co2 eq / kg water
process.makeup_water.set_CF(WU, 4.3e-2/1000) #m3/kg cooling water
process.makeup_water.set_CF(HTC, 3.6e-14)  #CTUh / kg cooling water
process.makeup_water.set_CF(HTNC, 1.4e-13) #CTUh / kg cooling water
process.makeup_water.set_CF(ETOX, 4.3e-7) #CTUe / kg cooling water
process.makeup_water.set_CF(ACD, 2.7918e-7) #
process.makeup_water.set_CF(OZD, 1.314e-12) #!!!!!!!
process.makeup_water.set_CF(POCP, 1.505e-7)#!!!!!!!
#%% CFs for cooling tower makeup water
#same values as makeup water...
process.cooling_tower_makeup_water.set_CF(FFC,5.09/1000, )
process.cooling_tower_makeup_water.set_CF(GWP,0.33/1000, )
process.cooling_tower_makeup_water.set_CF(WU, 4.3e-2/1000) #m3/kg cooling water
process.cooling_tower_makeup_water.set_CF(HTC, 3.6e-14)  #CTUh / kg cooling water
process.cooling_tower_makeup_water.set_CF(HTNC, 1.4e-13) #CTUh / kg cooling water
process.cooling_tower_makeup_water.set_CF(ETOX, 4.3e-7) #CTUe / kg cooling water
process.cooling_tower_makeup_water.set_CF(ACD, 2.7918e-7) #!!!!!!!
process.cooling_tower_makeup_water.set_CF(OZD, 1.314e-12) #!!!!!!!
process.cooling_tower_makeup_water.set_CF(POCP, 1.505e-7)#!!!!!!!

#%%CFs for adsorbent
#EF 2.0 activated silica production, production mix, at plant, technology mix, 100% active substance (activated carbon not available)
#adsorbent landfilling: EF 2.0 Landfill of polluted inorganic waste, production mix (region specific sites), at landfill site, landfill including leachate treatment and with transport without collection and pre-treatment
process.adsorption_column.ins[2].ID = 'Adsorbent'
process.adsorption_column.ins[2].set_CF(GWP, 1.78+0.02643,  ) 
process.adsorption_column.ins[2].set_CF(FFC, 23.99+0.3452,  ) 
process.adsorption_column.ins[2].set_CF(WU,1.0335 + .00207  )
process.adsorption_column.ins[2].set_CF(HTC, 5.247e-8+ 1.354e-8 )
process.adsorption_column.ins[2].set_CF(HTNC,3.247e-7 + 1.569e-8 )
process.adsorption_column.ins[2].set_CF(ETOX,1.471 + 0.2591 )
process.adsorption_column.ins[2].set_CF(ACD,0.01891 + .00015 )
process.adsorption_column.ins[2].set_CF(OZD,2.833e-9 + 4.311e-14 )
process.adsorption_column.ins[2].set_CF(POCP,0.00619 + .00207 )



#%%CFs for peracetic acid
process.system.feeds[2].ID = "peracetic_acid"
process.system.feeds[2].set_CF(GWP, 0.0)


#%%CFs for freshwater
process.U10.ins[1].ID = 'fresh_water'
process.U10.ins[1].set_CF(GWP,0.0)


#%%
process.tea.operating_hours = 309*8*2
process.system.operating_hours = process.tea.operating_hours

'''for feed in process.system.feeds:
    if feed.characterization_factors.get('GWP') is None:
        feed.characterization_factors['GWP'] = 0.0
'''
inventory_table = bst.report.lca_inventory_table(systems =[process.system],keys=GWP, items=process.system.outs)

#LCA
GWP = 'GWP'



for i in process.parameters:
    i.distribution = Uniform(*i.bounds)



#inventory_table = bst.report.lca_inventory_table([process.system], keys = GWP, items = [process.products])

GWP_table = bst.report.lca_displacement_allocation_table(
    [process.system],
    'GWP',
    products, 
)

revenue_allocation = bst.report.lca_property_allocation_factor_table(
    [process.system],
    property='revenue',
)

energy_allocation = bst.report.lca_property_allocation_factor_table(
    [process.system],
    property='energy',
)

gwp_allocation = bst.report.lca_displacement_allocation_factor_table(
    [process.system],
    items = process.products,
    key =   GWP
)

