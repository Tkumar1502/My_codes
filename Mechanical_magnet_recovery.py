# -*- coding: utf-8 -*-
"""
Created on Wed Jun  3 14:19:56 2026

@author: tlnu
"""

import biosteam as bst
from plastics.strap.property_package import STRAP_chemicals_outline, create_property_package, create_property_package_MSW
from biosteam import units
import numpy as np
from biosteam import Unit


#from biosteam import settings
custom_chemicals = create_property_package()

bst.settings.set_thermo(
    [custom_chemicals.NdFeB,
     custom_chemicals.HDPE],
    )

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

#%%setting the plant conditions
capacity = 1000*325     #kg per year
operating_days = 328.5
N_shifts = 3
operating_hours = operating_days*N_shifts*8

feed_per_hour = capacity/operating_hours

feed = bst.Stream(ID = 'Bioreactor Bags', Impeller = feed_per_hour, units = 'kg/hr')


#%%defining units
class Disinfection_Unit(bst.Unit):
    pass

class HandSorting_Unit (bst.Unit):
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
       

class Lathe_Machine (bst.Unit):
    _N_ins = 1
    _N_outs = 2
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, power_kw=5.5, purchase_cost=12000):
        super().__init__(ID, ins, outs, thermo)
        self.power_kw = power_kw          # kW rating of the industrial lathe
        self.custom_purchase_cost = purchase_cost # Capital cost ($) for TEA
        
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
        # 1. TEA: Tell BioSTEAM how much electricity this unit pulls while running
        # BioSTEAM automatically converts power (kW) * operating hours into utility costs
        self.power_utility.consumption = self.power_kw 
        
        # 2. TEA: Map out design parameters if scaling up (optional, defaults to 1)
        self.design_results['Power rating'] = self.power_kw

    def _cost(self):
        # 3. TEA: Define the purchase cost of the machine
        # BioSTEAM's TEA classes use this dictionary to calculate fixed capital investment (FCI)
        self.baseline_purchase_costs['Industrial Lathe'] = self.custom_purchase_cost

class Arbor_Press(bst.Unit):
    _N_ins = 1
    _N_outs = 2
    def __init__(self, ID='', ins=None, outs=(), thermo=None, power_kw=0, purchase_cost=5000):
        super().__init__(ID, ins, outs, thermo)
        self.power_kw = power_kw
        self.custom_purchase_cost = purchase_cost
        
    def _run(self):
        feed=self.ins[0]
        magnets = self.outs[0]
        waste_hdpe = self.outs[1]
        
        magnets.empty()
        waste_hdpe.empty()
        
        magnets.imass['NdFeB'] = feed.imass['NdFeB']
        waste_hdpe.imass['HDPE'] = feed.imass['HDPE']
        
    def _design(self):
        self.power_utility.consumption = self.power_kw

    def _cost(self):
        self.baseline_purchase_costs['Arbor Press'] = self.custom_purchase_cost



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
    
    def _DPI(self, installed_equipment_cost):
        return installed_equipment_cost
    def _TDC(self, DPI): return DPI
    def _FCI(self, TDC): return TDC
    def _FOC(self, FCI): 
        return (FCI*(self.property_tax + self.property_insurance
                     + self.maintenance + self.administration)
                + self.labor_cost*(1+self.fringe_benefits+self.supplies))  # Maintenance insurance fixed at 2% CapEx






#%%
U1 = Disinfection_Unit('U1', ins = feed)
U2 = HandSorting_Unit('U2', ins = U1-0, N_shifts =3, N_workers = 3)
U3 = Lathe_Machine('U3', ins = U2-0)
U4 = Arbor_Press('U4', U3-0)


U1.outs[0].ID = 'sterilized bags'
U2.outs[0].ID = 'impellers'
U3.outs[0].ID = 'trimmed Impellers'
U4.outs[0].ID = 'NdFeB Magnets'

mechanical_magnet_recovery = bst.System('mechanical_magnet_recovery', path = (U1, U2, U3, U4))



mechanical_magnet_recovery.simulate()


tea = MechanicalRecyclingTEA(
    system=mechanical_magnet_recovery,
    IRR=0.15,
    duration=(2026, 2046),
    depreciation='MACRS7',
    income_tax=0.21, # Previously 35% in published study
    operating_days=200,
    lang_factor=3,
    construction_schedule=(0.4, 0.6),
    WC_over_FCI=0.05,
    labor_cost=2.5e6,
    fringe_benefits=0.4,
    property_tax=0.001,
    property_insurance=0.005,
    supplies=0.20,
    maintenance=0.01,
    administration=0.005
)






