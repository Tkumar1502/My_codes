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
filterwarnings('ignore')

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
operating_days = 250
N_shifts = 2 
operating_hours = operating_days*N_shifts*8

feed_per_hour = capacity/operating_hours

feed = bst.Stream(ID = 'Bioreactor Bags', Impeller = feed_per_hour, units = 'kg/hr')


#%%defining units
class Disinfection_Unit(bst.Unit):
    _N_ins =1
    _N_outs = 1
    
    def _run (self):
        self.outs[0].copy_like(self.ins[0])

class HandSorting_Unit (bst.Unit):
    N_shifts = 2
    N_workers = 2 
    
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
    
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, power_kw=7.5, 
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
        # 1. TEA: Tell BioSTEAM how much electricity this unit pulls while running
        # BioSTEAM automatically converts power (kW) * operating hours into utility costs
        self.power_utility.consumption = self.power_kw 
        
        # 2. TEA: Map out design parameters if scaling up (optional, defaults to 1)
        self.design_results['Power rating'] = self.power_kw

    def _cost(self):
        # 3. TEA: Define the purchase cost of the machine
        # BioSTEAM's TEA classes use this dictionary to calculate fixed capital investment (FCI)
        self.baseline_purchase_costs['Industrial Lathe'] = self.custom_purchase_cost
        worker_salary = 5e4
        self.total_salay = worker_salary*self.N_shifts*self.N_workers*1.6
        

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
        self.power_utility.consumption = self.power_kw

    def _cost(self):
        worker_salary = 5e4
        self.total_salary = worker_salary*self.N_shifts*self.N_workers*1.6
        
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
    def _TDC(self, DPI): 
        isbl = DPI*0.9
        osbl = DPI*0.1
        warehouse = 0.04*isbl
        site_dev = 0.09*isbl
        piping=0.045*isbl
        return isbl + osbl + warehouse + site_dev + piping
    
    def _FCI(self, TDC): 
        proratable = 0.10*TDC
        field_exp = 0.10*TDC
        construction = 0.20*TDC
        contingency = 0.40*TDC
        other_startup = 0.10*TDC
        TIC = proratable + field_exp + construction + contingency + other_startup
        return TDC + TIC
    
    def _FOC(self, FCI): 
        return (FCI*(self.property_tax + self.property_insurance
                     + self.maintenance + self.administration)
                + self.labor_cost*(1+self.fringe_benefits+self.supplies))  # Maintenance insurance fixed at 2% CapEx
    
    
    
    def CAPEX_table(self):
        """Generates a standardized chemical engineering CAPEX ledger in Millions of USD."""
        
        fci = self.FCI/1e6
        tci = self.TCI/1e6
        wc = self.working_capital/1e6
        tdc = self.TDC/1e6
        DPI = self.DPI/1e6
        
        isbl = DPI*0.9
        osbl = DPI*0.1
        warehouse = 0.04*isbl
        site_dev = 0.09*isbl
        piping=0.045*isbl
        
        proratable = 0.10 * tdc
        field_exp = 0.10 * tdc
        construction = 0.20 * tdc
        contingency = 0.40 * tdc
        other_startup = 0.10 * tdc
        
        tic = proratable + field_exp + construction + contingency + other_startup
        # 6. Build the structural DataFrame
        index = [
            ('Direct costs', 'ISBL installed equipment cost'),
            ('Direct costs', 'OSBL installed equipment cost'),
            ('Direct costs', 'Warehouse'),
            ('Direct costs', 'Site development'),
            ('Direct costs', 'Additional piping'),
            ('Total direct cost (TDC)', ''),
            ('Indirect costs', 'Proratable costs'),
            ('Indirect costs', 'Field expenses'),
            ('Indirect costs', 'Construction'),
            ('Indirect costs', 'Contingency'),
            ('Indirect costs', 'Other (start-up, permits, etc.)'),
            ('Total indirect cost (TIC)', ''),
            ('Fixed capital investment (FCI)', ''),
            ('Working capital (WC)', ''),
            ('Total capital investment (TCI)', '')
        ]
        
        notes = [
            '', '', '4.0% of ISBL', '9.0% of ISBL', '4.5% of ISBL',
            '', '10.0% of TDC', '10.0% of TDC', '20.0% of TDC', '40.0% of TDC', '10.0% of TDC',
            '', 'TDC + TIC', f'{self.WC_over_FCI*100:.1f}% of FCI', 'FCI + WC'
        ]
        
        costs = [
            isbl, osbl, warehouse, site_dev, piping, tdc,
            proratable, field_exp, construction, contingency, other_startup, tic,
            fci, wc, tci
        ]
        
        df = pd.DataFrame(
            {'Notes': notes, 'Cost [MM$]': [round(c, 4) for c in costs]},
            index=pd.MultiIndex.from_tuples(index)
        )
        
        return df 




#%%
U1 = Disinfection_Unit('U1', ins = feed)
U2 = HandSorting_Unit('U2', ins = U1-0, N_shifts =2, N_workers = 2)
U3 = Lathe_Machine('U3', ins = U2-0)

S1 = bst.Splitter('S1', ins=U3.outs[0], outs = ('to Arbor press 1', 'to Arbor press B'), split=0.5)
U4 = Arbor_Press('U4', S1-0)
U5 = Arbor_Press('U5',ins=S1.outs[1])
M1 = bst.Mixer('M1', ins = (U4.outs[0], U5.outs[0]), outs='Total NdFeB Magnets')
M2 = bst.Mixer('M2', ins = (U4.outs[1], U5.outs[1]), outs = 'Total Waste HDPE casings')


U1.outs[0].ID = 'sterilized bags'
U2.outs[0].ID = 'impellers'
U3.outs[0].ID = 'trimmed Impellers'
U4.outs[0].ID = U5.outs[0].ID = 'NdFeB Magnets'
U4.outs[1].ID = U5.outs[1].ID = 'waste hdpe casings'
U3.outs[1].ID = 'waste hdpe shavings'


process = bst.System('mechanical_magnet_recovery', path = (U1, U2, U3,S1,U4,U5,M1,M2))

#products
NdFeB_magnets = M1.outs[0]
HDPE_casing = M2.outs[0] 
HDPE_shavings = U3.outs[1]

#setting up prices
NdFeB_magnets.price = 100
HDPE_casing.price= -0.072
HDPE_shavings.price = -0.072

#feed stock handling and transportation
feed.price = 0.25


process.simulate()
process.diagram()

#adjusted total days for 2 shifts and 250 days, to make 4000 operating_hours
adjusted_days = operating_hours/24

total_labor_cost = U2.total_salary + U4.total_salary + U5.total_salary

tea = MechanicalRecyclingTEA(
    system=process,
    IRR=0.15,
    duration=(2026, 2046),
    depreciation='MACRS7',
    income_tax=0.21, # Previously 35% in published study
    operating_days=adjusted_days,
    lang_factor=None,
    construction_schedule=(1.0,),
    WC_over_FCI=0.05,
    labor_cost=total_labor_cost,
    fringe_benefits=0.4,
    property_tax=0.001,
    property_insurance=0.005,
    supplies=0.05,
    maintenance=0.03,
    administration=0.005
)







