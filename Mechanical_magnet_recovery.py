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
filterwarnings('ignore')

#from biosteam import settings
custom_chemicals = create_property_package()


Films = bst.Chemical('Films',search_ID = 'Polyethylene',phase = 's')
FittingsFilters = bst.Chemical('FittingsFilters', search_ID = 'Polypropylene', phase = 's')
#BrownSupport = bst.Chemical('BrownSupport', search_ID = 'Polypropylene', phase = 's')
BrownSupport = FittingsFilters.copy('BrownSupport')
#SiliconeTubings = bst.Chemical('SiliconeTubings', search_ID = 'Polypropylene', phase = 's')
SiliconeTubings = FittingsFilters.copy('SiliconeTubings')

bst.settings.set_thermo(
    [custom_chemicals.NdFeB,
     custom_chemicals.HDPE,
    'Water',
    'AceticAcid',
    Films,
    FittingsFilters,
    BrownSupport,
    SiliconeTubings]
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

Impeller_composition = {
    HDPE.ID:0.14,
    NdFeB.ID: 0.86
    }




#%%setting the plant conditions
capacity = 1000*1600     #kg per year
operating_days = 250
N_shifts = 2 
operating_hours = operating_days*N_shifts*8

feed_per_hour = capacity/operating_hours

weight_per_bag = 4
water_per_bag = 50

fresh_water_utility = bst.Stream(ID = 'Fresh_water')
fresh_paa_utility = bst.Stream(ID = 'Peracetic Acid')


feed_composition = {'NdFeB': 0.15,
        'HDPE': 0.03,
        'Films': 0.10,
        'FittingsFilters': 0.36,
        'BrownSupport': 0.07,
        'SiliconeTubings': 0.29
        }

feed = bst.Stream(ID = 'Bioreactor Bags', units = 'kg/hr',
                  **{chem: feed_per_hour*frac for chem, frac in feed_composition.items()})
bst.settings.CEPCI = 836.9      #updated to 2025

#%%defining units
#@cost(basis='Volume',ID='Disinfection System',cost=31000.0,S=4.2,n=0.6,BM =2.0,CE=836.9)
class Disinfection_Unit(bst.Unit):
   
   _N_ins = 3
   _N_outs = 3
   
   _units = {'Volume':'m3'}
   
   def __init__(self,ID='',ins = None, outs = (), thermo=None, soak_time = 0.5, cycles_before_dump = 4, 
                N_workers=1,N_shifts=2):
       super().__init__(ID,ins,outs,thermo)
       self.soak_time = soak_time
       self.cycles_before_dump = cycles_before_dump
       self.N_workers = N_workers
       self.N_shifts = N_shifts
    
   def _run(self):
        bags_in, water_in, paa_in = self.ins
        clean_bags_out, evap_loss, wastewater = self.outs
        
        # 1. Mass Throughput Calculations
        plastic_flow = bags_in.F_mass 
        bags_per_hour = plastic_flow / 4.0  
        
        # 2. Fluid Volume Dynamics (70 L/kg water per bag)
        water_needed_per_hour = bags_per_hour * 70.0 
        vessel_water_pool = water_needed_per_hour * self.soak_time 
        
        # 3. Water Losses & Purge Calculations
        hourly_evap_water = 15.0 
        hourly_dragout_water = plastic_flow * 0.1 
        
        time_between_dumps = self.cycles_before_dump * self.soak_time 
        hourly_purge_water = vessel_water_pool / time_between_dumps
        
        total_water_makeup = hourly_purge_water + hourly_evap_water + hourly_dragout_water
        
        # 4. Disinfectant (PAA) Dosing
        paa_concentration = 0.002
        total_paa_makeup = total_water_makeup * paa_concentration
        dirt_flow = plastic_flow * 0.05
        
        # --- Output Streams ---
        wastewater.empty()
        wastewater.imass['Water'] = hourly_purge_water
        wastewater.imass[HDPE.ID] = dirt_flow 
        wastewater.imass['AceticAcid'] = total_paa_makeup 
        
        evap_loss.empty()
        evap_loss.imass['Water'] = hourly_evap_water
        
        clean_bags_out.copy_like(bags_in)
        for chem in ['HDPE','Films','FittingsFilters','BrownSupport','SiliconeTubings']:
            clean_bags_out.imass[chem] = bags_in.imass[chem]*0.95
        
        clean_bags_out.imass['Water'] = hourly_dragout_water
        
        # --- Input Utilities ---
        water_in.empty()
        water_in.imass['Water'] = total_water_makeup
        
        paa_in.empty()
        paa_in.imass['AceticAcid'] = total_paa_makeup

   def _design(self):
        bags_per_hour = self.ins[0].F_mass / 4.0
        batch_bags = bags_per_hour * self.soak_time
        
        working_fluid_m3 = (batch_bags * 70.0) / 1000.0  
        vessel_volume = working_fluid_m3 * 1.2 
        
        self.design_results['Vessels Required'] = 1
        self.design_results['Volume'] = vessel_volume  # Key is 'Volume'
        
        self.power_utility.consumption = 1.5 * vessel_volume
    
   def _cost(self):
        # Fetch using the exact same key: 'Volume'
        V = self.design_results['Volume']
        
        # Calculate Base Cost (baseline = $31,000 for 4.2 m3, scaled to exponent 0.6)
        base_cost = 31000.0 * (V / 4.2)**0.6
        
        # Adjust for Inflation
        
        
        # 1. Baseline Purchase Cost (The raw equipment cost)
        self.baseline_purchase_costs['Disinfection System'] = base_cost
        
        # 2. Final Purchase Cost 
        self.purchase_costs['Disinfection System'] = base_cost
        
        # 3. Installed Cost 
        self.installed_costs['Disinfection System'] = base_cost * 2.0
        
        worker_salary = 5e4
        self.total_salary = worker_salary*self.N_shifts*self.N_workers*1.6

class HandSorting_Unit (bst.Unit):
    _N_ins=1
    _N_outs=2
    
    
    def __init__(self, ID= '', ins =None, outs=(),thermo=None,N_workers=2,N_shifts=2):
        super().__init__(ID,ins,outs,thermo)
        self.N_workers = N_workers
        self.N_shifts = N_shifts
    
    def _run(self):
        feed = self.ins[0]
        other_components = self.outs[1]
        impellers = self.outs[0]
        
        #mass balance logic
        other_components.empty()
        impellers.empty()
        
        #other_components.copy_like(OtherComponents)
        
        impellers.imass['NdFeB'] = feed.imass['NdFeB']
        impellers.imass['HDPE']= feed.imass['HDPE']
        
        other_components.imass['Films'] = feed.imass['Films']
        other_components.imass['FittingsFilters'] = feed.imass['FittingsFilters']
        other_components.imass['BrownSupport'] = feed.imass['BrownSupport']
        other_components.imass['SiliconeTubings'] = feed.imass['SiliconeTubings']
        other_components.imass['Water'] = feed.imass['Water']
     
            
        

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
        self.total_salary = worker_salary*self.N_shifts*self.N_workers*1.6
        

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
        
    
    
    def _DPI(self, installed_equipment_cost):
        return installed_equipment_cost
    def _TDC(self, DPI): 
        isbl = DPI
        osbl = DPI*0.4
        warehouse = 0.04*isbl
        site_dev = 0.25*isbl
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
        
        isbl = DPI
        osbl = DPI*0.4
        warehouse = 0.04*isbl
        site_dev = 0.25*isbl
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
            '', '', '4.0% of ISBL', '25.0% of ISBL', '4.5% of ISBL',
            '', '10.0% of TDC', '10.0% of TDC', '20.0% of TDC', '40.0% of TDC', '10.0% of TDC',
            '', 'TDC + TIC', f'{self.WC_over_FCI*100:.1f}% of FCI', 'FCI + WC'
        ]
        
        costs = [
            isbl, osbl, warehouse, site_dev, piping, tdc,
            proratable, field_exp, construction, contingency, other_startup, tic,
            fci, wc, tci
        ]
        
        df = pd.DataFrame(
            {'Notes': notes, 'Cost [MM$]': [round(c, 3) for c in costs]},
            index=pd.MultiIndex.from_tuples(index)
        )
        
        return df 




#%%
#U1 = Disinfection_Unit('U1', ins = feed)
U1 = Disinfection_Unit('U1', ins=(feed, fresh_water_utility, fresh_paa_utility), 
                       outs=('clean_rinsed_bags', 'evaporative_loss', 'sewer_effluent'))
U2 = HandSorting_Unit('U2', ins = U1-0, N_shifts =2, N_workers = 2)
U3 = Lathe_Machine('U3', ins = U2-0)

S1 = bst.Splitter('S1', ins=U3.outs[0], outs = ('to Arbor press 1', 'to Arbor press B'), split=0.5)
U4 = Arbor_Press('U4', S1-0)
U5 = Arbor_Press('U5',ins=S1.outs[1])
M1 = bst.Mixer('M1', ins = (U4.outs[0], U5.outs[0]), outs='Total NdFeB Magnets')
M2 = bst.Mixer('M2', ins = (U4.outs[1], U5.outs[1]), outs = 'Total Waste HDPE casings')
S2 = bst.Splitter('S2',ins=U2.outs[1],split=1)
S2.isplit['Water']=0



U1.outs[0].ID = 'sterilized bags'
U2.outs[0].ID = 'impellers'
U2.outs[1].ID = 'Other components'
U3.outs[0].ID = 'trimmed Impellers'
U4.outs[0].ID = U5.outs[0].ID = 'NdFeB Magnets'
U4.outs[1].ID = U5.outs[1].ID = 'waste hdpe casings'
U3.outs[1].ID = 'waste hdpe shavings'
S2.outs[0].ID = 'other components'
S2.outs[1].ID = 'water'

process = bst.System('mechanical_magnet_recovery', path = (U1, U2, U3,S1,U4,U5,M1,M2,S2))

#products
NdFeB_magnets = M1.outs[0]
HDPE_casing = M2.outs[0] 
HDPE_shavings = U3.outs[1]
other_components = S2.outs[0]
waste_water = S2.outs[1]

#setting up prices
NdFeB_magnets.price = 100
HDPE_casing.price= -0.072
HDPE_shavings.price = -0.072
other_components.price = -0.072
waste_water.price = -0.002642                #considering 10usd for 1000 gallons

#feed stock handling and transportation
feed.price = 0.25


#raw materials cost
fresh_water_utility = 0.002
fresh_paa_utility = 2.50


process.simulate()
process.diagram()

#adjusted total days for 2 shifts and 250 days, to make 4000 operating_hours
adjusted_days = operating_hours/24



# handsorting=2 workers, 1 worker lathe machine, 2 workers on 2 arbor press
total_labor_cost = U1.total_salary +  U2.total_salary + U4.total_salary + U5.total_salary + U3.total_salary

tea = MechanicalRecyclingTEA(
    system=process,
    IRR=0.15,
    duration=(2026, 2046),
    depreciation='MACRS7',
    income_tax=0.21, # Previously 35% in published study
    operating_days=adjusted_days,
    lang_factor=None,
    construction_schedule=(1.0,),
    WC_over_FCI=0.12,
    labor_cost=total_labor_cost,
    fringe_benefits=0.4,
    property_tax=0.001,
    property_insurance=0.005,
    supplies=0.05,
    maintenance=0.03,
    administration=0.005
)







