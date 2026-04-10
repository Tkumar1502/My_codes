#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 14 13:31:25 2024

@author: charlesgranger
"""
import matplotlib.pyplot as plt
from plastics.strap import STRAPProcessPE
import math
import pandas as pd

#for quick tests##############################
"""feed_capacity =  5000 #metric tpy
hrly_feed_capacity = feed_capacity*1000/8760 #kg/hr

pm = STRAPProcessPE()
pm.set_processing_capacity(hrly_feed_capacity)
pm.system.simulate()
pm.Xylene.show('wt')
"""


def get_MSP(x):
    PE_fraction = 0.80  #this line doesn't do anything right now
    feed_capacity =  x #metric tpy
    hrly_feed_capacity = feed_capacity*1000/8760 #kg/hr
    
    pm = STRAPProcessPE()
    pm.set_processing_capacity(hrly_feed_capacity)
    
    pm.tea.IRR=0.20 #assumed company IRR
    pm.tea.duration=(2025, 2045)
    pm.tea.depreciation='MACRS7'  #7 year depreciation schedule
    pm.tea.income_tax=0.35
    pm.tea.operating_days=0.80 * 365  #88% time on stream
    pm.tea.lang_factor=None 
    pm.tea.construction_schedule= (0.08, 0.60, 0.32) #construction costs split across first 3 years
    pm.tea.startup_months=6 
    pm.tea.startup_FOCfrac=1 #100% fixed operating costs during startup
    pm.tea.startup_salesfrac=0.5 
    pm.tea.startup_VOCfrac=0.75 #75% variable operating costs during startup
    pm.tea.WC_over_FCI=0.05 # % of total equipment
    pm.tea.finance_interest=0.10  #10% interest
    pm.tea.finance_years=10
    pm.tea.finance_fraction=0.4 #remaining 60% is capital
    pm.tea.warehouse=0.04 # % of total equipment
    pm.tea.site_development=0.09 # % of total equipment
    pm.tea.additional_piping=0.045 # % of total equipment
    pm.tea.proratable_costs=0.10 # % of total equipment
    pm.tea.field_expenses=0.10 # % of total equipment
    pm.tea.construction=0.20 # % of total equipment
    pm.tea.contingency=0.5 #50% additional cost of capital for contingency
    pm.tea.other_indirect_costs=0.10 # % of total equipment
    pm.tea.labor_cost= math.ceil((2 + math.ceil(feed_capacity/ 12000))*4.8) * 80000
    pm.tea.labor_burden=0.90 #overhead + benefits = 90% of salaries
    pm.tea.property_insurance=0.007 # % of total equipment
    pm.tea.maintenance=0.03 # % of total equipment
    pm.tea.steam_power_depreciation='MACRS20'
    
    #pm.PE_resin.price = pm.tea.solve_price(pm.PE_resin)
    pm.system.simulate()
    #pm.system.diagram()
    #pm.show()
    operators = 2 + math.ceil(feed_capacity / 12000)
    #print('operators: ', operators, '/shift')
    #print(pm.MSP())
    #table = pm.tea.get_cashflow_table()
    return pm.MSP()


def get_DCFROR(x, y):
     
    feed_capacity =  x #metric tpy
    hrly_feed_capacity = feed_capacity*1000/8760 #kg/hr
    
    pm = STRAPProcessPE()
   
    pm.set_processing_capacity(hrly_feed_capacity)
    
    #set product price
    pm.PE_resin.price = y
    
    
    #pm.tea.IRR=0.20 #assumed company IRR
    pm.tea.duration=(2025, 2045)
    pm.tea.depreciation='MACRS7'  #7 year depreciation schedule
    pm.tea.income_tax=0.35
    pm.tea.operating_days=0.90 * 365  #90% time on stream
    pm.tea.lang_factor=None 
    pm.tea.construction_schedule= (0.08, 0.60, 0.32) #construction costs split across first 3 years
    pm.tea.startup_months=6 
    pm.tea.startup_FOCfrac=1 #100% fixed operating costs during startup
    pm.tea.startup_salesfrac=0.5 
    pm.tea.startup_VOCfrac=0.75 #75% variable operating costs during startup
    pm.tea.WC_over_FCI=0.05 # % of total equipment
    pm.tea.finance_interest=0.10  #10% interest
    pm.tea.finance_years=10
    pm.tea.finance_fraction=0.0 #remaining 60% is capital
    pm.tea.warehouse=0.04 # % of total equipment
    pm.tea.site_development=0.09 # % of total equipment
    pm.tea.additional_piping=0.045 # % of total equipment
    pm.tea.proratable_costs=0.10 # % of total equipment
    pm.tea.field_expenses=0.10 # % of total equipment
    pm.tea.construction=0.20 # % of total equipment
    pm.tea.contingency=0.5 #50% additional cost of capital for contingency
    pm.tea.other_indirect_costs=0.10 # % of total equipment
    pm.tea.labor_cost= math.ceil((2 + math.ceil(feed_capacity/ 12000))*4.8) * 80000
    pm.tea.labor_burden=0.90 #overhead + benefits = 90% of salaries
    pm.tea.property_insurance=0.007 # % of total equipment
    pm.tea.maintenance=0.03 # % of total equipment
    pm.tea.steam_power_depreciation='MACRS20'
    
    #pm.PE_resin.price = pm.tea.solve_price(pm.PE_resin)
    pm.system.simulate()
    #pm.system.diagram()
    #pm.show()
    operators = 2 + math.ceil(feed_capacity / 12000)
    #print('operators: ', operators, '/shift')
    #print(pm.MSP())
    #table = pm.tea.get_cashflow_table()
    return pm.tea.solve_IRR()
"""# #calculate MSP for a defined size range
size_range = [1000,2500,5000,10000,20000,30000,40000,50000]  #tpy
MSPs = []
TCIs = []
for i in range(len(size_range)):
    MSPs.append(get_MSP(size_range[i]))
    #TCIs.append()
print(MSPs)


#plot MSP vs capactiy
plt.plot(size_range, MSPs)
plt.xlabel('metric tons per year PIW/PCW')
plt.ylabel('MSP ($/kg)')

plt.savefig('./images/msp_capacity_cytiva.png',dpi=300)

plt.show()"""

#get MSP for defined size and IRR
#msp = get_MSP(5000)
#print('MSP = ', msp)

#get DCFROR for defined size and product price
DCFROR = get_DCFROR(50, 3)

#Run system
feed_capacity =  50 #metric tpy
hrly_feed_capacity = feed_capacity*1000/8760 #kg/hr
pm = STRAPProcessPE()
pm.set_processing_capacity(hrly_feed_capacity)
pm.system.simulate()

pm.PE_resin.price = 3
pm.tea.IRR = DCFROR

#get cashflow table based on MSP

cf_table = pm.tea.get_cashflow_table()
#pm.tea.save_report('./results/5ktkon_report_09202024.xlsx')

#cf_table.to_excel('./results/30kton_cashflow_09202024-2.xlsx', index=False)

