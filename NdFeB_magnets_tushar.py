# -*- coding: utf-8 -*-
"""
Created on Fri Feb 27 13:30:36 2026

@author: tlnu
"""

import biosteam as bst
from plastics import strap
from chaospy.distributions import Uniform
import pandas as pd
from warnings import filterwarnings
import numpy as np




filterwarnings('ignore')

bst.main_flowsheet.clear() #resets names for consistency


#define scenario

#total bags = 494,214; wt. of each bag = 657g
processing_capacity = 325 #tons/year                        #should this be in ktons or tons

process = strap.BaselineSTRAPProcess(
    scenario = 'PE/Xylene',
    target_plastic_percent = 0.137, #HDPE from the impellers (avg)
    
    processing_capacity = processing_capacity,
    sell_leftover_plastic = False,
    burn_leftover_plastic = False, 
    facilities = False,      #not sure what is this
    simulate = False,       #don't change this!
)


#create flowchart
process.plastic.ID = 'Bioreactor Bags'
#process.system.diagram(format = 'svg')


#now we have to remove shredding and add handsorting

#Create handsorting unit
class Handsorting (bst.Unit):
    pass


#connecting handsorting and removing shredding
units = process.system.units

HS = Handsorting() #creating handsorting as unit

#adjusting the streams going in and outs
HS.ins[0] = process.plastic #plastic here is the feed that you entered
process.T1.ins[0] = HS.outs[0]
process.U3.ins[0] = process.T1.outs[0]
units.remove(process.U1)

units.append(HS)            #have to make sure you append any new unit to see the changes in the flowsheet



#changing the names for the streams?
HS.outs[0].ID = 'Impellers'         #changing the sorted feed of interest to impellers
process.M2.outs[0].ID = 'NdFeB Magnets'
process.U9.outs[0].ID = 'HDPE resins'


process.system.update_configuration(units)  #we need update the system after adding handsorting


process.system.diagram(format = 'svg') #just to check if the changes are working or not



#TEA pareameters
process.tea.operating_days = 328.5
process.tea.labor_cost = 1350000 # 2 operators x $60K, 1 loader x $50K, 1 engineer x $100K


#define products and set sale prices
products = [process.PE_resin, process.leftover_plastic]

#price for PE
process.PE_resin.price = 1.20

#price for Neodymium (enter zero if not selling)
process.leftover_plastic.price = 100 


#bounds for sensitivity analysis
process.set_polymer_mass_fraction.bounds = (0.15, 0.35)

#probably need to change it, it is too diluted!
process.set_dissolution_capacity.bounds = (2.13, 2.25)        #percentage of plastics in solvent

process.set_solvent_loss.bounds = (0.01, 1)                 #in percent
process.set_dissolution_temperature.bounds = (363, 403)
process.set_precipitation_temperature.bounds = (313,323)     #not sure if it should be room temperature
process.set_centrifuged_plastic_solvent_content.bounds = (25,75)
process.set_feedstock_distance.bounds = (50,1000)
process.set_feedstock_price.bounds = (0.03, 0.08)           #can be updated when we include handsorting in it



#set baseline metrices
process.set_dissolution_temperature.baseline = 373
process.set_solvent_loss.baseline = 0.1
process.set_polymer_mass_fraction.baseline = 0.137
process.set_dissolution_capacity.baseline = 2.19            #needs to adjusted, too diluted     
process.set_precipitation_temperature.baseline = 45+273
process.set_feedstock_price.baseline = 0.25                 #will be adjusted with handsorting
process.set_IRR.baseline = 0.15


#why this loop is here
for i in process.parameters:
    i.distribution = Uniform(*i.bounds)                     #makes equal probability for points in between the bounds

assumptions,results = process.baseline()
assumption_table = pd.DataFrame(assumptions)


'''TEA analysis and other functions '''

def get_MSP(dt):
    #process.set_dissolution_temperature.baseline = dt
    process.dissolution_step.T = dt
    process.system.simulate()
    print('MSP for', process.set_dissolution_temperature.baseline, process.MSP())


def MSP_at_PE_mass_fraction_and_dissolution_capacity(mass_fraction, 
                                                     dissolution_capacity, process):
    process.set_polymer_mass_fraction(mass_fraction)
    process.set_dissolution_capacity(dissolution_capacity)
    process.system.simulate()
    return process.MSP()

def MSP(pc):                                    #MSP at any processing capacity (ktons)
     process.set_processing_capacity(pc)
     process.system.simulate()
     return process.MSP()
 
#creating different processing capacities
lb,ub = processing_capacity/2, processing_capacity*5
pcs = np.linspace(lb,ub,30)


#calculate unit operating cost for given capacity
def UOC(pc):
    process.set_processing_capacity(pc)
    process.system.similate()
    aoc = process.tea.AOC
    
    #not exactly sure what is this? and why it's here?
    landfill = process.leftover_plastic.price*process.leftover_plastic.F_mass*process.tea_operating_hours
    
    return (aoc-landfill)/(process.PE_resin.F_mass*process.tea.operating_hours)


def save_UOC():
    UOCs = [UOC(i) for i in pcs]
    
    #creating a dataframe with PCs, UOCs
    df = pd.DataFrame({'pcs': pcs, 'UOC': UOCs})
    
    return df


#calculte total capital investment (TCI) for given capacity
def TCI(pc):
    process.set_processing_capacity(pc)
    process.system.simulate()
    return process.tea.TCI

#calculate dROI for given capacity
def dROI(pc):
    process.set_processing_capacity(pc)
    process.system.simulate()
    droi = process.tea.NPV/process.tea.TCI
    return droi

