# -*- coding: utf-8 -*-
"""
Created on Thu Mar 19 15:32:58 2026

@author: tlnu
"""

#adding this on macbook 4:49pm-4/3/26

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

#creating a storage tankf for magnets
S_mag = bst.Splitter(split = 0.137)
S_mag.ins[0] = process.T4.outs[0]
S_mag.outs[0] = process.P3.ins[0]
units.append(S_mag)

T_mag = bst.StorageTank(ins= S_mag.outs[1])
units.append(T_mag)
T_mag.outs[0].ID = 'NdFeB Magnets'

#removing unnecessary units for magnets
process.U6.ins[0] = process.P3.outs[0]
units.remove(process.U4)
units.remove(process.U5)   
units.remove(process.M2)     


#removing the streams
process.M3.ins[:] = [process.M3.ins[i] for i in (0, 2, 3)]
#process.T5.ins.append(process.s22)

process.system.update_configuration(units)  #we need update the system after adding handsorting


process.system.diagram(format = 'svg') #just to check if the changes are working or not