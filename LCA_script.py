import biosteam as bst
import pandas as pd
from plastics import strap
import numpy as np
from matplotlib import pyplot as plt
from warnings import filterwarnings
from datetime import datetime

filterwarnings('ignore')

#define scenario
processing_capacity=5e+03
process = strap.BaselineSTRAPProcess(
    scenario='PE/Xylene',
    target_plastic_percent=90, # %
    processing_capacity=processing_capacity, # MT/yr
    sell_leftover_plastic=False,
    burn_leftover_plastic=False,
    facilities=False,
    simulate=False,
)

#create flowchart
process.system.diagram(format='png') # Note how the leftover plastic is sent to heat and power generation.


#bounds
process.set_polymer_mass_fraction.bounds=(.60, .95)
process.set_dissolution_capacity.bounds=(5,10)
process.set_solvent_loss.bounds=(0.1,0.3)
process.set_precipitation_temperature.bounds=(313,323)
process.set_dissolution_temperature.bounds=(363,393)
process.set_centrifuged_plastic_solvent_content.bounds=(25, 75)
process.set_feedstock_distance.bounds=(50,1000)
process.set_feedstock_price.bounds=(0,0.10)


#other process adjustments
process.plastic.ID='Feedstock_plastic'
process.set_dissolution_temperature.baseline = 373
process.set_solvent_loss.baseline = 0.2
process.set_polymer_mass_fraction.baseline = 0.8
process.set_dissolution_capacity.baseline = 7.5
process.set_precipitation_temperature.baseline = 45+273
process.set_feedstock_price.baseline=(0.05)

def MSP(pc):
    process.set_processing_capacity(pc)
    process.system.simulate()
    return process.MSP()

def msp_confidence():
    #remove capacity parameter to set as x in montecarlo simulation
    remove = {'Processing capacity'}
    to_remove = [p for p in process.parameters if p.name in remove]
    for param in to_remove:
        process.parameters.remove(param)
        
    # #run monte carlo
    N_samples = 1000
    rule = 'L' # For Latin-Hypercube sampling
    np.random.seed(1234) # For consistent results
    samples = process.model.sample(N_samples, rule)
    process.model.load_samples(samples)
    process.model.evaluate(
        notify=100 )
     
    # #evaluate across processing capacity
    capacities = np.array([2500,5000,10000,15000, ])
    now = datetime.now()
    filename = now.strftime("MC_%H-%M_%d_%m_%Y")
    process.model.evaluate_across_coordinate('Processing capacity', 
                                             MSP, 
                                             capacities, 
                                             xlfile=f'./{filename}.xlsx'
                                             )
    
    #plot across capacities
    df = pd.read_excel(f'./{filename}.xlsx', sheet_name='- MSP', index_col=0)
    plt.figure()
    bst.plots.plots.plot_montecarlo_across_coordinate(capacities, df)
    plt.ylabel('Cost to produce STRAP resin ($/kg PE)')
    plt.xlabel('Plant Capacity (mton/year)')

# def msp_confidence():
#     remove = {'Processing capacity'}
#     to_remove = [p for p in process.parameters if p.name in remove]
#     for param in to_remove:
#         process.parameters.remove(param)
        
#     #evaluate across processing capacity
#     capacities = np.array([1000,2000,3000,4000,5000])
#     now = datetime.now()
#     filename = now.strftime("MC_%H-%M_%d_%m_%Y")
#     process.model.evaluate_across_coordinate('Processing capacity', 
#                                              MSP, 
#                                              capacities, 
#                                              xlfile=f'./{filename}.csv'
#                                              )
    
#     df = pd.read_excel('./testMC.xlsx', sheet_name='- MSP', index_col=0)
#     bst.plots.plots.plot_montecarlo_across_coordinate(capacities, df)