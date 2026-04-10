#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 27 15:40:51 2025

@author: charlesgranger
"""
############### CHANGES TO BE MADE #####################################
# everything runs much slower
# MSP from MC is different
# parameters are not being excluded from single point sensitivity anymore


#4) can we add a method to remove/lose a fraction in the filters?
#5) add feedtsock distance to LCA parameters [lca done, price needed]

#8) need to account for burning of lost solvent emissions
#9) spearman's rank restets the parameters. is that because i am only defining bound 
#but not distributions, so it defaults to the original definition?

#add abands to msp plot --> evaluate across coordinate
#biosteam>plots>plots>evaluate_across...

#before sampling, remove processing capacity parameter, 
#monte carlo sampling 18.1.9

########################################################################

import biosteam as bst
import pandas as pd
from plastics import strap
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.patches as mpatches
from warnings import filterwarnings
from datetime import datetime
from chaospy.distributions import Uniform
filterwarnings('ignore')

#define scenario
processing_capacity=200
process = strap.BaselineSTRAPProcess(
    scenario='PE/Xylene',
    target_plastic_percent=50, # %
    processing_capacity=processing_capacity, # MT/yr
    sell_leftover_plastic=False,
    burn_leftover_plastic=False,
    facilities=False,
    simulate=False,
    precipitation_configuration = 'solvent mixing'
)

#create flowchart
process.system.diagram(format='png') # Note how the leftover plastic is sent to heat and power generation.


#bounds
process.set_polymer_mass_fraction.bounds=(.60, .95)
process.set_dissolution_capacity.bounds=(5,10)
process.set_solvent_loss.bounds=(0.01,1)
process.set_precipitation_temperature.bounds=(313,323)
process.set_dissolution_temperature.bounds=(363,393)
process.set_centrifuged_plastic_solvent_content.bounds=(25, 75)
process.set_feedstock_distance.bounds=(50,1000)
process.set_feedstock_price.bounds=(0,0.10)

#other process adjustments
process.plastic.ID='Feedstock_plastic'
process.set_dissolution_temperature.baseline = 373
process.set_solvent_loss.baseline = 0.1
process.set_polymer_mass_fraction.baseline = 0.8
process.set_dissolution_capacity.baseline = 7.5
process.set_precipitation_temperature.baseline = 45+273
process.set_feedstock_price.baseline=(0.05)
process.set_centrifuged_plastic_solvent_content.baseline = 30
process.set_IRR.baseline=0.15

for i in process.parameters:
    i.distribution = Uniform(*i.bounds)

#get assumptions for print
assumptions, results = process.baseline()

assumptions_table = pd.DataFrame(assumptions)


def get_MSP(dt):
    #process.set_dissolution_temperature.baseline = dt
    process.dissolution_step.T = dt
    process.system.simulate()
    print('MSP for ', process.dissolution_step.T, process.MSP())

def get_GWP(x):
    #process.set_dissolution_temperature.baseline = dt
    process.dissolution_step.T = x
    process.system.simulate()
    print('GWP for ', process.dissolution_step.T, process.GWP())

#####################contour sensitivity plot ########################
#plot MSP vs %PE and %in solution
def MSP_at_PE_mass_fraction_and_dissolution_capacity(mass_fraction, dissolution_capacity, process):
    process.set_polymer_mass_fraction(mass_fraction)
    process.set_dissolution_capacity(dissolution_capacity)
    process.system.simulate()
    return process.MSP()

def contour_plot():
    xlim = np.array([0.6, 0.9])
    ylim = np.array([4, 10])
    X, Y, Z = bst.plots.generate_contour_data(
        MSP_at_PE_mass_fraction_and_dissolution_capacity,
        xlim=xlim, ylim=ylim, args=(process,),
        n=6,
    )
    
    # Plot contours
    xlabel = "PE content [wt %]"
    ylabel = 'Dissolution capacity [wt %]'
    xticks = [60, 70, 80, 90]
    yticks = [4, 6, 8, 10]
    metric_bar = bst.plots.MetricBar(
        'MSP', 'USD$\cdot$kg$^{-1}$', plt.cm.get_cmap('viridis_r'), 
        [.80, .90, 1.00, 1.10, 1.20, 1.30, 1.40], 15, 2
    )
    fig, axes, CSs, CB, other_axes = bst.plots.plot_contour_single_metric(
        100 * X, Y, Z[:, :, None, None], xlabel, ylabel, xticks, yticks, metric_bar,  
        fillcolor=None, styleaxiskw=dict(xtick0=False), label=True,
    )

############### capacity vs MSP & CAPEX plot ############################
lb, ub = processing_capacity/2, processing_capacity*5
pcs = np.linspace(lb, ub, 30)

#calculate MSP for given capacity
def MSP(pc):
    process.set_processing_capacity(pc)
    process.system.simulate()
    return process.MSP()

#calculate total capital investment (TCI) for given capacity
def TCI(pc):
    process.set_processing_capacity(pc)
    process.system.simulate()
    return process.tea.TCI

#get data
def msp_plot():
    MSPs = [MSP(i) for i in pcs]
    TCIs = [TCI(i) for i in pcs]
    market = [1.17 for i in pcs]
    
    #plot data
    plt.figure()
    line1, = plt.plot(pcs, MSPs, label='STRAP resin cost') # plot MSP for different capacities
    line2, = plt.plot(pcs, market, 'r--', label='Market Value PE') #market value
    plt.ylabel('Cost to produce STRAP resin ($/kg PE)')
    plt.xlabel('Plant Capacity (mton/year)')
    
    # Create a secondary y-axis
    ax2 = plt.gca().twinx()
    line3, = ax2.plot(pcs, TCIs, label='TCI', color='green')
    ax2.set_ylabel('Total Capital Investment ($)')  # Label for the second y-axis
    
    # Combine legends for both axes
    lines = [line1, line2, line3]
    labels = [line.get_label() for line in lines]
    plt.legend(lines, labels, loc='upper center')

def msp_confidence():
    TCIs = [TCI(i) for i in pcs]
    market = [1.17 for i in pcs]
    
    #remove capacity parameter to set as x in montecarlo simulation
    # remove = {'Processing capacity'}
    # to_remove = [p for p in process.parameters if p.name in remove]
    # for param in to_remove:
    #     process.parameters.remove(param)
        
    # #run monte carlo
    N_samples = 100
    rule = 'L' # For Latin-Hypercube sampling
    np.random.seed(1234) # For consistent results
    samples = process.model.sample(N_samples, rule)
    process.model.load_samples(samples)
     
    # #evaluate across processing capacity
    now = datetime.now()
    filename = now.strftime("MC_%H-%M_%d_%m_%Y")
    process.model.evaluate_across_coordinate('Processing capacity', 
                                             process.set_processing_capacity, 
                                             pcs, 
                                             xlfile=f'./{filename}.xlsx'
                                             )
    
    #plot across capacities
    df = pd.read_excel(f'./{filename}.xlsx', sheet_name='- MSP', index_col=0)
    plt.figure(figsize=(10, 6))
    percentiles = bst.plots.plots.plot_montecarlo_across_coordinate(pcs, df)
   
   # Extract percentiles
    p5, p25, p50, p75, p95 = percentiles  # 5th, 25th, 50th (median), 75th, 95th percentiles
    
    # Add legend manually for Monte Carlo results
    line_median, = plt.plot(pcs, p50, 'k-', label='Median Cost')  # Center line
    fill = plt.fill_between(pcs, p25, p75, color='brown', alpha=0.3, label='Interquartile Range (25th-75th)')  # Shaded region
    line_ci, = plt.plot(pcs, p5, 'k--', label='Confidence Interval (5th-95th)')  # Lower bound
    plt.plot(pcs, p95, 'k--')  # Upper bound (same label as lower)
    
    # Add other plots
    line2, = plt.plot(pcs, market, 'r--', label='Market Value PE')  # Market value
    plt.ylabel('Cost to produce STRAP resin ($/kg PE)')
    plt.xlabel('Plant Capacity (mton/year)')
    
    # Create a secondary y-axis
    ax2 = plt.gca().twinx()
    line3, = ax2.plot(pcs, TCIs, label='TCI', color='green')
    ax2.set_ylabel('Total Capital Investment ($)')

    # Combine legends
    lines = [line_median, fill, line_ci, line2, line3]
    labels = [line.get_label() for line in lines]
    plt.legend(lines, labels, loc='upper center')
    # #plot data
   #  #line1, = plt.plot(pcs, MSPs, label='STRAP resin cost') # plot MSP for different capacities
 
   #  # Create a secondary y-axis
   #  ax2 = plt.gca().twinx()
   #  line3, = ax2.plot(pcs, TCIs, label='TCI', color='green')
   #  ax2.set_ylabel('Total Capital Investment ($)')  # Label for the second y-axis
    
   #  # Combine legends for both axes
   #  lines = [monte_carlo_patch,line2, line3]
   #  labels = [line.get_label() for line in lines]
   #  plt.legend(lines, labels, loc='upper center')
############### tornado plot ############################################

############################## LCA #####################################

######## SET CHARACTERIZATION FACTORS #########################
GWP = 'GWP'
FFC = 'FFC'

#natural gas CF
process.natural_gas.set_CF(GWP, 0.61146) # kg CO2 eq / kg NG

process.natural_gas.set_CF(
    FFC,
    51, # 51MJ/kg from EF 2.0
)

#add CFs for solvent
process.solvent.set_CF(
    FFC,
    52.33, #EF 2.0; Xylene production, production mix, at plant, technology mix, 100% active substance
)

process.solvent.set_CF(
    GWP,
    0.9383, # EF 2.0; Xylene production, production mix, at plant, technology mix, 100% active substance
)

#CFs for water
process.makeup_water.set_CF(
    FFC,
    5.09/1000, #MJ/kg
)

process.makeup_water.set_CF(
    GWP,
    0.33/1000,  #kg co2 eq / kg water
)

process.cooling_tower_makeup_water.set_CF(
    FFC,
    5.09/1000, #MJ/kg
)

process.cooling_tower_makeup_water.set_CF(
    GWP,
    0.33/1000,  #kg co2 eq / kg water
)

#CFs for adsorbent
process.adsorption_column.ins[2].ID = 'Adsorbent'
process.adsorption_column.ins[2].set_CF(GWP, 1.78,  ) #EF 2.0 activated silica production, production mix, at plant, technology mix, 100% active substance (activated carbon not available)
process.adsorption_column.ins[2].set_CF(FFC, 23.99,  ) 

#CFs for transportation
#values taken from aurora's paper, database EF 2.0 (EU)
#radius = 1000 #km plastic travels
#process.plastic.set_CF(GWP, process.feedstock_distance * 0.08 * 1/1000, ) #distance (km) * .08 kg co2 / (t*km) * t/1000kg
#process.plastic.set_CF(FFC, process.feedtsock_distance * 1.04 * 1/1000, ) #distance (km) * .08 kg co2 / (t*km) * t/1000kg

# ELECTRICITY / MJ "Electricity grid mix 1kV-60kV , consumption mix, to consumer, AC, technology mix, 1kV - 60kV"
#climate change - 0.60364 #kg Co2 / kWh
#FFC - 9.11521 MJ / kWh
#water use - 0.10916 m3 / kWh

#GWP table
products = [process.PE_resin]
GWP_table = bst.report.lca_displacement_allocation_table(
    [process.system],
    'GWP',
    products, # For including products without characterization factors
)


# #GWP bar chart
def gwp_barchart():
    production_methods = ("Virgin PE", "STRAP PE")
    
    # Initialize the weight_counts dictionary
    weight_counts = {"Virgin resin": np.array([2.091, 0])}  #add virgin resin GWP here --> EF 2.0: LDPE granulates, production mix, at plant, Polymerisation of ethylene, 0.91- 0.96 g/cm3, 28 g/mol per repeating unit
    
    # Populate weight_counts using the second value of the index and column 2 from GWP_table
    for factor_name, row in GWP_table.iterrows():
        # Check if the index is a tuple and exclude rows where the first entry starts with "Total"
        if isinstance(factor_name, tuple):
            if factor_name[0].startswith("Total"):
                continue  # Skip rows where the first entry of the tuple starts with "Total"
            factor_name = factor_name[1]  # Get the second element of the tuple
        else:
            factor_name = str(factor_name)  # Handle non-tuple cases (if any)
    
        value = row.iloc[1]  # Column 2: corresponding value
        # Explicitly check if value is zero
        if not isinstance(value, (int, float)) or np.isnan(value) or float(value) == 0.0:
            continue  # Skip rows with invalid or zero values
        weight_counts[factor_name] = np.array([0, value])
    
    # Sort weight_counts by total GWP contribution in descending order
    weight_counts = dict(sorted(weight_counts.items(), key=lambda x: x[1].sum(), reverse=True))
    
    # Set the bar width
    width = 0.5
    
    # Create the stacked bar chart
    fig, ax = plt.subplots()
    bottom = np.zeros(2)
    
    for label, weight_count in weight_counts.items():
        p = ax.bar(production_methods, weight_count, width, label=label, bottom=bottom)
        bottom += weight_count
    
    # Customize the plot
    ax.set_title("Global Warming Potential of PE Production")
    ax.set_ylabel("GWP (kg CO2 eq. / kg PE)")
    ax.legend(loc="upper right")
    plt.show()


#FFC table
products = [process.PE_resin]
FFC_table = bst.report.lca_displacement_allocation_table(
    [process.system],
    'FFC',
    products, # For including products without characterization factors
)

#FFC bar chart
def ffc_barchart():
    production_methods = ("Virgin PE", "STRAP PE")
    
    # Initialize the weight_counts dictionary
    weight_counts = {"Virgin resin": np.array([75.77, 0])}  # add virgin resin FFC here
    
    # Populate weight_counts using the second value of the index and column 2 from FFC_table
    for factor_name, row in FFC_table.iterrows():
        # Check if the index is a tuple and exclude rows where the first entry starts with "Total"
        if isinstance(factor_name, tuple):
            if factor_name[0].startswith("Total"):
                continue  # Skip rows where the first entry of the tuple starts with "Total"
            factor_name = factor_name[1]  # Get the second element of the tuple
        else:
            factor_name = str(factor_name)  # Handle non-tuple cases (if any)
    
        value = row.iloc[1]  # Column 2: corresponding value
        # Explicitly check if value is zero
        if not isinstance(value, (int, float)) or np.isnan(value) or float(value) == 0.0:
            continue  # Skip rows with invalid or zero values
        weight_counts[factor_name] = np.array([0, value])
    
    # Sort weight_counts by total FFC contribution in descending order
    weight_counts = dict(sorted(weight_counts.items(), key=lambda x: x[1].sum(), reverse=True))
    
    # Set the bar width
    width = 0.5
    
    # Create the stacked bar chart
    fig, ax = plt.subplots()
    ax.set_ylabel("FFC (MJ / kg PE)")
    bottom = np.zeros(2)
    
    for label, weight_count in weight_counts.items():
        p = ax.bar(production_methods, weight_count, width, label=label, bottom=bottom)
        bottom += weight_count
    
    # Customize the plot
    ax.set_title("FFC of PE Production")
    ax.legend(loc="upper right")
    plt.show()


#MSP single point sensitivity analysis 
def msp_tornado():
    model = process.model
    remove = {'Temperature', 'Feedstock distance', }
    parameters = [i for i in model.parameters if i.name not in remove]
    model.parameters = parameters
    baseline, lower, upper = model.single_point_sensitivity()
    metric_index = process.MSP.index
    index = [i.describe(distribution=False) # Instead of displaying distribution, it displays lower, baseline, and upper values
             for i in process.model.parameters]
    bst.plots.plot_single_point_sensitivity(baseline[metric_index],
                                            lower[metric_index],
                                            upper[metric_index],
                                            name='Cost to produce STRAP resin ($/kg PE)',
                                            w=1.0,
                                            index=index)


def get_FFC(x):
    #process.set_dissolution_temperature.baseline = dt
    process.precipitation_step.T = x
    process.system.simulate()
    print('FFC for ', process.dissolution_step.T, process.FFC())

# GWP single point sensitivity analysis 
def gwp_tornado():
    model = process.model
    remove = {'Price', 'IRR', 'Processing capacity'}
    parameters = [i for i in model.parameters if i.name not in remove]
    model.parameters = parameters
    baseline, lower, upper = model.single_point_sensitivity()
    metric_index = process.GWP.index
    index = [i.describe(distribution=False) # Instead of displaying distribution, it displays lower, baseline, and upper values
             for i in process.model.parameters]
    bst.plots.plot_single_point_sensitivity(baseline[metric_index],
                                            lower[metric_index],
                                            upper[metric_index],
                                            name='GWP [kg CO2 eq./kg PE]',
                                            w=1.0,
                                            index=index)

# FFC single point sensitivity analysis 
def ffc_tornado():
    from random import shuffle
    model = process.model
    remove = {'Price', 'IRR', 'Processing capacity'}
    parameters = [i for i in model.parameters if i.name not in remove]
    #shuffle(parameters)
    model.parameters = parameters
    baseline, lower, upper = model.single_point_sensitivity()
    metric_index = process.FFC.index
    index = [i.describe(distribution=False) # Instead of displaying distribution, it displays lower, baseline, and upper values
             for i in process.model.parameters]
    bst.plots.plot_single_point_sensitivity(baseline[metric_index],
                                            lower[metric_index],
                                            upper[metric_index],
                                            name='FFC [MJ/kg PE]',
                                            index=index)

#MSP spearman's rank sensitivity
# def msp_spearman():
#     process.model.evaluate()
#     df_rho, df_p = process.model.spearman_r()
#     bst.plots.plot_spearman_1d(df_rho['-', 'MSP'],index=[i.describe() for i in process.model.parameters], name='cost to produce STRAP PE')

def msp_spearman():
    # N_samples = 100  # Adjust the number of samples as needed
    # rule = 'L'  # Latin-Hypercube sampling
    # np.random.seed(1234)  # Ensure reproducibility
    
    # samples = process.model.sample(N_samples, rule)  # Generate samples
    # process.model.load_samples(samples)  # Load samples into the model
    # process.model.evaluate()  # Evaluate model over the samples
    process.system.simulate()
    df_rho, df_p = process.model.spearman_r()  # Perform Spearman analysis
    bst.plots.plot_spearman_1d(df_rho['-', 'MSP [USD/kg]'], 
                               index=[i.describe() for i in process.model.parameters], 
                               name='cost to produce STRAP PE')