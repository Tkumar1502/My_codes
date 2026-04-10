#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  1 13:32:58 2025

@author: charlesgranger
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 21 15:20:14 2025

@author: charlesgranger
"""
############### CHANGES TO BE MADE #####################################

#9) spearman's rank restets the parameters. is that because i am only defining bound 
#but not distributions, so it defaults to the original definition?
########################################################################
import os
from biorefineries.tea.cellulosic_ethanol_tea import foc_table, capex_table
from thermosteam.utils import array_roundsigfigs
import biosteam as bst
import pandas as pd
from plastics import strap
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.patches as mpatches
from warnings import filterwarnings
from datetime import datetime
from chaospy.distributions import Uniform
from matplotlib.colors import TwoSlopeNorm, Normalize
from matplotlib.ticker import PercentFormatter
filterwarnings('ignore')

#define scenario
processing_capacity=20e3 #tons / year feedstock

#single wash
process = strap.BaselineSTRAPProcess(
    scenario='PE/Xylene',
    target_plastic_percent=0.90, # %
    processing_capacity=processing_capacity, # MT/yr
    sell_leftover_plastic=True, #undissolved fraction
    burn_leftover_plastic=False, #waste to energy
    facilities=False,
    simulate=False, #don't change this!
    
)

#adjust tea parameters
process.tea.operating_hours = 7862.4 #50weeks, 2 shifts, 5 days a week 
process.tea.labor_cost = 2*420000 # 1 operators x $60K, 1 loader x $50k, 1 engineer X $100k per shift (X2)

#create flowchart
process.system.diagram(format='png') # Note how the leftover plastic is sent to heat and power generation.

#define products and set sale prices
products = [process.PE_resin, process.leftover_plastic]
process.PE_resin.price = 1.20 # $1.20 / kg
process.leftover_plastic.price = 0#8.00 # $8.00 / kg --> set to zero if not selling


#bounds 1-wash
process.set_polymer_mass_fraction.bounds=(.80, .95)
process.set_dissolution_capacity.bounds=(5, 10)
process.set_solvent_loss.bounds=(0.01,1)
process.set_precipitation_temperature.bounds=(313,323)
process.set_dissolution_temperature.bounds=(363,393)
process.set_centrifuged_plastic_solvent_content.bounds=(25, 75)
process.set_feedstock_distance.bounds=(50,1000)
process.set_feedstock_price.bounds=(0.00,0.08)

#set baseline 1-wash
process.plastic.ID='Feedstock_plastic'
process.set_dissolution_temperature.baseline = 373
process.set_solvent_loss.baseline = 0.1
process.set_polymer_mass_fraction.baseline = 0.9
process.set_dissolution_capacity.baseline = 7.5
process.set_precipitation_temperature.baseline = 45+273
process.set_feedstock_price.baseline=0.02
process.set_IRR.baseline=0.15

for i in process.parameters:
    i.distribution = Uniform(*i.bounds)

############################## LCA #####################################

######## SET CHARACTERIZATION FACTORS #########################
GWP = 'GWP'
FFC = 'FFC'
WU = 'WU'

#natural gas CF

process.natural_gas.set_CF(GWP, 0.61146) # kg CO2 eq / kg NG

process.natural_gas.set_CF(
    FFC,
    51, # 51MJ/kg from EF 2.0
)

process.natural_gas.set_CF(WU,1000)

#add CFs for solvent
process.solvent.set_CF(
    FFC,
    52.33 + 2.260, #EF 2.0; Xylene production, production mix, at plant, technology mix, 100% active substance
) # solvent vapor trap incineration: EF 2.0 Waste incineration of inert material, production mix, at consumer, waste-to-energy plant with dry flue gas treatment, including transport and pre-treatment, inert material waste

process.solvent.set_CF(
    GWP,
    0.9383 + 0.1563, # EF 2.0; Xylene production, production mix, at plant, technology mix, 100% active substance
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
process.adsorption_column.ins[2].set_CF(GWP, 1.78+0.02643,  ) #EF 2.0 activated silica production, production mix, at plant, technology mix, 100% active substance (activated carbon not available)
process.adsorption_column.ins[2].set_CF(FFC, 23.99+0.3452,  ) 
# adsorbent landfilling: EF 2.0 Landfill of polluted inorganic waste, production mix (region specific sites), at landfill site, landfill including leachate treatment and with transport without collection and pre-treatment


# ELECTRICITY / MJ "Electricity grid mix 1kV-60kV , consumption mix, to consumer, AC, technology mix, 1kV - 60kV"
#climate change - 0.60364 #kg Co2 / kWh
#FFC - 9.11521 MJ / kWh
#water use - 0.10916 m3 / kWh





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
lb, ub = 5000, 100000
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

#calculate dROI for given capacity
def dROI(pc):
    process.set_processing_capacity(pc)
    process.system.simulate()
    droi = process.tea.NPV / process.tea.TCI
    return droi

def cap_by_price(lb, ub, low_price, high_price):
    #define x y bounds
    capacities = np.linspace(lb, ub, 20)
    prices = np.linspace(low_price,high_price,20)
    
    #store dROI output
    # Create a 2D array to hold the function values
    IRR = np.zeros((len(prices), len(capacities)))
    
    for i, cap in enumerate(capacities): 
        for j, price in enumerate(prices):
            process.set_processing_capacity(cap)
            process.PE_resin.price = price
            process.system.simulate()
            IRR[i,j] =  process.tea.solve_IRR()
    
            
    # Create the normalization: center at 0
    norm = TwoSlopeNorm(vmin=np.min(IRR), vcenter=0, vmax=np.max(IRR))
    
    # Create the heat map
    img = plt.imshow(IRR, 
               extent=(capacities[0], capacities[-1], prices[0], prices[-1]), 
               origin='lower', 
               aspect='auto', 
               cmap='RdYlGn',  # Red (low) -> Yellow (0) -> Green (high)
               norm=norm)
    
    # Add and format the colorbar
    cbar = plt.colorbar(img)
    cbar.set_label('IRR (%)')
    cbar.ax.yaxis.set_major_formatter(PercentFormatter(xmax=1))
    
    plt.xlabel('Feedstock capacity')
    plt.ylabel('STRAP resin sale price ($/kg)')
    plt.title('Pallet Wrap')
    plt.show()

def irr_contour(lb, ub, low_price, high_price):
    # Define x and y bounds
    capacities = np.linspace(lb, ub, 20)
    prices = np.linspace(low_price, high_price, 20)

    # Create a 2D array to hold the IRR values
    IRR = np.zeros((len(prices), len(capacities)))

    for i, cap in enumerate(capacities): 
        for j, price in enumerate(prices):
            process.set_processing_capacity(cap)
            process.PE_resin.price = price
            process.system.simulate()
            IRR[j, i] = process.tea.solve_IRR()

    # Create meshgrid
    X, Y = np.meshgrid(capacities, prices)

    # Set figure size correctly using subplots
    fig, ax = plt.subplots(figsize=(6, 6))

    # Generate contour levels
    min_val = np.min(IRR)
    max_val = np.max(IRR)
    neg_levels = np.arange(np.floor(min_val * 100) / 100, 0, 1.0)
    pos_levels = np.arange(0.05, np.ceil(max_val * 100) / 100 + 0.001, 0.05)
    levels = np.concatenate((neg_levels, [0], pos_levels))

    norm = TwoSlopeNorm(vmin=min_val, vcenter=0, vmax=max_val)

    # Create filled contour plot
    contourf = ax.contourf(X, Y, IRR, levels=levels, cmap='RdYlGn', norm=norm)

    # Add contour lines and labels (only for positive levels)
    contour_lines = ax.contour(X, Y, IRR, levels=pos_levels, colors='black', linewidths=0.5)
    ax.clabel(contour_lines, fmt=lambda x: f"{x:.0%}", fontsize=8)

    # Add colorbar
    # cbar = fig.colorbar(contourf, ax=ax)
    # cbar.set_label('IRR (%)')
    # cbar.ax.yaxis.set_major_formatter(PercentFormatter(xmax=1))

    # Labels and title
    ax.set_xlabel('Feedstock capacity (metric tpy)')
    ax.set_ylabel('STRAP resin sale price ($/kg)')
    ax.set_title('Printed Film - IRR')

    # Add faint grid
    ax.grid(True, linestyle='--', color='gray', linewidth=0.5, alpha=0.4)

    plt.show()
    
def sales_minus_aoc_contour(lb, ub, low_price, high_price):
    # Define x and y bounds
    capacities = np.linspace(lb, ub, 20)
    prices = np.linspace(low_price, high_price, 20)

    # Create a 2D array to hold the Sales - AOC values
    sales_aoc = np.zeros((len(prices), len(capacities)))

    for i, cap in enumerate(capacities): 
        for j, price in enumerate(prices):
            process.set_processing_capacity(cap)
            process.PE_resin.price = price
            process.system.simulate()
            sales_aoc[j, i] = process.tea.sales - process.tea.AOC

    # Create meshgrid
    X, Y = np.meshgrid(capacities, prices)

    # Set figure size
    fig, ax = plt.subplots(figsize=(6, 6))

    # Determine value range
    min_val = np.min(sales_aoc)
    max_val = np.max(sales_aoc)

    # Define contour levels using linspace (limited to 15 levels)
    levels = np.linspace(min_val, max_val, 15)

    # Choose normalization based on range
    if min_val < 0 and max_val > 0:
        norm = TwoSlopeNorm(vmin=min_val, vcenter=0, vmax=max_val)
    else:
        norm = Normalize(vmin=min_val, vmax=max_val)

    # Create filled contour plot
    contourf = ax.contourf(X, Y, sales_aoc, levels=levels, cmap='RdYlGn', norm=norm)

    # Add contour lines and labels (for all levels)
    contour_lines = ax.contour(X, Y, sales_aoc, levels=levels, colors='black', linewidths=0.5)
    ax.clabel(contour_lines, fmt=lambda x: f"{x:.3g}", fontsize=8)

    # Labels and title
    ax.set_xlabel('Feedstock capacity (metric tpy)')
    ax.set_ylabel('STRAP resin sale price ($/kg)')
    ax.set_title('Cytiva PIW: Annual sales-OPEX')

    # Optional colorbar (uncomment if needed)
    # cbar = fig.colorbar(contourf, ax=ax)
    # cbar.set_label('Sales / AOC')

    # Add grid
    ax.grid(True, linestyle='--', color='gray', linewidth=0.5, alpha=0.4)

    plt.show()
    

def revenue(feed_price):
    process.set_feedstock_price(feed_price)
    process.system.simulate()
    sales = process.tea.sales
    opex = process.tea.AOC
    return sales - opex
    

def revenue_vs_feedstock(lb, ub):
    feed_price = np.linspace(lb, ub)
    rev = [revenue(i) for i in feed_price]
    
    plt.plot(feed_price, rev)
    plt.show()
    


#get data
def msp_plot():
    MSPs = [MSP(i) for i in pcs]
    TCIs = [TCI(i) for i in pcs]
    market = [1.20 for i in pcs]
    
   # Create a DataFrame with PCS, MSPs, and TCIs
    df = pd.DataFrame({'pcs': pcs, 'MSP': MSPs, 'TCI': TCIs})
    
   #plot data
    plt.figure()
    line1, = plt.plot(pcs, MSPs, label='STRAP resin cost') # plot MSP for different capacities
    #line2, = plt.plot(pcs, market, 'r--', label='Market Value PE') #market value
    plt.axhspan(1.00, 1.40, color='blue', alpha=0.15, label='Market Value Range')
    plt.ylabel('Cost to produce STRAP resin ($/kg PE)')
    plt.xlabel('Plant Capacity (mton/year)')
    
    # Create a secondary y-axis
    ax2 = plt.gca().twinx()
    line3, = ax2.plot(pcs, TCIs, label='TCI', color='green')
    ax2.set_ylabel('Total Capital Investment ($)')  # Label for the second y-axis
    ax2.ticklabel_format(axis='y', style='sci', scilimits=(6, 6)) 
    
    # Combine legends for both axes
    lines = [line1, line3]
    labels = [line.get_label() for line in lines]
    plt.legend(lines + [plt.Line2D([0], [0], color='blue', alpha=0.15, linewidth=10)], 
               labels + ['Vrigin PE - Market Range'], 
               loc='upper center')
    plt.show()
    return df

#returns diescounted return on investment vs plant capacity
def droi_plot():
    dROIs = [dROI(i) for i in pcs]
    TCIs = [TCI(i) for i in pcs]
    #zero = [0 for i in pcs]
    
    #plot data
    plt.figure()
    line1, = plt.plot(pcs, dROIs, label='discounted ROI') # plot MSP for different capacities
    #line2, = plt.plot(pcs, zero, 'r--', label=False) #market value
    plt.ylabel('discounted ROI (NPV/TCI)x100')
    plt.xlabel('Plant Capacity (mton/year)')
    
    # Create a secondary y-axis
    ax2 = plt.gca().twinx()
    line3, = ax2.plot(pcs, TCIs, label='TCI', color='green')
    ax2.set_ylabel('Total Capital Investment ($)')  # Label for the second y-axis
    
    # Combine legends for both axes
    lines = [line1, line3]
    labels = [line.get_label() for line in lines]
    plt.legend(lines, labels, loc='upper center')    

    # Save data to Excel
    output_folder = './results/cytiva/'
    os.makedirs(output_folder, exist_ok=True)
    df = pd.DataFrame({
        'Plant Capacity (mton/year)': pcs,
        'Discounted ROI': dROIs,
        'Total Capital Investment ($)': TCIs
    })
    df.to_excel(os.path.join(output_folder, 'droi_results_scenario2.xlsx'), index=False)
    

#make a plot of msp vs capacity with confidence bands
def msp_confidence():
    TCIs = [TCI(i) for i in pcs]
    #market = [1.17 for i in pcs]
        
    # #run monte carlo
    N_samples = 100
    rule = 'L' # For Latin-Hypercube sampling
    np.random.seed(1234) # For consistent results
    samples = process.model.sample(N_samples, rule)
    process.model.load_samples(samples)
     
    # #evaluate across processing capacity
    now = datetime.now()
    filename = now.strftime("MC")
    process.model.evaluate_across_coordinate('Processing capacity', 
                                             process.set_processing_capacity, 
                                             pcs, 
                                             xlfile=f'./{filename}.xlsx'
                                             )
    
    #plot across capacities
    df = pd.read_excel(f'./{filename}.xlsx', sheet_name='- MSP', index_col=0)
    plt.figure(figsize=(8, 6))
    percentiles = bst.plots.plots.plot_montecarlo_across_coordinate(pcs, df)
   
   # Extract percentiles
    p5, p25, p50, p75, p95 = percentiles  # 5th, 25th, 50th (median), 75th, 95th percentiles
    
    print(p50)
    
   # Add legend manually for Monte Carlo results
    line_median, = plt.plot(pcs, p50, 'k-', label='Median Cost')  # Center line
    fill = plt.fill_between(pcs, p25, p75, color='brown', alpha=0.3, label='Interquartile Range (25th-75th)')  # Shaded region
    line_ci, = plt.plot(pcs, p5, 'k--', label='Confidence Interval (5th-95th)')  # Lower bound
    plt.plot(pcs, p95, 'k--')  # Upper bound (same label as lower)
    
    # Add other plots
    #line2, = plt.plot(pcs, market, 'r--', label='Market Value PE')  # Market value
    plt.axhspan(1.00, 1.40, color='blue', alpha=0.15, label='Market Value Range')
    plt.ylabel('Cost to produce STRAP resin ($/kg)')
    plt.xlabel('Plant Feedstock Capacity (metric ton/year)')
    
    # Create a secondary y-axis
    ax2 = plt.gca().twinx()
    line3, = ax2.plot(pcs, TCIs, label='TCI', color='green')
    ax2.set_ylabel('Total Capital Investment ($)')
    ax2.ticklabel_format(axis='y', style='sci', scilimits=(6, 6))  # Forces 10^6 scale

    # Combine legends
    lines = [line_median, fill, line_ci, line3]
    labels = [line.get_label() for line in lines]
    plt.legend(lines + [plt.Line2D([0], [0], color='blue', alpha=0.15, linewidth=10)], 
               labels + ['Vrigin PE - Market Range'], 
               loc='upper center')

############### tornado plot ############################################



#GWP tables
GWP_table = bst.report.lca_displacement_allocation_table(
    [process.system],
    'GWP',
    products, 
)

inventory_table = bst.report.lca_inventory_table(
    [process.system], 
    'GWP', 
    products,)


# #GWP bar chart
def gwp_barchart(label_offsets=(-0.15, 0.05)):
    production_methods = ("Virgin PE", "STRAP PE")
    
    weight_counts = {"Virgin resin": np.array([2.091, 0])}
    
    total_gwp = 0
    for factor_name, row in GWP_table.iterrows():
        if isinstance(factor_name, tuple):
            if factor_name[0].startswith("Total"):
                continue
            factor_name = factor_name[1]
        else:
            factor_name = str(factor_name)
    
        value = row.iloc[1]
        if not isinstance(value, (int, float)) or np.isnan(value) or float(value) == 0.0:
            continue
        
        weight_counts[factor_name] = np.array([0, value])
        total_gwp += value
    
    weight_counts = {k: v for k, v in weight_counts.items() if v.sum() / total_gwp > 0.005}
    weight_counts = dict(sorted(weight_counts.items(), key=lambda x: x[1].sum(), reverse=True))
    
    width = 0.5
    fig, ax = plt.subplots()
    bottom = np.zeros(2)
    
    for label, weight_count in weight_counts.items():
        p = ax.bar(production_methods, weight_count, width, label=label, bottom=bottom)
        bottom += weight_count

    for i, total in enumerate(bottom):
        offset = label_offsets[i] if i < len(label_offsets) else 0.05
        ax.text(
            i, total + offset, f"{total:.2f}",
            ha='center', va='bottom',
            bbox=dict(facecolor='white', edgecolor='none', pad=1.0)
        )

    ax.set_title("Global Warming Potential of PE Production")
    ax.set_ylabel("GWP (kg CO2 eq. / kg PE)")
    ax.legend(loc="best")
    plt.show()



#FFC table
FFC_table = bst.report.lca_displacement_allocation_table(
    [process.system],
    'FFC',
    products, 
)

inventory_table = bst.report.lca_inventory_table(
    [process.system], 
    'FFC', 
    products,)

# LCA tables
def save_lca_tables_to_excel():
    # Generate the required tables
    GWP_table = bst.report.lca_displacement_allocation_table(
        [process.system], 'GWP', products
    )
    inventory_GWP_table = bst.report.lca_inventory_table(
        [process.system], 'GWP', products
    )
    
    FFC_table = bst.report.lca_displacement_allocation_table(
        [process.system], 'FFC', products
    )
    inventory_FFC_table = bst.report.lca_inventory_table(
        [process.system], 'FFC', products
    )
    
    # Save tables to an Excel file with two sheets
    filename="./results/cytiva_12tpm.xlsx"
    with pd.ExcelWriter(filename) as writer:
        GWP_table.to_excel(writer, sheet_name="GWP Tables")
        inventory_GWP_table.to_excel(writer, sheet_name="GWP Tables", startrow=GWP_table.shape[0] + 2)
        FFC_table.to_excel(writer, sheet_name="FFC Tables")
        inventory_FFC_table.to_excel(writer, sheet_name="FFC Tables", startrow=FFC_table.shape[0] + 2)
    
    print(f"LCA tables saved to {filename}")

#FFC bar chart
def ffc_barchart(label_offsets=(-5, 0.5)):
    production_methods = ("Virgin PE", "STRAP PE")
    
    # Initialize the weight_counts dictionary
    weight_counts = {"Virgin resin": np.array([75.77, 0])}  # add virgin resin FFC here

    # Populate weight_counts using the second value of the index and column 2 from FFC_table
    for factor_name, row in FFC_table.iterrows():
        if isinstance(factor_name, tuple):
            if factor_name[0].startswith("Total"):
                continue  # Skip rows where the first entry of the tuple starts with "Total"
            factor_name = factor_name[1]  # Get the second element of the tuple
        else:
            factor_name = str(factor_name)  # Handle non-tuple cases (if any)
    
        value = row.iloc[1]  # Column 2: corresponding value
        if not isinstance(value, (int, float)) or np.isnan(value) or float(value) == 0.0:
            continue  # Skip rows with invalid or zero values
        weight_counts[factor_name] = np.array([0, value])

    # Filter out components that contribute less than 1% to total FFC
    total_ffc = sum(val.sum() for val in weight_counts.values())
    weight_counts = {k: v for k, v in weight_counts.items() if v.sum() / total_ffc >= 0.0001}

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
        
    for i, total in enumerate(bottom):
        offset = label_offsets[i] if i < len(label_offsets) else 0.05
        ax.text(
            i, total + offset, f"{total:.2f}",
            ha='center', va='bottom',
            bbox=dict(facecolor='white', edgecolor='none', pad=1.0)
        )

    # Customize the plot
    ax.set_title("FFC of PE Production")
    ax.legend(loc="upper right")
    plt.show()


#MSP single point sensitivity analysis 
def msp_tornado():
    model = process.model
    remove = {'Temperature', 'Distance', 'Solvent content' }
    parameters = [i for i in model.parameters if i.name not in remove]
    model.parameters = parameters
    baseline, lower, upper = model.single_point_sensitivity()
    metric_index = process.MSP.index
    index = [i.describe(distribution=False) # Instead of displaying distribution, it displays lower, baseline, and upper values
             for i in process.model.parameters]
    bst.plots.plot_single_point_sensitivity(baseline[metric_index],
                                            lower[metric_index],
                                            upper[metric_index],
                                            name='Cost to produce STRAP resin ($/kg product)',
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
                                            name='GWP [kg CO2 eq./kg product]',
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
                                            name='FFC [MJ/kg product]',
                                            w=1.0,
                                            index=index)

#MSP spearman's rank sensitivity
# def msp_spearman():
#     process.model.evaluate()
#     df_rho, df_p = process.model.spearman_r()
#     bst.plots.plot_spearman_1d(df_rho['-', 'MSP'],index=[i.describe() for i in process.model.parameters], name='cost to produce STRAP PE')


#TO DO: FIX THIS!
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
    
############ REPORT TABLES #########################################
def save_system_reports():
    folder = './results/cytiva/'
    #folder = os.path.join(folder, 'results')
    for scenario in ('experimental', 'NREL'):
        pm = process
        filename = f'cytiva_12tpm_detailed_report.xlsx'
        file = os.path.join(folder, filename)
        pm.system.save_report(file)
        
def save_detailed_expenditure_tables(sigfigs=3):
    folder = './results/cytiva/'
    #folder = os.path.join(folder, 'results')
    filename = '12tpmexpenditures.xlsx'
    file = os.path.join(folder, filename)
    writer = pd.ExcelWriter(file)
    product = 'PE_resin'
    process_models = [
        process
    ]
    for pm in process_models: pm.PE_resin.price = pm.tea.solve_price(pm.PE_resin)
    systems = [i.system for i in process_models]
    names = ['baseline']
    teas = [i.tea for i in process_models]
    tables = {
        'VOC': bst.report.voc_table(systems, product, names, with_products=True),
        'FOC': foc_table(teas, names),
        'CAPEX': capex_table(teas, names),
    }
    for key, table in tables.items():
        values = array_roundsigfigs(table.values, sigfigs=3, inplace=True)
        if key == 'CAPEX': # Bug in pandas
            for i, col in enumerate(table):
                table[col] = values[:, i]
        table.to_excel(writer, key)
    writer.close()
    return tables


#FIX THIS TABLE!
def save_detailed_life_cycle_tables(sigfigs=3):
    # Define the process model and its resin ID
    process_model = process
    
    # Specify system and file details
    system = process_model.system
    folder = './results/cytiva'
    filename = 'life_cyclecytiva_12tpm.xlsx'
    file = os.path.join(folder, filename)
    
    # Create Excel writer
    writer = pd.ExcelWriter(file)
    
    # Define streams and system names
    streams = [getattr(process_model, 'PE_resin')]
    names = ['Baseline']
    
    # Generate tables
    tables = {
        'Inventory': bst.report.lca_inventory_table(
            [system], 'GWP', streams, system_names=names
        ),
        'Energy allocation factors': bst.report.lca_property_allocation_factor_table(
            [system], property='energy', basis='GGE', system_names=names
        ),
    }
    
    # Environmental impact table
    values = np.array([[process_model.GWP_electricity(),
                        process_model.GWP_PE_resin()]])
    index = ['Electricity [kg∙CO2e∙kWh-1]', 'Polymer resin [kg∙CO2e∙kg-1]']
    columns = ['Baseline']
    df_gwp = pd.DataFrame(values, index=index, columns=columns)
    tables['Estimated environmental impact'] = df_gwp

    # Save tables to Excel
    for key, table in tables.items():
        array_roundsigfigs(table.values, sigfigs=3, inplace=True)
        table.to_excel(writer, key)

    writer.close()
    return tables
