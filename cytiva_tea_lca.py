#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 14 11:05:08 2025

@author: charlesgranger
"""
############### CHANGES TO BE MADE #####################################

#currently "solvent loss" is being accounted for in allocation methods

# fix the detailed LCA table

#can we move centrifilter before surge tank and pump?
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
from matplotlib.ticker import PercentFormatter, FuncFormatter
filterwarnings('ignore')

# Added for inline SVG display in interactive environments
try:
    from IPython.display import SVG, display
except ImportError:
    # IPython not available; define dummy functions
    def display(obj):
        print('Display not available in this environment.')
    class SVG:
        def __init__(self, filename=None):
            self.filename = filename

def show_svg(path):
    """Display an SVG file inline if possible.

    Parameters
    ----------
    path : str or pathlib.Path
        Path to the SVG file.
    """
    display(SVG(filename=str(path)))

#define scenario
processing_capacity=325

#single wash
process = strap.BaselineSTRAPProcess(
    scenario='PE/Xylene',
    target_plastic_percent=0.9, # %
    processing_capacity=processing_capacity, # MT/yr
    sell_leftover_plastic=True,
    burn_leftover_plastic=False,
    facilities=False,
    simulate=False,
)

#adjust operating days
process.tea.operating_days = 328.5

#create flowchart
process.system.diagram(format='png') # Note how the leftover plastic is sent to heat and power generation.

#define products and set sale prices
products = [process.PE_resin, process.leftover_plastic]
process.PE_resin.price = 1.20 #+ 0.072 # $1.20 / kg
process.leftover_plastic.price = 100#-0.072#8.00 # $8.00 / kg

#bounds 1-wash
process.set_polymer_mass_fraction.bounds=(.6, .95)
process.set_dissolution_capacity.bounds=(5, 10)
process.set_solvent_loss.bounds=(0.01,1)
process.set_precipitation_temperature.bounds=(313,323)
process.set_dissolution_temperature.bounds=(363,393)
process.set_centrifuged_plastic_solvent_content.bounds=(25, 75)
process.set_feedstock_distance.bounds=(50,1000)
process.set_feedstock_price.bounds=(0.072,0.172)

#set baseline 1-wash
process.plastic.ID='Feedstock_plastic'
process.set_dissolution_temperature.baseline = 373
process.set_solvent_loss.baseline = 0.1
process.set_polymer_mass_fraction.baseline = 0.14
process.set_dissolution_capacity.baseline = 5.0
process.set_precipitation_temperature.baseline = 45+273
process.set_feedstock_price.baseline=(0.05)
process.set_IRR.baseline=0.15

for i in process.parameters:                                        #didn't get why it is here
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

#calculate unit operating cost for given capacity
def UOC(pc):
    process.set_processing_capacity(pc)
    process.system.simulate()
    aoc = process.tea.AOC 
    landfill = process.leftover_plastic.price*process.leftover_plastic.F_mass*process.tea.operating_hours
    return (aoc-landfill)/ (process.PE_resin.F_mass*process.tea.operating_hours)

def save_UOC():
    UOCs = [UOC(i) for i in pcs]
    
    # Create a DataFrame with PCS, UOCs
    df = pd.DataFrame({'pcs': pcs, 'UOC': UOCs})
    
    return df
    
        
    

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
    return df

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
    plt.figure(figsize=(10, 4))
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

#plot of IRR vs PE price and EVOH price at different capacities


def slice_irr():
    # Define your variable ranges
    scales = np.array([2500, 5000, 10000])        
    pe_prices = np.linspace(0.62, 2.5, 20)   
    evoh_prices = np.linspace(6, 9, 20)   
    
    # Initialize empty array for IRR results
    irr_values = np.empty((len(scales), len(pe_prices), len(evoh_prices)))
    
    # Loop and compute IRRs
    for i, scale in enumerate(scales):
        for j, p1 in enumerate(pe_prices):
            for k, p2 in enumerate(evoh_prices):
                try:
                    process.set_processing_capacity(scale)
                    process.PE_resin.price = p1 
                    process.leftover_plastic.price = p2 
                    process.system.simulate()
                    irr = process.tea.solve_IRR()
                except Exception:
                    irr = np.nan
                irr_values[i, j, k] = irr

    # Save results
    np.savez_compressed("./results/cytiva/PE_only/rr_grid.npz",
                        irr_values=irr_values,
                        scales=scales,
                        pe_prices=pe_prices,
                        evoh_prices=evoh_prices)
    
    records = []
    for i, scale in enumerate(scales):
        for j, p1 in enumerate(pe_prices):
            for k, p2 in enumerate(evoh_prices):
                records.append({
                    "Feedstock capacity (metric ton / year)": scale,
                    "PE price ($/kg) ": p1,
                    "EVOH price ($/kg)": p2,
                    "IRR": irr_values[i, j, k]
                })
    df = pd.DataFrame(records)
    df.to_excel("./results/cytiva/PE_only/irr_grid_data.xlsx", index=False)

    # ---- PLOTTING SECTION ----
    X, Y = np.meshgrid(pe_prices, evoh_prices)
    irr_percent = irr_values * 100  # Convert to percent
    levels = np.arange(0, 70, 5)    # 0% to 65% in 5% steps

    fig, axes = plt.subplots(3, 1, figsize=(6, 16), constrained_layout=True)

    # Compute consistent bounds for all axes
    x_min, x_max = pe_prices.min(), pe_prices.max()
    y_min, y_max = evoh_prices.min(), evoh_prices.max()
    x_range = x_max - x_min
    y_range = y_max - y_min
    aspect_ratio = x_range / y_range  # Match plot scaling

    for i, ax in enumerate(axes):
        Z = irr_percent[i].T  # Use percent values
        cf = ax.contourf(X, Y, Z, levels=levels, cmap='RdYlGn', alpha=0.9)
        cl = ax.contour(X, Y, Z, levels=levels, colors='k', linewidths=0.7)
        ax.clabel(cl, fmt='%.0f%%', fontsize=14)  # Format as percent

        ax.set_title(f'Feedstock capacity: {int(scales[i]):,} metric tpy', fontsize=14)
        ax.set_xlabel('STRAP PE Price ($/kg)', fontsize=14)
        ax.set_ylabel('STRAP EVOH Price ($/kg)', fontsize=14)
        ax.set_xlim(x_min, x_max)
        ax.set_ylim(y_min, y_max)
        ax.set_aspect(aspect_ratio)
        ax.tick_params(labelsize=14)

    # Custom colorbar placement (right side)
    cbar_ax = fig.add_axes([0.92, 0.15, 0.02, 0.7])
    cbar = fig.colorbar(cf, cax=cbar_ax)
    cbar.set_label('IRR', fontsize=14)
    cbar.ax.tick_params(labelsize=14)
    cbar.ax.yaxis.set_major_formatter(FuncFormatter(lambda x, _: f'{int(x)}%'))

    plt.suptitle('S2', fontsize=18, y=1.03)
    plt.savefig('./results/cytiva/PE_only/irr_contour_stack.pdf', bbox_inches='tight') 
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

    # Set figure size
    fig, ax = plt.subplots(figsize=(6, 5))

    # Generate contour levels
    min_val = np.nanmin(IRR)
    max_val = np.nanmax(IRR)

    neg_levels = np.arange(np.floor(min_val * 100) / 100, 0, 1.0) if min_val < 0 else np.array([])
    pos_levels = np.arange(0.05, np.ceil(max_val * 100) / 100 + 0.001, 0.05) if max_val > 0 else np.array([])
    levels = np.concatenate((neg_levels, [0] if min_val < 0 and max_val > 0 else [], pos_levels))

    # Choose appropriate normalization
    if min_val < 0 and max_val > 0:
        norm = TwoSlopeNorm(vmin=min_val, vcenter=0, vmax=max_val)
    else:
        norm = Normalize(vmin=min_val, vmax=max_val)

    # Create filled contour plot
    levels = np.arange(0, 0.65, 0.05)
    contourf = ax.contourf(X, Y, IRR, levels=levels, cmap='RdYlGn', alpha=0.9) #,norm=norm)

    # Add contour lines and labels (only for positive levels)
    if len(pos_levels) > 0:
        contour_lines = ax.contour(X, Y, IRR, levels=pos_levels, colors='black', linewidths=0.5)
        ax.clabel(contour_lines, fmt=lambda x: f"{x:.0%}", fontsize=8)

    # Add colorbar
    cbar = fig.colorbar(contourf, ax=ax)
    cbar.set_label('IRR (%)')
    cbar.ax.yaxis.set_major_formatter(PercentFormatter(xmax=1))

    # Labels and title
    ax.set_xlabel('Feedstock capacity (metric tpy)')
    ax.set_ylabel('STRAP PE sale price ($/kg)')
    ax.set_title('S1')

    # Add faint grid
    ax.grid(True, linestyle='--', color='gray', linewidth=0.5, alpha=0.4)
    
    plt.savefig('./results/cytiva/PE_only/irr_contour.pdf', bbox_inches='tight') 

    plt.show()
############### tornado plot ############################################

############################## LCA #####################################

######## SET CHARACTERIZATION FACTORS #########################
GWP = 'GWP'
FFC = 'FFC'
WU = 'WU'
HTC = 'HTC'
HTNC = 'HTNC'
ETOX = 'ETOX'
ACD = 'ACD'
OZD = 'OZD'
POCP = 'POCP'


#natural gas CF
process.natural_gas.set_CF(GWP, 0.61146) # kg CO2 eq / kg NG
process.natural_gas.set_CF(FFC,50.937, )# 51MJ/kg from EF 2.0
process.natural_gas.set_CF(WU,0.01633, )
process.natural_gas.set_CF(HTC,2.708e-10, )
process.natural_gas.set_CF(HTNC,3.5399e-9 )
process.natural_gas.set_CF(ETOX, 0.02579)
process.natural_gas.set_CF(ACD,0.00130)
process.natural_gas.set_CF(OZD, 2.8956e-11)
process.natural_gas.set_CF(POCP,.00124)                        

#add CFs for xylene
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

#CFs for water
process.makeup_water.set_CF(FFC,5.09/1000, )#MJ/kg
process.makeup_water.set_CF(GWP,0.33/1000, ) #kg co2 eq / kg water
process.makeup_water.set_CF(WU, 4.3e-2/1000) #m3/kg cooling water
process.makeup_water.set_CF(HTC, 3.6e-14)  #CTUh / kg cooling water
process.makeup_water.set_CF(HTNC, 1.4e-13) #CTUh / kg cooling water
process.makeup_water.set_CF(ETOX, 4.3e-7) #CTUe / kg cooling water
process.makeup_water.set_CF(ACD, 2.7918e-7) #
process.makeup_water.set_CF(OZD, 1.314e-12) #!!!!!!!
process.makeup_water.set_CF(POCP, 1.505e-7)#!!!!!!!

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

#CFs for adsorbent
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


#CFs for EVOH
# process.leftover_plastic.set_CF(GWP,7.30)
# process.leftover_plastic.set_CF(FFC,161.77565)
# process.leftover_plastic.set_CF(WU,1.32857)
# process.leftover_plastic.set_CF(HTC,2.52425e-8)
# process.leftover_plastic.set_CF(HTNC,4.23136e-7)
# process.leftover_plastic.set_CF(ETOX,1.06984)
# process.leftover_plastic.set_CF(ACD,.00868)
# process.leftover_plastic.set_CF(OZD,2.80388e-10)
# process.leftover_plastic.set_CF(POCP,.01236)




# ELECTRICITY / MJ "Electricity grid mix 1kV-60kV , consumption mix, to consumer, AC, technology mix, 1kV - 60kV"
#climate change - 0.60364 #kg Co2 / kWh
#FFC - 9.11521 MJ / kWh
#water use - 0.10916 m3 / kWh

#GWP tables
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
    items = products,
    key =   GWP
)

inventory_table = bst.report.lca_inventory_table(
    [process.system], 
    'GWP', 
    products,)


# #GWP bar chart
def gwp_barchart():
    production_methods = ("Virgin PE", "STRAP PE")
    
    # Initialize the weight_counts dictionary
    weight_counts = {"Virgin resin": np.array([2.091, 0])}  #add virgin resin GWP here --> EF 2.0: LDPE granulates, production mix, at plant, Polymerisation of ethylene, 0.91- 0.96 g/cm3, 28 g/mol per repeating unit
    
    # Populate weight_counts using the second value of the index and column 2 from GWP_table
    total_gwp = 0  # Track total GWP for percentage calculations
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
        total_gwp += value  # Add to total GWP
    
    # Filter out components contributing less than 1% of total GWP
    weight_counts = {k: v for k, v in weight_counts.items() if v.sum() / total_gwp > 0.001}
    
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
    # List of impact categories to export
    impact_categories = ['GWP', 'FFC', 'WU', 'HTC', 'HTNC', 'ETOX', 'ACD', 'OZD', 'POCP']
    
    # Excel output file path
    filename = "./results/cytiva/PE_only/LCA_results.xlsx"
    
    # Begin writing to the Excel file
    with pd.ExcelWriter(filename) as writer:
        for impact in impact_categories:
            # Generate both tables for each impact category
            displacement_table = bst.report.lca_displacement_allocation_table(
                [process.system], impact, products
            )
            inventory_table = bst.report.lca_inventory_table(
                [process.system], impact, products
            )

            # Save to a new sheet named e.g., "GWP Tables", "FFC Tables", etc.
            sheet_name = f"{impact} Tables"
            displacement_table.to_excel(writer, sheet_name=sheet_name)
            inventory_table.to_excel(
                writer,
                sheet_name=sheet_name,
                startrow=displacement_table.shape[0] + 2  # Leave a blank row in between
            )

    print(f"LCA tables saved to {filename}")

#FFC bar chart
def ffc_barchart():
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
    weight_counts = {k: v for k, v in weight_counts.items() if v.sum() / total_ffc >= 0.003}

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
    
    fig = plt.gcf()
    ax = fig.axes[0]  # get the first (and likely only) axes object
    xlim = ax.get_xlim()
    if xlim[1] < 2.2:
        ax.set_xlim(xlim[0], 2.2)
    
    plt.show()
    low = pd.DataFrame(lower[metric_index])
    print('baseline: ', baseline[metric_index])
    high = pd.DataFrame(upper[metric_index])
    
    return low, high



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
    plt.show()
    low = pd.DataFrame(lower[metric_index])
    print('baseline: ', baseline[metric_index])
    high = pd.DataFrame(upper[metric_index])
    return low, high

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
    folder = './results/cytiva/PE_only'
    #folder = os.path.join(folder, 'results')
    for scenario in ('experimental', 'NREL'):
        pm = process
        filename = f'TESTcytiva_{scenario}_detailed_report.xlsx'
        file = os.path.join(folder, filename)
        pm.system.save_report(file)
        
def save_detailed_expenditure_tables(sigfigs=3):
    folder = './results/cytiva/PE_only'
    #folder = os.path.join(folder, 'results')
    filename = 'DEBUGexpenditures.xlsx'
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
    folder = './results/cytiva/PE_only'
    filename = 'life_cycle.xlsx'
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
