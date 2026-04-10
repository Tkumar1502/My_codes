#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 14 11:05:08 2025

@author: charlesgranger
"""
############### CHANGES TO BE MADE #####################################

# how can we set bounds and baseline for a multiple wash system
# fix the detailed LCA table
#low pressure steam counts towards GWP but is not seen in the inventory?

#can we move centrifilter before surge tank and pump?

#8) need to account for burning of lost solvent emissions
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
processing_capacity=5e3

#single wash
process = strap.BaselineSTRAPProcess(
    scenario='PC/THF',
    target_plastic_percent=70, # %
    processing_capacity=processing_capacity, # MT/yr
    sell_leftover_plastic=False,
    burn_leftover_plastic=False,
    facilities=False,
    simulate=False,
)

#adjust operating days
process.tea.operating_days = 328.5
process.tea.labor_cost = 1350000 # 2 operators x $60K, 1 loader x $50k, 1 engineer X $100k

#create flowchart
process.system.diagram(format='png') # Note how the leftover plastic is sent to heat and power generation.

#define products and set sale prices
products = [process.PC_resin, process.leftover_plastic]
process.PC_resin.price = 2.50 # 
process.leftover_plastic.price = -0.062 # tipping fee: https://wasteadvantagemag.com/eref-report-shows-10-increase-in-u-s-landfill-tipping-fees-largest-increase-since-2022/#:~:text=EREF%20released%20their%202024%20Analysis,landfills%20charge%20about%2034%25%20more.

#bounds 1-wash
process.set_polymer_mass_fraction.bounds=(.6, .8)
process.set_dissolution_capacity.bounds=(5, 10)
process.set_solvent_loss.bounds=(0.01,1)
process.set_precipitation_temperature.bounds=(303,313)
process.set_dissolution_temperature.bounds=(330,339)
process.set_centrifuged_plastic_solvent_content.bounds=(25, 75)
process.set_feedstock_distance.bounds=(50,2000)
process.set_feedstock_price.bounds=(0,0.10)

process.set_solvent_price.bounds = (1.00, 3.00)

#set baseline 1-wash
process.plastic.ID='Feedstock_plastic'
process.set_dissolution_temperature.baseline = 336
process.set_solvent_loss.baseline = 0.1
process.set_polymer_mass_fraction.baseline = 0.7
process.set_dissolution_capacity.baseline = 7.5
process.set_precipitation_temperature.baseline = 35+273
process.set_feedstock_price.baseline=(0.05)
process.set_feedstock_distance.baseline=(2000)
process.set_IRR.baseline=0.15
process.set_solvent_price.baseline = 2.00

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
lb, ub = processing_capacity/5, processing_capacity*10
pcs = np.linspace(lb, ub, 20)

def irr_contour(lb, ub, low_price, high_price):
    # Define x and y bounds
    capacities = np.linspace(lb, ub, 50) #50
    prices = np.linspace(low_price, high_price, 20) #20

    # Create a 2D array to hold the IRR values
    IRR = np.zeros((len(prices), len(capacities)))

    for i, cap in enumerate(capacities): 
        for j, price in enumerate(prices):
            process.set_processing_capacity(cap)
            process.PC_resin.price = price
            process.system.simulate()
            IRR[j, i] = process.tea.solve_IRR()

    # Create meshgrid
    X, Y = np.meshgrid(capacities, prices)

    # Set figure size correctly using subplots
    fig, ax = plt.subplots(figsize=(5, 5))

    # Generate contour levels
    min_val = np.min(IRR)
    max_val = np.max(IRR)
    neg_levels = np.arange(np.floor(min_val * 100) / 100, 0, 1.0)
    pos_levels = np.arange(0.05, np.ceil(max_val * 100) / 100 + 0.001, 0.05)
    levels = np.concatenate((neg_levels, [0], pos_levels))

    norm = TwoSlopeNorm(vmin=min_val, vcenter=0, vmax=max_val)

    # Create filled contour plot
    contourf = ax.contourf(X, Y, IRR, levels=levels, cmap='RdYlGn', norm=norm, extend='both')

    # Add contour lines and labels (only for positive levels)
    contour_lines = ax.contour(X, Y, IRR, levels=np.concatenate(([0], pos_levels)), colors='black', linewidths=0.5)
    # # Manually chosen positions for labels along the contour lines
    # label_positions = [
    #     (2900, 2.3),   # capacity = 15000 tpy, price = $1.2/kg
    #     (3700, 2.5),   # capacity = 30000 tpy, price = $2.0/kg
    #     (4500, 2.65),    
    #     (6000, 2.75),
    #     (8700, 2.75),
    # ]
    
    ax.clabel(
        contour_lines,
        fmt=lambda x: f"IRR = {x:.0%}",  # Format labels
        inline=True,                     # Cuts line under label for clarity
        inline_spacing=5,                 # Space between label and line
        fontsize=11,                      # Font size
        manual=False            # Use these exact positions
    )

    # Add colorbar with adjusted position
    # cbar = fig.colorbar(contourf, ax=ax, pad=0.15)  # pad increases distance from plot
    # cbar.set_label('IRR (%)')
    # cbar.ax.yaxis.set_major_formatter(PercentFormatter(xmax=1))

    # Labels and title
    ax.set_xlabel('Feedstock capacity (metric tpy)')
    ax.set_ylabel('STRAP resin sale price ($/kg)')
    # ax.set_title('Printed Film - IRR')

    # Add faint grid
    ax.grid(True, linestyle='--', color='gray', linewidth=0.5, alpha=0.4)

    # # Calculate and plot TCI on secondary y-axis
    # TCIs = [TCI(i)/1000000 for i in capacities]
    # ax2 = ax.twinx()
    # ax2.plot(capacities, TCIs, color='blue', linewidth=2, label='TCI', linestyle='--')
    # ax2.set_ylabel('Total Capital Investment ($M)')
    # ax2.legend(loc='lower right')
    plt.savefig('./results/hoya_vision/irr_contour.TIF')
    plt.savefig('./results/hoya_vision/irr_contour.png', dpi=300)
    plt.show()

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

def dROI(pc):
    process.set_processing_capacity(pc)
    process.system.simulate()
    droi = process.tea.NPV / process.tea.TCI
    return droi

#get data
def msp_plot():
    MSPs = [MSP(i) for i in pcs]
    TCIs = [TCI(i) for i in pcs]
    market = [2.50 for i in pcs]
    margin = [process.PC_resin.price - MSP(i) for i in pcs]
    
    #plot data
    plt.figure()
    line1, = plt.plot(pcs, MSPs, label='STRAP resin production cost') # plot MSP for different capacities
    line2, = plt.plot(pcs, market, 'r--', label='Market Value PE') #market value
    line4, = plt.plot(pcs, margin, 'k.', label='Margin (sale price - cost') #market value
    plt.ylabel('Cost to produce STRAP resin ($/kg PE)')
    plt.xlabel('Plant Capacity (mton/year)')
    
    # Create a secondary y-axis
    ax2 = plt.gca().twinx()
    line3, = ax2.plot(pcs, TCIs, label='TCI', color='green')
    ax2.set_ylabel('Total Capital Investment ($)')  # Label for the second y-axis
    
    # Combine legends for both axes
    lines = [line1, line2, line3, line4]
    labels = [line.get_label() for line in lines]
    plt.legend(lines, labels, loc='upper center')
    
def margin_plot():
    TCIs = [TCI(i) for i in pcs]
    zero = [0 for i in pcs]
    margin = [process.PC_resin.price - MSP(i) for i in pcs]
    
    #plot data
    plt.figure()
    line1, = plt.plot(pcs, margin, label='Sale price - STRAP cost') # plot MSP for different capacities
    line2, = plt.plot(pcs, zero, 'r--', label=None) #market value
    plt.ylabel('Sale price - Cost to produce STRAP resin ($/kg PE)')
    plt.xlabel('Plant Capacity (mton/year)')
    
    # Create a secondary y-axis
    ax2 = plt.gca().twinx()
    line3, = ax2.plot(pcs, TCIs, label='TCI', color='green')
    ax2.set_ylabel('Total Capital Investment ($)')  # Label for the second y-axis
    
    # Combine legends for both axes
    lines = [line1, line2, line3]
    labels = [line.get_label() for line in lines]
    plt.legend(lines, labels, loc='upper center')
    
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
    output_folder = './results/hoya_vision/'
    os.makedirs(output_folder, exist_ok=True)
    df = pd.DataFrame({
        'Plant Capacity (mton/year)': pcs,
        'Discounted ROI': dROIs,
        'Total Capital Investment ($)': TCIs
    })
    df.to_excel(os.path.join(output_folder, 'droi_results.xlsx'), index=False)
    
    
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
    plt.axhspan(2.00, 3.00, color='blue', alpha=0.15, label='Market Value Range')
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
               labels + ['Vrigin PC - Market Range'], 
               loc='upper center')
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

#CFs for adsorbent
process.adsorption_column.ins[2].ID = 'Adsorbent'
process.adsorption_column.ins[2].set_CF(GWP, 1.78+0.02643,  ) #EF 2.0 activated silica production, production mix, at plant, technology mix, 100% active substance (activated carbon not available)
process.adsorption_column.ins[2].set_CF(FFC, 23.99+0.3452,  ) 
# adsorbent landfilling: EF 2.0 Landfill of polluted inorganic waste, production mix (region specific sites), at landfill site, landfill including leachate treatment and with transport without collection and pre-treatment

#CFs for transportation
#values taken from aurora's paper, database EF 2.0 (EU)
#radius = 1000 #km plastic travels
#process.plastic.set_CF(GWP, process.feedstock_distance * 0.08 * 1/1000, ) #distance (km) * .08 kg co2 / (t*km) * t/1000kg
#process.plastic.set_CF(FFC, process.feedtsock_distance * 1.04 * 1/1000, ) #distance (km) * .08 kg co2 / (t*km) * t/1000kg

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

inventory_table = bst.report.lca_inventory_table(
    [process.system], 
    'GWP', 
    products,)

#LCA barcharts
impact_categories = [GWP, FFC, WU, HTC, HTNC, ETOX, ACD, OZD, POCP]
LDPE_impact_factors = {'GWP': 2.091, 'FFC':75.77, 'WU':1.562, 'HTC':2.66E-8, 
                       'HTNC':2.40E-7, 'ETOX':.69707, 'ACD':.00554, 
                       'OZD':2.11E-10, 'POCP':.00751}

def lca_factor_barcharts(
    factors,
    virgin_values,
    label_offsets_percent=(-0.1, 0.15)  # (Virgin PE %, STRAP PE %)
):
    """
    Create stacked bar charts for multiple LCA impact factors.
    
    Parameters
    ----------
    factors : list of str
        Impact categories to plot (e.g., ["GWP", "AP", "EP"]).
    virgin_values : dict
        Dictionary of virgin resin values (e.g., {"GWP": 2.091, "AP": 0.34}).
    label_offsets_percent : tuple of float
        Relative offset for floating labels for each bar (as fraction of bar height),
        e.g., (0.05, 0.1) = 5% for Virgin PE, 10% for STRAP PE.
    """

    # Units for impact factors
    factor_units = {
        "GWP": "kg CO₂ eq",
        "FFC": "MJ",
        "WU": "m³",
        "HTC":"CTUh",
        "HTNC":"CTUh",
        "ETOX":"CTHe",
        "ACD":"mol H+ eq",
        "OZD": "kg CFC-11 eq",
        "POCP": "kg NMVOC eq",
        # add more factors as needed
    }
    
    for factor in factors:
        # Build LCA table for this factor
        table = bst.report.lca_displacement_allocation_table(
            process.system,
            factor,
            products,
        )
        
        production_methods = ("Virgin PE", "STRAP PE")
        weight_counts = {"Virgin resin": np.array([virgin_values.get(factor, 0), 0])}
        
        total_val = 0
        for factor_name, row in table.iterrows():
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
            total_val += value
        
        # Filter out negligible contributions (<0.5% of total)
        weight_counts = {
            k: v for k, v in weight_counts.items()
            if total_val == 0 or v.sum() / total_val > 0.005
        }
        weight_counts = dict(sorted(weight_counts.items(), key=lambda x: x[1].sum(), reverse=True))
        
        # Plot
        width = 0.5
        fig, ax = plt.subplots()
        bottom = np.zeros(2)
        
        for label, weight_count in weight_counts.items():
            ax.bar(production_methods, weight_count, width, label=label, bottom=bottom)
            bottom += weight_count
        
        # Add floating labels with individual relative offsets and 3 sig figs
        for i, total in enumerate(bottom):
            bar_height = total
            offset_pct = label_offsets_percent[i] if i < len(label_offsets_percent) else 0.05
            offset = bar_height * offset_pct if bar_height > 0 else 0.05
            
            # 3 sig figs, use scientific notation if < 0.01
            if abs(total) < 0.01 and total != 0:
                label_text = f"{total:.2e}"
            else:
                label_text = f"{total:.3g}"
            
            ax.text(
                i, bar_height + offset, label_text,
                ha='center', va='bottom',
                bbox=dict(facecolor='white', edgecolor='none', pad=1.0)
            )
        
        # Labels & legend
        unit = factor_units.get(factor, "arbitrary units")
        ax.set_ylabel(f"{factor} ({unit} per kg PE)")
        ax.legend(loc="best")
        #ax.set_title(f"{factor} of PE Production")
        plt.tight_layout()
        
        # Save each plot by factor name
        # plt.savefig("./results/pallet_wrap/" + f"{factor}_PE_barchart.png", dpi=300, bbox_inches="tight")
        plt.show()

# #GWP bar chart
def gwp_barchart():
    production_methods = ("Virgin PE", "STRAP PE")
    
    # Initialize the weight_counts dictionary
    weight_counts = {"Virgin resin": np.array([4.59, 0])}  #from EF 2.0 Polycarbonate (PC) granulate, production mix, at plant, Technology mix, dipenyl carbonate route and phosgene route, 1.20–1.22 g/cm3 - world w/o EU -28+EFTA 
    
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
    weight_counts = {k: v for k, v in weight_counts.items() if v.sum() / total_gwp > 0.005}
    
    # Sort weight_counts by total GWP contribution in descending order
    weight_counts = dict(sorted(weight_counts.items(), key=lambda x: x[1].sum(), reverse=True))
    
    # Set the bar width
    width = 0.5
    
    # Create the stacked bar chart
    fig, ax = plt.subplots(figsize=(5,5))
    bottom = np.zeros(2)
    
    for label, weight_count in weight_counts.items():
        p = ax.bar(production_methods, weight_count, width, label=label, bottom=bottom)
        bottom += weight_count
    
    # Customize the plot
    # ax.set_title("Global Warming Potential of PE Production")
    ax.set_ylabel("GWP (kg CO2 eq. / kg PE)")
    ax.legend(loc="upper right")
    #plt.show()



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
    filename="./results/hoya_vision/LCA_results.xlsx"
    with pd.ExcelWriter(filename) as writer:
        GWP_table.to_excel(writer, sheet_name="GWP Tables")
        inventory_GWP_table.to_excel(writer, sheet_name="GWP Tables", startrow=GWP_table.shape[0] + 2)
        FFC_table.to_excel(writer, sheet_name="FFC Tables")
        inventory_FFC_table.to_excel(writer, sheet_name="FFC Tables", startrow=FFC_table.shape[0] + 2)
    
    print(f"LCA tables saved to {filename}")

#FFC bar chart
def ffc_barchart():
    production_methods = ("Virgin PE", "STRAP PE")
    
    # Initialize the weight_counts dictionary
    weight_counts = {"Virgin resin": np.array([101.57, 0])}  # add virgin resin FFC here

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
    folder = './results/hoya_vision'
    #folder = os.path.join(folder, 'results')
    for scenario in ('experimental', 'NREL'):
        pm = process
        filename = f'hoya_vision_{scenario}_detailed_report.xlsx'
        file = os.path.join(folder, filename)
        pm.system.save_report(file)
        
def save_detailed_expenditure_tables(sigfigs=3):
    folder = './results/hoya_vision'
    #folder = os.path.join(folder, 'results')
    filename = 'expenditures.xlsx'
    file = os.path.join(folder, filename)
    writer = pd.ExcelWriter(file)
    product = 'PC_resin'
    process_models = [
        process
        for i in ('best-experimental', 'best-NREL')
    ]
    for pm in process_models: pm.PC_resin.price = pm.tea.solve_price(pm.PC_resin)
    systems = [i.system for i in process_models]
    names = ['baseline', 'potential']
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
    folder = './results/hoya_vision'
    filename = 'life_cycle.xlsx'
    file = os.path.join(folder, filename)
    
    # Create Excel writer
    writer = pd.ExcelWriter(file)
    
    # Define streams and system names
    streams = [getattr(process_model, 'PC_resin')]
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
                        process_model.GWP_PC_resin()]])
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
