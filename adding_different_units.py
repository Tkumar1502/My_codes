# -*- coding: utf-8 -*-
"""
Created on Thu Mar 19 15:32:58 2026

@author: tlnu


git commit -am "your messsage" && git push
"""

#adding this on macbook 4:49pm-4/3/26

#adding back from work pc 5:16-4/3/26
import biosteam as bst
#bst.settings.ID_magic = False
from plastics import strap
from chaospy.distributions import Uniform
import pandas as pd
from warnings import filterwarnings
import numpy as np
import matplotlib as plt
import os
from datetime import datetime
from matplotlib.ticker import PercentFormatter, FuncFormatter
from matplotlib.colors import TwoSlopeNorm, Normalize


filterwarnings('ignore')

bst.main_flowsheet.clear() #resets names for consistency


#define scenario

#total bags = 494,214; wt. of each bag = 657g
processing_capacity = 325 #tons/year                        #should this be in ktons or tons

process = strap.BaselineSTRAPProcess(
    scenario = 'PE/Xylene',
    target_plastic_percent = 13.7, #HDPE from the impellers (avg)
    
    processing_capacity = processing_capacity,
    sell_leftover_plastic = True,
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
HS.ins[0] = process.plastic        #plastic here is the feed that you entered
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

T_mag = bst.StorageTank(ins= S_mag.outs[1], outs='NdFeB Magnets')
# T_mag.outs[0] = process.leftover_plastic

units.append(T_mag)
process.NdFeB_Magnets = T_mag.outs[0]

#removing unnecessary units for magnets
process.U6.ins[0] = process.P3.outs[0]
units.remove(process.U4)
units.remove(process.U5)   
units.remove(process.M2)     


process.U6.outs[0].ID = 'Impurities'

#removing the streams
process.M3.ins[:] = [process.M3.ins[i] for i in (0, 2, 3)]
#process.T5.ins.append(process.s22)

process.system.update_configuration(units)  #we need update the system after adding handsorting


process.system.diagram(format = 'svg') #just to check if the changes are working or not



#adjust tea parameters
process.tea.labor_cost = 1350000    #2 operators X $60k, 1 loader x $50k, 1 engineer x $100k


#define products and set sale prices
products = [process.PE_resin, process.NdFeB_Magnets]             #shouldn't this be process.products
process.PE_resin.price = 1.20
process.NdFeB_Magnets.price = 100


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

'''
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
    
    
'''