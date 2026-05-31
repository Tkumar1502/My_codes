# -*- coding: utf-8 -*-
"""
Created on Tue Apr 21 15:27:00 2026

@author: tlnu
"""

import biosteam as bst
import importlib

bst.settings.ID_magic = False
from plastics import strap
from chaospy.distributions import Uniform
import pandas as pd
from warnings import filterwarnings
import numpy as np
import matplotlib.pyplot as plt
import os
from datetime import datetime
from matplotlib.ticker import PercentFormatter, FuncFormatter
from matplotlib.colors import TwoSlopeNorm, Normalize

importlib.reload(strap.process_model)

filterwarnings('ignore')

bst.main_flowsheet.clear()

if hasattr(strap.MagnetRecovery, 'cache'):
    strap.MagnetRecovery.cache.clear()


process = strap.MagnetRecovery(
    processing_capacity = 325,      #tons
    sell_leftover_plastic = True,
    simulate=False)

process.tea.labor_cost = 580000 + process.HS.total_salary
process.S2.outs[1].ID = 'NdFeB Magnets'
process.s22.ID = 'Impurities'

products = [process.NdFeB_Magnets, process.HDPE_resins]
process.products[:] = [process.HDPE_resins, process.NdFeB_Magnets]

process.HDPE_resins.price = 1.20
process.NdFeB_Magnets.price = 100

#update the processing conditions if model starts behaving funky
process.tea.operating_days = 328.5
process.set_processing_capacity(325)
process.system.simulate()

process.system.diagram()

#TeA parameters
#process.set_IRR(0.15)



#%% functions for TEA
def profit(conc, capacity,solvent_loss,plastic_conc,solvent_price,feedstock_price, NdFeB_price, HDPE_price):
    conc = 1-conc
    process.set_processing_capacity(capacity)
    process.set_polymer_mass_fraction(conc)
    process.set_solvent_loss(solvent_loss)
    process.set_dissolution_capacity(plastic_conc)
    process.set_solvent_price(solvent_price)
    process.set_feedstock_price(feedstock_price)
    process.NdFeB_Magnets.price = NdFeB_price
    process.HDPE_resins.price = HDPE_price
    
    
    process.system.simulate()
    return process.tea.net_earnings 



def my_irr(conc, capacity, solvent_loss, plastic_conc,solvent_price, feedstock_price, NdFeB_price, HDPE_price):
    conc = 1-conc
    process.set_processing_capacity(capacity)
    process.set_polymer_mass_fraction(conc)
    process.set_solvent_loss(solvent_loss)
    process.set_dissolution_capacity(plastic_conc)
    process.set_solvent_price(solvent_price)
    process.set_feedstock_price(feedstock_price)
    process.NdFeB_Magnets.price = NdFeB_price
    process.HDPE_resins.price = HDPE_price
    
    process.system.simulate()
    return process.tea.solve_IRR()*100
 
    

#%% making a tornado plot (sensitivity plot) Here we are showing different parameters vs net earnnings
'''
first let's have our baseline scenario
processing_capacity = 325
NdFeB concentration = 0.86 (86% of the feed)
Solvent loss = 0.1%
Dissolution capacity (HDPE conc. in Xylenes) = 5.0 wt.%>
Solvent price = 2.19 usd/kg
Feedstock price = 0.1 usd/kg 


'''

def tornado_plot():
    base_params = {
        'conc': 0.86,          # Changed 'Conc' to 'conc'
        'capacity': 325,      # Changed 'Capacity' to 'capacity'
        'solvent_loss': 0.1, 
        'plastic_conc': 5.0, 
        'solvent_price': 2.19, 
        'feedstock_price': 0.1,
        'NdFeB_price': 100,
        'HDPE_price': 1.20
       
        
    }
    
    # 2. Update these keys to match the dictionary above
    ranges = {
        'conc': (0.80, 0.90),
        'capacity': (300, 350),
        'solvent_loss': (0.01, 1.0),
        'plastic_conc': (2.0, 7.5),
        'solvent_price': (2.19, 3.25),
        'feedstock_price': (0.1, 0.25),
        'NdFeB_price': (90.0,150.0),
        'HDPE_price': (1.00, 1.50)
        
    }
    
    label_map = {
        'conc': 'NdFeB Concentration',
        'capacity': 'Processing Capacity (tons)',
        'solvent_loss': 'Solvent Loss Rate',
        'plastic_conc': 'Feed-to-Solvent (wt%)',
        'solvent_price': 'Solvent Price ($/kg)',
        'feedstock_price': 'Feedstock Price ($/kg)',
        'NdFeB_price': 'Neodymium Magnet Price ($/kg)',
        'HDPE_price': 'HDPE Price ($/kg)'
        
    }
    
   
    
    results = []
    results_irr = []  
    
    
    base_profit_m = profit(**base_params) / 1e6
    base_irr_m = my_irr(**base_params)
  
    
    for param, (low, high) in ranges.items():
        # Test Low
        
        p_low = base_params.copy(); p_low[param] = low
        p_high = base_params.copy(); p_high[param] = high
        
        
       

        # Now call profit with the FILTERED dictionaries
        low_profit_m = profit(**p_low) / 1e6
        high_profit_m = profit(**p_high) / 1e6
        
        low_irr = my_irr(**p_low)
        high_irr = my_irr(**p_high)
        
        
       
        results.append(
            (label_map[param], low_profit_m, high_profit_m, abs(high_profit_m - low_profit_m), low, high)
            )
        results_irr.append((label_map[param], low_irr, high_irr, abs(high_irr - low_irr), low,high)
                           )
        
    plot_it(results, base_profit_m, "Sensitivity: Net Earnings", "Profit ($ Millions)", "skyblue")
    plot_it(results_irr, base_irr_m, "Sensitivity: IRR", "IRR (%)", "teal")



#%%plot function with plot design
def plot_it(data, base_val, title, x_label, color):
        data.sort(key=lambda x: x[3], reverse=True)
        names, lows, highs = [d[0] for d in data], [d[1] for d in data], [d[2] for d in data]
        p_lows = [d[4] for d in data]
        p_highs = [d[5] for d in data]
        
        
       
        plt.figure(figsize=(10, 5))
        
        widths = [h-l for h,l in zip (highs,lows)]
        
        bars = plt.barh(
                names,
                widths,
                left=lows,
                color=color,
                edgecolor='black',
                alpha=0.8
            )
 # low value to baseline = red
        plt.barh(
                names,
                [base_val - l for l in lows],
                left=lows,
                color='red',
                edgecolor='black',
                alpha=0.8
            )

        # baseline to high value = green
        plt.barh(
                names,
                [h - base_val for h in highs],
                left=[base_val] * len(names),
                color='green',
                edgecolor='black',
                alpha=0.8
            )        
        
        #plt.barh(names, [h - l for h, l in zip(highs, lows)], left=lows, color=color, edgecolor='black', alpha=0.8)
        
        
       
        
        

        
        plt.axvline(base_val, color='white', linestyle='--', label=f'Baseline: {base_val:.2f}')
        
        
       
        x_min, x_max = plt.xlim()
        x_gap = 0.03 * (x_max - x_min)   # spacing from bar
        
        for bar, lo, hi in zip(bars, p_lows, p_highs):
            y = bar.get_y() + bar.get_height() / 2
        
            # Low value (left of bar)
            plt.text(
                bar.get_x() - x_gap,
                y,
                f'{lo}',
                ha='right',
                va='center',
                fontsize=12,
                fontweight = 'bold',
                color = 'grey',
                alpha = 0.85
            )
        
            # High value (right of bar)
            plt.text(
                bar.get_x() + bar.get_width() + x_gap,
                y,
                f'{hi}',
                ha='left',
                va='center',
                fontsize=12,
                fontweight = 'bold',
                color = 'grey',
                alpha = 0.85
            )
 
       


            
        # Add padding to x-axis
        all_points = lows + highs + [base_val]
        padding = (max(all_points) - min(all_points)) * 0.25
        plt.xlim(min(all_points) - padding, max(all_points) + padding)
    
        
        #plt.title(title); plt.xlabel(x_label); plt.legend(); plt.grid(axis='x', alpha=0.3)
        
        ax = plt.gca()

        plt.xlabel(x_label, fontsize=14, fontweight='bold')
        plt.xticks(fontsize=12, fontweight='bold')
        plt.yticks(fontsize=12, fontweight='bold')
        plt.legend(prop={'weight': 'bold', 'size': 12})
        
        plt.grid(False)
        
        for spine in ax.spines.values():
            spine.set_linewidth(2)
        
        leg = plt.legend(prop = {'weight': 'bold', 'size':12})
        leg.get_frame().set_edgecolor('black')
        leg.get_frame().set_linewidth(2)
        
        plt.tight_layout(); plt.show()   
        




