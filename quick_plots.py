# -*- coding: utf-8 -*-
"""
Created on Sun Sep 15 18:34:20 2024

@author: yoelr
"""
from warnings import filterwarnings
from .process_model import BaselineSTRAPProcess, STRAPMSWProcess
import matplotlib.pyplot as plt
import biosteam as bst
import os
import numpy as np

__all__ = (
    'STRAP_coproduct_tornado_plot',
    'STRAP_incineration_tornado_plot',
    'STRAP_MSW_tornado_plot',
    'plot_MSP_across_capacity_composition',
)

images_folder = os.path.join(os.path.dirname(__file__), 'images')
results_folder = os.path.join(os.path.dirname(__file__), 'results')

def MSP_at_capacity_composition(mass_fraction, capacity, process_model):
    process_model.set_mass_fraction.setter(mass_fraction / 100)
    process_model.set_processing_capacity.setter(capacity * 1e3)
    process_model.system.simulate()
    MSP = process_model.MSP()
    return np.array([MSP])

def plot_MSP_across_capacity_composition(load=True):
    bst.plots.set_font(size=12, family='sans-serif', font='Arial')
    pm = BaselineSTRAPProcess(simulate=True, burn_coproducts=False)
    pm.last_capacity = None
    xlim = np.array(pm.set_mass_fraction.bounds) * 100
    ylim = np.array(pm.set_processing_capacity.bounds) / 1e3
    X, Y, Z = bst.plots.generate_contour_data(
        MSP_at_capacity_composition,
        file=os.path.join(results_folder, 'MSP_capacity_composition.npy'),
        load=load, save=True,
        xlim=xlim, ylim=ylim,
        args=(pm,),
        n=10,
    )
    
    # Plot contours
    ylabel = "Processing capacity [$\mathrm{10}^{3} \cdot \mathrm{MT} \cdot \mathrm{yr}^{\mathrm{-1}}$]"
    xlabel = 'Polymer fraction [wt %]'
    yticks = [2.5, 5, 7.5, 10]
    xticks = [10, 30, 50, 70, 90]
    metric_bar = bst.plots.MetricBar(
        'MSP', '$\mathrm{USD} \cdot \mathrm{kg}^{\mathrm{-1}}$', plt.cm.get_cmap('copper_r'), 
        bst.plots.rounded_tickmarks_from_data(Z, 5, 1, expand=0, p=0.5), 
        15, 1
    )
    fig, axes, CSs, CB, other_axes = bst.plots.plot_contour_single_metric(
        X, Y, Z[:, :, None], xlabel, ylabel, xticks, yticks, metric_bar,  
        fillcolor=None, styleaxiskw=dict(xtick0=False), label=True,
    )
    plt.subplots_adjust(left=0.2, right=0.9, wspace=0.15, hspace=0.2, top=0.85, bottom=0.15)
    for i in ('svg', 'png'):
        file = os.path.join(images_folder, f'MSP_capacity_composition.{i}')
        plt.savefig(file, dpi=900, transparent=True)

def STRAP_coproduct_tornado_plot():
    scenario = BaselineSTRAPProcess.Scenario('Toluene', 'PE', 50, 5000, True, False, True)
    pm = BaselineSTRAPProcess(scenario=scenario)
    pm.system.simulate()
    baseline, lower, upper = pm.model.single_point_sensitivity(etol=0.01)
    GWP, MSP = pm.model.metrics
    metric_index = MSP.index
    index = [i.describe(distribution=False).replace('(', '\n(')
             for i in pm.model.parameters]
    MSP_plot, _ = bst.plots.plot_single_point_sensitivity(
        baseline[metric_index],
        lower[metric_index],
        upper[metric_index],
        name='MSP [USD/kg]',
        sort=True,
        top=5,
        index=index
    )
    plt.subplots_adjust(left=0.45, right=0.975, top=0.98, bottom=0.1)
    for i in ('svg', 'png'):
        name = f'STRAP_coproduct_tornado_plot_MSP.{i}'
        file = os.path.join(images_folder, name)
        plt.savefig(file, dpi=900, transparent=True)
    metric_index = GWP.index
    GWP_plot, _ = bst.plots.plot_single_point_sensitivity(
        baseline[metric_index],
        lower[metric_index],
        upper[metric_index],
        name='GWP [kg-CO2e/kg]',
        sort=True,
        top=5,
        index=index
    )
    plt.subplots_adjust(left=0.50, right=0.975, top=0.98, bottom=0.1)
    for i in ('svg', 'png'):
        name = f'STRAP_coproduct_tornado_plot_GWP.{i}'
        file = os.path.join(images_folder, name)
        plt.savefig(file, dpi=900, transparent=True)
        
def STRAP_incineration_tornado_plot():
    scenario = BaselineSTRAPProcess.Scenario('Toluene', 'PE', 50, 5000, False, True, True)
    pm = BaselineSTRAPProcess(scenario=scenario)
    pm.system.simulate()
    baseline, lower, upper = pm.model.single_point_sensitivity(etol=0.01)
    GWP, MSP = pm.model.metrics
    metric_index = MSP.index
    index = [i.describe(distribution=False).replace('(', '\n(')
             for i in pm.model.parameters]
    MSP_plot, _ = bst.plots.plot_single_point_sensitivity(
        baseline[metric_index],
        lower[metric_index],
        upper[metric_index],
        name='MSP [USD/kg]',
        sort=True,
        top=5,
        index=index
    )
    plt.subplots_adjust(left=0.45, right=0.975, top=0.98, bottom=0.1)
    for i in ('svg', 'png'):
        name = f'STRAP_incineration_tornado_plot_MSP.{i}'
        file = os.path.join(images_folder, name)
        plt.savefig(file, dpi=900, transparent=True)
    metric_index = GWP.index
    GWP_plot, _ = bst.plots.plot_single_point_sensitivity(
        baseline[metric_index],
        lower[metric_index],
        upper[metric_index],
        name='GWP [kg-CO2e/kg]',
        sort=True,
        top=5,
        index=index
    )
    plt.subplots_adjust(left=0.50, right=0.975, top=0.98, bottom=0.1)
    for i in ('svg', 'png'):
        name = f'STRAP_incineration_tornado_plot_GWP.{i}'
        file = os.path.join(images_folder, name)
        plt.savefig(file, dpi=900, transparent=True)
        
def STRAP_MSW_tornado_plot():
    filterwarnings('ignore')
    pm = STRAPMSWProcess()
    pm.system.simulate()
    baseline, lower, upper = pm.model.single_point_sensitivity(etol=0.01)
    GWP, MSP = pm.model.metrics
    metric_index = MSP.index
    index = [i.describe(distribution=False).replace('(', '\n(')
             for i in pm.model.parameters]
    MSP_plot, _ = bst.plots.plot_single_point_sensitivity(
        baseline[metric_index],
        lower[metric_index],
        upper[metric_index],
        name='MSP [USD/kg]',
        sort=True,
        top=8,
        index=index
    )
    plt.subplots_adjust(left=0.45, right=0.975, top=0.98, bottom=0.1)
    for i in ('svg', 'png'):
        name = f'STRAP_MSW_tornado_plot_MSP.{i}'
        file = os.path.join(images_folder, name)
        plt.savefig(file, dpi=900, transparent=True)
    metric_index = GWP.index
    GWP_plot, _ = bst.plots.plot_single_point_sensitivity(
        baseline[metric_index],
        lower[metric_index],
        upper[metric_index],
        name='GWP [kg-CO2e/kg]',
        sort=True,
        top=8,
        index=index
    )
    plt.subplots_adjust(left=0.50, right=0.975, top=0.98, bottom=0.1)
    for i in ('svg', 'png'):
        name = f'STRAP_MSW_tornado_plot_GWP.{i}'
        file = os.path.join(images_folder, name)
        plt.savefig(file, dpi=900, transparent=True)