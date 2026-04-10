# -*- coding: utf-8 -*-
"""
"""
from plastics import strap as sp
from warnings import filterwarnings, warn
from biosteam.utils import GG_colors, CABBI_colors, colors
import biosteam as bst
import numpy as np
import pandas as pd
import os
from thermosteam.units_of_measure import format_units
from thermosteam.utils import roundsigfigs
from matplotlib import pyplot as plt
import matplotlib.patches as mpatches
from colorpalette import Color
import yaml
import matplotlib.colors as clr

__all__ = (
    'plot_monte_carlo',
    'run_monte_carlo_tipping_fee',
    'plot_monte_carlo_across_tipping_fee',
    'plot_IRR_across_processing_capacity_tipping_fee',
    'plot_monte_carlo_across_processing_capacity',
    'run_monte_carlo_across_process_capacity',
    'run_monte_carlo_STRAPMSW',
    'run_monte_carlo', 'plot_spearman', 'plot_kde',
    'plot_spearman_both',
    'sobol_analysis',
    'plot_sobol',
    'plot_breakdowns',
    'plot_MSP_GWP_across_titer_glucose_yield',
    'montecarlo_results',
    'get_distributions',
    'run_monte_carlo_STRAP_CUWP_single_wash_conceptual_configurations',
    'plot_kde_STRAP_CUWP_single_wash_conceptual_configurations',
    'plot_kde_TCI_MSP',
    'plot_kde_CI_MSP',
    'plot_monte_carlo_across_processing_capacity_IRR',
)

letter_color = colors.neutral.shade(25).RGBn
results_folder = os.path.join(os.path.dirname(__file__), 'results')
images_folder = os.path.join(os.path.dirname(__file__), 'images')
line_color = Color(fg='#8E9BB3').RGBn
vertical_line_color = GG_colors.red.shade(40).tint(20).RGBn
INST_EQ_COST = 'Inst. eq. cost'
ELEC_CONS = 'Elec. cons.'
ELEC_PROD = 'Elec. prod.'
COOLING = 'Cooling'
HEATING = 'Heating'
CAPITAL_UNITS = 'MM$'
ELEC_UNITS = 'kW/kg'
DUTY_UNITS = 'MJ/kg'

def functional_unit(f, value):
    return lambda: f() / value

def sobol_file(name, extension='xlsx'):
    filename = name + '_sobol'
    filename += '.' + extension
    return os.path.join(results_folder, filename)

def monte_carlo_file(name, extension='xlsx'):
    filename = name + '_monte_carlo'
    filename += '.' + extension
    return os.path.join(results_folder, filename)

def spearman_file(name):
    filename = name + '_spearman'
    filename += '.xlsx'
    return os.path.join(results_folder, filename)

def autoload_file_name(name):
    filename = name
    return os.path.join(results_folder, filename)

def set_font(size=8, family='sans-serif', font='Arial'):
    import matplotlib
    fontkwargs = {'size': size}
    matplotlib.rc('font', **fontkwargs)
    params = matplotlib.rcParams
    params['font.' + family] = font
    params['font.family'] = family

def set_figure_size(width=None, aspect_ratio=None, units=None): 
    # units default to inch
    # width defaults 6.614 inches
    # aspect ratio defaults to 0.65
    if aspect_ratio is None:
        aspect_ratio = 0.65
    if width is None or width == 'full':
        width = 6.6142
    elif width == 'half':
        width = 6.6142 / 2
    else:
        if units is not None:
            from thermosteam.units_of_measure import convert
            width = convert(width, units, 'inch')
    import matplotlib
    params = matplotlib.rcParams
    params['figure.figsize'] = (width, width * aspect_ratio)

def get_spearman_names(parameters):
    from plastics.strap import STRAPMSWProcess
    pm = STRAPMSWProcess(simulate=False)
    name = 'name'
    full_name = 'full_name'
    spearman_labels = {
        i: full_name for i in parameters
    }
    # spearman_labels[pm.set_IRR] = 'IRR'
    
    def with_units(f, name, units=None):
        d = f.distribution
        dname = type(d).__name__
        if units is None: units = f.units
        if dname == 'Triangle':
            distribution = ', '.join([format(j, '.3g')
                                      for j in d._repr.values()])
        elif dname == 'Uniform':
            distribution = ' $-$ '.join([format(j, '.3g')
                                         for j in d._repr.values()])
        elif dname == 'Trunc': # Must be truncated gaussian
            normal, *_ = d._repr.values()
            distribution = r' $\pm$ '.join([format(j, '.3g')
                                         for j in normal._repr.values()])
        if units is None:
            return f"{name}\n[{distribution}]"
        else:
            return f"{name}\n[{distribution} {format_units(units)}]"
        
    def get_full_name(f):
        a = f.element_name
        if a == 'Cofermentation':
            a = 'Co-Fermentation'
        b = f.name
        if b == 'GWP': 
            return f"{a} {b}"
        else:
            return f"{a} {b.lower()}"
        
    for i, j in tuple(spearman_labels.items()):
        if j == name:
            spearman_labels[i.index] = with_units(i, i.name)
        elif j == full_name:
            spearman_labels[i.index] = with_units(i, get_full_name(i))
        elif isinstance(j, tuple):
            spearman_labels[i.index] = with_units(i, *j)
        elif isinstance(j, str):
            spearman_labels[i.index] = with_units(i, j)
        else:
            raise TypeError(str(j))
        del spearman_labels[i]
    
    return spearman_labels

def plot_monte_carlo():
    set_font(size=8)
    color_wheel = CABBI_colors.wheel([
        'blue_light', 'green_dirty', 'orange', 'green', 'grey',
        'orange', 'orange', 'green', 'orange', 'green',
    ])
    metrics_TEA = (
        'resin_production', 
        'ethanol_production', 
        'electricity_production',
        'TCI', 
        'IRR',
        'MSP', 
    )
    metrics_LCA = (
        'GWP_ethanol', 
        'GWP_electricity', 
        'GWP_polymer_resin',
        'FFC_ethanol',
        'FFC_electricity',
        'FFC_polymer_resin',
        'WU_ethanol',
        'WU_electricity',
        'WU_polymer_resin',
    )
    
    # TEA
    set_figure_size(aspect_ratio=0.65)
    fig, axes = _plot_monte_carlo(
        ncols=3,
        expand=0.15, 
        labels=['Full',
                'Baseline',
                'Potential'],
        xrot=30,
        metrics=metrics_TEA,
        color_wheel=color_wheel,
    )
    for ax, letter in zip(axes, 'ABCDEFGHIJKLMNO'):
        plt.sca(ax)
        ylb, yub = plt.ylim()
        letter = f'({letter.lower()})'
        plt.text(-0.25, ylb + (yub - ylb) * 0.92, letter, color=letter_color,
                 horizontalalignment='center',verticalalignment='center',
                 fontsize=10, fontweight='bold')
    plt.subplots_adjust(left=0.12, right=0.95, wspace=0.55, top=0.98, bottom=0.2)
    for i in ('svg', 'png'):
        file = os.path.join(images_folder, f'MC_TEA.{i}')
        plt.savefig(file, transparent=True)
    
    # LCA
    set_figure_size(aspect_ratio=1.05)
    fig, axes = _plot_monte_carlo(
        ncols=3,
        expand=0.15, 
        labels=['Full',
                'Baseline',
                'Potential'],
        xrot=30,
        metrics=metrics_LCA,
        color_wheel=color_wheel,
    )
    for ax, letter in zip(axes, 'ABCDEFGHIJKLMNO'):
        plt.sca(ax)
        ylb, yub = plt.ylim()
        letter = f'({letter.lower()})'
        plt.text(-0.25, ylb + (yub - ylb) * 0.92, letter, color=letter_color,
                 horizontalalignment='center',verticalalignment='center',
                 fontsize=10, fontweight='bold')
    plt.subplots_adjust(left=0.12, right=0.95, wspace=0.55, top=0.98, bottom=0.2)
    for i in ('svg', 'png'):
        file = os.path.join(images_folder, f'MC_LCA.{i}')
        plt.savefig(file, transparent=True)

def _plot_monte_carlo(scenarios=None, metrics=None, labels=None, 
                     ncols=2, color_wheel=None, tickmarks=None,
                     ylabels=None, xrot=None, expand=None,):
    if scenarios is None: scenarios = ['all', 'baseline', 'potential']
    if color_wheel is None: color_wheel = CABBI_colors.wheel()
    pm = sp.STRAPMSWProcess(simulate=False)
    if ylabels is None: 
        subnames = ('Ethanol', 'Electricity', 'Polymer resin')
        names = ('Carbon intensity', 'Water usage', 'Fossil fuel consumption')
        shorthands = ('GWP', 'WU', 'FFC')
        replacements = {
            f'{j}\n{i}': f'{k}$_\mathrm' + '{' + j.replace('Polymer r', 'R').replace(' ', r'\ ') + '}$'
            for i, k in zip(names, shorthands)
            for j in subnames
        }
        ylabels = []
        for i in metrics:
            label = getattr(pm, i).label() 
            for i, j in replacements.items():
                if i in label:
                    label = label.replace(i, j)
                    break
            ylabels.append(label)
    # Data array
    rows = metrics
    N_rows = len(rows)
    columns = scenarios
    N_cols = len(columns)
    # Subplots
    nrows = int(round(N_rows / ncols))
    fig, axes_box = plt.subplots(ncols=ncols, nrows=nrows)
    plt.subplots_adjust(wspace=0.45)
    axes = axes_box.transpose()
    axes = axes.flatten()
    xtext = labels or scenarios
    N_marks = len(xtext)
    xticks = tuple(range(N_marks))
    
    def get_data(metric, name):
        df = get_monte_carlo('STRAPMSW_' + name, getattr(pm, metric).index)
        values = df.values
        return values
    
    color = color_wheel.next()
    
    def plot(arr, position):
        if arr.ndim == 2:
            N = arr.shape[1]
            width = 0.618 / N
            boxwidth = 0.618 / (N + 2)
            plots = []
            for i in range(N):
                boxplot = bst.plots.plot_montecarlo(
                    data=arr[:, i], positions=[position + (i-(N-1)/2)*width], 
                    light_color=color.RGBn, 
                    dark_color=color.shade(60).RGBn,
                    width=boxwidth,
                    hatch=getattr(color, 'hatch', None),
                )
                plots.append(boxplot)
            return plots
        else:
            return bst.plots.plot_montecarlo(
                data=arr, positions=[position], 
                light_color=color.RGBn, 
                dark_color=color.shade(60).RGBn,
                width=0.618,
            )
    
    data = np.zeros([N_rows, N_cols], dtype=object)
    data[:] = [[get_data(i, j) for j in columns] for i in rows]
    if tickmarks is None: 
        tickmarks = [
            bst.plots.rounded_tickmarks_from_data(
                i, N_ticks=6, lb_max=None if any([(j < 0).any() for j in i]) else 0,
                f=lambda x: roundsigfigs(x, sigfigs=2, index=1),
                f_min=lambda x: np.min(x),
                f_max=lambda x: np.max(x),
                expand=expand,
            ) 
            for i in data
        ]
        tickmarks = [roundsigfigs(i, sigfigs=3, index=1) for i in tickmarks]
        

    xf = len(columns) - 0.5
    for i in range(N_rows):
        ax = axes[i]
        plt.sca(ax)
        plt.xlim(-0.5, xf)
        
    for j in range(N_cols):
        color_wheel.restart()
        for i in range(N_rows):
            ax = axes[i]
            plt.sca(ax)
            plot(data[i, j], j)
            plt.ylabel(ylabels[i])
    
    for i in range(N_rows):
        ax = axes[i]
        plt.sca(ax)
        yticks = tickmarks[i]
        plt.ylim([yticks[0], yticks[1]])
        if yticks[0] < 0.:
            bst.plots.plot_horizontal_line(0, color=CABBI_colors.black.RGBn, lw=0.8, linestyle='--')
        try:
            xticklabels = xtext if (ax in axes_box[-1] or i == N_rows - 1) else []
        except:
            xticklabels = xtext if i == N_rows - 1 else []
        bst.plots.style_axis(ax,  
            xticks = xticks,
            yticks = yticks,
            xticklabels= xticklabels, 
            ytick0=True,
            ytickf=False,
            offset_xticks=False,
            xrot=xrot,
        )
    for i in range(N_rows, nrows * ncols):
        ax = axes[i]
        plt.sca(ax)
        plt.axis('off')
    if fig is None:
        fig = plt.gcf()
    else:
        plt.subplots_adjust(hspace=0)
    fig.align_ylabels(axes)
    return fig, axes

def plot_spearman_both(scenario, ProcessModel=None, **kwargs):
    set_font(size=11)
    set_figure_size(aspect_ratio=0.8)
    labels = ['TEA', 'LCA']
    if ProcessModel is None: ProcessModel = sp.STRAPMSWProcess
    pm = ProcessModel(simulate=False, scenario=scenario, **kwargs)
    rhos = []
    file = spearman_file('STRAPMSW_' + scenario)
    df = pd.read_excel(file, header=[0, 1], index_col=[0, 1])
    names = get_spearman_names(pm.model.parameters)
    names = [names[i].replace('Polymer', 'PEPP') for i in df.index]
    metric_names = []
    for label in labels:
        if label == 'TEA':
            metric = pm.MSP
            metric_name = metric.name
            values = df[metric.index]
            for i in pm.model.parameters:
                if i not in pm.general_parameters: 
                    values[i.index] = 0
        elif label == 'LCA':
            metric = pm.GWP_ethanol
            metric_name = r'GWP$_{\mathrm{energy}}$'
            values = df[metric.index]
            for i in pm.model.parameters:
                if i not in pm.general_parameters: 
                    values[i.index] = 0
        else:
            raise ValueError(f"invalid label '{label}'")
        rhos.append(values)
        metric_names.append(metric_name)
    color_wheel = [Color(fg='#0c72b9'), Color(fg='#d34249')]
    fig, ax = bst.plots.plot_spearman_2d(rhos, index=names,
                                         color_wheel=color_wheel,
                                         name=metric_name,
                                         cutoff=0.02,
                                         xlabel="Spearman's rank correlation coefficient",
                                         **kwargs)
    legend_kwargs = {'loc': 'lower left'}
    plt.legend(
        handles=[
            mpatches.Patch(
                color=color_wheel[i].RGBn, 
                label=metric_names[i],
            )
            for i in range(len(labels))
        ], 
        **legend_kwargs,
    )
    plt.subplots_adjust(
        hspace=0.05, wspace=0.05,
        top=0.98, bottom=0.15,
        left=0.45, right=0.9,
    )
    for i in ('svg', 'png'):
        name = f'spearman.{i}'
        file = os.path.join(images_folder, name)
        plt.savefig(file, dpi=900, transparent=True)
    return fig, ax

def plot_spearman(kind=None, ProcessModel=None, scenario='all', 
                  MSP=True,
                  **kwargs):
    set_font(size=10)
    if kind is None: kind = 'TEA'
    if ProcessModel is None: ProcessModel = sp.STRAPMSWProcess
    pm = ProcessModel(simulate=False, scenario=scenario)
    if kind == 'TEA':
        if MSP:
            metric = pm.MSP
            set_figure_size(aspect_ratio=1, width=6.6142)
            top = 12
        else:
            metric = pm.IRR
            set_figure_size(aspect_ratio=1.5, width=6.6142 * 2.5/4)
            top = 12
        metric_name = metric.name
    elif kind == 'LCA':
        metric = pm.GWP_ethanol
        metric_name = r'carbon intensity'
        set_figure_size(aspect_ratio=0.42)
        top = 6
    else:
        raise ValueError(f"invalid kind '{kind}'")
    rhos = []
    file = spearman_file('STRAPMSW_' + scenario)
    df = pd.read_excel(file, header=[0, 1], index_col=[0, 1])
    rhos = df[metric.index]
    names = get_spearman_names(pm.model.parameters)
    bad_name = 'extracted polymer fraction'
    index = [names[i] for i in rhos.index]
    index = [
        ('Extracted polymer fraction\n[66.6 - 86.6 wt % total plastic]' if bad_name in i else i.replace('cellulase cellulase', 'Enzyme').replace('cellulase', 'Enzyme').replace('Cofermenation s', 'S').replace('resin', 'Resin').replace('MSW', 'MRF residue'))
        for i in index
    ]
    color_wheel = [GG_colors.blue, GG_colors.orange]
    fig, ax = bst.plots.plot_spearman_2d([rhos], index=index,
                                         color_wheel=color_wheel,
                                         w=1.0,
                                         xlabel=f"Spearman's $\\rho$ with {metric_name}",
                                         top=top,
                                         **kwargs)
    left = 0.4
    bottom = 0.08
    if not MSP and kind == 'TEA':
        left = 0.6
    if kind == 'LCA':
        bottom = 0.15
    plt.subplots_adjust(
        hspace=0.05, wspace=0.05,
        top=0.98, bottom=bottom,
        left=left, right=0.95,
    )
    for i in ('svg', 'png'):
        name = f'spearman_{kind}.{i}'
        file = os.path.join(images_folder, name)
        plt.savefig(file, dpi=900, transparent=True)
    return fig, ax

def plot_kde_STRAP_CUWP_single_wash_conceptual_configurations():
    set_font(size=10)
    set_figure_size(width='half', aspect_ratio=1.1)
    scenarios = (
        sp.Scenario(
            solvent='Toluene',
            target_plastic='PE',
            target_plastic_percent=50, # %
            processing_capacity=5e+03, # MT-plastic/yr
            sell_leftover_plastic=True,
            burn_leftover_plastic=False,
            facilities=True,
        ),
        sp.Scenario(
            solvent='Toluene',
            target_plastic='PE',
            target_plastic_percent=50, # %
            processing_capacity=5e+03, # MT-plastic/yr
            sell_leftover_plastic=False,
            burn_leftover_plastic=True,
            facilities=True,
        ),
        sp.Scenario(
            solvent='Toluene',
            target_plastic='PE',
            target_plastic_percent=50, # %
            processing_capacity=5e+03, # MT-plastic/yr
            sell_leftover_plastic=False,
            burn_leftover_plastic=False,
            facilities=True,
        ),
    )
    # axes = None
    # unlabeled = True
    for scenario in scenarios:
        pm = sp.BaselineSTRAPProcess(simulate=False, scenario=scenario)
        metrics = [pm.GWP_ethanol.index, pm.MSP.index]
        Xi, Yi = metrics
        df = get_monte_carlo(pm.name, metrics)
        y = df[Yi].values
        x = df[Xi].values
        # TODO: figure out why GWP is so high
        yticks = [0, 2, 4, 6, 8]
        xticks = [0, 1, 2, 3, 4]
        fig, ax, axes = bst.plots.plot_kde(
            y=y, x=x, xticks=xticks, yticks=yticks,
            xticklabels=True, yticklabels=True,
            xbox_kwargs=dict(light=CABBI_colors.orange.RGBn, dark=CABBI_colors.orange.shade(60).RGBn),
            ybox_kwargs=dict(light=CABBI_colors.blue.RGBn, dark=CABBI_colors.blue.shade(60).RGBn),
            xbox_width=400,
            aspect_ratio=1.1,
            # axes=axes,
        )
        # if unlabeled:
        plt.sca(ax)
        plt.ylabel(r'MSP $[\mathrm{USD} \cdot \mathrm{kg}^{\mathrm{-1}}]$')
        plt.xlabel(r'GWP $[\mathrm{kgCO2e} \cdot \mathrm{L}^{\mathrm{-1}}]$')
        # unlabeled = False
        plt.subplots_adjust(
            hspace=0.05, wspace=0.05,
            top=0.98, bottom=0.15,
            left=0.15, right=0.98,
        )
        for i in ('svg', 'png'):
            file = os.path.join(images_folder, f'{pm.name}_kde.{i}')
            plt.savefig(file, dpi=900, transparent=True)

def run_monte_carlo_STRAPMSW(scenarios=None):
    if isinstance(scenarios, str): scenarios = [scenarios]
    if scenarios is None: scenarios = ['all', 'baseline', 'potential'] 
    for i in scenarios:
        run_monte_carlo(
            10000, ProcessModel=sp.STRAPMSWProcess,
            scenario=i
        )

def run_monte_carlo_STRAP_CUWP_single_wash_conceptual_configurations():
    scenarios = (
        sp.Scenario(
            solvent='Toluene',
            target_plastic='PE',
            target_plastic_percent=50, # %
            processing_capacity=5e+03, # MT-plastic/yr
            sell_leftover_plastic=True,
            burn_leftover_plastic=False,
            facilities=True,
        ),
        sp.Scenario(
            solvent='Toluene',
            target_plastic='PE',
            target_plastic_percent=50, # %
            processing_capacity=5e+03, # MT-plastic/yr
            sell_leftover_plastic=False,
            burn_leftover_plastic=True,
            facilities=True,
        ),
        sp.Scenario(
            solvent='Toluene',
            target_plastic='PE',
            target_plastic_percent=50, # %
            processing_capacity=5e+03, # MT-plastic/yr
            sell_leftover_plastic=False,
            burn_leftover_plastic=False,
            facilities=True,
        ),
    )
    for i in scenarios:
        run_monte_carlo(
            1000, ProcessModel=sp.BaselineSTRAPProcess,
            scenario=i
        )

def run_monte_carlo(
        N, rule='L',
        sample_cache={},
        autosave=True,
        autoload=True,
        ProcessModel=None,
        scenario=None,
        **kwargs
    ):
    filterwarnings('ignore', category=bst.exceptions.DesignWarning)
    filterwarnings('ignore', category=bst.exceptions.CostWarning)
    if ProcessModel is None: ProcessModel = sp.STRAPMSWProcess
    pm = ProcessModel(simulate=False, scenario=scenario, **kwargs)
    pm.model.exception_hook = 'warn'
    filename = 'STRAPMSW_' + scenario
    N_notify = min(int(N/10), 20)
    autosave = N_notify if autosave else False
    autoload_file = autoload_file_name(filename)
    np.random.seed(1)
    samples = pm.model.sample(N, rule)
    pm.model.load_samples(samples)
    pm.model.evaluate(
        notify=int(N/10),
        autosave=autosave,
        autoload=autoload,
        file=autoload_file,
    )
    file = monte_carlo_file(filename)
    pm.model.table.to_excel(file)
    pm.model.table = pm.model.table.dropna(how='all', axis=1)
    for i in pm.model.metrics:
        if i.index not in pm.model.table: pm.model._metrics.remove(i)
    pm.model.table = pm.model.table.dropna(how='any', axis=0)
    rho, p = pm.model.spearman_r(filter='omit nan')
    file = spearman_file(filename)
    rho.to_excel(file)

def run_monte_carlo_across_process_capacity(
        scenario='baseline', N=500, rule='L',
    ):
    filterwarnings('ignore')
    ProcessModel = sp.STRAPMSWProcess
    pm = ProcessModel(
        simulate=False, 
        scenario=scenario,
        preprocessing=False,
    )
    pm.model.exception_hook = 'raise'
    filename = f'STRAPMSW_{scenario}_monte_carlo_across_processing_capacity.xlsx'
    coordinate = np.linspace(
        30,
        150, 
        12
    )
    pm.set_processing_capacity.units = '10^3 MT/yr'
    def setter(processing_capacity):
        feedstock = pm.feedstock
        pm.system.rescale(
            feedstock, 
            processing_capacity * 1e6 / pm.tea.operating_hours / feedstock.F_mass
        )
        
    pm.set_processing_capacity.setter = setter
    file = os.path.join(results_folder, filename)
    np.random.seed(1)
    samples = pm.model.sample(N, rule)
    pm.model.load_samples(samples)
    pm.model.evaluate_across_coordinate(
        name='Processing capacity',
        notify=int(N/10),
        f_coordinate=pm.set_processing_capacity,
        coordinate=coordinate,
        notify_coordinate=True,
        xlfile=file,
    )

def plot_monte_carlo_across_processing_capacity_IRR(uncertainty=True, full=False):
    bst.plots.set_font(size=10, family='sans-serif', font='Arial')
    if full:
        bst.plots.set_figure_size(aspect_ratio=0.55, width='full')
    else:
        bst.plots.set_figure_size(aspect_ratio=1.1, width='half')
    # colors = (Color(fg='#8172B3'), Color(fg='#55A868'))
    colors = (GG_colors.purple, GG_colors.green)
    scenarios = ('baseline', 'potential')
    for scenario, color in zip(scenarios, colors):
        filename = f'STRAPMSW_{scenario}_monte_carlo_across_processing_capacity.xlsx'
        file = os.path.join(results_folder, filename)
        ProcessModel = sp.STRAPMSWProcess
        pm = ProcessModel(
            simulate=False, 
            scenario=scenario,
        )
        df = pd.read_excel(file, sheet_name=pm.IRR.short_description, index_col=0)
        df = df.dropna()
        processing_capacity = np.array(df.columns)
        if uncertainty:
            IRR = bst.plots.plot_montecarlo_across_coordinate(
                processing_capacity, df, 
                fill_color=[*color.RGBn, 0.5],
                median_color=color.shade(20).RGBn,
                p5_color=[*color.shade(20).RGBn, 0.80],
                smooth=2,
            )
        else:
            def f(processing_capacity):
                feedstock = pm.feedstock
                pm.system.rescale(
                    feedstock, 
                    processing_capacity * 1e6 / pm.tea.operating_hours / feedstock.F_mass
                )
                pm.system.simulate()
                return pm.IRR()
            IRR = [f(i) for i in processing_capacity]
            plt.plot(processing_capacity, IRR,
                c=color.RGBn, ls='--' if scenario == 'potential' else '-',
            )
    if full:
        xticks = [30, 45, 60, 75, 90, 105, 120, 135, 150]
    else:
        xticks = [30, 60, 90, 120, 150]
    plt.xlim([30, 150])
    plt.ylim([0, 20])
    plt.xlabel(f"Processing capacity [{format_units('10^3*MT/yr')}]", fontsize=10)
    plt.ylabel("IRR [%]", fontsize=10)
    bst.plots.style_axis(
        fs=10,
        xticks=xticks,
        ytick0=False,
        ytickf=True,
        # right=False,
        # top=False
    )
    # ax = plt.gca()
    # ax.spines[['right', 'top']].set_visible(False)
    plt.subplots_adjust(hspace=0.05, left=0.2, right=0.92, bottom=0.2, top=0.90)
    for i in ('svg', 'png'):
        file = os.path.join(images_folder, f'IRR_monte_carlo_across_processing_capacity.{i}')
        plt.savefig(file, dpi=900, transparent=True)

def plot_monte_carlo_across_processing_capacity(
        scenario='baseline', verification=False,
    ):
    bst.plots.set_font(size=10, family='sans-serif', font='Arial')
    bst.plots.set_figure_size(aspect_ratio=0.5)
    if verification:
        filename = f'verification_STRAPMSW_{scenario}_monte_carlo_across_processing_capacity.xlsx'
    else:
        filename = f'STRAPMSW_{scenario}_monte_carlo_across_processing_capacity.xlsx'
    file = os.path.join(results_folder, filename)
    ProcessModel = sp.STRAPMSWProcess
    pm = ProcessModel(
        simulate=False, 
        scenario=scenario,
    )
    df = pd.read_excel(file, sheet_name=pm.MSP.short_description, index_col=0)
    df = df.dropna()
    if verification:
        processing_capacity = np.array(df.columns)
    else:
        processing_capacity = 907.185e-6 * np.array(df.columns)
    median_purple = GG_colors.purple.shade(40).RGBn
    median_green = GG_colors.green.shade(40).RGBn
    MSPs = bst.plots.plot_montecarlo_across_coordinate(
        processing_capacity, df, 
        fill_color=[*GG_colors.purple.RGBn, 0.5],
        median_color=median_purple,
        p5_color=[*GG_colors.purple.shade(20).RGBn, 0.80],
        smooth=2,
    )
    market_price_range = np.array(pm.resin_price_range) / 907.185
    if verification:
        xticks = [50, 75, 100, 125, 150]
        plt.xlim([50, 150])
        xtickmarks = [str(i) for i in xticks]
        if scenario == 'NREL':
            lb = -0.5
            ub = 2
        else:
            lb = 0
            ub = 3.5
        plt.ylim([lb, ub])
    else:
        xticks = [28.500, 100, 200, 300, 400, 500, 600, 699]
        plt.xlim([907.185e-6 * pm.landfill_sizes['small'], 907.185e-6 * pm.landfill_sizes['large']])
        xtickmarks = [str(i) for i in xticks]
        xtickmarks[0] = ''
        xtickmarks[-1] = ''
        lb = -1
        ub = 2.5
        plt.ylim([lb, ub])
    plt.xlabel(f"Processing capacity [{format_units('10^3*MT/yr')}]")
    plt.ylabel(f"MSP [{format_units('USD/kg')}]")
    # bst.plots.style_axis(
    #     xticks=[0, 3, 6, 9, 12, 15],
    #     ytick0=False,
    #     ytickf=True,
    # )
    plt.subplots_adjust(hspace=0.05, left=0.2, right=0.96, bottom=0.15, top=0.95)
    if verification:
        TCI_ub = 350
        plt.xlabel(f"Processing capacity [{format_units('10^3*MT/yr')}]")
    else:
        TCI_ub = 1600
        plt.xlabel(f"Processing capacity [{format_units('10^3*ton/yr')}]")
    plt.ylabel(f"MSP [{format_units('USD/kg')}]", color=median_purple)
    ax = plt.gca()
    ax.tick_params(axis='y', color=median_purple, labelcolor=median_purple)
    ax.spines['left'].set_color(median_purple)
    ax.spines['right'].set_color(median_green)
    new_axis = ax.twinx() 
    plt.sca(new_axis)
    df = pd.read_excel(file, sheet_name=pm.TCI.short_description, index_col=0)
    df = df.dropna()
    TCI = bst.plots.plot_montecarlo_across_coordinate(
        processing_capacity, df, 
        fill_color=[*GG_colors.green.RGBn, 0.5],
        median_color=median_green,
        p5_color=[*GG_colors.green.shade(20).RGBn, 0.80],
        smooth=0,
    )
    
    price_range = (market_price_range - lb) / (ub - lb) * TCI_ub
    
    plt.fill_between(plt.xlim(), *price_range,
                     color=[*line_color, 0.5],
                     linewidth=2,
                     zorder=-1)
    
    # bst.plots.plot_horizontal_line(pos, lw=2, ls='-', color=line_color, zorder=-1)
    if not verification:
        bst.plots.plot_vertical_line(pm.landfill_sizes['medium'] * 907.185e-6, lw=2, ls='-', color=vertical_line_color, zorder=-1)
    new_axis.tick_params(axis='y', color=median_green, labelcolor=median_green)
    plt.ylim([0, TCI_ub])
    plt.ylabel(f"TCI [{format_units('10^6*USD')}]", color=median_green)
    plt.gcf().tight_layout()
    new_axis.spines['left'].set_color(median_purple)
    new_axis.spines['right'].set_color(median_green)
    plt.sca(ax)
    ax.tick_params(axis='y', left=True, direction="inout", length=4)
    ax.tick_params(axis='x', bottom=True, direction="inout", length=4)
    new_axis.tick_params(axis='y', right=True, direction="inout", length=4)
    y_twin = ax.twiny()
    y_twin.tick_params(axis='x', top=True, direction="in", length=4)
    y_twin.zorder = 2
    # breakpoint()
    plt.xticks(xticks, len(xticks) * [''])
    plt.xlim([xticks[0], xticks[-1]])
    y_twin.spines['left'].set_color(median_purple)
    y_twin.spines['right'].set_color(median_green)
    plt.subplots_adjust(hspace=0.05, left=0.1, right=0.90, bottom=0.15, top=0.95)

    for i in ('svg', 'png'):
        if verification:
            file = os.path.join(images_folder, f'verification_{scenario}_monte_carlo_across_processing_capacity.{i}')
        else:
            file = os.path.join(images_folder, f'monte_carlo_{scenario}_across_processing_capacity.{i}')
        plt.savefig(file, dpi=900, transparent=True)

def run_monte_carlo_tipping_fee(
        scenarios=['experimental', 'NREL'], N=500, rule='L',
    ):
    filterwarnings('ignore')
    ProcessModel = sp.STRAPMSWProcess
    for scenario in scenarios:
        pm = ProcessModel(
            simulate=False, 
            scenario=scenario,
        )
        pm.model.exception_hook = 'raise'
        filename = f'STRAPMSW_{scenario}_monte_carlo_across_tipping_fee.xlsx'
        file = os.path.join(results_folder, filename)
        np.random.seed(1)
        samples = pm.model.sample(N, rule)
        pm.model.load_samples(samples)
        lb = pm.get_tipping_fee('MT')[0]
        ub = pm.get_tipping_fee('MA')[0]
        pm.model.evaluate_across_coordinate(
            name='Tipping fee',
            notify=int(N/10),
            f_coordinate=pm.set_MSW_tipping_fee,
            coordinate=np.linspace(lb, ub, 12),
            notify_coordinate=True,
            simulation_independent_coordinate=True,
            xlfile=file,
        )

def plot_monte_carlo_across_tipping_fee(
        scenarios=['experimental', 'NREL'],
    ):
    bst.plots.set_font(size=10, family='sans-serif', font='Arial')
    bst.plots.set_figure_size(aspect_ratio=0.5)
    colors = [GG_colors.blue, GG_colors.orange]
    for scenario, color in zip(scenarios, colors):
        filename = f'STRAPMSW_{scenario}_monte_carlo_across_tipping_fee.xlsx'
        file = os.path.join(results_folder, filename)
        ProcessModel = sp.STRAPMSWProcess
        pm = ProcessModel(
            simulate=False, 
            scenario=scenario,
        )
        df = pd.read_excel(file, sheet_name=pm.MSP.short_description, index_col=0)
        df = df.dropna()
        tipping_fees = np.array(df.columns)
        MSPs = bst.plots.plot_montecarlo_across_coordinate(
            tipping_fees, df, 
            fill_color=[*color.RGBn, 0.5],
            median_color=color.shade(40).RGBn,
            p5_color=[*color.shade(20).RGBn, 0.80],
            smooth=1,
        )
    lb = pm.get_tipping_fee('MT')[0]
    ub = pm.get_tipping_fee('MA')[0]
    xticks = [lb, 40, 50, 60, 70, 80, 90, ub]
    market_price = pm.baseline_resin_price / 907.185
    bst.plots.plot_horizontal_line(market_price, lw=2, ls='-', color=line_color, zorder=-1)
    bst.plots.plot_vertical_line(pm.get_tipping_fee('nationwide')[0], lw=2, ls='-', color=vertical_line_color, zorder=-1)
    plt.ylim([-1.5, 2.5])
    plt.xlabel(f"Tipping fee [{format_units('USD/ton')}]")
    plt.ylabel(f"MSP [{format_units('USD/kg')}]")
    plt.xlim([xticks[0], xticks[-1]])
    bst.plots.style_axis(
        xticks=xticks,
        xtick0=False,
        xtickf=False,
    )
    plt.subplots_adjust(hspace=0.05, left=0.1, right=0.90, bottom=0.15, top=0.95)
    for i in ('svg', 'png'):
        file = os.path.join(images_folder, f'monte_carlo_across_tipping_fee.{i}')
        plt.savefig(file, dpi=900, transparent=True)

def sobol_analysis(N=256):
    from SALib.analyze import sobol
    filterwarnings('ignore', category=bst.exceptions.DesignWarning)
    filterwarnings('ignore', category=bst.exceptions.CostWarning)
    pm = sp.STRAPMSWProcess(simulate=False)
    pm.model.exception_hook = 'raise'
    samples = pm.Sobol_model.sample(N=N, rule='sobol', seed=0)
    pm.Sobol_model.load_samples(samples)
    pm.Sobol_model.evaluate(
        notify=int(len(samples)/10),
    )
    problem = pm.Sobol_model.problem()
    for metric in (pm.MSP, pm.GWP):
        Y = pm.Sobol_model.table[metric.index].values
        results = sobol.analyze(problem, Y)
        results = {i: j.tolist() for i, j in results.items()}
        file = sobol_file('_'.join([pm.name, metric.name]), extension='yaml')
        with open(file, 'w') as file:
            yaml.dump(results, file)
    
def plot_sobol(metric=None):
    # TODO: plot as stacked bar plot
    # plot_stacked_bar
    filterwarnings('ignore', category=bst.exceptions.DesignWarning)
    filterwarnings('ignore', category=bst.exceptions.CostWarning)
    pm = sp.STRAPMSWProcess(simulate=False)
    if metric is None: metric = 'MSP'
    metric = getattr(pm, metric)
    file = sobol_file('_'.join([pm.name, metric.name]), extension='yaml')
    color = GG_colors.blue.RGBn
    with open(file, 'r') as stream: 
        data = yaml.full_load(stream)
        S1 = 100 * np.array(data['S1'])
        S1_conf = 100 * np.array(data['S1_conf'])
    index, = np.where(S1 > 1)
    parameters = [pm.Sobol_model.parameters[i] for i in index]
    S1 = S1[index]
    S1_conf = S1_conf[index]
    y = [i for i in range(S1.size)]
    plt.errorbar(S1, y, xerr=S1_conf, c=color, marker='s', linestyle='')
    plt.yticks(y, [i.short_description for i in parameters])
    plt.xlim(0, 100)
    plt.xlabel('Share of total MSP variance [%]')
    plt.show()

def MSP_at_processing_capacity_tipping_fee(
        processing_capacity, tipping_fee, process_models,
    ):
    values = []
    for process_model in process_models:
        process_model.set_MSW_tipping_fee.setter(tipping_fee * 0.907185)
        process_model.set_processing_capacity.setter(processing_capacity * 1e3 / 0.907185)
        process_model.system.simulate()
        values.append(process_model.IRR())
    return np.array(values)

def plot_IRR_across_processing_capacity_tipping_fee(load=True):
    bst.plots.set_font(size=10, family='sans-serif', font='Arial')
    bst.plots.set_figure_size(aspect_ratio=0.55)
    xlim = np.array([30, 130])
    ylim = np.array([25, 150])
    scenarios = ['baseline', 'potential']
    process_models = [sp.STRAPMSWProcess(simulate=True, scenario=i) for i in scenarios]
    X, Y, Z = bst.plots.generate_contour_data(
        MSP_at_processing_capacity_tipping_fee,
        file=os.path.join(results_folder, 'IRR_across_processing_capacity_tipping_fee.npy'),
        load=load, save=True,
        xlim=xlim, ylim=ylim,
        args=(process_models,),
        n=12,
    )
    # Plot contours
    ylabel = r"Tipping fee [USD $\cdot \mathrm{MT}^{\mathrm{-1}}$]"
    yticks = [30, 50, 70, 90, 110, 130, 150]
    xlabel = r"Processing capacity [$10^3 \cdot$ MT $\cdot \mathrm{yr}^{\mathrm{-1}}$]"
    xticks = [30, 55, 80, 105, 130]
    metric_bars = [
        bst.plots.MetricBar(
            'IRR [%]', ' ', plt.cm.get_cmap('viridis_r'), 
            bst.plots.rounded_tickmarks_from_data(Z, 5, 1, expand=0, p=0.1), 
            10, 1, shrink=1, pad=0.05, forced_size=1,
        ),
    ]
    fig, axes, CSs, CB, other_axes = bst.plots.plot_contour_2d(
        X, Y, Z, '', '', xticks, yticks, metric_bars, 
        contour_label_interval=1,
        fillcolor=None, styleaxiskw=dict(xtick0=True), label=True,
    )
    plt.subplots_adjust(left=0.1, right=0.9, wspace=0.18, hspace=0.2, top=0.9, bottom=0.16)
    fig.text(0, 0.5, ylabel, va='center', rotation='vertical')
    fig.text(0.5, 0.04, xlabel, ha='center')
    for i in ('svg', 'png'):
        file = os.path.join(images_folder, f'IRR_across_processing_capacity_tipping_fee.{i}')
        plt.savefig(file, dpi=900, transparent=True)

def MSP_GWP_at_titer_glucose_yield(
        titer, glucose_yield, process_model, convergence_model, productivities,
    ):
    process_model.set_glucose_yield.setter(glucose_yield)
    process_model.set_cofermentation_ethanol_titer.setter(titer)
    values = []
    for productivity in productivities:
        process_model.set_cofermentation_ethanol_productivity.setter(productivity)
        with convergence_model.practice([titer, glucose_yield, productivity]):
            process_model.system.simulate()
        values.append([process_model.MSP(), process_model.GWP()])
    return np.array(values).T

def plot_MSP_GWP_across_titer_glucose_yield(load=True):
    bst.plots.set_font(size=10, family='sans-serif', font='Arial')
    bst.plots.set_figure_size(aspect_ratio=0.68)
    pm = sp.STRAPMSWProcess(simulate=False)
    xlim = np.array([9, 54])
    ylim = np.array([55, 90])
    productivities = np.array([0.1722, 1.5])
    params = (pm.set_cofermentation_ethanol_titer, pm.set_glucose_yield, pm.set_cofermentation_ethanol_productivity)
    convergence_model = bst.ConvergenceModel(predictors=params)
    X, Y, Z = bst.plots.generate_contour_data(
        MSP_GWP_at_titer_glucose_yield,
        file=os.path.join(results_folder, 'MSP_dissolution_capacity_boiling_point.npy'),
        load=load, save=True,
        xlim=xlim, ylim=ylim,
        args=(pm, convergence_model, productivities),
        n=12,
    )
    Z[:, :, 1] *= 1000
    # Plot contours
    ylabel = "Hydrolysis glucose yield [% theoretical]"
    xlabel = r"Ethanol titer [ g $\cdot \mathrm{L}^{\mathrm{-1}}$]"
    yticks = [55, 60, 65, 70, 75, 80, 85, 90]
    xticks = [9, 18, 27, 36, 45, 54]
    metric_bars = [
        bst.plots.MetricBar(
            'MSP', r'$[\mathrm{USD} \cdot \mathrm{kg}^{\mathrm{-1}}]$', plt.cm.get_cmap('viridis_r'), 
            bst.plots.rounded_tickmarks_from_data(Z[:, :, 0], 5, 0.02, expand=0, p=0.02), 
            13, 2
        ),
        bst.plots.MetricBar(
            'GWP', r'$[\mathrm{gCO2e} \cdot \mathrm{L}^{\mathrm{-1}}]$', plt.cm.get_cmap('inferno_r'), 
            bst.plots.rounded_tickmarks_from_data(Z[:, :, 1], 5, 1, expand=0, p=1), 
            13, 0
        )
    ]
    fig, axes, CSs, CB, other_axes = bst.plots.plot_contour_2d(
        X, Y, Z, xlabel, '', xticks, yticks, metric_bars,  
        fillcolor=None, styleaxiskw=dict(xtick0=False), label=True,
        titles=['productivity [g/L/h]: 0.1722', '1.5'],
    )
    plt.subplots_adjust(left=0.1, right=0.9, wspace=0.15, hspace=0.2, top=0.9, bottom=0.12)
    fig.text(0, 0.5, ylabel, va='center', rotation='vertical')
    for i in ('svg', 'png'):
        file = os.path.join(images_folder, f'MSW_titer_yield_contours.{i}')
        plt.savefig(file, dpi=900, transparent=True)

def get_monte_carlo(name, index, cache={}):
    key = name
    if key in cache:
        df = cache[key]
    else:
        file = monte_carlo_file(key)
        cache[key] = df = pd.read_excel(file, header=[0, 1], index_col=[0])
    df = df[index]
    mc = df.dropna(how='all', axis=0)
    return mc
    
def plot_breakdowns(scenario='baseline', **kwargs):
    # plot_stacked_bar
    bst.plots.set_font(size=10)
    bst.plots.set_figure_size(aspect_ratio=0.68)
    biorefinery = sp.STRAPMSWProcess(scenario=scenario, **kwargs, simulate=True)
    unit_groups = biorefinery.unit_groups
    colors = [biorefinery.area_colors[i.name].RGBn for i in unit_groups]
    hatches = [biorefinery.area_hatches[i.name] for i in unit_groups]
    units = sum([i.units for i in unit_groups], [])
    joint_group = bst.UnitGroup(None, units)
    production = biorefinery.product.F_mass / 1000
    for i in (*unit_groups, joint_group):
        if i.metrics: continue
        i.metric(i.get_installed_cost,
                 INST_EQ_COST, CAPITAL_UNITS)
        i.metric(functional_unit(i.get_cooling_duty, production),
                 COOLING, DUTY_UNITS)
        i.metric(functional_unit(i.get_heating_duty, production),
                 HEATING, DUTY_UNITS)
        i.metric(functional_unit(i.get_electricity_consumption, production),
                 ELEC_CONS, ELEC_UNITS)
    
    def format_total(x):
        if x < 1e3:
            return format(x, '.3g')
        else:
            x = int(x)
            n = 10 ** (len(str(x)) - 3)
            value = int(round(x / n) * n)
            return format(value, ',')
    
    fig, axes = bst.plots.plot_unit_groups(
        unit_groups,
        colors=colors,
        hatches=hatches,
        format_total=format_total,
        fraction=True,
        joint_group=joint_group,
        legend_kwargs=dict(
            loc='upper right',
            ncol=1,
            bbox_to_anchor=(1.5, 1),
            labelspacing=1.5, handlelength=2.8,
            handleheight=1, scale=0.8,
        ),
    )
    plt.subplots_adjust(left=0.12, right=0.7, wspace=0.1, top=0.9, bottom=0.10)
    for i in ('svg', 'png'):
        name = f'breakdowns_{scenario}.{i}'
        file = os.path.join(images_folder, name)
        plt.savefig(file, dpi=900, transparent=True)

def plot_kde(scenario):
    set_font(size=10)
    set_figure_size(width='half', aspect_ratio=1.1)
    pm = sp.STRAPMSWProcess(simulate=False, scenario=scenario)
    metrics = [pm.GWP_ethanol.index, pm.IRR.index]
    # metrics = [pm.GWP.index, pm.IRR.index]
    Xi, Yi = metrics
    df = get_monte_carlo('STRAPMSW_' + scenario, metrics)
    y = df[Yi].values
    x = 1000 * df[Xi].values
    yticks = [0, 5, 10, 15, 20, 25, 30]
    xticks = [0, 60, 120, 180, 240, 300, 360]
    fig, ax, axes = bst.plots.plot_kde(
        y=y, x=x, xticks=xticks, yticks=yticks,
        xticklabels=True, yticklabels=True,
        xbox_kwargs=dict(light=CABBI_colors.orange.RGBn, dark=CABBI_colors.orange.shade(60).RGBn),
        ybox_kwargs=dict(light=CABBI_colors.blue.RGBn, dark=CABBI_colors.blue.shade(60).RGBn),
        xbox_width=400,
        aspect_ratio=1.1,
    )
    ax.set_clip_on(False)
    plt.sca(ax)
    xlb, xub = plt.xlim()
    ylb, yub = plt.ylim()
    xpos = lambda x: xlb + (xub - xlb) * x
    ypos = lambda y: ylb + (yub - ylb) * y
    # xleft = 0.02
    xright = 0.86
    # ytop = 0.94
    ybottom = 0.02
    first = True
    if first: description = ' under market price'
    else: description = ''
    plt.sca(ax)
    # price_range = np.array(pm.resin_price_range) / 907.185
    # plt.fill_between(plt.xlim(), *price_range,
    #                  color=[*line_color, 0.5],
    #                  linewidth=2,
    #                  zorder=-1)
    cellulosic_ethanol_GWP = 297.2 # GREET 2023
    bst.plots.plot_vertical_line(cellulosic_ethanol_GWP, lw=2, ls='--', color=line_color, zorder=-1)
    plt.ylabel(r'MSP $[\mathrm{USD} \cdot \mathrm{kg}^{\mathrm{-1}}]$')
    plt.xlabel(r'Carbon intensity $[\mathrm{gCO2e} \cdot \mathrm{L}^{\mathrm{-1}}]$')
    plt.subplots_adjust(
        hspace=0.05, wspace=0.05,
        top=0.98, bottom=0.15,
        left=0.2, right=0.98,
    )
    for i in ('svg', 'png'):
        file = os.path.join(images_folder, f'{scenario}_kde.{i}')
        plt.savefig(file, dpi=900, transparent=True)
    
def plot_kde_CI_MSP(scenarios=None):
    from warnings import filterwarnings
    filterwarnings('ignore')
    set_font(size=10)
    set_figure_size(width='half', aspect_ratio=1.1)
    scenarios = ['baseline', 'potential']
    pm = sp.STRAPMSWProcess(simulate=False, scenario=scenarios[0])
    price_range = np.array(pm.resin_price_range) / 907.185
    scenarios_list = []
    metrics = [pm.GWP_ethanol.index, pm.MSP.index]
    Xi, Yi = metrics
    y1d = []
    x1d = []
    opportunity_space_1d = []
    # opportunity_space_2d.append(opportunity_space_1d)
    for scenario in scenarios:
        pm = sp.STRAPMSWProcess(simulate=False, scenario=scenario)
        df = get_monte_carlo('STRAPMSW_' + scenario, metrics)
        y = df[Yi].values
        x = df[Xi].values * 1e3
        # print(scenario, 'y', y.min(), y.max())
        y1d.append(y)
        x1d.append(x)
        # print(scenario, 'x', x.min(), x.max())
        under_within = (y <= price_range[1]).sum() / y.size
        opportunity_space_1d.append(under_within)
        scenarios_list.append(scenario)
    yticks = [0, 1, 2, 3, 4, 5, 6, 7]
    xticks = [0, 85, 170, 255, 340]
    # yticks = [-2, -1, 0, 1, 2, 3, 4]
    # xticks = [0, 400, 800, 1200, 1600]
    # breakpoint()
    # xticks = [0, 12, 24, 36, 48]
    fig, ax, axes = bst.plots.plot_kde(
        y=y1d, x=x1d, xticks=xticks, yticks=yticks,
        xticklabels=True, yticklabels=True,
        colors = [
            clr.LinearSegmentedColormap.from_list(
                c.ID,
                [*[c.shade(60 - 20 * j).RGBn for j in range(3)],
                 *[c.tint(30 * j).RGBn for j in range(3)]],
                N=256
            )
            for c in (GG_colors.purple, GG_colors.green)
        ],
        xbox_kwargs=[
            dict(light=GG_colors.purple.tint(50).RGBn, dark=GG_colors.purple.shade(60).RGBn),
            dict(light=GG_colors.green.tint(50).RGBn, dark=GG_colors.green.shade(60).RGBn)
        ],
        ybox_kwargs=[
            dict(light=GG_colors.purple.tint(50).RGBn, dark=GG_colors.purple.shade(60).RGBn),
            dict(light=GG_colors.green.tint(50).RGBn, dark=GG_colors.green.shade(60).RGBn)
        ],
        xlabel=r'Carbon intensity$_{\mathrm{EtOH}}$ $[\mathrm{gCO2e} \cdot \mathrm{L}^{\mathrm{-1}}]$',
        # xlabel='TCI $[10^6 \cdot \mathrm{USD}]$',
        ylabel=r'MSP$_{\mathrm{Resin}}$ $[\mathrm{USD} \cdot \mathrm{kg}^{\mathrm{-1}}]$',
        xbox_width=400,
        aspect_ratio=1.1,
    )
    plt.sca(ax)
    opportunity_space_1d = np.array(opportunity_space_1d)
    print(opportunity_space_1d)
    plt.fill_between(plt.xlim(), *price_range,
                     color=[*line_color, 0.5],
                     linewidth=2,
                     zorder=-1)
    cellulosic_ethanol_GWP = 297.2 # GREET 2023
    bst.plots.plot_vertical_line(cellulosic_ethanol_GWP, lw=2, ls='--', color=line_color, zorder=-1)
    plt.subplots_adjust(
        hspace=0, wspace=0,
        top=0.9, bottom=0.12,
        left=0.2, right=0.98,
    )
    for i in ('svg', 'png'):
        file = os.path.join(images_folder, f"{'_'.join(scenarios_list)}_kde.{i}")
        plt.savefig(file, dpi=900, transparent=True)
    
# def plot_kde_CI_MSP(scenarios=None):
#     from warnings import filterwarnings
#     filterwarnings('ignore')
#     set_font(size=10)
#     set_figure_size(width='full', aspect_ratio=0.6)
#     opportunity_space_2d = []
#     scenarios = ['all', 'baseline', 'NREL']
#     pm = sp.STRAPMSWProcess(simulate=False, scenario=scenarios[0])
#     price_range = np.array(pm.resin_price_range) / 907.185
#     scenarios_list = []
#     metrics = [pm.GWP_ethanol.index, pm.MSP.index]
#     Xi, Yi = metrics
#     y1d = []
#     x1d = []
#     opportunity_space_1d = []
#     # opportunity_space_2d.append(opportunity_space_1d)
#     for scenario in scenarios:
#         pm = sp.STRAPMSWProcess(simulate=False, scenario=scenario)
#         df = get_monte_carlo('STRAPMSW_' + scenario, metrics)
#         y = df[Yi].values
#         x = df[Xi].values * 1e3
#         # print(scenario, 'y', y.min(), y.max())
#         y1d.append(y)
#         x1d.append(x)
#         # print(scenario, 'x', x.min(), x.max())
#         under = (y <= price_range[0]).sum() / y.size
#         within = (y <= price_range[1]).sum() / y.size - under
#         opportunity_space_1d.append(
#             [under, within]
#         )
#         scenarios_list.append(scenario)
#     yticks = [-2, -1, 0, 1, 2, 3, 4]
#     xticks = [0, 75, 150, 225, 300, 375]
#     # yticks = [-2, -1, 0, 1, 2, 3, 4]
#     # xticks = [0, 400, 800, 1200, 1600]
#     # breakpoint()
#     # xticks = [0, 12, 24, 36, 48]
#     fig, axes = bst.plots.plot_kde_1d(
#         ys=y1d, xs=x1d, xticks=xticks, yticks=yticks,
#         xticklabels=True, yticklabels=True,
#         xbox_kwargs=dict(light=CABBI_colors.orange.RGBn, dark=CABBI_colors.orange.shade(60).RGBn),
#         ybox_kwargs=dict(light=CABBI_colors.blue.RGBn, dark=CABBI_colors.blue.shade(60).RGBn),
#         xlabel=r'Carbon intensity $[\mathrm{gCO2e} \cdot \mathrm{L}^{\mathrm{-1}}]$',
#         # xlabel='TCI $[10^6 \cdot \mathrm{USD}]$',
#         ylabel=r'MSP $[\mathrm{USD} \cdot \mathrm{kg}^{\mathrm{-1}}]$',
#         xbox_width=400,
#         aspect_ratio=1.1,
#     )
#     xlb, xub = plt.xlim()
#     ylb, yub = plt.ylim()
#     xpos = lambda x: xlb + (xub - xlb) * x
#     ypos = lambda y: ylb + (yub - ylb) * y
#     xleft = 0.02
#     # xright = 0.98
#     # ytop = 0.94
#     ybottom = 0.02
#     first = True
#     opportunity_space_1d = np.array(opportunity_space_1d)
#     for ax, ps in zip(axes[0], opportunity_space_1d):
#         plt.sca(ax)
#         # bst.plots.plot_horizontal_line(market_price, lw=2, ls='-', color=line_color)
#         plt.fill_between(plt.xlim(), *price_range,
#                          color=[*line_color, 0.5],
#                          linewidth=2,
#                          zorder=-1)
#         cellulosic_ethanol_GWP = 297.2 # GREET 2023
#         bst.plots.plot_vertical_line(cellulosic_ethanol_GWP, lw=2, ls='--', color=line_color, zorder=-1)
#         if first: 
#             values = f"{ps[0]:.0%} under, {ps[1]:.0%} within\nmarket price"
#             first = False
#         else:
#             values = f"{ps[0]:.0%}, {ps[1]:.0%}"
#         plt.text(xpos(xleft), ypos(ybottom), values, color=line_color,
#                  horizontalalignment='left', verticalalignment='bottom',
#                  fontsize=10, fontweight='bold', zorder=10)
#     plt.subplots_adjust(
#         hspace=0, wspace=0,
#         top=0.9, bottom=0.12,
#         left=0.1, right=0.98,
#     )
#     for i in ('svg', 'png'):
#         file = os.path.join(images_folder, f"{'_'.join(scenarios_list)}_kde.{i}")
#         plt.savefig(file, dpi=900, transparent=True)
    
def plot_kde_TCI_MSP(scenarios_2d=None):
    from warnings import filterwarnings
    filterwarnings('ignore')
    set_font(size=10)
    set_figure_size(width='full', aspect_ratio=0.62)
    y2d = []
    x2d = []
    opportunity_space_2d = []
    pm = sp.STRAPMSWProcess(simulate=False, scenario=scenarios_2d[0][0])
    market_price = pm.baseline_resin_price / 907.185
    scenarios_list = []
    metrics = [pm.TCI.index, pm.MSP.index]
    Xi, Yi = metrics
    for scenarios_1d in scenarios_2d:
        y1d = []
        x1d = []
        opportunity_space_1d = []
        y2d.append(y1d)
        x2d.append(x1d)
        opportunity_space_2d.append(opportunity_space_1d)
        for scenario in scenarios_1d:
            pm = sp.STRAPMSWProcess(simulate=False, scenario=scenario)
            df = get_monte_carlo('STRAPMSW_' + scenario, metrics)
            y = df[Yi].values
            y1d.append(y)
            x1d.append(df[Xi].values / 1e6)
            opportunity_space_1d.append(
                (y <= market_price).sum() / y.size
            )
            scenarios_list.append(scenario)
    yticks = [-3, -2, -1, 0, 1, 2]
    # yticks = [-2, -1, 0, 1, 2, 3, 4]
    xticks = [0, 400, 800, 1200, 1600]
    # breakpoint()
    # xticks = [0, 1, 2, 3, 4, 5, 6]
    fig, axes = bst.plots.plot_kde_1d(
        ys=y2d, xs=x2d, xticks=xticks, yticks=yticks,
        xticklabels=True, yticklabels=True,
        xbox_kwargs=dict(light=CABBI_colors.orange.RGBn, dark=CABBI_colors.orange.shade(60).RGBn),
        ybox_kwargs=dict(light=CABBI_colors.blue.RGBn, dark=CABBI_colors.blue.shade(60).RGBn),
        # xlabel='Carbon intensity $[\mathrm{gCO2e} \cdot \mathrm{L}^{\mathrm{-1}}]$',
        xlabel=r'TCI $[10^6 \cdot \mathrm{USD}]$',
        ylabel=r'MSP $[\mathrm{USD} \cdot \mathrm{kg}^{\mathrm{-1}}]$',
        xbox_width=400,
        aspect_ratio=1.1,
    )
    xlb, xub = plt.xlim()
    ylb, yub = plt.ylim()
    xpos = lambda x: xlb + (xub - xlb) * x
    ypos = lambda y: ylb + (yub - ylb) * y
    # xleft = 0.02
    xright = 0.98
    # ytop = 0.94
    ybottom = 0.02
    first = True
    for ax, ps in zip(axes[0], np.transpose(opportunity_space_2d)):
        if first: 
            description = ' under market price'
            first = False
        else:
            description = ''
        plt.sca(ax)
        values = ' and '.join([f"{p:.0%}" for p in ps])
        bst.plots.plot_horizontal_line(market_price, lw=2, ls='-', color=line_color)
        plt.text(xpos(xright), ypos(ybottom), values + description, color=line_color,
                 horizontalalignment='right', verticalalignment='bottom',
                 fontsize=10, fontweight='bold', zorder=10)
    plt.subplots_adjust(
        hspace=0, wspace=0,
        top=0.9, bottom=0.1,
        left=0.1, right=0.98,
    )
    for i in ('svg', 'png'):
        file = os.path.join(images_folder, f"{'_'.join(scenarios_list)}_kde.{i}")
        plt.savefig(file, dpi=900, transparent=True)
    
def get_monte_carlo_key(index, dct, with_units=False):
    key = index[1] if with_units else index[1].split(' [')[0]
    if key in dct: key = f'{key}, {index[0]}'
    return key
    
def montecarlo_results(scenario='all'):
    f = sp.STRAPMSWProcess(scenario=scenario)
    metrics = [
        f.GWP_ethanol.index, f.MSP.index
    ]
    results = {}
    df = get_monte_carlo('STRAPMSW_' + scenario, metrics)
    results[f.name] = dct = {}
    for index in metrics:
        data = df[index].values
        q05, q50, q95 = roundsigfigs(np.percentile(data, [5, 50, 95], axis=0), 3)
        key = get_monte_carlo_key(index, dct, False)
        if q50 < 0:
            dct[key] = f"{-q50} [{-q95}, {-q05}] -negative-"
        else:
            dct[key] = f"{q50} [{q05}, {q95}]"
    return results

def _get_full_name(f):
    a = f.element_name
    b = f.name
    if a == '-':
        name = b.capitalize().replace('gwp', 'GWP')
    elif b == 'GWP': 
        name = f"{a} {b}"
    else:
        name = f"{a} {b.lower()}"
    return name.replace('co2', 'CO2')

def get_distributions(parameters, save=True, table=True, **kwargs):
    name = 'name'
    full_name = 'full_name'
    parameters = {
        i: full_name for i in parameters
    }
    
    def with_units(f, name, units=None):
        if units is None: units = f.units
        if units is None: units = '-'
        return name, units
        
    def nested_roundsigfigs(x, n):
        try:
            return roundsigfigs(float(x), n)
        except:
            return [roundsigfigs(i, n) for i in x._repr.values()]
        
    def get_distribution(f):
        d = f.distribution
        dname = type(d).__name__
        if dname == 'Trunc': dname = type(d._repr['dist']).__name__
        values = tuple([nested_roundsigfigs(j, 3) for j in d._repr.values()])
        if dname == 'Triangle':
            return f"Triangular{values}"
            
        elif dname == 'Uniform':
            return f'Uniform{values}'
        elif dname == 'Normal':
            mu, sigma = [str(i) for i in values[0]]
            return f"Normal(μ={mu}, σ={sigma})"
        else: 
            raise RuntimeError('unknown distribution')
        
    rows = []
    for i, j in parameters.items():
        if j == name:
            parameter_name, units = with_units(i, i.name)
        elif j == full_name:
            parameter_name, units = with_units(i, _get_full_name(i))
        elif isinstance(j, tuple):
            parameter_name, units = with_units(i, *j)
        elif isinstance(j, str):
            parameter_name, units = with_units(i, j)
        else:
            raise TypeError(str(j))
        rows.append({
            'Parameter': parameter_name,
            'Units': units,
            'Baseline': roundsigfigs(i.baseline, 3),
            'Distribution': get_distribution(i),
        })
    if table:
        table = pd.DataFrame(
            rows, 
            index=list(
                range(1, len(parameters) + 1)
            )
        )
        # table.drop(columns=['Shape', 'Parameter(s)'], inplace=True)
        if save:
            file = os.path.join(results_folder, 'parameter_distributions.xlsx')
            table.to_excel(file)
            table.index.name = '#'
        return table
    else:
        return rows