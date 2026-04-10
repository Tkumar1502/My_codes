# -*- coding: utf-8 -*-
"""
Created on Fri Nov  5 01:57:46 2021

@author: yrc2
"""
import os
import numpy as np
import pandas as pd
import biosteam as bst
from plastics import strap
from biorefineries.tea.cellulosic_ethanol_tea import foc_table, capex_table
from thermosteam.utils import array_roundsigfigs
import matplotlib.pyplot as plt
import os

__all__ = (
    'save_all_tables',
    'save_system_reports',
    'save_feedstock_composition',
    'save_detailed_expenditure_tables',
    'save_detailed_life_cycle_tables',  
    'save_water_material_balance_table',
    'LCA_verification',
    'TEA_verification',
    'plot_bars',
    'plot_utilities',
    'print_relative_flows',
)

folder = os.path.dirname(__file__)
folder = os.path.join(folder, 'results')
images_folder = os.path.join(os.path.dirname(__file__), 'images')

def save_all_tables():
    save_system_reports()
    save_feedstock_composition()
    save_detailed_expenditure_tables()
    save_detailed_life_cycle_tables()
    save_water_material_balance_table()

def save_feedstock_composition():
    pm = strap.STRAPMSWProcess(scenario='baseline')
    df = pm.feedstock._info(
        'cwt', T=None, P=None, flow=None, 
        composition=None, N=100, IDs=None,
        sort=False, df=True
    )
    composition = df.loc['Composition [%]']
    composition.columns = ['MRF residue\ncomposition [%]']
    file = os.path.join(folder, 'feedstock_composition.xlsx')
    composition.to_excel(file)

def save_system_reports():
    for scenario in ('baseline', 'NREL'):
        pm = strap.STRAPMSWProcess(scenario=scenario)
        filename = f'STRAPMSW_{scenario}_detailed_report.xlsx'
        file = os.path.join(folder, filename)
        pm.system.save_report(file)

def save_detailed_expenditure_tables(sigfigs=3):
    filename = 'expenditures.xlsx'
    file = os.path.join(folder, filename)
    writer = pd.ExcelWriter(file)
    product = 'resin'
    process_models = [
        strap.STRAPMSWProcess(scenario=i)
        for i in ('baseline', 'NREL')
    ]
    for pm in process_models: pm.resin.price = pm.tea.solve_price(pm.resin)
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
    
def save_detailed_life_cycle_tables(sigfigs=3):
    process_models = [strap.STRAPMSWProcess(scenario=i) for i in ('baseline', 'potential')]
    for i in process_models: i.resin.ID = 'Polymer_resin'
    systems = [i.system for i in process_models]
    filename = 'life_cycle.xlsx'
    file = os.path.join(folder, filename)
    writer = pd.ExcelWriter(file)
    streams = [getattr(i, 'resin') for i in process_models] + [getattr(i, 'ethanol') for i in process_models]
    names = ['Baseline', 'Potential']
    tables = {
        'Inventory': bst.report.lca_inventory_table(
            systems, items=streams, system_names=names
        ),
        'Energy allocation factors': bst.report.lca_property_allocation_factor_table(
            systems, property='energy', basis='GGE', system_names=names, groups=('ethanol',),
        ),
        'Characterization factors': bst.report.lca_characterization_factor_table(
            systems
        ),
    }
    index = [
        ('Electricity', 'GWP [kg∙CO2e∙kWh-1]'), 
        ('Electricity', 'FFC [MJ∙kWh-1]'), 
        ('Electricity', 'WU [kg-water∙kWh-1]'), 
        ('Ethanol', 'GWP [kg∙CO2e∙L-1]'), 
        ('Ethanol', 'FFC [MJ∙L-1]'), 
        ('Ethanol', 'WU [kg-water∙L-1]'), 
        ('Polymer resin', 'GWP [kg∙CO2e∙kg-1]'),
        ('Polymer resin', 'FFC [MJ∙kg-1]'),
        ('Polymer resin', 'WU [kg-water∙kg-1]')
    ]
    values = np.zeros([len(index), 2])
    for j, pm in enumerate(process_models):
        indicators = [
            pm.GWP_electricity,  pm.FFC_electricity, pm.WU_electricity, 
            pm.GWP_ethanol, pm.FFC_ethanol, pm.WU_ethanol, 
            pm.GWP_polymer_resin, pm.FFC_polymer_resin, pm.WU_polymer_resin,
        ]
        for i, indicator in enumerate(indicators): values[i, j] = indicator()
    columns = ['Baseline', 'Potential']
    df_gwp = pd.DataFrame(values, index=pd.MultiIndex.from_tuples(index), columns=columns)
    tables['Estimated environmental impact'] = df_gwp
    for key, table in tables.items(): 
        array_roundsigfigs(table.values, sigfigs=3, inplace=True)
        table.to_excel(writer, key) 
    writer.close()
    return tables

def save_water_material_balance_table():
    process_models = [strap.STRAPMSWProcess(scenario=i) for i in ('baseline', 'potential')]
    for i in process_models: i.resin.ID = 'Polymer_resin'
    baseline, potential = [i.system for i in process_models]
    filename = 'water_mass_balance.xlsx'
    file = os.path.join(folder, filename)
    writer = pd.ExcelWriter(file)
    names = ['Baseline', 'Potential']
    tables = {
        'Baseline water balance': bst.report.water_mass_balance_table(
            baseline
        ),
        'Baseline water balance': bst.report.water_mass_balance_table(
            potential
        ),
    }
    for key, table in tables.items(): 
        array_roundsigfigs(table.values, sigfigs=3, inplace=True)
        table.to_excel(writer, key) 
    writer.close()
    return tables

def LCA_verification(sigfigs=3):
    process_models = [strap.STRAPMSWProcess(scenario=i, preprocessing=True) for i in ('baseline', 'NREL')]
    for i in process_models: i.resin.ID = 'Polymer_resin'
    systems = [i.system for i in process_models]
    names = ['Baseline', 'Potential']
    streams = [getattr(i, 'resin') for i in process_models] + [getattr(i, 'ethanol') for i in process_models]
    tables = {
        'Inventory': bst.report.lca_inventory_table(
            systems, 'GWP', streams, system_names=names
        ),
        'Energy allocation factors': bst.report.lca_property_allocation_factor_table(
            systems, property='energy', basis='GGE', system_names=names, groups=('ethanol',),
        ),
    }
    values = np.zeros([3, 2])
    index = ['Electricity [kg∙CO2e∙kWh-1]', 
             'Ethanol [kg∙CO2e∙L-1]', 
             'Polymer resin [kg∙CO2e∙kg-1]']
    for j, pm in enumerate(process_models):
        GWPs = [pm.GWP_electricity, pm.GWP_ethanol, pm.GWP_polymer_resin]
        for i, GWP in enumerate(GWPs):
            values[i, j] = GWP()
    columns = ['Baseline', 'Potential']
    df_gwp = pd.DataFrame(values, index=index, columns=columns)
    tables['Estimated environmental impact'] = df_gwp
    for key, table in tables.items(): 
        array_roundsigfigs(table.values, sigfigs=3, inplace=True)
    
    inventory = tables['Inventory'].loc['Inputs']
    inventory = inventory.drop(['Natural gas', 'Bisulfite', 'Citric acid', 'NaOCl'])
    inventory.loc['Direct non-biogenic emissions'] = [pm.direct_nonbiogenic_emissions() * pm.tea.operating_hours for i in process_models]
    CFs = [0.525, 0.82, 0.404, 0.84, 1.81, 1.]
    impact = (np.array(CFs).reshape([len(CFs), 1]) * inventory.iloc[:, [0, 1]])
    impact /= impact.sum(axis=0)
    impact *= 100
    impact.columns = ['Baseline', 'Potential']
    impact = impact.sort_values(by=['Baseline'], ascending=False)
    impact.T.plot.bar(stacked=True, rot=0)
    plt.ylabel('Contribution to carbon intensity [%]')
    for i in ('svg', 'png'):
        file = os.path.join(images_folder, f'verification_impact_breakdown.{i}')
        plt.savefig(file, dpi=900, transparent=True)
    
    allocation = tables['Energy allocation factors'] * 100
    allocation.columns = ['Baseline', 'Potential']
    allocation.T.plot.bar(stacked=True, rot=0)
    plt.ylabel('Allocation factors [%]')
    for i in ('svg', 'png'):
        file = os.path.join(images_folder, f'verification_allocation.{i}')
        plt.savefig(file, dpi=900, transparent=True)
    
    allocated_impact = tables['Estimated environmental impact']
    allocated_impact.index = [i.replace('-1', '$^{-1}$') for i in allocated_impact.index]
    allocated_impact['Grid, Cellulosic, Virgin'] = [0.603, 0.319, 1.89]
    allocated_impact.T.plot.bar(stacked=False, rot=0)
    plt.ylabel('Carbon intensity')
    for i in ('svg', 'png'):
        file = os.path.join(images_folder, f'verification_allocated_carbon_intensity.{i}')
        plt.savefig(file, dpi=900, transparent=True)
        
    # allocated_impact.columns = [i.replace('-1', '$^{-1}$') for i in allocated_impact.columns]
    # electricity, ethanol, resin = [allocated_impact[i] for i in allocated_impact]
    # electricity['Grid'] = 0.603
    # ethanol['Cellulosic'] = 0.319
    # resin['Virgin'] = 1.89
    # for i, j in zip([electricity, ethanol, resin], ['electricity', 'ethanol', 'resin']):
    #     plt.figure()
    #     i.plot.bar(stacked=False, rot=0)
    #     units = i.name[i.name.index('['):]
    #     plt.ylabel(f'Carbon intensity {units}')
    #     for i in ('svg', 'png'):
    #         file = os.path.join(images_folder, f'verification_allocated_{j}_carbon_intensity.{i}')
    #         plt.savefig(file, dpi=900, transparent=True)
        
    return tables

def TEA_verification():
    from plastics import strap
    from warnings import filterwarnings
    filterwarnings('ignore')
    process = strap.STRAPMSWProcess(scenario='baseline', preprocessing=True, simulate=False) # You may also want scenario='all' and scenario='NREL'.
    capex, opex = process.get_production_cost_contribution()
    capex = pd.DataFrame([*capex.values( )], index=[*capex], columns=['CAPEX'])
    capex.T.plot.bar(stacked=True)
    opex = pd.DataFrame([*opex.values()], index=[*opex], columns=['OPEX']) 
    opex.T.plot.bar(stacked=True)
    plt.show()
    
# Corn ethanol is 1.10 kg CO2e / L
# Gasoline 0.6267 + 2.15 kg CO2e / L  basen on greet and assuming octane to CO2.

# from thermosteam import Chemical
# octane = Chemical('octane')
# rho = octane.get_property('rho', 'kg/L', T=298.15, P=101325, phase='l')
# rho_kmol_per_L = rho / octane.MW
# MW_CO2 = 12.01 + 16 * 2
# kg_CO2_per_octane = 8 * MW_CO2
# emissions = kg_CO2_per_octane * rho_kmol_per_L

# def plot_LCA(sigfigs=3):
#     import numpy as np
#     import pandas as pd
#     import biosteam as bst
#     import seaborn as sns
#     import matplotlib.pyplot as plt
#     from plastics import strap
#     sns.set(style='ticks')
#     bst.set_figure_size(aspect_ratio=0.6, width='full')
#     bst.set_font(size=10)
#     fig, (ethanol_ax, resin_ax, electricity_ax) = plt.subplots(1, 3)
#     scenarios = ('baseline', 'potential')
#     process_models = [strap.STRAPMSWProcess(scenario=i, preprocessing=False) for i in scenarios]
#     for i in process_models: i.resin.ID = 'Polymer_resin'
#     systems = [i.system for i in process_models]
#     streams = [getattr(i, 'resin') for i in process_models] + [getattr(i, 'ethanol') for i in process_models]
#     inventory_table = bst.report.lca_inventory_table(
#         systems, 'GWP', streams, system_names=scenarios
#     )
#     CFs = bst.report.lca_displacement_allocation_table(systems, 'GWP', streams, system_names=scenarios
#     )['Characterization factor [kg*CO2e/kg]'].loc['Inputs']
#     impact = inventory_table.loc['Inputs'] * CFs.values[:, None]
#     other = ['Natural gas', 'Bisulfite', 'Citric acid', 'NaOCl', 'Urea', 'Denaturant', 'FGD lime', 'Xylene']
#     impact.loc['Other'] = impact.loc[other].sum()
#     impact = impact.drop(other)
#     impact.loc['Direct non-biogenic emissions'] = [i.direct_nonbiogenic_emissions() * i.tea.operating_hours for i in process_models]
#     impact.columns = ['baseline', 'potential']
#     impact = impact.sort_values(by=['baseline'], ascending=False)
#     products = [
#         (ethanol_ax, 'GWP_ethanol', 'EtOH', 'L', 319),
#         (resin_ax, 'GWP_resin', 'Resin', 'kg', 603),
#         (electricity_ax, 'GWP_electricity', 'Electricity', 'kWh',  1890),
#     ]
#     for ax, method, name, units, alternative in products:
#         impact *= [getattr(i, method)() for i in process_models] / impact.sum(axis=0) * 1000
#         impact.T.plot.bar(stacked=True, rot=30, ax=ax, fontsize=10)
#         plt.sca(ax)
#         ax.tick_params(axis='x', which='major', length=6,
#                        direction="inout")
#         ax.tick_params(axis='y', which='major', length=6,
#                        direction="inout")
#         ax.get_legend().remove()
#         ax.spines[['right', 'top']].set_visible(False)
#         plt.ylabel(r'name carbon intensity [g∙CO$_2$e∙units$^{-1}$]'.replace('name', name).replace('units', units), fontsize=10)
#         plt.axhline(y=alternative, color='darkgray', linestyle='--')
#     plt.subplots_adjust(wspace=0.9, hspace=0.9, right=0.9)
#     for i in ('svg', 'png'):
#         file = os.path.join(images_folder, f'LCA_scenarios.{i}')
#         plt.savefig(file, dpi=900, transparent=True)
    
def plot_IRR_bars():
    import numpy as np
    import pandas as pd
    import biosteam as bst
    import seaborn as sns
    import matplotlib.pyplot as plt
    from plastics import strap
    from biosteam.utils import GG_colors
    sns.set(style='ticks')
    bst.set_figure_size(aspect_ratio=0.55, width='full')
    bst.set_font(size=10)
    fig, ax = plt.subplots(1, 1)
    scenario_names = ('baseline', 'potential')
    process_models = [
        strap.STRAPMSWProcess(scenario=name, simulate=True) 
        for name in scenario_names
    ]
    representative_stages = (
        ('WA', 'pacific'),
        ('MA', 'northeast'),
        ('WY', 'mountains'),
        ('MO', 'midwest'),
        ('KY', 'southeast'),
        ('AR', 'southcentral'),
    )
    data = []
    for state, region in representative_stages:
        for scenario, pm in zip(scenario_names, process_models):
            mean, _ = pm.get_tipping_fee(state)
            pm.set_MSW_tipping_fee(mean)
            data.append([f'{state}\n{region}', scenario, pm.IRR()])
    df = pd.DataFrame(data, columns=['State', 'Scenario', 'IRR [%]'])
    df = df.sort_values('IRR [%]')
    sns.barplot(
        df, x="State", y="IRR [%]", hue="Scenario",
        palette={
            'baseline': GG_colors.purple.RGBn,
            'potential': GG_colors.green.RGBn,
        }
    )
    plt.xlabel('')
    ax.get_legend().remove()

def plot_utilities():
    import numpy as np
    import pandas as pd
    import biosteam as bst
    import seaborn as sns
    import matplotlib.pyplot as plt
    from plastics import strap
    from biosteam.utils import GG_colors
    from warnings import filterwarnings
    filterwarnings('ignore')
    sns.set(style='ticks')
    bst.set_figure_size(aspect_ratio=0.6, width='full')
    bst.set_font(size=10)
    fig, (heating_ax, cooling_ax, electricity_ax) = plt.subplots(1, 3)
    scenario_names = ('baseline', 'potential')
    process_models = [
        strap.STRAPMSWProcess(scenario=name, simulate=True) 
        for name in scenario_names
    ]
    for pm in process_models:
        for facility in pm.system.facilities:
            facility.heat_utilities = [hu for hu in facility.heat_utilities if hu.flow > 0]
    def get(name):
        values = [
            getattr(pm.system, f'get_{name}')() / pm.MRF_residues.F_mass / pm.system.operating_hours
            for pm in process_models
        ]
        if 'duty' in name:
            values = [i / 1000 for i in values]
        return values
    index = ('heating_duty', 'cooling_duty', 'electricity_consumption')
    data = [get(i) for i in index]
    df = pd.DataFrame(
        data, 
        columns=['baseline', 'potential'], 
        index=[i.replace('_', ' ') for i in index],
    )
    options = [
        (heating_ax, 'heating duty', r'Heating duty [MJ$\cdot$kg$^{-1}$]'),
        (cooling_ax, 'cooling duty', r'Cooling duty [MJ$\cdot$kg$^{-1}$]'),
        (electricity_ax, 'electricity consumption', r'Electricity [kW$\cdot$kg$^{-1}$]'),
    ]
    for ax, name, label in options:
        plt.sca(ax)
        df.loc[name].plot.bar(
            ax=ax, fontsize=10,
            color=[GG_colors.purple.RGBn, GG_colors.green.RGBn],
            rot=30,
        )
        plt.ylabel(label, fontsize=10)
        ax = plt.gca()
        ax.tick_params(axis='x', which='major', length=6,
                       direction="inout")
        ax.tick_params(axis='y', which='major', length=6,
                       direction="inout")
        ax.spines[['right', 'top']].set_visible(False)
    
    plt.subplots_adjust(wspace=0.8, hspace=0.8, right=0.95, bottom=0.2, top=0.95)
    for i in ('svg', 'png'):
        file = os.path.join(images_folder, f'energy_bars.{i}')
        plt.savefig(file, dpi=900, transparent=True)

def print_relative_flows():
    scenario_names = ('baseline', 'potential')
    process_models = [
        strap.STRAPMSWProcess(scenario=name, simulate=True) 
        for name in scenario_names
    ]
    for name, pm in zip(scenario_names, process_models):
        print('baseline')
        print('--------')
        print('solvent', 
            pm.U105.ins[1].F_mass 
            / pm.MRF_residues.F_mass
        )
        print('wastewater', 
            pm.U301.ins[0].F_mass 
            / pm.MRF_residues.F_mass
        )
        print('resin', 
            pm.resin.F_mass
            / pm.MRF_residues.F_mass
        )
        print('ethanol', 
            pm.ethanol.F_mass
            / pm.MRF_residues.F_mass
        )
    
def plot_bars(LCA=True):
    import numpy as np
    import pandas as pd
    import biosteam as bst
    import seaborn as sns
    import matplotlib.pyplot as plt
    from plastics import strap
    sns.set(style='ticks')
    bst.set_figure_size(aspect_ratio=0.6, width='full')
    bst.set_font(size=10)
    fig, (CAPEX_ax, OPEX_ax, LCA_ax) = plt.subplots(1, 3)
    scenario_names = ('baseline', 'potential')
    process_models = [
        strap.STRAPMSWProcess(scenario=name, simulate=True) 
        for name in scenario_names
    ]
    for pm in process_models:
        for group in pm.unit_groups: 
            try:
                group.autofill_metrics(
                    shorthand=False, 
                    installed_cost=True,
                    cooling_duty=False,
                    heating_duty=False,
                    electricity_consumption=False,
                    electricity_production=False,
                    material_cost=False
                )
            except:
                pass
    CAPEX = pd.concat([
        bst.UnitGroup.df_from_groups(pm.unit_groups)
        for pm in process_models
    ], axis=1)
    CAPEX.loc['Facilities'] = CAPEX.loc['WWT'] + CAPEX.loc['Facilities']
    CAPEX.loc['Indirect costs'] = [(i.tea.TCI - i.tea.DPI) /1e6 for i in process_models]
    CAPEX = CAPEX.drop(['WWT'])
    CAPEX.loc['Other'] = [i.tea.TCI / 1e6 for i in process_models] - CAPEX.sum()
    CAPEX = CAPEX.sort_index(key=lambda x:[-abs(sum(CAPEX.loc[i])) for i in x])
    CAPEX.columns = ['baseline', 'potential']
    plt.sca(CAPEX_ax)
    CAPEX.T.plot.bar(stacked=True, rot=30, ax=CAPEX_ax, fontsize=10)
    plt.ylabel(r'CAPEX [$10^6\cdot$USD]', fontsize=10)
    ax = plt.gca()
    ax.tick_params(axis='x', which='major', length=6,
                   direction="inout")
    ax.tick_params(axis='y', which='major', length=6,
                   direction="inout")
    ax.get_legend().remove()
    ax.spines[['right', 'top']].set_visible(False)
    # ax.legend(bbox_to_anchor=(1.05, 1.05))
    # ax.legend(
    #     loc='upper center', bbox_to_anchor=(0.5, 1.05),
    #     ncol=3, fancybox=True
    # )
    VOC_table = bst.report.voc_table(
        [i.system for i in process_models], 
        system_names=scenario_names,
        product_IDs=[]
    )
    VOC_table = VOC_table.drop('Price [$/MT]', axis=1)
    materials = VOC_table.loc['Raw materials']
    other_cost = VOC_table.loc['Other utilities & fees'].sum()
    key_raw_materials = ['MRF residues', 'Cellulase']
    OPEX = materials.loc[key_raw_materials]
    OPEX.loc['Other'] = other_cost + materials.sum() - OPEX.sum() + [i.tea.FOC/1e6 for i in process_models]
    products = VOC_table.loc['Co-products & credits']
    products.loc['Ethanol/RIN D3'] = products.loc['Ethanol RIN D3'] + products.loc['Ethanol']
    products.loc['Electricity prod.'] = products.loc['Electricity production']
    products = products.drop(['Ethanol RIN D3', 'Ethanol', 'Electricity production'])
    OPEX_and_revenue = pd.concat([-OPEX, products])
    OPEX_and_revenue = OPEX_and_revenue.sort_index(key=lambda x:[-sum(OPEX_and_revenue.loc[i]) for i in x])
    OPEX_and_revenue.columns = ['baseline', 'potential']
    print(OPEX_and_revenue)
    print(OPEX_and_revenue.sum())
    plt.sca(OPEX_ax)
    OPEX_and_revenue.T.plot.bar(stacked=True, rot=30, ax=OPEX_ax, fontsize=10)
    plt.axhline(y=0, color='darkgray', linestyle='--')
    plt.ylabel(r'OPEX & Revenue [$10^6\cdot$USD$\cdot$yr$^{-1}$]', fontsize=10)
    ax = plt.gca()
    ax.tick_params(axis='x', which='major', length=6,
                   direction="inout")
    ax.tick_params(axis='y', which='major', length=6,
                   direction="inout")
    ax.get_legend().remove()
    ax.spines[['right', 'top']].set_visible(False)
    # ax.legend(bbox_to_anchor=(1.05, 1.05))
    # ax.legend(
    #     loc='upper center', bbox_to_anchor=(0.5, 1.05),
    #     ncol=3, fancybox=True
    # )
    if LCA:
        systems = [i.system for i in process_models]
        streams = [getattr(i, 'resin') for i in process_models] + [getattr(i, 'ethanol') for i in process_models]
        inventory_table = bst.report.lca_inventory_table(
            systems, 'GWP', streams, system_names=scenario_names
        )
        CFs = bst.report.lca_displacement_allocation_table(systems, 'GWP', streams, system_names=scenario_names
        )['Characterization factor [kg*CO2e/kg]'].loc['Inputs']
        impact = inventory_table.loc['Inputs'] * CFs.values[:, None]
        other = ['Bisulfite', 'Citric acid', 'NaOCl', 'Urea', 'Denaturant', 'FGD lime', 'Xylene', 'Bioreactor cleaning chemicals']
        impact.loc['Other'] = impact.loc[other].sum()
        impact = impact.drop(other)
        impact.loc['Direct non-biogenic emissions'] = [i.direct_nonbiogenic_emissions() * i.tea.operating_hours for i in process_models]
        impact.columns = ['baseline', 'potential']
        impact = impact.sort_values(by=['baseline'], ascending=False)
        impact *= [i.GWP_ethanol() for i in process_models] / impact.sum(axis=0) * 1000
        ax = LCA_ax
        plt.sca(ax)
        impact.T.plot.bar(stacked=True, rot=30, ax=ax, fontsize=10)
        ax.tick_params(axis='x', which='major', length=6,
                       direction="inout")
        ax.tick_params(axis='y', which='major', length=6,
                       direction="inout")
        ax.get_legend().remove()
        ax.spines[['right', 'top']].set_visible(False)
        plt.ylabel(r'Carbon intensity$_{\mathrm{EtOH}}$ [g∙CO$_2$e∙L$^{-1}$]', fontsize=10)
        plt.axhline(y=297.2, color='gray', linestyle='--')
        plt.axhline(y=400, color='gray', linestyle='--')
    else:
        IRR = pd.DataFrame([i.IRR() for i in process_models], index=['baseline', 'potential'])
        ax = LCA_ax
        plt.sca(ax)
        IRR.plot.bar(stacked=True, rot=30, ax=ax, fontsize=10)
        ax.tick_params(axis='x', which='major', length=6,
                       direction="inout")
        ax.tick_params(axis='y', which='major', length=6,
                       direction="inout")
        ax.get_legend().remove()
        ax.spines[['right', 'top']].set_visible(False)
        plt.ylabel(r'IRR [%]', fontsize=10)
    plt.subplots_adjust(wspace=1.5, hspace=1.5, right=0.95, bottom=0.2, top=0.95)
    for i in ('svg', 'png'):
        file = os.path.join(images_folder, f'bars.{i}')
        plt.savefig(file, dpi=900, transparent=True)