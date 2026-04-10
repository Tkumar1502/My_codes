# -*- coding: utf-8 -*-
"""
"""
import biosteam as bst

__all__ = (
    'load_process_settings',
    'load_STRAP_MSW_process_settings',
)

# TODO: Potentially add Boiler-Turbogenerator and Cooling Tower

GWP = 'GWP'
FFC = 'FFC'
WU = 'WU'
HTC = 'HTC'
HTNC = 'HTNC'
ETOX = 'ETOX'
ACD = 'ACD'
OZD = 'OZD'
POCP = 'POCP'

def load_process_settings():
    settings = bst.settings
    chilled_water_agent = settings.get_cooling_agent('chilled_water')
    chilled_water_agent.T = 273.15 + 5 # Optimistic assumption
    settings.define_impact_indicator(GWP, 'kg*CO2e')
    settings.define_impact_indicator(FFC, 'MJ')
    settings.define_impact_indicator(WU, 'kg')#'m3')
    settings.define_impact_indicator(HTC,'kg') #'CTUh')
    settings.define_impact_indicator(HTNC,'kg') #'CTUh')
    settings.define_impact_indicator(ETOX,'kg') #'CTUe')
    settings.define_impact_indicator(ACD, 'kg')#'mol H+ eq')
    settings.define_impact_indicator(OZD, 'kg')#'CFC11 eq')
    settings.define_impact_indicator(POCP,'kg') #'NMVOC eq')
    for i in settings.heating_agents[:-1]:
        #from EF2.0 Process steam from natural gas, production mix, at heat plant, technology mix regarding firing and flue gas cleaning, MJ, 90% efficiency - RNA
        settings.set_utility_agent_CF(
            i.ID, GWP, 0.084797, basis ='MJ' # 
        )
        settings.set_utility_agent_CF(
            i.ID, FFC, 1.3747, basis='MJ' # [MJ / kg NG] From Open-LCA Environmental Footprint 2.0
        )
        settings.set_utility_agent_CF(i.ID, WU, .00053, basis='MJ')
        settings.set_utility_agent_CF(i.ID, HTC, 1.6938e-8, basis='MJ')
        settings.set_utility_agent_CF(i.ID, HTNC, 4.00298E-10, basis='MJ')
        settings.set_utility_agent_CF(i.ID, ETOX, 0.17130, basis='MJ')
        settings.set_utility_agent_CF(i.ID, ACD, 6.7669E-5, basis='MJ')
        settings.set_utility_agent_CF(i.ID, OZD, 2.1524E-15, basis='MJ')
        settings.set_utility_agent_CF(i.ID, POCP, 8.07768E-5, basis='MJ')  
        
        
        
    settings.set_electricity_CF(
        GWP, 0.60364 #from EF 2.0 #0.38 # [kg*CO2*eq / kWhr] From GREET 2020; NG-Fired Simple-Cycle Gas Turbine CHP Plant
    )
    settings.set_electricity_CF(
        FFC, 9.11521 # [MJ / kWhr] From EF 2.0
    )
    settings.set_electricity_CF(
        WU, 0.10916 # [M3 / kWhr] From EF 2.0
    )
    settings.set_electricity_CF(
        HTC, 3.75277e-8 # [ / kWhr] From EF 2.0
    )
    settings.set_electricity_CF(
        HTNC, 8.3305e-9 # [ / kWhr] From EF 2.0
    )
    settings.set_electricity_CF(
        ETOX, 0.39034 # [ / kWhr] From EF 2.0
    )
    settings.set_electricity_CF(
        ACD, 0.00211 # [ / kWhr] From EF 2.0
    )
    settings.set_electricity_CF(
        OZD, 1.15189e-10 # [ / kWhr] From EF 2.0
    )
    settings.set_electricity_CF(
        POCP, 0.00078 # [MJ / kWhr] From EF 2.0
    )
    
    
    settings.CEPCI = 816.0 # 2022
    settings.electricity_price = 0.07
    # settings.set_utility_agent_CF(
    #     # Assuming cooling tower uses 0.00454 kWh / kg of recirculated water
    #     'cooling_water', GWP, 0.00454 * 0.38, 'kg'
    # )
    # settings.set_utility_agent_CF(
    #     # Assuming chilled water uses cooling water plus electricity (3400*0.7457 / 14 / 418400 kWh / kJ).
    #     'chilled_water', GWP, -0.38 * (3400*0.7457 / 14*4184000), 'kJ'
    # )
    hps = settings.get_heating_agent("high_pressure_steam")
    hps.heat_transfer_efficiency = 0.85
    hps.regeneration_price = 0.08064
    hps.T = 529.2
    hps.P = 44e5
    mps = settings.get_heating_agent("medium_pressure_steam")
    mps.heat_transfer_efficiency = 0.90
    mps.regeneration_price = 0.07974
    mps.T = 480.3
    mps.P = 18e5
    lps = settings.get_heating_agent("low_pressure_steam")
    lps.heat_transfer_efficiency = 0.95
    lps.regeneration_price = 0.06768
    lps.T = 428.6
    lps.P = 55e4
    
def load_STRAP_MSW_process_settings():
    settings = bst.settings
    settings.define_impact_indicator(GWP, 'kg*CO2e')
    settings.define_impact_indicator(FFC, 'MJ')
    settings.define_impact_indicator(WU, 'L')
    settings.CEPCI = 816.0 # 2022
    settings.electricity_price = 0.07
    hps = settings.get_heating_agent("high_pressure_steam")
    hps.heat_transfer_efficiency = 0.85
    hps.regeneration_price = 0.08064
    hps.T = 529.2
    hps.P = 44e5
    mps = settings.get_heating_agent("medium_pressure_steam")
    mps.heat_transfer_efficiency = 0.90
    mps.regeneration_price = 0.07974
    mps.T = 480.3
    mps.P = 18e5
    lps = settings.get_heating_agent("low_pressure_steam")
    lps.heat_transfer_efficiency = 0.95
    lps.regeneration_price = 0.06768
    lps.T = 428.6
    lps.P = 55e4