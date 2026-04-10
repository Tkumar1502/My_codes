# -*- coding: utf-8 -*-
import biosteam as bst
from plastics import strap
import numpy as np

categories = (
    'GWP',
    'HTC',
    'HTNC',
    'ETOX',
)
GWP_INDEX, HTC_INDEX, HTNC_INDEX, ETOX_INDEX = range(4)

def create_process(facilities=False):
    # STRAP should define impact indicators internally
    # bst.settings.define_impact_indicator('GWP', units='kg*CO2e')
    # ...
    process = strap.BaselineSTRAPProcess(
        precipitation_temperature_format ='constant',
        target_plastic=('PE'),
        target_plastic_percent=(80), # %
        solvent=('Toluene'),
        processing_capacity=8e+03, # MT/yr
        sell_leftover_plastic=False,
        burn_leftover_plastic=False,
        facilities=facilities,
        save=True, # Use cache
    )
    process.get_indicators = lambda: np.array([process.GWP(), process.HTC(), process.HTNC(), process.ETOX()])
    process.set_electricity_CF = lambda value: [bst.settings.set_electricity_CF(i, value, basis='MJ') for i in categories]
    process.set_natural_gas_CF = lambda value: [process.natural_gas.set_CF(i, value) for i in categories]
    return process
    

def test_electricity():
    process = create_process()
    process.set_electricity_CF(0.5)
    indicators_low = process.get_indicators()
    process.set_electricity_CF(50)
    indicators_high = process.get_indicators()
    assert (indicators_high >= 1.1 * indicators_low).all()

def test_natural_gas():
    process = create_process(facilities=True)
    process.set_natural_gas_CF(0.5)
    indicators_low = process.get_indicators()
    process.set_natural_gas_CF(50)
    indicators_high = process.get_indicators()
    assert (indicators_high >= 1.1 * indicators_low).all()

def test_cost():
    process = create_process(facilities=False)
    TCI_low = process.tea.TCI
    process = create_process(facilities=True)
    TCI_high = process.tea.TCI
    assert (TCI_high >= 1.1 * TCI_low).all()
    
if __name__ == '__main__':
    test_electricity()
    test_natural_gas()
    test_cost()