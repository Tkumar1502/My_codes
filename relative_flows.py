# -*- coding: utf-8 -*-
"""
Created on Fri Oct 24 17:18:51 2025

@author: yoelr
"""
import numpy as np
import pandas as pd
import biosteam as bst
import seaborn as sns
import matplotlib.pyplot as plt
from plastics import strap
from warnings import filterwarnings
filterwarnings('ignore')

def relative_flows(scenario=None):
    pm = strap.STRAPMSWProcess(scenario=scenario)
    flows = dict(
        MRF_residue = pm.feedstock.F_mass,
        ethanol = pm.ethanol.F_mass,
        solvent = pm.M101.outs[0].F_mass,
        resin = pm.resin.F_mass,
        wastewater = pm.M301.outs[0].F_mass,
    )
    MRF_residue = flows['MRF_residue']
    for i, j in flows.items():
        flows[i] = round(j / MRF_residue, 3)
    return pd.DataFrame([*flows.values()], [*flows.keys()])
        
if __name__ == '__main__':
    print('baseline')
    print('--------')
    print(relative_flows('baseline'))
    print()
    print('potential')
    print('---------')
    print(relative_flows('potential'))