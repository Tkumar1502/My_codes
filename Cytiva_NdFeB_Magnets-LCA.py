# -*- coding: utf-8 -*-
"""
Created on Tue Apr 21 15:31:58 2026

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

bst.settings.CEPCI = 836.9

process = strap.MagnetRecovery(
    processing_capacity = 325,      #tons
    sell_leftover_plastic = True,
    simulate=True)

process.tea.labor_cost = 580000 + process.HS.total_salary
process.S2.outs[1].ID = 'NdFeB Magnets'
process.s22.ID = 'Impurities'

products = [process.NdFeB_Magnets, process.HDPE_resins]
process.products[:] = [process.HDPE_resins, process.NdFeB_Magnets]

process.HDPE_resins.price = 1.20
process.NdFeB_Magnets.price = 100

process.system.diagram()

#LCA
inventory_table = bst.report.lca_inventory_table()