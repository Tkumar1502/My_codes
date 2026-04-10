# -*- coding: utf-8 -*-
"""
Created on Wed Feb 12 11:22:53 2025

@author: yoelr
"""
from pint import Quantity as Q_
from thermosteam import *
from math import pi

Tc0 = Q_(45 + 273.15, 'K')
TH = Q_(273.15 + 190, 'K')
Cp = Q_(Chemical('Toluene').Cp('l', Tc0.magnitude), 'kJ / kg / K')
U = Q_(0.5, 'kW / m2 / K')
rho = Q_(Chemical('Toluene').rho('l', Tc0.magnitude, 101325), 'kg / m3')
F = rho * Q_(60 / 20, 'L / s')
A = pi * Q_(8, 'in') ** 2 / 4
T_cf = U * A * (TH - Tc0) / F / Cp