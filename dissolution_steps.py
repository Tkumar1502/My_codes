# -*- coding: utf-8 -*-
"""
"""
import thermosteam as tmo
import numpy as np


class DissolutionStep:
    __slots__ = (
        'plastic', # Plastic ID
        'dissolved_plastic', # Dissolved plastic ID
        'solvent', # Solvent ID
        'reaction', # Dissolution reaction
        'capacity', # Maximum plastic concentration in solvent before viscosity limitations [wt / vol]
        'solvent_content', # Solvent content after centrifugation
        'T', # Dissolution temperature
        'tau', # Dissolution time
    )
    
    def __init__(self,
            plastic: str,
            dissolved_plastic: str,
            solvent: str,
            reaction: tmo.Reaction,
            capacity: float,
            solvent_content: float,
            T: float,
            tau: float,
        ):
        self.plastic = plastic
        self.dissolved_plastic = dissolved_plastic
        self.solvent = solvent
        self.reaction = reaction
        self.capacity = capacity
        self.solvent_content = solvent_content
        self.T = T
        self.tau = tau
    
    def __eq__(self, other):
        return type(self) is type(other) and (self.solvent, self.plastic) == (other.solvent, other.plastic)
    
    def __hash__(self):
        return hash((self.solvent, self.plastic))
        
    def __repr__(self):
        return f"{type(self).__name__}({self.plastic}, {self.dissolved_plastic}, {self.solvent}, {self.reaction}, {self.capacity}, {self.T}, {self.tau})"
    

def EVOH_DMSOWater_dissolution():
    chemicals = tmo.settings.chemicals
    rho_DMSO = chemicals.DMSO.rho(T=298.15, P=101325, phase='l')
    rho_water = chemicals.Water.rho(T=298.15, P=101325, phase='l')
    mass_DMSO = 0.6 * rho_DMSO
    mass_water = 0.4 * rho_water
    total_mass = mass_DMSO + mass_water
    x_DMSO = mass_DMSO / total_mass
    x_water = mass_water / total_mass
    chemicals.define_group('DMSOWater', ('DMSO', 'Water'), np.array([x_DMSO, x_water]), wt=True)
    return DissolutionStep(
        'EVOH', 'EVOHoligomer', 'DMSOWater', tmo.Reaction('EVOH -> EVOHoligomer', 'EVOH', X=1.0, basis='wt'), 0.03, 
        0.4, 383.15, 4
    )

def PC_THF_dissolution():
    return DissolutionStep(
        'PC', 'PColigomer', 'THF', tmo.Reaction('PC -> PColigomer', 'PC', X=1.0, basis='wt'), 0.03, 
        0.5, 336.15, 0.5
    )

def PE_Toluene_dissolution():
    return DissolutionStep(
        'PE', 'PEoligomer', 'Toluene', tmo.Reaction('PE -> PEoligomer', 'PE', X=1.0, basis='wt'), 0.03, 
        0.5, 368.15, 0.5
    )

def PEPP_Xylene_dissolution():
    return DissolutionStep(
        'PEPP', 'PEPPoligomer', 'Xylene', 
        tmo.Reaction('PEPP -> PEPPoligomer', 'PEPP', X=1.0, basis='wt'), 0.03, 
        0.5, 130 + 273.15, 0.5
    )

#testing this
def PVDF_DMF_dissolution():
    return DissolutionStep(
        'PVDF', 'PVDFoligomer', 'DMF', 
        tmo.Reaction('PVDF -> PVDFoligomer', 'PVDF', X=1.0, basis='wt'), 0.03, 
        0.5, 130 + 273.15, 0.5
    )
    