# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2023, Yoel Cortes-Pena <yoelcortes@gmail.com>, Priscilla Lee <thescillalee@icloud.com>
#
# Adsorption modeling was contributed by Pricilla Lee
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
import biosteam as bst
import thermosteam as tmo
import numpy as np
from scipy.integrate import solve_ivp
from math import sqrt, pi, ceil
from biosteam.units.design_tools import PressureVessel
from biosteam.units.decorators import cost
from thermosteam import separations as sep
from biosteam.units.design_tools import (
    CEPCI_by_year,
    cylinder_diameter_from_volume,
    cylinder_area,
    size_batch,
)
from math import exp, log
from numpy.testing import assert_allclose
from biosteam.units import design_tools as design
from biorefineries.cellulosic.units import (
    SeedTrain,
    CoFermentation
)
from flexsolve import secant, wegstein
import matplotlib.pyplot as plt
from numba import njit

JacketedSurgeTank = bst.StorageTank  # TODO: Add cost of jacket

# TODO: Find cost for microfilter, the following only works for nanofilter.
# Membrane separation processes. Perry's Chemical Engineer's Handbook 7th Edition.

CEPCI2022 = 816.0  # 2022
Europe_investment_site_factor = 1.2
euro_to_dollar = 1.04
installation_cost = 115000  # euro
volume_treated = 10 * 24 * 30 * 18  # m3
membrane_area = 27.5  # m2
membrane_cost = 95  # euro / m2
cleaning_cost = 50  # euro / m2
maintenance_cost = installation_cost * 0.02
yearly_operating_hours = 8000
operating_cost = membrane_area * \
    (membrane_cost + cleaning_cost) + maintenance_cost * \
    yearly_operating_hours / (18 * 30 * 24)
operating_cost_per_volume_treated = operating_cost / volume_treated
US_operating_cost = euro_to_dollar * \
    operating_cost_per_volume_treated / Europe_investment_site_factor
US_installation_cost = euro_to_dollar * \
    installation_cost / Europe_investment_site_factor
CEPCI2013 = bst.units.design_tools.CEPCI_by_year[2013]
electricity_cost = 2e3  # euro / yr
electricity_price = 0.04  # euro / kWh
electricity_demand = electricity_cost / \
    yearly_operating_hours / electricity_price

# %% Preprocessing

# TODO: Bare module factor for all equipment is 2, update later.
# TODO: Check power with Srikar
@cost('Flow rate', units='ton/hr', cost=0, CE=567.3, n=1, S=1, kW=0, BM=1, lifetime=30)
class HandSorting(bst.Splitter): pass

@cost('Flow rate', units='ton/hr', cost=87000, CE= 800.8, n=0.6, S=2, kW=0.87*2, BM=2, lifetime=15) # Assumed eddy current kw/ton is same as Vecoplan, 800.8 is the CEPCI for 2023
class EddyCurrent(bst.Splitter): pass

@cost('Flow rate', units='ton/hr', cost=38000, CE=576.1, n=0.6, S=7.8, kW=30.1*7.8, BM=2, lifetime=10) # Vecoplan Kw is 30kwh/dryton , lifetime is 10 hours, CEPCI is 576.1 in 2014
class Vecoplan(bst.Unit): pass

@cost('Flow rate', units='ton/hr', cost=60000, CE=800.8, n=1, S=7.8, kW=0.38*7.8, BM=2, lifetime=15) # Got these from Pralhad's values. Not sure about CEPCI, have used 2023 values
class Magnet(bst.Splitter): pass

@cost('Flow rate', units='ton/hr', cost=2250000, CE=816, n=0.6, S=11.71, kW=35.5*11.71, BM=2, lifetime=30) # Vecoplan Kw is 30kwh/dryton , lifetime is 10 hours, CEPCI is 576.1 in 2014
class Crumble(bst.Unit):
    _N_ins=1
    _N_outs=3
    def _run(self):
        overs,unders,residue = self.outs
        feed = self.ins[0]
        residue.copy_flow(feed)
        residue.imol['Unders','Overs'] = 0
        unders.imol['Unders'] = feed.imol['Unders']
        overs.imol['Overs'] = feed.imol['Overs']
        unders.phase = 's'
        overs.phase = 's'
        residue.phase = 's'

# %% Ethanol production


class SeedTrain(SeedTrain):

    def _init(self, reactions=None, saccharification=False):
        self.saccharification = saccharification
        chemicals = self.chemicals
        self.reactions = reactions or bst.ParallelReaction([
            #   Reaction definition                   Reactant    Conversion
            bst.Reaction('Glucose -> 2 Ethanol + 2 CO2',
                         'Glucose',   0.9000, chemicals),
            bst.Reaction('3 Xylose -> 5 Ethanol + 5 CO2',
                         'Xylose',    0.8000, chemicals),
            bst.Reaction('Glucose -> Cellmass',                'Glucose',
                         0.0473, chemicals, correct_mass_balance=True),
            bst.Reaction('Xylose -> Cellmass',                 'Xylose',
                         0.0421, chemicals, correct_mass_balance=True),
        ])

    def _setup(self):
        super()._setup()
        self.outs[0].phase = 'g'

    def _run(self):
        vent, effluent = self.outs
        effluent.mix_from(self.ins, energy_balance=False)
        self.reactions.force_reaction(effluent)
        effluent.mol.remove_negatives()
        effluent.T = self.T
        vent.empty()
        vent.copy_flow(effluent, ('CO2', 'O2', 'N2'), remove=True)


class CoFermentation(CoFermentation):

    def _init(self, tau=36, N=None, V=3785.4118, T=305.15, P=101325,
              Nmin=2, Nmax=36, cofermentation=None):
        bst.BatchBioreactor._init(self, tau, N, V, T, P, Nmin, Nmax)
        self.P = P
        chemicals = self.chemicals
        self.loss = None
        self.cofermentation = cofermentation or bst.ParallelReaction([
            #   Reaction definition                   Reactant    Conversion
            bst.Reaction('Glucose -> 2 Ethanol + 2 CO2',
                         'Glucose',   0.9500, chemicals),
            bst.Reaction('3 Xylose -> 5 Ethanol + 5 CO2',
                         'Xylose',    0.8500, chemicals),
            bst.Reaction('Glucose -> Cellmass',                'Glucose',
                         0.05, chemicals, correct_mass_balance=True),
            bst.Reaction('Xylose -> Cellmass',                 'Xylose',
                         0.05, chemicals, correct_mass_balance=True),
        ])

        if 'CSL' in chemicals:
            self.CSL_to_constituents = bst.Reaction(
                'CSL -> 0.5 H2O + 0.25 LacticAcid + 0.25 Protein', 'CSL', 1.0000, chemicals, basis='wt',
            )
            self.CSL_to_constituents.basis = 'mol'
        else:
            self.CSL_to_constituents = None
        self.lipid_reaction = None


# %% STRAP

@cost('Flow rate', 'Screw feeder', units='ft3/hr',
      lb=400, ub=10e4, CE=567, cost=1096, n=0.22)
class VacuumDryer(design.PressureVessel, bst.Unit):
    """
    Create a vacuum dryer that dries solids by vacuum, and heat.

    Parameters
    ----------
    ins : 
        * [0] Wet solids.
    outs : 
        * [0] Dried solids
        * [1] Hot gas
    split : dict[str, float]
        Component splits to hot gas (stream [1]).
    H : float, optional
        Specific evaporation rate [kg/hr/ft3]. Defaults to 0.566. 
    length_to_diameter : float, optional
        Note that the drum is horizontal. Defaults to 25.
    T : float, optional
        Operating temperature [K]. Defaults to 343.15.
    moisture_content : float
        Moisture content of solids [wt / wt]. Defaults to 0.10.

    Notes
    -----
    The default parameter values are based on heuristics for drying 
    dried distillers grains with solubles (DDGS). However, values should 
    be updated to match design for vacuum drier.

    """
    auxiliary_unit_names = (
        'vacuum_system', 'condenser', 'pump'
    )
    _units = {
        'Evaporation': 'kg/hr',
        'Volume': 'ft3',
    }
    _N_ins = 1
    _N_outs = 2

    @property
    def isplit(self):
        """[ChemicalIndexer] Componentwise split of feed to 0th outlet stream."""
        return self._isplit

    @property
    def split(self):
        """[Array] Componentwise split of feed to 0th outlet stream."""
        return self._isplit.data

    def _init(self, split, H=20., length_to_diameter=25, T=380.15, P=101325 * 0.1,
              moisture_content=1e-6, moisture_ID=None):
        self._isplit = self.chemicals.isplit(split)
        self.P = P
        self.T = T
        self.H = H
        self.length_to_diameter = length_to_diameter
        self.moisture_content = moisture_content
        self.moisture_ID = moisture_ID
        self.vessel_material = 'Carbon steel'
        condenser = self.auxiliary(
            'condenser', bst.HXutility, ins=[''],
            V=0,
        )
        self.auxiliary(
            'pump', bst.HXutility, ins=condenser-0, outs=[self.outs[1]],
            V=0,
        )

    def _run(self):
        wet_solids, = self.ins
        dry_solids, condensate, = self.outs
        hot_gas, = self.condenser.ins
        wet_solids.split_to(hot_gas, dry_solids, self.split)
        sep.adjust_moisture_content(
            dry_solids, hot_gas, self.moisture_content, self.moisture_ID)
        hot_gas.P = self.P
        hot_gas.phase = 'g'
        design_results = self.design_results
        design_results['Evaporation'] = hot_gas.F_mass
        dry_solids.T = hot_gas.T = self.T
        self.condenser.run()
        self.pump.run()

    def _design(self):
        self._decorated_design()
        feed = self.ins[0]
        F_mass = 0.0006124 * feed.F_mass  # lb/s
        length_to_diameter = self.length_to_diameter
        design_results = self.design_results
        design_results['Volume'] = volume = design_results['Evaporation'] / self.H
        design_results['Diameter'] = diameter = cylinder_diameter_from_volume(
            volume, length_to_diameter)
        design_results['Length'] = length = diameter * length_to_diameter
        design_results.update(
            self._horizontal_vessel_design(self.P, diameter, length)
        )
        self.power_utility(
            0.0146 * F_mass**0.85 * length * 0.7457
        )
        self.add_heat_utility(self.H_out - self.H_in, self.T)
        self.vacuum_system = bst.VacuumSystem(self)
        self.condenser._design()
        self.pump._design()

    def _cost(self):
        self._decorated_cost()
        design_results = self.design_results
        dct = self.baseline_purchase_costs
        dct.update(
            self._horizontal_vessel_purchase_cost(
                design_results['Weight'],
                design_results['Diameter'],
                design_results['Length']
            )
        )
        dct['Heating jacket'] = dct.pop('Horizontal pressure vessel')
        self.condenser._cost()
        self.pump._cost()


class MeltDegasser(bst.Flash):
    pass  # TODO: Make this into a system with 3-stages, a condenser, and a single vacuum system

############# my attempt ###############
# @cost('flowrate', 'stage 1 flash', units = 'kg/hr',CE=CEPCI2022,
#       cost=30000, S=250, n=0.8, BM=4)
# @cost('flowrate', 'stage 2 flash', units = 'kg/hr',CE=CEPCI2022,
#       cost=30000, S=250, n=0.8, BM=4)
# @cost('flowrate', 'stage 3 flash', units = 'kg/hr',CE=CEPCI2022,
#       cost=30000, S=250, n=0.8, BM=4)
# @cost('flowrate', 'stage 1 melt pump', units = 'kg/hr',CE=CEPCI2022,
#       cost=20000, S=250, n=0.8, BM=4)
# @cost('flowrate', 'stage 2 melt pump', units = 'kg/hr',CE=CEPCI2022,
#       cost=20000, S=250, n=0.8, BM=4)
# @cost('flowrate', 'stage 3 melt pump', units = 'kg/hr',CE=CEPCI2022,
#       cost=20000, S=250, n=0.8, BM=4)
# @cost('flowrate', 'stage 1 hx', units = 'kg/hr',CE=CEPCI2022,
#       cost=10000, S=250, n=0.8, BM=4)
# @cost('flowrate', 'stage 2 hx', units = 'kg/hr',CE=CEPCI2022,
#       cost=10000, S=250, n=0.8, BM=4)
# @cost('flowrate', 'stage 2 hx', units = 'kg/hr',CE=CEPCI2022,
#       cost=10000, S=250, n=0.8, BM=4)
# class MeltDegasser(bst.Flash):
#     pass  # TODO: Make this into a system with 3-stages, a condenser, and a single vacuum system
########################################


@cost('Flow rate', 'additional melt devolitalization', units='kg/hr', CE=CEPCI2022,
      cost=400e3 ,BM=3, n=0.6, S=500, kW=100) 
@cost('Flow rate', 'melt pump', units='kg/hr', CE=CEPCI2022,
      cost=50e3 ,BM=3, n=0.6, S=500, kW=15)
@cost('Flow rate', units='kg/hr', CE=CEPCI2022,
      cost=800e3 ,BM=2.4, n=0.6, S=500, kW=294.2)  # TODO: Double check electricity consumption.
class ScrewPressDegasser(bst.Unit):
    """
    Create a degasser that dries solids by expression, vacuum, and heat.

    Parameters
    ----------
    ins : 
        * [0] Wet solids.
    outs : 
        * [0] Dried solids
        * [1] Hot gas
    split : dict[str, float]
        Component splits to hot gas (stream [1]).
    T : float, optional
        Operating temperature [K]. Defaults to 343.15.
    moisture_content : float
        Moisture content of solids [wt / wt]. Defaults to 0.10.

    """
    _units = {
        'Evaporation': 'kg/hr',
        'Volume': 'ft3',
    }
    _N_ins = 1
    _N_outs = 2
    _vessel_material = '-'

    @property
    def isplit(self):
        """[ChemicalIndexer] Componentwise split of feed to 0th outlet stream."""
        return self._isplit

    @property
    def split(self):
        """[Array] Componentwise split of feed to 0th outlet stream."""
        return self._isplit.data

    def _init(self, split, T=380.15, P=101325 * 0.05,
              moisture_content=1e-6, moisture_ID=None):
        self._isplit = self.chemicals.isplit(split)
        self.P = P
        self.T = T
        self.moisture_content = moisture_content
        self.moisture_ID = moisture_ID

    def _run(self):
        wet_solids, = self.ins
        dry_solids, condensate, = self.outs
        wet_solids.split_to(condensate, dry_solids, self.split)
        sep.adjust_moisture_content(
            dry_solids, condensate, self.moisture_content, self.moisture_ID
        )
        condensate.P = self.P
        condensate.phase = 'l'
        dry_solids.T = condensate.T = self.T


#: TODO: Remove old cost correlation (or maybe compare)
# @cost('Flow rate', units='lb/hr', CE=567, lb=150, ub=12000, BM=1.39,
#       f=lambda S: exp((11.0991 - 0.3580*log(S) + 0.05853*log(S)**2)),
#       kW=0.0372) # kWh per biomass
# class ScrewPressDegasser(design.PressureVessel, bst.Unit):
#     """
#     Create a degasser that dries solids by expression, vacuum, and heat.

#     Parameters
#     ----------
#     ins :
#         * [0] Wet solids.
#     outs :
#         * [0] Dried solids
#         * [1] Hot gas
#     split : dict[str, float]
#         Component splits to hot gas (stream [1]).
#     H : float, optional
#         Specific evaporation rate [kg/hr/ft3]. Defaults to 0.566.
#     length_to_diameter : float, optional
#         Note that the drum is horizontal. Defaults to 25.
#     T : float, optional
#         Operating temperature [K]. Defaults to 343.15.
#     moisture_content : float
#         Moisture content of solids [wt / wt]. Defaults to 0.10.

#     Notes
#     -----
#     The default parameter values are based on heuristics for drying
#     dried distillers grains with solubles (DDGS). However, values should
#     be updated to match design for vacuum drier.

#     """
#     auxiliary_unit_names = (
#         'vacuum_system', 'condenser', 'pump'
#     )
#     _units = {
#         'Evaporation': 'kg/hr',
#         'Volume': 'ft3',
#     }
#     _N_ins = 1
#     _N_outs = 2

#     @property
#     def isplit(self):
#         """[ChemicalIndexer] Componentwise split of feed to 0th outlet stream."""
#         return self._isplit
#     @property
#     def split(self):
#         """[Array] Componentwise split of feed to 0th outlet stream."""
#         return self._isplit.data

#     def _init(self, split, H=20., length_to_diameter=25, T=380.15, P=101325 * 0.05,
#               moisture_content=1e-6, moisture_ID=None):
#         self._isplit = self.chemicals.isplit(split)
#         self.P = P
#         self.T = T
#         self.H = H
#         self.vessel_material = 'Carbon steel'
#         self.length_to_diameter = length_to_diameter
#         self.moisture_content = moisture_content
#         self.moisture_ID = moisture_ID
#         condenser = self.auxiliary(
#             'condenser', bst.HXutility, ins=[''], V=0,
#         )
#         self.auxiliary(
#             'pump', bst.HXutility, ins=condenser-0, outs=[self.outs[1]],
#             V=0,
#         )

#     def _run(self):
#         wet_solids, = self.ins
#         dry_solids, condensate, = self.outs
#         hot_gas, = self.condenser.ins
#         wet_solids.split_to(hot_gas, dry_solids, self.split)
#         sep.adjust_moisture_content(dry_solids, hot_gas, self.moisture_content, self.moisture_ID)
#         hot_gas.P = self.P
#         hot_gas.phase = 'g'
#         design_results = self.design_results
#         design_results['Evaporation'] = hot_gas.F_mass
#         dry_solids.T = hot_gas.T = self.T
#         self.condenser.run()
#         self.pump.run()

#     def _design(self):
#         self._decorated_design()
#         length_to_diameter = self.length_to_diameter
#         design_results = self.design_results
#         design_results['Volume'] = volume = design_results['Evaporation'] / self.H
#         design_results['Diameter'] = diameter = cylinder_diameter_from_volume(volume, length_to_diameter)
#         design_results['Length'] = length = diameter * length_to_diameter
#         design_results.update(
#             self._horizontal_vessel_design(self.P, diameter, length)
#         )
#         self.add_heat_utility(self.H_out - self.H_in, self.T)
#         self.vacuum_system = bst.VacuumSystem(self)
#         self.condenser._design()
#         self.pump._design()

#     def _cost(self):
#         self._decorated_cost()
#         design_results = self.design_results
#         dct = self.baseline_purchase_costs
#         dct.update(
#             self._horizontal_vessel_purchase_cost(
#                 design_results['Weight'],
#                 design_results['Diameter'],
#                 design_results['Length']
#             )
#         )
#         dct['Heating jacket'] = dct.pop('Horizontal pressure vessel')
#         self.condenser._cost()
#         self.pump._cost()


class BatchPlasticDissolution(bst.STR):
    _N_ins = 2
    _N_outs = 1
    tau_0_default = 0.1

    def _init(self, dissolution_step, minimum_solvent_to_feed=4, **kwargs):
        kwargs['tau'] = dissolution_step.tau
        super()._init(**kwargs)
        self.dissolution_step = dissolution_step
        self.minimum_solvent_to_feed = minimum_solvent_to_feed

    def _run(self):
        plastic, solvent = self.ins
        effluent, = self.outs
        ds = self.dissolution_step
        F_plastic = sum(
            [i.imass[ds.plastic, ds.dissolved_plastic].sum() for i in self.ins]
        )
        solvent.imass[ds.solvent] = max(
            F_plastic / ds.capacity, self.minimum_solvent_to_feed * plastic.F_mass
        )
        effluent.mix_from([plastic, solvent], energy_balance=True)
        ds.reaction(effluent)


class PrecipitationTank(bst.BatchCrystallizer):
    tau_0 = 0.1

    def _init(self, precipitation_step, **kwargs):
        kwargs['tau'] = precipitation_step.tau
        self.precipitation_step = precipitation_step
        super()._init(**kwargs)
        self.V = 1e3  # TODO: Base this on common sizes

    @property
    def T(self):
        return self.precipitation_step.T

    @T.setter
    def T(self, T):
        self.precipitation_step.T = T

    def _run(self):
        feed, = self.ins
        outlet, = self.outs
        outlet.copy_like(feed)
        ps = self.precipitation_step
        outlet.T = ps.T
        ps.precipitate(outlet)


class PseudoContinuousPlasticDissolutionTank(PressureVessel, bst.Unit):
    _N_ins = 3
    _N_outs = 3
    _ins_size_is_fixed = False
    _outs_size_is_fixed = False
    refill_time = 0.5
    auxiliary_unit_names = (
        'compressor', 'air_heat_exchanger', 'heaters',
    )

    def _init(self,
              dissolution_steps,  # Iterable[DissolutionStep]
              # m / hr; typical velocities are 4 to 14.4 m /hr for liquids; Adsorption basics Alan Gabelman (2017) Adsorption basics Part 1. AICHE
              solvent_superficial_velocity=7.2,
              # Mid point in velocity range for gasses, m / hr; Alan Gabelman (2017) Adsorption basics Part 1. AICHE
              air_superficial_velocity=1332,
              vessel_material='Stainless steel 316',
              vessel_type='Vertical',
              # Optional[float] Necessary for sizing length
              void_fraction=0.47,
              # Optional[float] Additional length of a column to account for mass transfer limitations (due to unused bed). Defaults to +2 ft per column.
              length_unused=1.219,
              # Optional[float] Time for drying after regeneration
              drying_time=0.05,
              # Optional[float] Defaults to maximum dissolution temperature +10 degC
              T_air=None,
        ):
        inlets = self.ins
        n_inlets = len(inlets)
        n_solvents = n_inlets - 2
        if n_solvents < 1:
            raise ValueError('no solvent inlet given; number of inlets must be greater than 2 '
                             '(the first inlet is the plastic, the last is air, and everything in between are the solvents)')
        outlets = self.outs
        n_outlets = len(outlets)
        n_missing = n_outlets - n_inlets
        if n_missing:
            outlets.extend([bst.Stream(thermo=self.thermo)
                           for i in range(n_missing)])
        if T_air is None:
            T_air = max([i.T for i in dissolution_steps]) + 50
        if len(dissolution_steps) != n_solvents:
            raise ValueError(
                'the number of dissolution groups does not match the number of solvents')
        self.solvent_superficial_velocity = solvent_superficial_velocity
        self.air_superficial_velocity = air_superficial_velocity
        self.vessel_material = vessel_material
        self.vessel_type = vessel_type
        self.length_unused = length_unused
        self.void_fraction = void_fraction
        self.T_air = T_air
        self.drying_time = drying_time
        self.dissolution_steps = dissolution_steps
        self.heaters = []
        for inlet, ds in zip(inlets[1:-1], dissolution_steps):
            self.auxiliary(
                'heaters', bst.HXutility, inlet, T=ds.T
            )
        compressor = self.auxiliary(
            'compressor', bst.IsentropicCompressor, self.air, eta=0.85, P=10 * 101325
        )
        self.auxiliary(
            'air_heat_exchanger', bst.HXutility, compressor-0,
            T=T_air, rigorous=False,
        )

    @property
    def air(self):
        return self._ins[-1]

    @property
    def vent(self):
        return self._outs[-1]

    @property
    def hot_compressed_air(self):
        return self.air_heat_exchanger.outs[0]

    def _run(self):
        plastic, *solvents, air = self.ins
        air.phase = 'g'
        plastic_outlet, *spent_solvents, spent_air = self.outs
        solvent_superficial_velocity = self.solvent_superficial_velocity
        air_superficial_velocity = self.air_superficial_velocity
        dissolution_steps = self.dissolution_steps
        drying_time = self.drying_time
        refill_time = self.refill_time
        dissolution_time = sum([i.tau for i in dissolution_steps])
        n_dissolutions = len(dissolution_steps)
        total_drying_time = n_dissolutions * drying_time
        self.cycle_time = cycle_time = dissolution_time + total_drying_time + refill_time
        plastic_outlet.copy_flow(plastic)
        plastic_outlet.T = dissolution_steps[-1].T
        heaters = self.heaters
        volumetric_solvent_flows = []
        for i, ds in enumerate(dissolution_steps):
            heater = heaters[i]
            solvent = solvents[i]
            spent_solvent = spent_solvents[i]
            F_plastic = plastic_outlet.imass[ds.plastic]
            plastic_outlet.imass[ds.plastic] = 0.  # Remove plastic film
            solvent.imass[ds.solvent] = F_plastic / ds.capacity
            heater.run()
            heated_solvent = heater.outlet
            spent_solvent.copy_like(heated_solvent)
            spent_solvent.imass[ds.plastic] = F_plastic
            ds.reaction(spent_solvent)
            volumetric_solvent_flows.append(
                heated_solvent.F_vol * cycle_time / ds.tau)
        required_diameters = [
            2 * sqrt(i / (solvent_superficial_velocity * pi)) for i in volumetric_solvent_flows]
        self.diameter = diameter = max(required_diameters)
        self.area = area = pi * diameter * diameter / 4
        self.length = length = (
            cycle_time * plastic.F_vol / (self.void_fraction * diameter)
        )
        self.vessel_volume = length * area
        hot_compressed_air = self.hot_compressed_air
        hot_compressed_air.T = self.air_heat_exchanger.T = self.T_air
        mean_air_flow = total_drying_time / cycle_time * \
            air_superficial_velocity * diameter
        hot_compressed_air.reset_flow(
            N2=0.78, O2=0.32, phase='g', total_flow=mean_air_flow, units='m3/hr'
        )
        air.copy_flow(hot_compressed_air)
        self.compressor._run()
        self.air_heat_exchanger._run()
        self.vent.copy_like(hot_compressed_air)

    def _design(self):
        design_results = self.design_results
        diameter = self.diameter
        length = self.length
        design_results['Number of vessels'] = 1
        design_results.update(
            self._vessel_design(
                self.hot_compressed_air.P * 0.000145038,  # Pa to psi
                diameter * 3.28084,  # m to ft
                length * 3.28084,  # m to ft
            )
        )
        for i in self.auxiliary_units:
            i._design()

    def _cost(self):
        design_results = self.design_results
        baseline_purchase_costs = self.baseline_purchase_costs
        baseline_purchase_costs.update(
            self._vessel_purchase_cost(
                design_results['Weight'],
                design_results['Diameter'],
                design_results['Length']
            )
        )
        N_reactors = design_results['Number of vessels']
        for i, j in baseline_purchase_costs.items():
            baseline_purchase_costs[i] *= N_reactors
        for i in self.auxiliary_units:
            i._cost()

# %% Pilot data from Charles


@cost('Flow rate', 'Shredder', CE=CEPCI2022, units='kg/hr',
      S=500, cost=217e3, BM=1.39, n=0.6, kW=90)
@cost('Flow rate', 'Feed belt', CE=CEPCI2022, units='kg/hr',
      S=500, cost=30e3, BM=1.74, n=0.38, kW=0.746)
@cost('Flow rate', 'Bucket elevator', CE=CEPCI2022, units='kg/hr',
      S=500, cost=15e3, BM=1.74, n=0.38, kW=0.746)
@cost('Flow rate', 'Hopper', CE=CEPCI2022, units='kg/hr',
      S=500, cost=35e3, BM=1.61, n=0.38)
@cost('Flow rate', 'Air lock', CE=CEPCI2022, units='kg/hr',
      S=30, cost=30e3, BM=4.16, n=0.5)
@cost('Flow rate', 'Ball valve', CE=CEPCI2022, units='kg/hr',
      S=30, cost=3e3, n=0.5)  # TODO: Consider removing this
class Shredding(bst.Unit):
    pass

# @cost('Solids loading', 'Centrifilter', CE=CEPCI2022, units='kg/hr',
#       S=91, cost=65e3, BM=2.03, n=0.5, kW=0.55)
# class Centrifilter(bst.SolidsSeparator):
#     _units = {'Solids loading': 'kg/hr'}

#     def _init(self, split, order=None, solids=None, moisture_content=0.40,
#               moisture_ID=None, strict_moisture_content=None):
#         bst.SolidsSeparator._init(
#             self, moisture_content=moisture_content,
#             split=split, order=order, moisture_ID=moisture_ID,
#             strict_moisture_content=strict_moisture_content
#         )
#         if solids is None:
#             self.solids = [i.ID for i in self.chemicals if i.locked_state == 's']
#         else:
#             self.solids = solids

#     def _design(self):
#         solids = self.solids
#         self.design_results['Solids loading'] = sum([s.imass[solids].sum() for s in self.ins if not s.isempty()]) # Total solids


@cost('Flow rate', 'Candle filter', CE=CEPCI2022, units='L/hr',
      # n=0.71 -> Pressure lead (Warren)
      S=400, cost=250e3, BM=2.32, n=0.5, kW=0.55)
class CandleFilter(bst.SolidsSeparator):

    def _init(self, split, order=None, solids=None, moisture_content=0.40,
              moisture_ID=None, strict_moisture_content=None):
        bst.SolidsSeparator._init(
            self, moisture_content=moisture_content,
            split=split, order=order, moisture_ID=moisture_ID,
            strict_moisture_content=strict_moisture_content
        )
        if solids is None:
            self.solids = [
                i.ID for i in self.chemicals if i.locked_state == 's']
        else:
            self.solids = solids

# %% Based on cost of pressure vessels


class Precipitator(PressureVessel, bst.Unit):

    @property
    def T(self):
        return self.precipitation_step.T

    @T.setter
    def T(self, T):
        self.precipitation_step.T = T

    @property
    def tau(self):
        return self.precipitation_step.tau

    @tau.setter
    def tau(self, tau):
        self.precipitation_step.tau = tau

    def _run(self):
        feed, = self.ins
        outlet, = self.outs
        outlet.copy_like(feed)
        ps = self.precipitation_step
        outlet.T = ps.T
        ps.precipitate(outlet)

    _units = {'Vessel volume': 'm3',
              'Batch time': 'hr',
              'Loading time': 'hr'}

    _N_ins = _N_outs = 1

    #: [float] Cleaning and unloading time (hr).
    tau_0 = 0.1

    #: [float] Fraction of filled tank to total tank volume.
    V_wf = 0.99

    def _get_design_info(self):
        return (('Cleaning and unloading time', self.tau_0, 'hr'),
                ('Working volume fraction', self.V_wf, ''))

    def _init(self,
              vessel_material='Carbon steel',
              kW=0, aspect_ratio=20, diameter=0.25, precipitation_step=None,
              ):
        self.precipitation_step = precipitation_step
        self.vessel_material = vessel_material
        self.vessel_type = 'Vertical'
        self.aspect_ratio = aspect_ratio
        self.diameter = diameter

        #: [float] Electricity usage per volume in kW/gal
        self.kW = kW

    @property
    def V(self):
        """[float] Crystallizer volume."""
        d = self.diameter
        h = self.aspect_ratio * d
        return h * d**2 / 4 * pi

    @property
    def height(self):
        return self.diameter * self.aspect_ratio

    def _design(self):
        effluent = self.effluent
        v_0 = effluent.F_vol
        tau = self.tau
        tau_0 = self.tau_0
        V_wf = self.V_wf
        design_results = self.design_results
        N = v_0 / self.V / V_wf * (tau+tau_0) + 1
        if N < 2:
            N = 2
        else:
            N = ceil(N)
        dct = size_batch(v_0, tau, tau_0, N, V_wf)
        design_results['Vessel volume'] = volume = dct.pop('Reactor volume')
        design_results.update(dct)
        design_results['Number of vessels'] = N
        self.parallel['Vertical pressure vessel'] = N
        self.parallel['Platform and ladders'] = N
        self.add_heat_utility(self.H_out - self.H_in, self.T)
        self.add_power_utility(self.kW * V_wf * volume * N)
        outlet = self.outs[0]
        P = outlet.P
        design_results.update(
            self._vessel_design(P, self.diameter * 3.28, self.height * 3.28)
        )

    def _cost(self):
        self.baseline_purchase_costs.update(
            self._vessel_purchase_cost(
                self.design_results['Weight'],
                self.design_results['Diameter'],
                self.design_results['Length'],
            )
        )

# %% More rigorous modeling

class Filter(bst.Unit):
    _N_outs = 2
    _units = {'Solids loading': 'kg/hr'}

    def _init(self,
              solutes=None,
              mesh_size=600,  # um
              moisture_content=0.50,
              solvent='Solvent',
              solid_polymer='SolidPolymer',
        ):
        if solutes is None:
            solutes = self.chemicals.chemical_group_members('Solutes')
        self.mesh_size = mesh_size
        self.moisture_content = moisture_content
        self.solvent = solvent
        self.solid_polymer = solid_polymer
        self.solutes = solutes

    def get_solids_split(self):
        mesh = self.mesh_size
        fallthrough = 1.23e-7 * mesh * mesh + 9.47e-5 * mesh - 0.078
        if fallthrough > 1: fallthrough = 1
        elif fallthrough < 0: fallthrough = 0
        split = {self.solid_polymer: 1 - fallthrough}
        if 'BiogenicMaterial' in self.chemicals:
            split['BiogenicMaterial'] = 1 - fallthrough
        return self.chemicals.kwsplit(split)

    def _run(self):
        feed = self.ins[0]
        solvent = self.solvent
        retentate, permeate = self.outs
        moisture_content = self.moisture_content
        split = self.get_solids_split()
        feed_solvent = feed.imass[solvent]
        feed.split_to(retentate, permeate, split)
        retentate_solvent = (
            moisture_content / (1 - moisture_content)
            * (retentate.F_mass - retentate.imass[solvent])
        )
        retentate.imass[solvent] = min(retentate_solvent, 0.99 * feed_solvent)
        solutes = self.solutes
        permeate.imass[solvent] = feed.imass[solvent] - retentate.imass[solvent]
        x = feed.imass[solutes] / feed.imass[solvent]
        retentate.imass[solutes] = retentate.imass[solvent] * x
        permeate.imass[solutes] = feed.imass[solutes] - retentate.imass[solutes]

    def _design(self):
        self.design_results['Solids loading'] = self.outs[0].F_mass

@cost('Solids loading', 'Centrifilter', CE=CEPCI2022, units='kg/hr',
      S=91, cost=65e3, BM=2.03, n=0.5, kW=0.55)
class Centrifilter(Filter):
    pass

@cost('Solids loading', 'Pressure filter',
      cost=3294700, CE=551, S=31815, n=0.8, BM=1.7)
@cost('Solids loading', 'Pressing air compressor receiver tank',
      cost=8e3, CE=551, S=31815, n=0.7, BM=3.1)
@cost('Solids loading', 'Dry air compressor receiver tank',
      cost=17e3, CE=551, S=31815, n=0.7, BM=3.1)
@cost('Solids loading', 'Pressing air pressure filter',
      cost=75200, CE=521.9, S=31815, n=0.6, kW=112, BM=1.6)
@cost('Solids loading', 'Dry air pressure filter (2)',
      cost=405000, CE=521.9, S=31815, n=0.6, kW=1044, BM=1.6)
class PrecipitatorFilter(Filter):
   pass

def test_centrifilter():
    Solvent = bst.Chemical(
        'Solvent', search_ID='Dodecane'
    )

    DissolvedPolymer = bst.Chemical(
        'DissolvedPolymer',
        search_ID='Dodecane', phase='l'
    )

    SolidPolymer = bst.Chemical(
        'SolidPolymer',
        search_ID='Dodecane', phase='s'
    )
    chemicals = [DissolvedPolymer, SolidPolymer, Solvent]
    for i, c in enumerate(chemicals):
        c._CAS = i
    bst.settings.set_thermo(chemicals)
    feed = bst.Stream('feed',
                      SolidPolymer=0.03, DissolvedPolymer=0.07, Solvent=0.9,
                      total_flow=1000, units='kg/hr'
                      )

    CF = Centrifilter('CF', ins=feed, outs=['permeate', 'retentate'])
    CF.simulate()
    CF.show('cwt')
    CF.diagram()

