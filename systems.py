# -*- coding: utf-8 -*-
"""
"""
import biosteam as bst
import thermo as tm
import thermosteam as tmo
import numpy as np
from biorefineries import cellulosic
from warnings import warn
import biorefineries as br
from .property_package import create_property_package
from .units import (
    PseudoContinuousPlasticDissolutionTank, BatchPlasticDissolution, Shredding,
    Precipitator, JacketedSurgeTank, PrecipitatorFilter, CandleFilter, VacuumDryer, ScrewPressDegasser,
    SeedTrain, CoFermentation, HandSorting, EddyCurrent, Vecoplan, Magnet, Crumble,
)
s = bst.stream_kwargs

__all__ = (
    # 'create_psuedo_continuous_system',
    'create_single_layer_batch_separation_system',
    'create_multilayer_batch_separation_system',
)

# TODO: 
# [-] Make sure feed and product storage is well accounted for

@bst.SystemFactory(
    ID='preprocessing_sys',
    ins=[s('MSW', Plastic=2.5443, BiogenicMaterial=11.6401, Solubles=0.01,
           total_flow=2e5, units="kg/hr")],
    outs=['MSW_residue'],
)
def create_preprocessing_system(ins, outs):
    feed, = ins
    MSW_residue, = outs
    
    HE_Loss = bst.Splitter('hand_sorting_loss', ins=feed, outs=('','sorting_loss'), split=0.98)
    S1 = HandSorting('hand_sorting', ins=HE_Loss-0, outs=('sorted_material','EC_feed'), split={'HandSorted':1})
    S2 = EddyCurrent('eddy_current', ins=S1-1, outs=('nonferrous_metal','vecoplan_feed'), split={'Metal':0.82})
    VS_Loss = bst.Splitter('vecoplan_loss', ins=S2-1, outs=('','vecoplan_loss'), split=0.993)
    S3 = Vecoplan('vecoplan', ins=VS_Loss-0, outs=('magnet_feed'))
    MC_Loss = bst.Splitter('MC_Loss', ins=S3-0, outs=('','crumble_loss'), split=0.98)
    S4 = Magnet('magnet', ins=MC_Loss-0, outs=('ferrous_metal','crumble_feed'), split={'Metal':1})
    S5 = Crumble('crumble', ins=S4-1, outs=('overs', 'unders', MSW_residue))
    
    # TODO: Update prices
    S1.outs[0].price = 1e-16
    S2.outs[0].price = 1e-16
    S4.outs[0].price = 1e-16
    S5.outs[0].price = 1e-16
    S5.outs[1].price = 1e-16

@bst.SystemFactory(
    ID='batch_sys',
    ins=[s("plastic", Plastic=300, bst="kg/hr")],
    outs=[s("leftover_plastic", price=2.4023)],
    fthermo=create_property_package,
)
def create_multilayer_batch_separation_system(
        ins, outs,
        dissolution_steps, 
        precipitation_steps,
        shred=True,
        facilities=False,
        core_facilities=True,
        precipitation_configuration=None,
        turbogenerator=True,
    ):
    # TODO: Work on solvent and plastic prices
    solvents = [bst.Stream(i.solvent, price=2.17) for i in dissolution_steps]
    products = [bst.Stream(i.plastic + '_resin', price=2.4023) for i in precipitation_steps]
    ins.extend(solvents)
    outs.extend(products)
    first = True
    plastic = ins[0]
    for ds, ps, solvent, product in zip(dissolution_steps, precipitation_steps, solvents, products):
        intermediate = bst.Stream()
        sys = create_single_layer_batch_separation_system(
            ins=[plastic, solvent], 
            outs=[intermediate, product], 
            shred=first, mockup=True,
            dissolution_step=ds,
            precipitation_step=ps,
            burn_leftover_plastic=False,
            facilities=first and facilities,
            core_facilities=first and core_facilities,
            turbogenerator=first and turbogenerator,
            precipitation_configuration=precipitation_configuration,
        )
        plastic = intermediate
        first = False
    sys.outs[0] = outs[0]

@bst.SystemFactory(
    ID='STRAPMSW_sys',
    ins=[s('MRF_residues', Plastic=2.5443, BiogenicMaterial=11.6401, Solubles=0.01,
           total_flow=2e5, units="kg/hr")],
    outs=[s('resin'), s('ethanol')],
    fixed_ins_size=False,
    fixed_outs_size=False,
)
def create_STRAPMSW_system(
        ins, outs, 
        dissolution_step, precipitation_step, precipitation_configuration,
        preprocessing=True
    ):
    MSW, = ins
    PEPP, ethanol = outs
    PEPP.register_alias('product')
    PEPP.register_alias(dissolution_step.plastic)
    if preprocessing:
        sys = create_preprocessing_system(
            ins=MSW,
            area=100,
        )
        # Add preprocessing non-residues
        residue = MSW.F_mass
        total = residue / 0.78
        MSW.imass['HandSorted', 'Metal', 'Overs', 'Unders'] = total * np.array([10.10, 2.382, 2.823, 6.382]) / 100
        MSW = sys.outs[-1] # Actual MSW-Residue
    n = 100 if preprocessing else 0
    sys = create_single_layer_batch_separation_system(
        ins=MSW,
        outs=['biogenics_and_plastic', PEPP],
        dissolution_step=dissolution_step, 
        precipitation_step=precipitation_step,
        precipitation_configuration=precipitation_configuration,
        mockup=True,
        area=100 + n,
        core_facilities=False,
    )
    biogenics_and_plastic, PEPP = sys.outs
    PRxn = bst.PRxn
    Rxn = bst.Rxn
    saccharification_reactions = PRxn([
        Rxn('Glucan + H2O -> Glucose',       'Glucan', 0.9000),
        Rxn('Glucan -> GlucoseOligomer',     'Glucan', 0.0400),
        Rxn('Xylan + H2O -> Xylose',         'Xylan',  0.9000),
        Rxn('Xylan + H2O -> XyloseOligomer', 'Xylan',  0.0240),
    ])
    
    splitter = bst.Splitter(100 + n, ins=biogenics_and_plastic, split=1)
    splitter.register_alias('incineration_splitter')
    splitter.outs[1].ID = 'to_boiler' # This tell's BioSTEAM to send it to CH&P
    cellulosic_fermentation_sys, cfdct = cellulosic.create_cellulosic_fermentation_system(
        ins=(splitter-0,),
        outs=['vent', 'beer', 'lignin'],
        mockup=True,
        udct=True,
        kind='Saccharification and Co-Fermentation',
        saccharification_reactions=saccharification_reactions,
        SeedTrain=SeedTrain,
        CoFermentation=CoFermentation,
        insoluble_solids_loading=0.23,
        add_nutrients=False,
        solids_loading=0.23, # 30 wt/vol % solids content in saccharification
        include_scrubber=True,
        area=200 + n,
    )
    seed_splitter = cfdct['S304']
    pressure_filter = cfdct['S303']
    pressure_filter.outs[1] = bst.Stream()
    pressure_filter.strict_moisture_content = False
    mixer = bst.Mixer(200 + n, ins=[pressure_filter.outs[1], 'dilution_water'], outs=seed_splitter.ins[0])
    cellulase_mixer = cfdct['M301']
    hydrolysis = cfdct['R301']
    seed_train = cfdct['R302']
    fermentor = cfdct['R303'] # Cofermentation
    cellulosic.add_urea_nutrient(fermentor, seed_train, area=200 + n)
    seed_train.register_alias('seed_train')
    fermentor.register_alias('fermentor')
    hydrolysis.register_alias('hydrolysis')
    seed_splitter.register_alias('seed_splitter')
    cellulase_mixer.register_alias('cellulase_mixer')
    
    def get_titer():
        s = fermentor.outs[1]
        return (s.imass['Ethanol']) / s.F_vol # kg / m3
    fermentor.get_titer = get_titer
    
    def get_additional_dilution_water():
        target = fermentor.titer
        current = get_titer()
        effluent = fermentor.outs[1]
        T = effluent.T
        P = effluent.P
        rho = fermentor.chemicals.Water.rho('l', T, P)
        return (1./target - 1./current) * effluent.imass['Ethanol'] * rho
    
    @fermentor.add_specification(run=True, impacted_units=[mixer])
    def adjust_titer():
        fermentor.run()
        dilution_stream = mixer.ins[1]
        new_dilution_water = dilution_stream.imass['Water'] + get_additional_dilution_water()
        if new_dilution_water < 0 :
            # TODO: Add evaporation if needed
            # raise RuntimeError('need evaporation')
            new_dilution_water = 0
            warn('need evaporation')
        dilution_stream.imass['Water'] = new_dilution_water
        productivity = fermentor.productivity
        fermentor.tau = fermentor.titer / productivity
    
    vent, beer, lignin = cellulosic_fermentation_sys.outs
    cellulosic_beer_distillation_sys = br.ethanol.create_beer_distillation_system(
        ins=beer,
        outs=[''],
        mockup=True,
        area=200 + n,
    )
    stripper_process_water = bst.Stream('')
    distilled_beer, stillage = cellulosic_beer_distillation_sys.outs
    ethanol_purification_sys, ep_dct = br.ethanol.create_ethanol_purification_system_after_beer_column(
        ins=distilled_beer,
        outs=[ethanol, stripper_process_water],
        mockup=True,
        udct=True,
        area=200 + n,
    )
    bst.SolidsCentrifuge(200 + n, stillage, 
                         split={'Cellmass': 0.99}, solids=('Cellmass',))
    bst.create_all_facilities(
        feedstock=biogenics_and_plastic,
        WWT_kwargs=dict(
            kind='high-rate',
            area=300 + n,
        ),
        HXN=False,
        # HXN_kwargs=dict(
        #     ID=500 + n, 
        #     Qmin=1e3,
        #     force_ideal_thermo=True, 
        #     cache_network=True,
        #     acceptable_energy_balance_error=0.01,
        # ),
        area=400 + n,
    )
    BT = bst.F(bst.BoilerTurbogenerator)
    BT.register_alias('BT')
    for i in bst.F(bst.SolidsSeparator):
        i.strict_moisture_content = False

@bst.SystemFactory(
    ID='batch_sys',
    ins=[s("plastic", Plastic=300, Solubles=0.0003, units="kg/hr")],
    outs=[s("leftover_plastic")],
    fixed_ins_size=False,
    fixed_outs_size=False,
    fthermo=create_property_package,
)
def create_single_layer_batch_separation_system(
        ins, outs,
        dissolution_step=None, 
        precipitation_step=None,
        shred=True,
        vent_sink=None,
        facilities=False,
        core_facilities=True,
        solvent_loss=1e-4,  
        burn_leftover_plastic=True,
        precipitation_configuration=None,
        turbogenerator=True,
    ):
    """This configuration follows recommendations from the Michigan 
    Technological University (MTU). Because centrifugation and ultrafiltrations
    is performed between every plastic dissolution, it can handle extreemly 
    small plastic sizes."""
    if precipitation_configuration is None:
        precipitation_configuration = 'solvent mixing'
    # TODO: Make prices part of the dissolution and precipitation step variables
    if len(ins) == 1:
        plastic, = ins
        solvent = bst.Stream(dissolution_step.solvent, price=2.17)
        solvent.register_alias('solvent', override=False, safe=False)
        ins.append(solvent)
    else:
        plastic, solvent = ins
    plastic.register_alias('feedstock', override=False, safe=False)
    if len(outs) == 1:
        undissolved_plastic, = outs
        product = bst.Stream(precipitation_step.plastic + '_resin', price=2.4023)
        outs.append(product)
    else:
        undissolved_plastic, product = outs
    product.register_alias('product', override=False, safe=False)
    plastic_storage = bst.StorageTank(
        ins=plastic, tau=7 * 24, vessel_material="Carbon steel"
    )
    plastic = plastic_storage-0
    if shred:
        shredder = Shredding(ins=plastic_storage-0)
        plastic = shredder-0
    storage = bst.StorageTank(
        ins=solvent, tau=7 * 24, vessel_material="Carbon steel",
        vessel_type='Floating roof', 
    )
    pump = bst.Pump(ins=storage-0, P=2 * 101325)
    recycle = bst.Stream()
    solvent_mixer = bst.Mixer(ins=[pump-0, recycle])
    surge_tank = bst.StorageTank(
        ins=solvent_mixer-0, vessel_type='Floating roof', 
        vessel_material="Stainless steel", tau=0.15,
    )
    surge_tank.line = 'Surge tank'
    splitter = bst.Splitter(ins=surge_tank-0, outs=('solvent_purge', ''), split=solvent_loss)
    splitter.register_alias('solvent_loss')
    splitter.register_alias(dissolution_step.solvent + '_loss')
    adsorption_column = bst.AdsorptionColumn(
        ins=[splitter-1, '', 'activated_carbon'],
        cycle_time=168, # h
        superficial_velocity=14.4, # m / h
        isotherm_model='Freundlich',
        isotherm_args = (48.4, 2.58), # K, n
        # k = 0.00931, # [1 / hr]
        void_fraction=1 - 0.38 / 0.8, 
        C_final_scaled=0.05, # Breakthrough condition
        adsorbate='Solubles',
        LUB_forced=0.6096, # For faster simulation
    )
    adsorption_column.register_alias('adsorption_column')
    adsorption_column.outs[2].ID = 'spent_activated_carbon'
    solvent_pump = bst.Pump(ins=adsorption_column-0, P=10 * 101325)
    solvent2dissolution = solvent_pump-0
    solvent_dissolution_heater = bst.HXutility(
        ins=solvent2dissolution, heat_only=True, T=dissolution_step.T
    )
    
    @solvent_dissolution_heater.add_specification(run=True)
    def adjust_temperature():
        solvent_dissolution_heater.T = dissolution_step.T
    
    dissolution = BatchPlasticDissolution(
        ins=(plastic, solvent_dissolution_heater-0), 
        vessel_material="Stainless steel 304", 
        dissolution_step=dissolution_step,
        minimum_solvent_to_feed=4,
        adiabatic=True
    )
    dissolution.register_alias('dissolution_vessel')
    # dissolution.breakpoint = False
    
    @dissolution.add_specification(
        args=(solvent_mixer, dissolution, dissolution_step, solvent_dissolution_heater),
        impacted_units=[storage], 
    )
    def adjust_fresh_solvent_flow(mixer, dissolution, step, hx):
        feed, solvent = dissolution.ins
        solvent_key = step.solvent
        # if dissolution.breakpoint: breakpoint()
        solvent_recycled = solvent.imol[solvent_key]
        dissolution.run()
        solvent_out = dissolution.outs[0].imol[solvent_key]
        solvent_loss = solvent_out - solvent_recycled
        if solvent_loss < 0:
            recycle.F_mass *= solvent.imol[solvent_key] / solvent_recycled
            solvent_loss = 0
        storage.ins[0].imol[solvent_key] = solvent_out * splitter.isplit[solvent_key]
        storage.run_until(dissolution, inclusive=True)
        mixed = solvent + feed
        H_actual = mixed.H
        mixed.T = dissolution_step.T
        H_target = mixed.H
        solvent.H += H_target - H_actual
        hx.T = solvent.T
        mixer.run_until(dissolution, inclusive=True)
        print('MAKING SURE')
    
    effluent, = dissolution.outs
    surge_tank = JacketedSurgeTank(
        ins=effluent, vessel_type='Floating roof', 
        vessel_material="Carbon steel", tau=dissolution.tau,
    )
    # TODO: Add option for design specification
    # @surge_tank.add_design_specification(design=True)
    # def adjust_surge_time():
    #     surge_tank.tau = dissolution.design_results['Batch time']
        
    surge_tank.line = 'Jacketed surge tank'
    pump = bst.Pump(ins=surge_tank-0, P=2 * 101325)
    centrifuge = bst.Centrifilter(
        ins=pump-0, 
        moisture_content=dissolution_step.solvent_content,
        solvent=dissolution_step.solvent,
        solid_polymer='Plastic',
        solutes=('Solutes', dissolution_step.dissolved_plastic),
    )
    @centrifuge.add_specification(run=True, args=(centrifuge,))
    def adjust_moisture(centrifuge):
        centrifuge.moisture_content = dissolution_step.solvent_content
       
    centrifuge.register_alias('filter_centrifuge')
    solids = [i.ID for i in centrifuge.chemicals if i.locked_state == 's']
    # centrifuge.isplit[solids] = 0.999
    dummy = bst.Stream(None)
    dummy.imol[dissolution_step.solvent] = 1
    P = 0.1 * 101325
    bp = dummy.bubble_point_at_P(P)
    dryer0 = VacuumDryer(
        ins=centrifuge-0,
        split=0,
        T=50 + bp.T, # TODO: expose assumption as a parameter
        moisture_content=0,
        moisture_ID=dissolution_step.solvent,
        P=P
    )
    condensate = dryer0-1
        
    # @flash.add_specification(args=(flash,))
    # def use_CEOS(flash):
    #     inlet = flash.ins[0]
    #     chemicals = inlet.available_chemicals
    #     IDs = [i.ID for i in chemicals]
    #     flows = inlet.imol[IDs]
    #     F_mol = flows.sum()
    #     zs = flows / F_mol
    #     zs = zs.tolist()
    #     T = flash.T
    #     P = flash.P
    #     # Use equations of state for the gas and liquid phases.
    #     flashpkg = tmo.equilibrium.FlashPackage(
    #         chemicals=chemicals,
    #         G=tm.CEOSGas, L=tm.CEOSLiquid, S=tm.GibbsExcessSolid,
    #         GE=tm.UNIFAC, GEkw=dict(version=1), Gkw=dict(eos_class=tm.PRMIX),
    #     )
    #     flasher = flashpkg.flasher(N_liquid=1)
    #     PT = flasher.flash(zs=zs, T=T, P=P)
    #     V = PT.VF
    #     gas_zs = PT.gas.zs
    #     vapor, liquid = flash.outs
    #     for i in (vapor, liquid): 
    #         i.T = T
    #         i.P = P
    #     vapor.imol[IDs] = vapor_flows = F_mol * V * np.asarray(gas_zs)
    #     liquid.imol[IDs] = flows - vapor_flows
    #     vapor.imol['N2', 'O2'] += liquid.imol['N2', 'O2']
    #     liquid.imol['N2', 'O2'] = 0.
    
    microfilter = CandleFilter(
        ins=centrifuge-1, 
        split=0,
        moisture_ID=dissolution_step.solvent,
        moisture_content=0,
    )
    mixer = bst.Mixer(ins=(microfilter-0, dryer0-0), outs=undissolved_plastic)
    precipitation_tank = Precipitator(
        ins=microfilter-1, 
        precipitation_step=precipitation_step,
    )
    microfilter.isplit[solids] = 1.0
    
    # TODO: Cost this as a filter press instead
    centrifuge = PrecipitatorFilter(
        ins=precipitation_tank-0,
        mesh_size=10, # um
        moisture_content=precipitation_step.precipitate_solvent_content,
        solvent=precipitation_step.solvent,
        solid_polymer='Plastic',
    )
    @centrifuge.add_specification(run=True, args=(centrifuge,))
    def adjust_moisture(centrifuge):
        centrifuge.moisture_content = precipitation_step.precipitate_solvent_content
        
    melt_vacuum = bst.Flash(
        ins=centrifuge-0, P=0.2 * 101325, 
        V=0.999, # TODO: Make this multiple stages and adjust model to vaporize set fraction at a T and P
        thermo=centrifuge.thermo.ideal(), # Avoid VLLE
        has_vapor_condenser=True,
    )
    melt_vacuum.line = 'Melt degasser'
    degasser = ScrewPressDegasser(
        ins=melt_vacuum-1,
        outs=[product, ''],
        split=0,
        T=50 + bp.T, # TODO: expose assumption as a parameter
        moisture_content=1e-9,
        moisture_ID=precipitation_step.solvent,
        P=101325 * 0.05
    )
    condenser = bst.HXutility(ins=degasser-1, V=0, rigorous=True)
    liq_mixer = bst.Mixer(
        ins=[centrifuge-1, condensate, melt_vacuum-0, condenser-0],
        outs=[recycle],
    )
    CHP = bst.BoilerTurbogenerator if turbogenerator else bst.Boiler
    if facilities:
        bst.create_all_facilities(
            HXN=False, CIP=False, WWT=False, PWC=False,
            CHP_kwargs=dict(autopopulate=burn_leftover_plastic,
                            gas_feed=False,
                            cls=CHP),
        )
    elif core_facilities:
        bst.create_all_facilities(
            HXN=False, CIP=False, WWT=False, PWC=False, CHP=False, 
        )
    match precipitation_configuration:
        case 'solvent mixing':
            splitter = bst.Splitter(
                ins=solvent2dissolution,
                split=0.5
            )
            hx = bst.HXutilities(
                ins=splitter-0,
                T=273.15 + 10,
                rigorous=False,
            )
            mixer = bst.Mixer(
                ins=[hx-0, microfilter-1],
            )
            solvent_dissolution_heater.ins[0] = splitter-1
            precipitation_tank.ins[0] = mixer-0
            
            @mixer.add_specification(args=(mixer, hx))
            def adjust_split_to_meet_precipitation_temperature(mixer, hx):
                # sum(Hin) = Hout
                # F_0 * h_in_0 + F_1 * h_in_1 = (F_0 * h_out_0) + (F_1 * h_out_1)
                # F_0 * (h_in_0 - h_out_0)  = F_1 * h_out_1 - F_1 * h_in_1
                # F_0 = (F_1 * h_out_1 - F_1 * h_in_1) / (h_in_0 - h_out_0)
                s0, s1 = mixer.ins
                s0.T = hx.T
                s0.empty()
                s0.imol[precipitation_step.solvent] = 1
                s0_out = s0.copy()
                s1_out = s1.copy()
                s0_out.T = s1_out.T = precipitation_step.T
                h_in_0 = s0.h
                h_in_1 = s1.h
                h_out_0 = s0_out.h
                h_out_1 = s1_out.h
                F_1 = s1.F_mol
                F_0 = (F_1 * h_out_1 - F_1 * h_in_1) / (h_in_0 - h_out_0)
                s0.F_mol = F_0
                mixer.simulate()
                s_out = mixer.outs[0]
                assert (
                    abs(s_out.F_mol - (F_1 + F_0)) < 1e-3 and 
                    abs(s_out.T - precipitation_step.T) < 1e-3
                )
                F_mol_recycle = splitter.ins[0].F_mol
                if F_mol_recycle < F_0:
                    splitter.split[:] = 1
                else:
                    splitter.split[:] = F_0 / F_mol_recycle
                
        case 'integrated heat transfer':
            hot_solvent = microfilter-1
            cooling = hot_solvent.sink
            cold_solvent = solvent2dissolution
            heating = cold_solvent.sink
            hx_integration = bst.HXprocess(ins=(hot_solvent, cold_solvent), outs=('cooled_solvent', 'heated_solvent'), phase0='l', phase1='l')
            cooling.ins[0], heating.ins[0] = hx_integration.outs
        case _:
            raise ValueError("precipitation_configuration must be either 'integrated heat transfer' or 'solvent mixing'")
    
    
    # for i in recycles: bst.mark_disjunction(i) # Flow is completely determined by dissolution tank

    
# @bst.SystemFactory(
#     ID='psuedo_continuous_sys',
#     ins=[s("plastic", Plastic=300, units="kg/hr")],
#     outs=[s("product", price=2.4023)],
#     fthermo=create_property_package,
# )
# def create_psuedo_continuous_system(
#         ins, outs,
#         dissolution_steps=None, 
#         precipitation_steps=None,
#     ):
#     """This configuration may only work with big plastic packing. Otherwise,
#     best to use the batch system with ultrafiltration."""
#     #: TODO: Cost filter.
#     plastic, = ins
#     undissolved_plastic, = outs
#     if dissolution_steps is None:
#         dissolution_steps = [EVOH_dissolution(), PE_dissolution()]
#     if precipitation_steps is None:
#         precipitation_steps = [EVOH_precipitation(), PE_precipitation()]
#     # TODO: Make prices part of the dissolution and precipitation step variables
#     solvents = [bst.Stream(i.solvent, price=2.17) for i in dissolution_steps]
#     products = [bst.Stream(i.plastic, price=2.4023) for i in precipitation_steps]
#     ins.extend(solvents)
#     outs.extend(products)
    
#     plastic_storage = bst.StorageTank(
#         ins=plastic, tau=7 * 24, vessel_material="Carbon steel"
#     )
#     conveyor = bst.ConveyingBelt(ins=plastic_storage-0)
#     mixers = []
#     recycles = []
#     solvents2dissolution = []
#     dissolution_times = [i.tau for i in dissolution_steps]
#     drying_time = 0.05
#     surge_time = sum(dissolution_times) + len(dissolution_times) * drying_time # This is also the cycle time
#     for solvent in solvents: 
#         storage = bst.StorageTank(
#             ins=solvent, tau=7 * 24, vessel_material="Carbon steel",
#             vessel_type='Floating roof', 
#         )
#         pump = bst.Pump(ins=storage-0, P=2 * 101325)
#         recycle = bst.Stream()
#         recycles.append(recycle)
#         mixer = bst.Mixer(ins=[pump-0, recycle])
#         mixers.append(mixer)
#         surge_tank = bst.StorageTank(
#             ins=mixer-0, vessel_type='Floating roof', 
#             vessel_material="Stainless steel", tau=surge_time,
#         )
#         surge_tank.line = 'Surge tank'
#         pump = bst.Pump(ins=surge_tank-0, P=10 * 101325)
#         solvents2dissolution.append(pump-0)
        
#     dissolution_tanks = PseudoContinuousPlasticDissolutionTank(
#         ins=(conveyor-0, *solvents2dissolution, 'air'), 
#         outs=(undissolved_plastic, *len(solvents2dissolution) * [''], 'vent'),
#         vessel_material="Stainless steel 304", 
#         dissolution_steps=dissolution_steps, drying_time=drying_time
#     )
#     for mixer, solvent, step in zip(mixers, solvents2dissolution, dissolution_steps):
#         @mixer.add_specification(run=True, args=(mixer, solvent, step))
#         def adjust_fresh_solvent_flow(mixer, solvent, step):
#             fresh, recycle = mixer.ins
#             fresh.imol[step.solvent] = solvent.imol[step.solvent] - recycle.imol[step.solvent]
    
#     plastic, *dissolved_plastics, vent = dissolution_tanks.outs
#     vents = []
#     for dissolved_plastic, precipitation_step, recycle, product in zip(dissolved_plastics, precipitation_steps, recycles, products):
#         surge_tank = bst.StorageTank(
#             ins=dissolved_plastic, vessel_type='Floating roof', 
#             vessel_material="Carbon steel", tau=surge_time,
#         )
#         surge_tank.line = 'Surge tank'
#         pump = bst.Pump(ins=surge_tank-0, P=2 * 101325)
#         cooler = bst.HXutility(ins=pump-0, T=precipitation_step.T)
#         @cooler.add_specification(args=(cooler, precipitation_step,))
#         def precipitate(cooler, precipitation_step):
#             cooler.T = precipitation_step.T
#             cooler.run()
#             precipitation_step.precipitate(cooler.outs[0])
        
#         centrifuge = PrecipitatorFilter(
#             ins=cooler-0,
#             moisture_content=precipitation_step.solvent_content,
#             solvent=precipitation_step.solvent,
#             solid_polymer='Plastic',
#         )
#         screw_press = bst.ScrewPress(
#             ins=centrifuge-0,
#             split=1,
#             moisture_content=precipitation_step.solvent_contents[1],
#             moisture_ID=precipitation_step.solvent
#         )
#         pellet_mill = bst.PelletMill(ins=screw_press-0)
#         dummy = bst.Stream(None)
#         dummy.imol[precipitation_step.solvent] = 1
#         bp = dummy.bubble_point_at_P()
#         dryer = bst.DrumDryer(
#             ins=pellet_mill-0,
#             outs=(product, '', ''),
#             split=0,
#             R=1.4,
#             H=20,
#             length_to_diameter=25,
#             T=10 + bp.T,
#             moisture_content=precipitation_step.solvent_contents[2],
#             moisture_ID=precipitation_step.solvent,
#             gas_composition=[('N2', 1.0)],
#             utility_agent="Steam",
#         )
        
#         flash = bst.Flash(
#             ins=dryer-1, 
#             T=precipitation_step.T_condensation,
#             P=101325
#         )
        
#         @flash.add_specification(args=(flash,))
#         def use_CEOS(flash):
#             inlet = flash.ins[0]
#             chemicals = inlet.available_chemicals
#             IDs = [i.ID for i in chemicals]
#             flows = inlet.imol[IDs]
#             F_mol = flows.sum()
#             zs = flows / F_mol
#             zs = zs.tolist()
#             T = flash.T
#             P = flash.P
#             # Use activity coefficients for the liquid phase 
#             # and equations of state for the gas phase.
#             flashpkg = tmo.equilibrium.FlashPackage(
#                 chemicals=chemicals,
#                 G=tm.CEOSGas, L=tm.GibbsExcessLiquid, S=tm.GibbsExcessSolid,
#                 GE=tm.UNIFAC, GEkw=dict(version=1), Gkw=dict(eos_class=tm.PRMIX),
#             )
#             flasher = flashpkg.flasher(N_liquid=1)
#             PT = flasher.flash(zs=zs, T=T, P=P)
#             V = PT.VF
#             gas_zs = PT.gas.zs
#             vapor, liquid = flash.outs
#             for i in (vapor, liquid): 
#                 i.T = T
#                 i.P = P
#             vapor.imol[IDs] = vapor_flows = F_mol * V * np.asarray(gas_zs)
#             liquid.imol[IDs] = flows - vapor_flows
            
#         vents.append(flash-0)
#         mixer = bst.Mixer(
#             ins=[centrifuge-1, screw_press-1, flash-1],
#             outs=recycle
#         )
#     mixer = bst.Mixer(ins=vents)
#     bst.ThermalOxidizer(ins=mixer-0)
#     # for i in recycles: bst.mark_disjunction(i) # Flow is completely determined by dissolution tank

    

    
