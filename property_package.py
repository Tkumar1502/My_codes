# -*- coding: utf-8 -*-
"""
"""
import biosteam as bst
from chemicals.elements import (
    molecular_weight as compute_molecular_weight,
    get_atoms,
)

MW = lambda formula: compute_molecular_weight(get_atoms(formula))

__all__ = (
    'STRAP_chemicals_outline',
    'create_property_package',
)

STRAP_chemicals_outline = bst.ChemicalsOutline([
    "Toluene", 
    "DMSO", 
    bst.ChemicalDraft("THF", search_ID='109-99-9'),
    bst.ChemicalDraft(
        'PE', 
        aliases=set(['Polyethylene']), 
        formula='C2H4',
        search_db=False,
        phase='s',
        rho=0.5 * (880 + 960), # kg / m3
        Cp=0.5 * (1.330 + 2.400), # J / g
        Tm=0.5 * (115 + 135) + 273.15, # K
        LHV=261.0 * MW('C2H4'), # https://nvlpubs.nist.gov/nistpubs/jres/78A/jresv78An5p611_A1b.pdf
        default=True,
    ),
    bst.ChemicalDraft(
        'PEoligomer', 
        search_ID='1-Hexene',
    ),
    bst.ChemicalDraft('ActivatedCarbon', search_db=False, rho=1540, phase='s', default=True, MW=1),
    bst.ChemicalDraft(
        'PC', 
        aliases=set(['Polycarbonate']), 
        formula='C2H4',
        search_db=False,
        phase='s',
        rho=0.5 * (880 + 960), # kg / m3
        Cp=0.5 * (1.330 + 2.400), # J / g
        Tm=0.5 * (115 + 135) + 273.15, # K
        default=True,
    ),
    bst.ChemicalDraft( # TODO: is this a good proxy?
        'PColigomer', 
        search_ID='1-Heptene',
    ),
    bst.ChemicalDraft(
        "PET", 
        aliases=set(['Polyethylene terephthalate']),
        formula='C10H8O4',
        search_db=False,
        phase='s',
        rho=1380, # kg / m3,
        Tm=523,
        Tb=623,
        Cp=1,
        default=True,
        LHV= 21285 * 192.16812,
    ),
    bst.ChemicalDraft( # Modeled as PET
        "BulkPlastic", 
        formula='C10H8O4',
        search_db=False,
        phase='s',
        rho=1380, # kg / m3,
        Tm=523,
        Tb=623,
        Cp=1,
        default=True,
        LHV= 21285 * 192.16812,
    ),
    
#NdFeB magnets
     bst.ChemicalDraft(
        "NdFeB",
        aliases  = set(['Neodymium Magnet']),
        search_db = False,
        phase = 's',
        rho = 0.5*(7300+7700),
        Tm = 273 + 1024,
        Cp = 0.5*(0.12+0.15),
        default = True,
        
        
        ),
 
 
#HDPE
     bst.ChemicalDraft(
         'HDPE', 
         aliases=set(['High Density Polyethylene']), 
         formula='C2H4',
         search_db=False,
         phase='s',
         rho=0.5 * (930 + 970), # kg / m3
         Cp=0.5 * (1.330 + 2.400), # J / g
         Tm=0.5 * (120 + 135) + 273.15, # K
         LHV=261.0 * MW('C2H4'), # https://nvlpubs.nist.gov/nistpubs/jres/78A/jresv78An5p611_A1b.pdf
         default=True,
     ),


    bst.ChemicalDraft(
    'HDPEoligomer', 
    search_ID='1-Hexene',
    CAS = 'HDPEoligomer',
    ),
    
    bst.ChemicalDraft(
        'EVOH',
        aliases=set(['Ethylene vinyl alcohol']),
        formula='C2H4OC2H4',
        search_db=False,
        phase='s',
        rho=0.5 * (1120 + 1140),
        Cp=2.4,
        default=True,
        LHV= 21285 * MW('C2H4OC2H4'), # TODO: Adjust lower heating value and add reference.
    ),
    bst.ChemicalDraft(
        'EVOHoligomer', 
        search_ID='3-buten-2-ol', # Approximate monomer
    ),
    "Water",
    bst.ChemicalDraft("N2", phase='g'), 
    bst.ChemicalDraft("O2", phase='g'), 
    "CH4", 
    "CO2",
    "SO2",
    "Xylene", 
    bst.ChemicalDraft("Ash", search_db=False, rho=1540, phase='s', default=True, MW=1),
    bst.ChemicalDraft("Minerals", search_db=False, rho=1540, phase='l', default=True, MW=1),
    bst.ChemicalDraft("Solubles", search_db=False, rho=1540, phase='l', default=True, MW=1),
])

def create_property_package():
    chemicals = STRAP_chemicals_outline.to_chemicals()
    chemicals.compile()
    return chemicals

def create_property_package_MSW():
    from biorefineries import cane
    chemicals = bst.Chemicals([
        "Xylene", 
        # bst.Chemical(
        #     'PE', 
        #     aliases=set(['Polyethylene']), 
        #     formula='C2H4',
        #     search_db=False,
        #     phase='s',
        #     rho=0.5 * (880 + 960), # kg / m3
        #     Cp=0.5 * (1.330 + 2.400), # J / g
        #     Tm=0.5 * (115 + 135) + 273.15, # K
        #     default=True,
        # ),
        # bst.Chemical(
        #     'PEoligomer', 
        #     search_ID='1-Hexene',
        # ),
        # bst.Chemical(
        #     'PP', 
        #     aliases=set(['Polypropylene']), 
        #     formula='C3H6',
        #     search_db=False,
        #     phase='s',
        #     rho=0.5 * (880 + 960), # kg / m3
        #     Cp=0.5 * (1.330 + 2.400), # J / g
        #     Tm=0.5 * (115 + 135) + 273.15, # K
        #     default=True,
        # ),
        # bst.Chemical(
        #     'PPoligomer', 
        #     search_ID='1-Hexene',
        # ),
        bst.Chemical(
            'PEPP', 
            aliases=set(['Polyolefins']), 
            formula='C3H6',
            search_db=False,
            phase='s',
            rho=0.5 * (880 + 960), # kg / m3
            Cp=0.5 * (1.330 + 2.400), # J / g
            Tm=0.5 * (115 + 135) + 273.15, # K
            default=True,
            LHV=0.5 * (45600 + 261.0) * MW('C3H6'),
        ),
        bst.Chemical(
            'PEPPoligomer', 
            search_ID='1-Hexene',
        ),
        bst.Chemical(
            "BulkPlastic", # Modeled as PET
            formula='C10H8O4',
            search_db=False,
            phase='s',
            rho=1380, # kg / m3,
            Tm=523,
            Tb=623,
            Cp=1,
            default=True,
        ),
        bst.Chemical('ActivatedCarbon', search_db=False, rho=1540, phase='s', default=True, MW=1),
        bst.Chemical("Minerals", search_db=False, rho=1540, phase='l', default=True, MW=1),
        bst.Chemical("Solubles", search_db=False, rho=1540, phase='l', default=True, MW=1),
        *cane.create_cellulosic_oilcane_chemicals(),
        bst.Chemical('HandSorted', search_db=False, default=True, phase='s'),
        bst.Chemical('Metal', search_db=False, default=True, phase='s'),
        bst.Chemical('Unders', search_db=False, default=True, phase='s'),
        bst.Chemical('Overs', search_db=False, default=True, phase='s'),
    ])
    chemicals.BulkPlastic.LHV = 23.22e3 * chemicals.BulkPlastic.MW # https://www.fire.tc.faa.gov/pdf/tn97-8.pdf
    chemicals.compile()
    chemicals.set_alias('Yeast', 'DryYeast')
    chemicals.set_alias('Yeast', 'Cellmass')
    chemicals.set_alias('Cellmass', 'Cells')
    chemicals.set_alias('Cellmass', 'cellmass')
    return chemicals