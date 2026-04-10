# -*- coding: utf-8 -*-
"""
"""
import biosteam as bst
import thermosteam as tmo
from biorefineries.tea import CellulosicEthanolTEA as STRAPTEA
__all__ = (
    'STRAPTEA',
    'create_baseline_tea',
)

# class STRAPTEA(bst.TEA):
    
#     __slots__ = ('OSBL_units', 'warehouse', 'site_development',
#                  'additional_piping', 'proratable_costs', 'field_expenses',
#                  'construction', 'contingency', 'other_indirect_costs', 
#                  'labor_cost', 'labor_burden', 'property_insurance',
#                  'maintenance', '_ISBL_DPI_cached', '_FCI_cached',
#                  '_utility_cost_cached', '_steam_power_depreciation',
#                  '_steam_power_depreciation_array',
#                  'boiler_turbogenerator', '_DPI_cached')
    
#     def __init__(self, system, IRR, duration, depreciation, income_tax,
#                  operating_days, lang_factor, construction_schedule,
#                  startup_months, startup_FOCfrac, startup_VOCfrac,
#                  startup_salesfrac, WC_over_FCI,  finance_interest,
#                  finance_years, finance_fraction, OSBL_units, warehouse,
#                  site_development, additional_piping, proratable_costs,
#                  field_expenses, construction, contingency,
#                  other_indirect_costs, labor_cost, labor_burden,
#                  property_insurance, maintenance, steam_power_depreciation,
#                  boiler_turbogenerator):
#         super().__init__(system, IRR, duration, depreciation, income_tax,
#                          operating_days, lang_factor, construction_schedule,
#                          startup_months, startup_FOCfrac, startup_VOCfrac,
#                          startup_salesfrac, WC_over_FCI,  finance_interest,
#                          finance_years, finance_fraction)
#         self.OSBL_units = OSBL_units
#         self.warehouse = warehouse
#         self.site_development = site_development
#         self.additional_piping = additional_piping
#         self.proratable_costs = proratable_costs
#         self.field_expenses = field_expenses
#         self.construction = construction
#         self.contingency = contingency
#         self.other_indirect_costs = other_indirect_costs
#         self.labor_cost = labor_cost
#         self.labor_burden = labor_burden
#         self.property_insurance = property_insurance
#         self.maintenance = maintenance
#         self.steam_power_depreciation = steam_power_depreciation
#         self.boiler_turbogenerator = boiler_turbogenerator
        
#     @property
#     def steam_power_depreciation(self):
#         """[str] 'MACRS' + number of years (e.g. 'MACRS7')."""
#         return self._steam_power_depreciation
#     @steam_power_depreciation.setter
#     def steam_power_depreciation(self, depreciation):
#         self._steam_power_depreciation_array = self._depreciation_array_from_key(
#             self._depreciation_key_from_name(depreciation)
#         )
#         self._steam_power_depreciation = depreciation
    
#     @property
#     def ISBL_installed_equipment_cost(self):
#         return self._ISBL_DPI(self.installed_equipment_cost)
    
#     @property
#     def OSBL_installed_equipment_cost(self):
#         if self.lang_factor:
#             raise NotImplementedError('lang factor cannot yet be used')
#         elif isinstance(self.system, bst.AgileSystem):
#             unit_capital_costs = self.system.unit_capital_costs
#             OSBL_units = self.OSBL_units
#             return sum([unit_capital_costs[i].installed_cost for i in OSBL_units])
#         else:
#             return sum([i.installed_cost for i in self.OSBL_units])
    
#     def _fill_depreciation_array(self, D, start, years, TDC):
#         depreciation_array = self._get_depreciation_array()
#         N_depreciation_years = depreciation_array.size
#         if N_depreciation_years > years:
#             raise RuntimeError('depreciation schedule is longer than plant lifetime')
#         system = self.system
#         BT = self.boiler_turbogenerator
#         if BT is None:
#             D[start:start + N_depreciation_years] = TDC * depreciation_array
#         else:
#             if isinstance(system, bst.AgileSystem): BT = system.unit_capital_costs[BT]
#             BT_TDC = BT.installed_cost 
#             D[start:start + N_depreciation_years] = (TDC - BT_TDC) * depreciation_array
            
#             depreciation_array = self._steam_power_depreciation_array
#             N_depreciation_years = depreciation_array.size
#             if N_depreciation_years > years:
#                 raise RuntimeError('steam power depreciation schedule is longer than plant lifetime')
#             D[start:start + N_depreciation_years] += BT_TDC * depreciation_array
    
#     def _ISBL_DPI(self, installed_equipment_cost):
#         """Direct permanent investment of units inside battery limits."""
#         if self.lang_factor:
#             raise NotImplementedError('lang factor cannot yet be used')
#         else:
#             self._ISBL_DPI_cached = installed_equipment_cost - self.OSBL_installed_equipment_cost
#         return self._ISBL_DPI_cached
        
#     def _DPI(self, installed_equipment_cost): # Direct Permanent Investment
#         factors = self.warehouse + self.site_development + self.additional_piping
#         self._DPI_cached = DPI = installed_equipment_cost + self._ISBL_DPI(installed_equipment_cost) * factors
#         return DPI
    
#     def _TDC(self, DPI): # Total Depreciable Capital
#         return DPI + self._depreciable_indirect_costs(DPI)
    
#     def _nondepreciable_indirect_costs(self, DPI):
#         return DPI * self.other_indirect_costs
    
#     def _depreciable_indirect_costs(self, DPI):
#         return DPI * (self.proratable_costs + self.field_expenses
#                       + self.construction + self.contingency)
    
#     def _FCI(self, TDC): # Fixed Capital Investment
#         self._FCI_cached = FCI = TDC + self._nondepreciable_indirect_costs(self._DPI_cached)
#         return FCI
    
#     def _FOC(self, FCI): # Fixed Operating Costs
#         return (FCI * self.property_insurance
#                 + self._ISBL_DPI_cached * self.maintenance
#                 + self.labor_cost * (1 + self.labor_burden))

def create_baseline_tea(sys, OSBL_units=None, cls=None, years=30):
    # Default parameters are from NREL 2011 report for cornstover ethanol
    if OSBL_units is None: OSBL_units = bst.get_OSBL(sys.cost_units)
    try:
        BT = tmo.utils.get_instance(OSBL_units, (bst.BoilerTurbogenerator, bst.Boiler))
    except:
        BT = None
    if cls is None: cls = STRAPTEA
    tea = cls(
        system=sys, 
        IRR=0.20, 
        duration=(2025, 2025 + years),
        depreciation='MACRS7', 
        income_tax=0.21, #switch from 35 on 8/1/2025
        operating_days=350.4,
        lang_factor=None, 
        construction_schedule=(0.08, 0.60, 0.32),
        startup_months=3, 
        startup_FOCfrac=1,
        startup_salesfrac=0.5,
        startup_VOCfrac=0.75,
        WC_over_FCI=0.05,
        finance_interest=0.08,
        finance_years=10,
        finance_fraction= 0.0,        #0.6,
        OSBL_units=OSBL_units,
        warehouse=0.04, 
        site_development=0.09, 
        additional_piping=0.045,
        proratable_costs=0.10,
        field_expenses=0.10,
        construction=0.20,
        contingency=0.4,
        other_indirect_costs=0.10, 
        labor_cost=10 * 60000,
        labor_burden=0.90,
        property_insurance=0.007, 
        maintenance=0.03,
        steam_power_depreciation='MACRS20',
        boiler_turbogenerator=BT)
    return tea



# def create_baseline_tea(system):
#     return STRAPTEA(
#         system=system,
#         IRR=0.10,
#         duration=[2022, 2042],
#         depreciation="MACRS7",
#         income_tax=0.25,
#         operating_days=8000/24,
#         construction_schedule=(0.4, 0.6),
#         WC_over_FCI=0.10,
#         labor_cost=208050,
#         fringe_benefits=0.5,
#         property_tax=0,
#         property_insurance=0.007,
#         supplies=0,
#         maintenance=0.03,
#         administration=0.00,
#         contingency=1.0,
#     )

# class STRAPTEA(bst.TEA):
#     """
#     Create a STRAPTEA object for techno-economic analysis.

#     Parameters
#     ----------
#     system : System
#         Should contain feed and product streams.
#     IRR : float
#         Internal rate of return (fraction).
#     duration : tuple[int, int]
#         Start and end year of venture (e.g. (2018, 2038)).
#     depreciation : str
#         'MACRS' + number of years (e.g. 'MACRS7').
#     operating_days : float
#         Number of operating days per year.
#     income_tax : float
#         Combined federal and state income tax rate (fraction).
#     lang_factor : float
#         Lang factor for getting fixed capital investment from
#         total purchase cost. If no lang factor, estimate capital investment
#         using bare module factors.
#     startup_schedule : tuple[float]
#         Startup investment fractions per year
#         (e.g. (0.5, 0.5) for 50% capital investment in the first year and 50%
#         investment in the second).
#     WC_over_FCI : float
#         Working capital as a fraction of fixed capital investment.
#     labor_cost : float
#         Total labor cost (USD/yr).
#     fringe_benefits : float
#         Cost of fringe benefits as a fraction of labor cost.
#     property_tax : float
#         Fee as a fraction of fixed capital investment.
#     property_insurance : float
#         Fee as a fraction of fixed capital investment.
#     supplies : float
#         Yearly fee as a fraction of labor cost.
#     maintenance : float
#         Yearly fee as a fraction of fixed capital investment.
#     administration : float
#         Yearly fee as a fraction of fixed capital investment.

#     References
#     ----------
#     .. [1] Huang, H., Long, S., & Singh, V. (2016). Techno-economic analysis of biodiesel
#         and ethanol co-production from lipid-producing sugarcane. Biofuels, Bioproducts
#         and Biorefining, 10(3), 299–315. https://doi.org/10.1002/bbb.1640
    
#     """

#     __slots__ = (
#         "labor_cost",
#         "fringe_benefits",
#         "maintenance",
#         "property_tax",
#         "property_insurance",
#         "_FCI_cached",
#         "supplies",
#         "maintanance",
#         "administration",
#         'OSBL_factor',
#         'engineering_factor',
#         'contingency'
#     )

#     def __init__(
#         self,
#         system,
#         IRR,
#         duration,
#         depreciation,
#         income_tax,
#         operating_days,
#         lang_factor,
#         construction_schedule,
#         WC_over_FCI,
#         labor_cost,
#         fringe_benefits,
#         property_tax,
#         property_insurance,
#         supplies,
#         maintenance,
#         administration,
#         contingency,
#     ):
#         super().__init__(
#             system,
#             IRR,
#             duration,
#             depreciation,
#             income_tax,
#             operating_days,
#             lang_factor,
#             construction_schedule,
#             startup_months=0,
#             startup_FOCfrac=0,
#             startup_VOCfrac=0,
#             startup_salesfrac=0,
#             finance_interest=0,
#             finance_years=0,
#             finance_fraction=0,
#             WC_over_FCI=WC_over_FCI,
#         )
#         self.labor_cost = labor_cost
#         self.fringe_benefits = fringe_benefits
#         self.property_tax = property_tax
#         self.property_insurance = property_insurance
#         self.supplies = supplies
#         self.maintenance = maintenance
#         self.administration = administration
#         self.contingency = contingency

#     def _DPI(self, installed_equipment_cost):
#         return installed_equipment_cost * (1 + self.contingency)
    
#     def _FOC(self, FCI):
#         return FCI * (
#             self.property_tax
#             + self.property_insurance
#             + self.maintenance
#             + self.administration
#         ) + self.labor_cost + (self.fringe_benefits * self.labor_cost + self.supplies)