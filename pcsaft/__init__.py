"""PC-SAFT and, in time, PC-SAFT with a two-state association term.

    import pcsaft
    hexane = pcsaft.get("hexane")
    P = pcsaft.pressure(300.0, rho, hexane)
    Psat, rho_L, rho_V = pcsaft.saturation(300.0, hexane)

Reference:
    Gross, J.; Sadowski, G. Perturbed-Chain SAFT: An Equation of State Based on
    a Perturbation Theory for Chain Molecules.
    Ind. Eng. Chem. Res. 2001, 40, 1244-1260.
"""

from .constants import kB, NA, R, ETA_MAX
from .parameters import Component, get, LIBRARY, N_ALKANES
from .eos import (segment_diameter, zetas, packing_fraction, eta_to_rho,
                  radial_distribution,
                  a_hs, a_hc, a_disp, a_assoc, a_res,
                  compressibility, compressibility_numeric,
                  pressure, ln_fugacity_coefficient)
from .equilibrium import (density, spinodal, saturation, saturation_curve,
                          critical_point)
from .data import (load_nist_saturation, check_saturation_data, save_results,
                   DATA_DIR, RESULTS_DIR)
from .validation import (check_consistency, check_ideal_gas_limit,
                         check_diameter_limits, compare_saturation, aad, report)
from .caloric import (u_res, h_res, s_res, cv_res, cp_res, cp_minus_cv,
                      enthalpy_of_vaporisation, compare_caloric)
from . import plotting, caloric

__version__ = "0.1.0"
