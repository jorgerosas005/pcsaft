"""PC-SAFT for pure components: density and saturation pressure.

    import pcsaft
    ethane = pcsaft.get("ethane")
    rho = pcsaft.density(273.15, 1e5, ethane, "vapour")
    Psat, rho_L, rho_V = pcsaft.saturation(250.0, ethane)

Reference:
    Gross, J.; Sadowski, G. Perturbed-Chain SAFT: An Equation of State Based on
    a Perturbation Theory for Chain Molecules.
    Ind. Eng. Chem. Res. 2001, 40, 1244-1260.
"""

from .constants import kB, NA, R, ETA_MAX
from .parameters import Component, get, LIBRARY, N_ALKANES, WAGNER25
from .eos import (segment_diameter, zetas, packing_fraction, eta_to_rho,
                  radial_distribution,
                  a_hs, a_hc, a_disp, a_assoc, a_res,
                  compressibility, compressibility_numeric,
                  pressure, ln_fugacity_coefficient)
from .equilibrium import density, wagner25, saturation

__version__ = "0.2.0"
