"""Density and Psat two ways: this package vs Milton's (Diego's) method.

    py scripts\\05_two_methods.py

OURS
    density : brentq on  P_EoS(eta) - P,  brackets from the spinodal / a scan
    Psat    : brentq on  ln phi_L - ln phi_V,  pressure window from the spinodal

MILTON
    density : Newton on a Z residual, starting from eta = 0.5 (liquid)
              or 1e-10 (vapour), the values in Gross & Sadowski 2001
    Psat    : successive substitution  P_new = P * phi_L / phi_V
              starting from a guessed pressure P0

Change only the CONFIG block.
"""

import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pcsaft
from pcsaft.eos import Z_hc, Z_disp
from pcsaft.constants import kB

# ============================== CONFIG ==================================
COMPOUND = "ethane"
T = 273.15        # K
P = 1.0e5         # Pa, for the density
P0 = 1.0e5        # Pa, Milton's starting guess for Psat
# ========================================================================

c = pcsaft.get(COMPOUND)


# ----------------------------------------------------------- Milton's way

def residual(eta, T, P):
    """Zero when the EoS and PV/RT give the same Z."""
    rho = pcsaft.eta_to_rho(eta, T, c)
    Z_eos = 1 + Z_hc(T, rho, c) + Z_disp(T, rho, c)
    Z_pv = P / (rho * 1e30 * kB * T)
    return Z_eos - Z_pv
    # Milton writes it as  1 - Z_pv / Z_eos  (same root).
    # Swap the line above for that one and Newton never finds the liquid
    # from eta = 0.5: it is almost flat there, and the first step flies off.


def density_milton(T, P, phase):
    """Newton from the paper's starting values. Returns (eta, iterations)."""
    eta = 0.5 if phase == "liquid" else 1e-10
    for k in range(100):
        h = 1e-6 * eta
        slope = (residual(eta + h, T, P) - residual(eta - h, T, P)) / (2 * h)
        new = eta - residual(eta, T, P) / slope
        if not 0 < new < 0.7405:
            return None, k
        if abs(new - eta) < 1e-12 * eta:
            return new, k + 1
        eta = new
    return None, 100


def psat_milton(T, P0):
    """Successive substitution. Returns (Psat, eta_L, eta_V, iterations)."""
    P = P0
    for k in range(500):
        eL, _ = density_milton(T, P, "liquid")
        eV, _ = density_milton(T, P, "vapour")
        if eL is None or eV is None or abs(eL - eV) < 1e-8:
            return None, None, None, k    # no liquid root at this P
        lnphiL = pcsaft.ln_fugacity_coefficient(T, pcsaft.eta_to_rho(eL, T, c), c)
        lnphiV = pcsaft.ln_fugacity_coefficient(T, pcsaft.eta_to_rho(eV, T, c), c)
        P_new = P * np.exp(lnphiL - lnphiV)
        if abs(P_new / P - 1) < 1e-10:
            return P_new, eL, eV, k + 1
        P = P_new
    return None, None, None, 500


# ------------------------------------------------------------------ run

print(f"{COMPOUND}   T = {T} K   P = {P/1e5:g} bar   P0 = {P0/1e5:g} bar\n")

print("DENSITY             eta            rho [kg/m3]   Newton its")
for phase in ("liquid", "vapour"):
    r = pcsaft.density(T, P, c, phase)
    print(f"  ours    {phase:7s} {pcsaft.packing_fraction(T, r, c):.10f}  {c.mass_density(r):12.4f}")
    e, k = density_milton(T, P, phase)
    if e is None:
        print(f"  Milton  {phase:7s} did not converge")
    else:
        print(f"  Milton  {phase:7s} {e:.10f}  {c.mass_density(pcsaft.eta_to_rho(e, T, c)):12.4f}   {k:4d}")

print("\nPSAT                Psat [bar]     eta_L          eta_V          its")
sat = pcsaft.saturation(T, c)
if sat is None:
    print("  ours    above Tc")
else:
    Ps, rl, rv = sat
    print(f"  ours    {Ps/1e5:14.8f}  {pcsaft.packing_fraction(T, rl, c):.10f}  "
          f"{pcsaft.packing_fraction(T, rv, c):.6e}")
Ps, eL, eV, k = psat_milton(T, P0)
if Ps is None:
    print(f"  Milton  failed: no liquid root at P0 = {P0/1e5:g} bar. Try a P0 closer to Psat.")
else:
    print(f"  Milton  {Ps/1e5:14.8f}  {eL:.10f}  {eV:.6e}  {k:4d}")
