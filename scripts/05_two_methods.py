"""Density and Psat two ways: this package vs Milton's residual, no spinodal.

    py scripts\\05_two_methods.py

Both solve the same two equations with brentq. They differ in the residual
used for the density and in where the brackets come from.

OURS
    density : P_EoS(eta) - P = 0
    Psat    : ln phi_L(T, P) - ln phi_V(T, P) = 0
    brackets: from the spinodal

MILTON
    density : 1 - Z_PV / Z_EoS = 0
    Psat    : ln phi_L(T, P) - ln phi_V(T, P) = 0
    brackets: from scanning, eta between 1e-10 and 0.7405, P between 1 Pa
              and 1000 bar. No spinodal.

Change only the CONFIG block.
"""

import warnings
warnings.filterwarnings("ignore")

import time
import numpy as np
from scipy.optimize import brentq
import pcsaft
from pcsaft.eos import Z_hc, Z_disp
from pcsaft.constants import kB

# ============================== CONFIG ==================================
COMPOUND = "ethane"
T = 273.15        # K
P = 1.0e5         # Pa, for the density
# ========================================================================

c = pcsaft.get(COMPOUND)


# ----------------------------------------------------------- Milton's way

def Z_eos(eta, T):
    rho = pcsaft.eta_to_rho(eta, T, c)
    return 1 + Z_hc(T, rho, c) + Z_disp(T, rho, c)


def residual(eta, T, P):
    """Milton's density residual. Zero when the EoS and PV/RT give the same Z."""
    rho = pcsaft.eta_to_rho(eta, T, c)
    Z_pv = P / (rho * 1e30 * kB * T)
    return 1 - Z_pv / Z_eos(eta, T)


def roots_milton(T, P):
    """Every eta where the residual is zero, smallest first.

    Scan eta from the paper's vapour start (1e-10) to closest packing (0.7405)
    and solve each sign change with brentq.

    The residual has poles wherever Z_EoS = 0, and one of them sits just below
    the liquid root. If the pole and the root fall between the same two grid
    points, the residual changes sign twice there and the scan sees nothing:
    the liquid root is lost. So first locate the poles and split the grid at
    each one, then never bracket across a pole.
    """
    etas = list(np.logspace(-10, np.log10(0.7405), 600))
    z = [Z_eos(e, T) for e in etas]
    for k in range(len(etas) - 1):
        if z[k] * z[k + 1] < 0:
            pole = brentq(Z_eos, etas[k], etas[k + 1], args=(T,), xtol=1e-16)
            etas += [pole * (1 - 1e-9), pole * (1 + 1e-9)]
    etas = np.sort(etas)

    z = [Z_eos(e, T) for e in etas]
    r = [residual(e, T, P) for e in etas]
    roots = []
    for k in range(len(etas) - 1):
        crosses_pole = z[k] * z[k + 1] < 0
        if r[k] * r[k + 1] < 0 and not crosses_pole:
            roots.append(brentq(residual, etas[k], etas[k + 1], args=(T, P), xtol=1e-16))
    return roots


def density_milton(T, P, phase):
    roots = roots_milton(T, P)
    if not roots:
        return None
    return roots[-1] if phase == "liquid" else roots[0]


def f(P, T):
    """ln phi_L - ln phi_V.  None if the liquid and the vapour root are the same."""
    roots = roots_milton(T, P)
    if len(roots) < 2:
        return None
    rho_l = pcsaft.eta_to_rho(roots[-1], T, c)
    rho_v = pcsaft.eta_to_rho(roots[0], T, c)
    return pcsaft.ln_fugacity_coefficient(T, rho_l, c) - pcsaft.ln_fugacity_coefficient(T, rho_v, c)


def psat_milton(T):
    """Walk up in P until f changes from + to -, then brentq on that interval.

    Below Psat the liquid is the metastable phase, so f > 0; above it f < 0.
    Where only one root exists, f is meaningless and the point is skipped.
    """
    Ps = np.logspace(0, 8, 401)          # 1 Pa to 1000 bar
    last_P, last_f = None, None
    for p in Ps:
        val = f(p, T)
        if val is None or not np.isfinite(val):
            continue
        if last_f is not None and last_f > 0 and val < 0:
            return brentq(lambda x: f(x, T), last_P, p, xtol=1e-6, rtol=1e-12)
        last_P, last_f = p, val
    return None


# ------------------------------------------------------------------ run

print(f"{COMPOUND}   T = {T} K   P = {P/1e5:g} bar\n")

print("DENSITY             eta            rho [kg/m3]")
for phase in ("liquid", "vapour"):
    r = pcsaft.density(T, P, c, phase)
    print(f"  ours    {phase:7s} {pcsaft.packing_fraction(T, r, c):.10f}  {c.mass_density(r):12.4f}")
    e = density_milton(T, P, phase)
    if e is None:
        print(f"  Milton  {phase:7s} no root found")
    else:
        print(f"  Milton  {phase:7s} {e:.10f}  {c.mass_density(pcsaft.eta_to_rho(e, T, c)):12.4f}")

print("\nPSAT                Psat [bar]        time")
t0 = time.time()
sat = pcsaft.saturation(T, c)
t_ours = time.time() - t0
print("  ours    above Tc" if sat is None else f"  ours    {sat[0]/1e5:16.10f}   {t_ours:6.2f} s")

t0 = time.time()
Ps = psat_milton(T)
t_milton = time.time() - t0
print("  Milton  not found" if Ps is None else f"  Milton  {Ps/1e5:16.10f}   {t_milton:6.2f} s")
