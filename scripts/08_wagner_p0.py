"""Psat by successive substitution, starting from Wagner 25.

    py scripts\\08_wagner_p0.py

    P0      = Wagner 25 (NIST form, as in Aspen's WAGNER25 parameter)
    density = Newton from eta = 0.5 / 1e-10 on the paper's residual
    Psat    = P <- P * phi_L / phi_V   until it stops changing

If P0 falls outside the loop (no liquid root there), it says so, and the
spinodal-based pcsaft.saturation() is used instead.

Wagner 25 (Aspen documentation, "NIST Wagner 25 Liquid Vapor Pressure"):

    ln P = ln Pc + (C1 tau + C2 tau^1.5 + C3 tau^2.5 + C4 tau^5) / Tr
    Tr = T / Tc,  tau = 1 - Tr
"""

import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pcsaft
from pcsaft.eos import Z_hc, Z_disp
from pcsaft.constants import kB

# ============================== CONFIG ==================================
COMPOUND = "ethane"

# WAGNER25/1..8 for this compound, copied from Aspen.
# C1, C2, C3, C4 are dimensionless. Write ln_Pc for Pc in PASCAL and Tc in K
# (if Aspen shows ln Pc for another pressure unit, convert it: ln(Pc in Pa)).
C1, C2, C3, C4 = None, None, None, None
ln_Pc = None        # ln(Pc / Pa)
Tc = None           # K

TEMPERATURES = [150, 200, 250, 280, 295, 300, 302, 304, 305, 307, 308]   # K
# ========================================================================

c = pcsaft.get(COMPOUND)
count = [0]


def wagner(T):
    Tr = T / Tc
    if Tr >= 1:
        return None
    tau = 1 - Tr
    return np.exp(ln_Pc + (C1 * tau + C2 * tau**1.5 + C3 * tau**2.5 + C4 * tau**5) / Tr)


def residual(eta, T, P):
    count[0] += 1
    rho = pcsaft.eta_to_rho(eta, T, c)
    Z_eos = 1 + Z_hc(T, rho, c) + Z_disp(T, rho, c)
    Z_pv = P / (rho * 1e30 * kB * T)
    return Z_eos / Z_pv - 1                       # P_EoS / P - 1, the paper's


def density(T, P, phase):
    eta = 0.5 if phase == "liquid" else 1e-10
    for k in range(200):
        h = 1e-6 * eta
        slope = (residual(eta + h, T, P) - residual(eta - h, T, P)) / (2 * h)
        step = -residual(eta, T, P) / slope
        while not 0 < eta + step < 0.7405:
            step /= 2
        if abs(step) < 1e-12 * eta:
            return eta + step
        eta += step
    return None


def psat_ss(T, P0):
    P = P0
    for k in range(300):
        eL, eV = density(T, P, "liquid"), density(T, P, "vapour")
        if eL is None or eV is None or abs(eL - eV) < 1e-8:
            return None, k                         # P is outside the loop
        lnL = pcsaft.ln_fugacity_coefficient(T, pcsaft.eta_to_rho(eL, T, c), c)
        lnV = pcsaft.ln_fugacity_coefficient(T, pcsaft.eta_to_rho(eV, T, c), c)
        P_new = P * np.exp(lnL - lnV)
        if abs(P_new / P - 1) < 1e-8:
            return P_new, k + 1
        P = P_new
    return None, 300


if None in (C1, C2, C3, C4, ln_Pc, Tc):
    raise SystemExit("Fill in the WAGNER25 constants in the CONFIG block first.")

print(f"{COMPOUND}: Wagner Tc = {Tc} K,  PC-SAFT Tc = {pcsaft.critical_point(c)[0]:.2f} K\n")
print("  T [K]   P0 Wagner   Psat SS      its  evals   Psat spinodal   note")
for T in TEMPERATURES:
    P0 = wagner(T)
    sat = pcsaft.saturation(T, c)
    ref = "  above Tc" if sat is None else f"{sat[0]/1e5:14.6f}"
    if P0 is None:
        print(f"  {T:6.1f}   {'-':>9s}   {'-':>10s}   {'':>4s} {'':>6s}  {ref}   T above Wagner Tc")
        continue
    count[0] = 0
    P, its = psat_ss(T, P0)
    if P is None:
        print(f"  {T:6.1f}   {P0/1e5:9.4f}   {'-':>10s}   {its:4d} {count[0]:6d}  {ref}   P0 outside the loop")
    else:
        print(f"  {T:6.1f}   {P0/1e5:9.4f}   {P/1e5:10.6f}   {its:4d} {count[0]:6d}  {ref}")
