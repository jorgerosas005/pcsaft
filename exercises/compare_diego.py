"""
The residual is the one Diego solves:

    Z_PV  = P / (rho_number * kB * T)        from the proposed eta
    Z_EoS = 1 + Z_hc + Z_disp                from PC-SAFT
    residual(eta) = 1 - Z_PV / Z_EoS


CAUTION on the ratio form. Z_EoS passes through zero just below the saturated
liquid root, so 1 - Z_PV/Z_EoS has a POLE there and changes sign for a reason
that has nothing to do with the root. A bracketing solver can latch onto the
pole instead of the solution. The difference form, Z_EoS - Z_PV, has the same
root and no pole, which is why the solver below uses it.
"""

import warnings
warnings.filterwarnings("ignore")

import numpy as np
from scipy.optimize import brentq

import pcsaft
from pcsaft.constants import kB, ETA_MAX
from pcsaft.eos import Z_hc, Z_disp

# ======================== CONFIG - edit these ============================

COMPOUND = "ethane"

T = 273.15            # temperature, K
P = 100000           # pressure, Pa   (1.0e6 Pa = 10 bar)
ETA_TRIAL = 0.50     # the eta you want to test by hand

SHOW_TABLE = True    # also print a reference table to send Diego
TABLE_STATES = [     # (T [K], P [Pa])
    (273.15, 1.0e5),
    (273.15, 1.0e6),
    (273.15, 2.0e6),
    (273.15, 2.0e6),
    (273.15, 2.0e6),
    (273.15, 4.0e6),
]
TABLE_TSAT = [180.0, 200.0, 220.0, 240.0, 260.0, 280.0, 300.0]   # for Psat

# =========================================================================

comp = pcsaft.get(COMPOUND)


def z_from_PV(T, P, eta):
    """Z = P V / (R T), evaluated at the proposed eta.

    rho comes out of eta in molecules/A^3, so multiply by 1e30 to reach
    molecules/m^3:
        Pa / (m^-3 * J K^-1 * K) = Pa / (J m^-3) = Pa / Pa = dimensionless
    """
    rho = pcsaft.eta_to_rho(eta, T, comp)          # molecules / A^3
    return P / (rho * 1e30 * kB * T)


def z_from_eos(T, eta):
    """Z from PC-SAFT at the same eta."""
    rho = pcsaft.eta_to_rho(eta, T, comp)
    return 1.0 + Z_hc(T, rho, comp) + Z_disp(T, rho, comp)


def residual(eta, T, P):
    """Diego's objective:  1 - Z_PV / Z_EoS.  Zero when the two agree."""
    return 1.0 - z_from_PV(T, P, eta) / z_from_eos(T, eta)


def residual_abs(eta, T, P):
    """Same root, no pole:  Z_EoS - Z_PV."""
    return z_from_eos(T, eta) - z_from_PV(T, P, eta)


def all_roots(T, P, n=4000):
    """Every eta satisfying the residual at this (T, P), by bracketing."""
    rho_ideal = P / (kB * T) / 1e30
    lo = 0.05 * np.pi / 6 * rho_ideal * comp.m * pcsaft.segment_diameter(T, comp)**3
    grid = np.logspace(np.log10(max(lo, 1e-20)), np.log10(ETA_MAX), n)
    f = np.array([residual_abs(e, T, P) for e in grid])
    out = []
    for k in range(len(grid) - 1):
        if np.isfinite(f[k]) and np.isfinite(f[k+1]) and f[k] * f[k+1] < 0:
            out.append(brentq(residual_abs, grid[k], grid[k+1], args=(T, P), xtol=1e-16))
    return out


def label_roots(roots):
    return ["vapour", "unstable", "liquid"] if len(roots) == 3 else ["single"] * len(roots)


def line(title):
    print("\n" + "=" * 72)
    print(title)
    print("=" * 72)


# ---------------------------------------------------------------- header

print(f"{comp.name}:  m = {comp.m}   sigma = {comp.sigma} A   eps/k = {comp.eps_k} K")
print(f"M = {comp.M} g/mol")
print(f"source: {comp.source}")

# ------------------------------------------- 1. residual at the trial eta

line(f"1. residual at the proposed eta = {ETA_TRIAL}")
rho_t = pcsaft.eta_to_rho(ETA_TRIAL, T, comp)
zpv_t, zeos_t = z_from_PV(T, P, ETA_TRIAL), z_from_eos(T, ETA_TRIAL)
print(f"  T            = {T:.4f} K")
print(f"  P            = {P:.6g} Pa = {P/1e5:.6g} bar")
print(f"  d(T)         = {pcsaft.segment_diameter(T, comp):.8f} A")
print(f"  eta          = {ETA_TRIAL:.10f}")
print(f"  rho          = {rho_t:.8e} molecules/A^3")
print(f"  rho_mass     = {comp.mass_density(rho_t):.6f} kg/m3")
print(f"  Z_PV         = {zpv_t:.10f}")
print(f"  Z_EoS        = {zeos_t:.10f}")
print(f"  1 - Zpv/Zeos = {1 - zpv_t/zeos_t:+.8e}")
print(f"  Zeos - Zpv   = {zeos_t - zpv_t:+.8e}")

# ------------------------------------------------- 2. roots at this (T,P)

line(f"2. roots at T = {T} K, P = {P/1e5:g} bar")
roots = all_roots(T, P)
print("  root        eta                Z_PV            Z_EoS        rho[kg/m3]")
for lab, e in zip(label_roots(roots), roots):
    rho = pcsaft.eta_to_rho(e, T, comp)
    print(f"  {lab:9s} {e:.12f}  {z_from_PV(T,P,e):14.10f} {z_from_eos(T,e):14.10f}  "
          f"{comp.mass_density(rho):10.4f}")
print("\n  package density(), for comparison:")
for phase in ("vapour", "liquid"):
    r = pcsaft.density(T, P, comp, phase)
    if r is not None:
        print(f"    {phase:7s} eta = {pcsaft.packing_fraction(T, r, comp):.12f}")

# ------------------------------- 3. the four Z at the saturation point

line(f"3. phase equilibrium at T = {T} K  (the four Z)")
sat = pcsaft.saturation(T, comp)
if sat is None:
    print("  no phase split at this T (above Tc)")
else:
    Psat, rho_l, rho_v = sat
    eta_l = pcsaft.packing_fraction(T, rho_l, comp)
    eta_v = pcsaft.packing_fraction(T, rho_v, comp)
    print(f"  Psat = {Psat:.8f} Pa = {Psat/1e5:.10f} bar\n")
    print("  phase      eta                 Z_PV            Z_EoS       "
          "1 - Zpv/Zeos       ln phi")
    for lab, e in (("liquid", eta_l), ("vapour", eta_v)):
        zpv = z_from_PV(T, Psat, e)
        zeos = z_from_eos(T, e)
        rho = pcsaft.eta_to_rho(e, T, comp)
        lnphi = pcsaft.ln_fugacity_coefficient(T, rho, comp)
        print(f"  {lab:8s} {e:.12f}  {zpv:14.10f} {zeos:14.10f}  "
              f"{1 - zpv/zeos:+13.4e}  {lnphi:+.10f}")

    lnphi_l = pcsaft.ln_fugacity_coefficient(T, rho_l, comp)
    lnphi_v = pcsaft.ln_fugacity_coefficient(T, rho_v, comp)
    print(f"\n  equilibrium criterion:  ln phi_L - ln phi_V = {lnphi_l - lnphi_v:+.3e}")
    print(f"  densities:  rho_L = {comp.mass_density(rho_l):.6f} kg/m3   "
          f"rho_V = {comp.mass_density(rho_v):.6f} kg/m3")

# ---------------------------------------------------------- 4. the table

if SHOW_TABLE:
    line("4. reference table for Diego")
    print("  (a) roots at fixed (T, P)\n")
    print("      T[K]      P[bar]  root       eta                Z_EoS")
    for Tx, Px in TABLE_STATES:
        for lab, e in zip(*(lambda r: (label_roots(r), r))(all_roots(Tx, Px))):
            print(f"  {Tx:8.2f} {Px/1e5:10.4g}  {lab:9s} {e:.12f} {z_from_eos(Tx, e):15.10f}")

    print("\n  (b) saturation\n")
    print("      T[K]       Psat[bar]        eta_L            eta_V"
          "          Z_L           Z_V")
    for Tx in TABLE_TSAT:
        s = pcsaft.saturation(Tx, comp)
        if s is None:
            print(f"  {Tx:8.2f}   no phase split")
            continue
        Ps, rl, rv = s
        el = pcsaft.packing_fraction(Tx, rl, comp)
        ev = pcsaft.packing_fraction(Tx, rv, comp)
        print(f"  {Tx:8.2f} {Ps/1e5:15.10f}  {el:.12f}  {ev:.12e}  "
              f"{z_from_eos(Tx, el):11.8f}  {z_from_eos(Tx, ev):11.8f}")

    print("\n  'single' = only one root: the isotherm is above the spinodal maximum,")
    print("  so vapour and liquid are not distinguishable at that (T, P).")
