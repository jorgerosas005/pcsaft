"""Which inner density residual works best inside saturation()?

    py scripts\\10_residual_comparison.py

Three residuals for the density at (T, P), all with the same Newton:
start at eta = 0.5 (liquid) or 1e-10 (vapour), central-difference derivative,
a step that leaves 0 < eta < ETA_MAX is halved, stop when |d eta| < 1e-12 eta.

    paper  :  P_EoS / P - 1
    milton :  1 - Z_calc / Z_EoS,     Z_calc = P / (rho k T)
    main   :  P_EoS - P               (what main's density() solves, with brentq)

milton = paper / (1 + paper): same roots, but large values are squashed to ~1
(flat, so Newton takes a huge step) and there is a pole where P_EoS = 0.

Saturation: P0 from Wagner 25, then Newton on g = ln phi_L - ln phi_V with a
forward difference dg/dP (2 extra density solves per iteration). Only the inner
residual changes between rows. "paper, analytic" uses dg/dP = (Z_L - Z_V)/P,
as in the core branch, for reference.

Reference values: main's density() (grid + brentq) and saturation() (spinodal).
Cost = evaluations of Z_EoS, counted everywhere (residuals and ln phi).
"""

import time
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pcsaft
import pcsaft.eos
from pcsaft.constants import kB, ETA_MAX

# ============================== CONFIG ==================================
COMPOUNDS = ["methane", "ethane", "propane", "butane", "pentane", "hexane",
             "heptane", "octane", "nonane", "decane", "dodecane"]
N_T = 25                 # temperatures from the lowest NIST T to 0.98 Tc, plus 3 near Tc
DP_REL = 1e-4            # forward-difference step for dg/dP, relative to P
FIG_COMPOUND, FIG_T = "ethane", 250.0
# ========================================================================

# lowest NIST temperature of each compound, K
T_MIN = {"methane": 90.69, "ethane": 90.37, "propane": 184.42, "butane": 134.89,
         "pentane": 143.47, "hexane": 180.0, "heptane": 182.55, "octane": 216.37,
         "nonane": 219.70, "decane": 243.50, "dodecane": 263.60}

# count every evaluation of Z_EoS (pressure() and ln phi call it too)
calls = [0]
Z_original = pcsaft.eos.compressibility


def Z_counted(T, rho, comp):
    calls[0] += 1
    return Z_original(T, rho, comp)


pcsaft.eos.compressibility = Z_counted


def residual(kind, eta, T, P, c):
    rho = pcsaft.eta_to_rho(eta, T, c)
    Z_eos = pcsaft.eos.compressibility(T, rho, c)
    Z_calc = P / (rho * 1e30 * kB * T)
    if kind == "paper":
        return Z_eos / Z_calc - 1
    if kind == "milton":
        return 1 - Z_calc / Z_eos
    return Z_eos * rho * 1e30 * kB * T - P          # main: P_EoS - P, in Pa


def newton_density(kind, T, P, c, phase):
    """Returns (eta or None, iterations)."""
    eta = 0.5 if phase == "liquid" else 1e-10
    for k in range(100):
        h = 1e-6 * eta
        f = residual(kind, eta, T, P, c)
        slope = (residual(kind, eta + h, T, P, c) - residual(kind, eta - h, T, P, c)) / (2 * h)
        if slope == 0 or not np.isfinite(slope) or not np.isfinite(f):
            return None, k + 1
        step = -f / slope
        while not 0 < eta + step < ETA_MAX:
            step = step / 2
        eta = eta + step
        if abs(step) < 1e-12 * eta:
            return eta, k + 1
    return None, 100


def wagner25(T, c):
    C1, C2, C3, C4, ln_Pc, Tc = c.wagner25
    Tr = T / Tc
    if Tr >= 1:
        return None
    tau = 1 - Tr
    return np.exp(ln_Pc + (C1 * tau + C2 * tau**1.5 + C3 * tau**2.5 + C4 * tau**5) / Tr)


def g_at(kind, T, P, c, its):
    """ln phi_L - ln phi_V at P, plus Z_L, Z_V; None if a density fails."""
    eL, nL = newton_density(kind, T, P, c, "liquid")
    eV, nV = newton_density(kind, T, P, c, "vapour")
    its[1] += nL + nV
    if eL is None or eV is None or abs(eL / eV - 1) < 1e-6:
        return None
    rL, rV = pcsaft.eta_to_rho(eL, T, c), pcsaft.eta_to_rho(eV, T, c)
    g = pcsaft.eos.ln_fugacity_coefficient(T, rL, c) - pcsaft.eos.ln_fugacity_coefficient(T, rV, c)
    return g, pcsaft.eos.compressibility(T, rL, c), pcsaft.eos.compressibility(T, rV, c)


def newton_saturation(kind, T, c, analytic=False):
    """Returns (P or None, why, [outer its, inner its])."""
    its = [0, 0]
    P = wagner25(T, c)
    if P is None:
        return None, "above Wagner Tc", its
    step_old = 1.0
    for k in range(50):
        its[0] += 1
        r = g_at(kind, T, P, c, its)
        if r is None:
            return None, "density failed", its
        g, Z_l, Z_v = r
        if analytic:
            step = g / (Z_v - Z_l)                          # relative step dP/P
        else:
            r2 = g_at(kind, T, P * (1 + DP_REL), c, its)
            if r2 is None:
                return None, "density failed", its
            dgdP = (r2[0] - g) / (P * DP_REL)
            step = -g / dgdP / P
        if abs(step) < 1e-10:
            return P, "ok", its
        if abs(step) >= abs(step_old) and abs(step) < 1e-5:  # round-off floor
            return P, "ok", its
        step_old = step
        P = P * (1 + step)
    return None, "no convergence", its


# ================================================== 1. density, on its own
print("1. Density at (T, P), 6 compounds x 6 T x 6 P, both phases\n")
rows = []
for name in ["methane", "ethane", "propane", "hexane", "decane", "dodecane"]:
    c = pcsaft.get(name)
    for T in np.linspace(c.T_range[0] + 5, c.T_range[1], 6):
        sp = pcsaft.spinodal(T, c)                       # (eta_v, eta_l, P_max, P_min) or None
        for P in [1e2, 1e4, 1e5, 1e6, 5e6, 2e7]:
            for phase in ("liquid", "vapour"):
                # does the asked phase exist at this (T, P)? (above Tc: one root, counts as both)
                if sp is None:
                    exists = True
                elif phase == "vapour":
                    exists = P < sp[2]
                else:
                    exists = P > sp[3]
                ref = pcsaft.density(T, P, c, phase)
                e_ref = pcsaft.packing_fraction(T, ref, c)
                for kind in ("paper", "milton", "main"):
                    calls[0] = 0
                    t0 = time.perf_counter()
                    e, its = newton_density(kind, T, P, c, phase)
                    dt = time.perf_counter() - t0
                    if not exists:
                        why = "absent, returned None" if e is None else "absent, returned a root"
                    elif e is None:
                        why = "no convergence"
                    elif abs(e / e_ref - 1) < 1e-6:
                        why = "ok"
                    else:
                        why = "other root"
                    rows.append({"compound": name, "T": T, "P": P, "phase": phase, "exists": exists,
                                 "residual": kind,
                                 "result": why, "its": its, "Z_calls": calls[0], "time_ms": 1e3 * dt})
dens = pd.DataFrame(rows)
dens.to_csv("results/residual_comparison_density.csv", index=False)

for phase in ("liquid", "vapour"):
    n_exist = len(dens[(dens.residual == "paper") & (dens.phase == phase) & dens.exists])
    n_absent = len(dens[(dens.residual == "paper") & (dens.phase == phase) & ~dens.exists])
    print(f"  {phase}: phase exists at {n_exist} states, absent at {n_absent}")
    print("  residual     ok %   other root   no conv.   its (ok)   Z calls (ok)   ms (ok)   absent -> None")
    for kind in ("paper", "milton", "main"):
        d = dens[(dens.residual == kind) & (dens.phase == phase) & dens.exists]
        a = dens[(dens.residual == kind) & (dens.phase == phase) & ~dens.exists]
        ok = d[d.result == "ok"]
        print(f"  {kind:8s}  {100 * len(ok) / len(d):7.1f}   {(d.result == 'other root').sum():10d}   "
              f"{(d.result == 'no convergence').sum():8d}   {ok.its.mean():8.1f}   {ok.Z_calls.mean():12.1f}   {ok.time_ms.mean():7.2f}"
              f"   {(a.result == 'absent, returned None').sum():5d} / {len(a)}")
    print()

# ================================================== 2. inside saturation
print(f"2. Saturation, {len(COMPOUNDS)} n-alkanes x {N_T + 3} T, P0 = Wagner 25\n")
rows = []
variants = [("paper", False), ("milton", False), ("main", False), ("paper", True)]
for name in COMPOUNDS:
    c = pcsaft.get(name)
    Tc = c.wagner25[5]
    Ts = list(np.linspace(T_MIN[name], 0.98 * Tc, N_T)) + [0.99 * Tc, 0.995 * Tc, 0.999 * Tc]
    for T in Ts:
        ref = pcsaft.saturation(T, c)
        for kind, analytic in variants:
            calls[0] = 0
            t0 = time.perf_counter()
            P, why, its = newton_saturation(kind, T, c, analytic)
            dt = time.perf_counter() - t0
            if why == "ok" and ref is not None and abs(P / ref[0] - 1) > 1e-5:
                why = "wrong Psat"
            rows.append({"compound": name, "T": T, "Tr": T / Tc,
                         "residual": kind + (", analytic" if analytic else ""),
                         "result": why, "outer_its": its[0], "inner_its": its[1],
                         "Z_calls": calls[0], "time_ms": 1e3 * dt})
    print(f"  done {name}")
sat = pd.DataFrame(rows)
sat.to_csv("results/residual_comparison_saturation.csv", index=False)

print("\n  residual            ok %   density failed   other   outer its   inner its   Z calls   ms   (means over ok)")
for kind in sat.residual.unique():
    d = sat[sat.residual == kind]
    ok = d[d.result == "ok"]
    other = len(d) - len(ok) - (d.result == "density failed").sum()
    print(f"  {kind:17s} {100 * len(ok) / len(d):6.1f}   {(d.result == 'density failed').sum():14d}   {other:5d}   "
          f"{ok.outer_its.mean():9.1f}   {ok.inner_its.mean():9.1f}   {ok.Z_calls.mean():7.0f}   {ok.time_ms.mean():5.1f}")

for kind in sat.residual.unique():
    bad = sat[(sat.residual == kind) & (sat.result != "ok")]
    if len(bad):
        print(f"\n  {kind} fails at:")
        for name, d in bad.groupby("compound", sort=False):
            print(f"    {name:9s} " + ", ".join(f"{tr:.3f}" for tr in d.Tr) + "   (T/Tc)   " + ", ".join(sorted(set(d.result))))

# ================================================== 3. the figure for Milton
c = pcsaft.get(FIG_COMPOUND)
T = FIG_T
P, rl, rv = pcsaft.saturation(T, c)
eL, eV = pcsaft.packing_fraction(T, rl, c), pcsaft.packing_fraction(T, rv, c)
etas = np.linspace(0.002, 0.6, 3000)
Z = np.array([Z_original(T, pcsaft.eta_to_rho(e, T, c), c) for e in etas])
poles = etas[np.where(np.diff(np.sign(Z)) != 0)[0]]          # P_EoS = 0: Milton's poles

fig, ax = plt.subplots(1, 2, figsize=(10, 4))
for a, kind, title in [(ax[0], "paper", r"paper: $P_{EoS}/P - 1$"),
                       (ax[1], "milton", r"Milton: $1 - Z_{calc}/Z_{EoS}$")]:
    f = np.array([residual(kind, e, T, P, c) for e in etas])
    f[np.abs(f) > 1e3] = np.nan
    a.plot(etas, f, color="black", lw=1)
    a.axhline(0, color="gray", lw=0.5)
    a.axvline(eL, color="tab:blue", ls="--", lw=1, label=r"liquid root")
    a.axvline(eV, color="tab:green", ls="--", lw=1, label=r"vapour root")

    path = [0.5]
    eta = 0.5
    for k in range(8):
        h = 1e-6 * eta
        s = (residual(kind, eta + h, T, P, c) - residual(kind, eta - h, T, P, c)) / (2 * h)
        step = -residual(kind, eta, T, P, c) / s
        while not 0 < eta + step < ETA_MAX:
            step = step / 2
        eta = eta + step
        path.append(eta)
    path = np.array(path)
    a.plot(path, [residual(kind, e, T, P, c) for e in path], "o", color="tab:red", ms=4,
           label="Newton from 0.5")
    for k, e0 in enumerate(poles):
        a.axvline(e0, color="tab:red", ls=":", lw=1, label=r"$P_{EoS} = 0$" if k == 0 else None)
    a.set_title(title)
    a.set_xlabel(r"$\eta$")
    if kind == "paper":
        a.set_yscale("asinh")
    else:
        a.set_ylim(-4, 3)
    a.legend(frameon=False, fontsize=8, loc="lower right")
fig.suptitle(f"{FIG_COMPOUND}, T = {T} K, P = Psat = {P/1e5:.3f} bar, liquid asked", fontsize=10)
fig.tight_layout()
fig.savefig("figures/residual_milton.png", dpi=200)
print("\nfigures/residual_milton.png")
