"""How much does the choice of density residual matter?

    py scripts\\07_residual_analysis.py

Every density solve starts from the paper's values, eta = 0.5 (liquid) or
1e-10 (vapour), and uses the same Newton iteration. Only the residual changes:

    paper   :  P_EoS / P - 1          (Pcalc = Psys, written dimensionless)
    ratio   :  1 - Z_PV / Z_EoS       (Milton / Diego)
    diff    :  Z_EoS - Z_PV           (same Z's, no division)

The only safeguard is the obvious one: a step that leaves 0 < eta < 0.7405 is
halved until it lands inside.

Psat uses successive substitution, P <- P * phi_L / phi_V, from P0 = 0.8 Psat.
That P0 is generous (from Psat/2 every residual fails at T/Tc = 0.95,
because there is no liquid root at P0): it is the same for the three residuals, so whatever
differs comes from the density solver alone.

Cost is counted in EoS evaluations (one evaluation = one Z_EoS), which is what
time is proportional to.
"""

import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pcsaft
from pcsaft.eos import Z_hc, Z_disp
from pcsaft.constants import kB
from pcsaft.plotting import figure, finish, save

COMPOUNDS = ["methane", "ethane", "propane", "butane", "hexane", "octane", "decane"]
TR = np.linspace(0.45, 0.95, 11)          # T / Tc of the EoS
RESIDUALS = ["paper", "ratio", "diff"]

count = [0]


def residual(kind, eta, T, P, c):
    count[0] += 1
    rho = pcsaft.eta_to_rho(eta, T, c)
    Z_eos = 1 + Z_hc(T, rho, c) + Z_disp(T, rho, c)
    Z_pv = P / (rho * 1e30 * kB * T)
    if kind == "paper":
        return Z_eos / Z_pv - 1          # = P_EoS / P - 1
    if kind == "ratio":
        return 1 - Z_pv / Z_eos
    return Z_eos - Z_pv


def newton(kind, T, P, c, phase):
    eta = 0.5 if phase == "liquid" else 1e-10
    for k in range(200):
        h = 1e-6 * eta
        r = residual(kind, eta, T, P, c)
        slope = (residual(kind, eta + h, T, P, c) - residual(kind, eta - h, T, P, c)) / (2 * h)
        if slope == 0 or not np.isfinite(slope):
            return None
        step = -r / slope
        while not 0 < eta + step < 0.7405:
            step /= 2
        if abs(step) < 1e-12 * eta:
            return eta + step
        eta += step
    return None


def psat_ss(kind, T, c, P0):
    P = P0
    for k in range(300):
        eL = newton(kind, T, P, c, "liquid")
        eV = newton(kind, T, P, c, "vapour")
        if eL is None or eV is None or abs(eL - eV) < 1e-8:
            return None
        lnL = pcsaft.ln_fugacity_coefficient(T, pcsaft.eta_to_rho(eL, T, c), c)
        lnV = pcsaft.ln_fugacity_coefficient(T, pcsaft.eta_to_rho(eV, T, c), c)
        P_new = P * np.exp(lnL - lnV)
        if abs(P_new / P - 1) < 1e-8:      # ln Z_L limits the digits at low P
            return P_new
        P = P_new
    return None


rows = []
for name in COMPOUNDS:
    c = pcsaft.get(name)
    Tc, _ = pcsaft.critical_point(c)
    for tr in TR:
        T = tr * Tc
        Ps, rl, rv = pcsaft.saturation(T, c)
        ref = {"liquid": pcsaft.packing_fraction(T, rl, c),
               "vapour": pcsaft.packing_fraction(T, rv, c)}
        for kind in RESIDUALS:
            row = {"compound": name, "m": c.m, "Tr": tr, "residual": kind}
            for phase in ("liquid", "vapour"):
                count[0] = 0
                e = newton(kind, T, Ps, c, phase)
                ok = e is not None and abs(e / ref[phase] - 1) < 1e-6
                row[f"{phase}_ok"] = ok
                row[f"{phase}_evals"] = count[0]
            count[0] = 0
            P = psat_ss(kind, T, c, 0.8 * Ps)
            row["psat_ok"] = P is not None and abs(P / Ps - 1) < 1e-6
            row["psat_evals"] = count[0]
            rows.append(row)
    print(f"done {name}")

df = pd.DataFrame(rows)
pcsaft.save_results(df, "residual_analysis")

# ---------------------------------------------------------------- summary
print("\nShare of cases solved correctly, and mean EoS evaluations when solved\n")
summary = []
for kind in RESIDUALS:
    d = df[df.residual == kind]
    line = {"residual": kind}
    for what in ("liquid", "vapour", "psat"):
        ok = d[f"{what}_ok"]
        line[f"{what} ok %"] = 100 * ok.mean()
        line[f"{what} evals"] = d.loc[ok, f"{what}_evals"].mean() if ok.any() else np.nan
    summary.append(line)
print(pd.DataFrame(summary).to_string(index=False, float_format=lambda v: f"{v:.0f}"))

# ---------------------------------------------------------------- figure
# rows: liquid density, vapour density, Psat.  columns: the three residuals.
# colour: chain length (one sequential colour map, light = short).
colors = plt.cm.Blues(np.linspace(0.35, 0.95, len(COMPOUNDS)))
labels = {"paper": r"paper: $P_{EoS}/P - 1$", "ratio": r"Milton: $1 - Z_{PV}/Z_{EoS}$",
          "diff": r"$Z_{EoS} - Z_{PV}$"}
what = [("liquid", "liquid density"), ("vapour", "vapour density"), ("psat", r"$P^{sat}$")]

fig, ax = figure(3, 3, width=4.0, height=3.0, sharex=True, sharey="row")
for j, kind in enumerate(RESIDUALS):
    for i, (key, title) in enumerate(what):
        a = ax[i][j]
        for name, col in zip(COMPOUNDS, colors):
            d = df[(df.residual == kind) & (df.compound == name)]
            good = d[d[f"{key}_ok"]]
            bad = d[~d[f"{key}_ok"]]
            a.plot(good["Tr"], good[f"{key}_evals"], marker="o", markersize=3,
                   color=col, label=name)
            a.plot(bad["Tr"], [0.04] * len(bad), "x", color="crimson", markersize=5,
                   transform=a.get_xaxis_transform())
        a.set_yscale("log")
        ticks = {"liquid": [15, 20, 30, 40], "vapour": [10, 20, 50, 100], "psat": [100, 300, 1000, 3000]}[key]
        a.set_yticks(ticks)
        a.set_yticklabels([str(t) for t in ticks])
        a.minorticks_off()
        if i == 0:
            a.set_title(labels[kind], fontsize=10)
        finish(a, r"$T/T_c$" if i == 2 else None,
               f"{title}\nEoS evaluations" if j == 0 else None,
               legend=(i == 0 and j == 0))
ax[0][0].legend(frameon=False, fontsize=7)
print(save(fig, "residual_analysis"))
print("red x along the bottom = failed or converged to the wrong root")
