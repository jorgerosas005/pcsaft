"""Our residual vs Diego's, for ethane at the state used in the comparison.

    py exercises\\plot_residuals.py

    f(eta) = P_EoS - P              (ours, in bar)
    r(eta) = 1 - Zpv / Zeos         (Diego's)  =  1 - P / P_EoS

Same roots. r also blows up wherever P_EoS = 0.
"""

import warnings
warnings.filterwarnings("ignore")

import numpy as np
from scipy.optimize import brentq

import pcsaft
from pcsaft.equilibrium import _P_of_eta
from pcsaft.plotting import figure, finish, save

COMPOUND = "ethane"
T = 273.15      # K
P = 1.0e5       # Pa

c = pcsaft.get(COMPOUND)

eta = np.logspace(-6, np.log10(0.7), 4000)
PE = np.array([_P_of_eta(e, T, c) for e in eta])
f = (PE - P) / 1e5
r = 1 - P / PE


def zeros(g):
    """Every sign change of g on the grid, refined with brentq."""
    vals = np.array([g(e) for e in eta])
    out = []
    for k in range(len(eta) - 1):
        if vals[k] * vals[k + 1] < 0:
            out.append(brentq(g, eta[k], eta[k + 1], xtol=1e-16))
    return out


roots = zeros(lambda e: _P_of_eta(e, T, c) - P)
poles = zeros(lambda e: _P_of_eta(e, T, c))

print(f"{COMPOUND} at T = {T} K, P = {P/1e5} bar")
print("roots (f = 0 and r = 0):", [f"{e:.6f}" for e in roots])
print("poles of r (P_EoS = 0): ", [f"{e:.6f}" for e in poles])

r_signs = np.sum(r[:-1] * r[1:] < 0)
print(f"sign changes of r on the grid: {r_signs}  (3 roots + {len(poles)} poles)")

fig, ax = figure(1, 3, width=4.4, height=3.8)

# (a) ours
ax[0].semilogx(eta, f, color="black")
ax[0].axhline(0, color="gray", linewidth=0.6)
ax[0].plot(roots, [0] * len(roots), "o", color="crimson", label="roots")
ax[0].set_ylim(-25, 40)
finish(ax[0], r"$\eta$", r"$f = P_{EoS} - P$  / bar", legend=True)

# (b) Diego's
ax[1].semilogx(eta, r, color="black")
ax[1].axhline(0, color="gray", linewidth=0.6)
for p in poles:
    ax[1].axvline(p, color="tab:blue", linestyle="--", linewidth=1)
ax[1].plot(roots, [0] * len(roots), "o", color="crimson", label="roots")
ax[1].plot([], [], "--", color="tab:blue", label=r"poles ($P_{EoS}=0$)")
ax[1].set_ylim(-5, 5)
finish(ax[1], r"$\eta$", r"$r = 1 - Z_{PV}/Z_{EoS}$", legend=True)

# (c) zoom on the liquid root
ez = np.linspace(poles[-1] - 0.0015, roots[-1] + 0.0015, 2000)
rz = np.array([1 - P / _P_of_eta(e, T, c) for e in ez])
ax[2].plot(ez, rz, color="black")
ax[2].axhline(0, color="gray", linewidth=0.6)
ax[2].axvline(poles[-1], color="tab:blue", linestyle="--", linewidth=1, label="pole")
ax[2].plot(roots[-1], 0, "o", color="crimson", label="liquid root")
ax[2].set_ylim(-5, 5)
finish(ax[2], r"$\eta$", r"$r$  (zoom)", legend=True)

print(save(fig, "residuals_ethane_273K"))
