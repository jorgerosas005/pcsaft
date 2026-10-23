"""Aspen Plus vs this package where the liquid root does not exist.

    py scripts\\04_aspen_pseudoroot.py

Ethane, PC-SAFT, same parameters in both (m = 1.6069, sigma = 3.5206 A,
eps/k = 191.42 K), T = 300 K. The Aspen numbers are typed in from a Pure
analysis (Valid phases Liquid-Only / Vapor-Only), Aspen Plus, 22-oct-2026.

Below P_min the isotherm has no liquid root. Aspen reports a pseudo-root
anyway; this package returns the only root there is, the vapour.
"""

import warnings
warnings.filterwarnings("ignore")

import numpy as np

import pcsaft
from pcsaft.plotting import figure, finish, save

T = 300.0
c = pcsaft.get("ethane")

# Aspen Plus, 300 K:  P [bar], vapour rho, liquid rho [kg/m3]
aspen = np.array([
    [1.0,   1.21458, 155.976],
    [10.0, 13.0753,  160.474],
    [20.0, 28.9233,  170.756],
    [30.0, 49.6993,  194.928],
    [38.0, 74.3348,  244.841],
    [39.0, 78.5076,  255.734],
    [39.3, 79.8411,  266.944],
    [40.0, 83.1271,  279.195],
])
P_a, rhoV_a, rhoL_a = aspen[:, 0], aspen[:, 1], aspen[:, 2]

# the isotherm
eta_v, eta_l, P_max, P_min = pcsaft.spinodal(T, c, n=4000)
eta = np.linspace(1e-4, 0.24, 1500)
rho = np.array([c.mass_density(pcsaft.eta_to_rho(e, T, c)) for e in eta])
P = np.array([pcsaft.pressure(T, pcsaft.eta_to_rho(e, T, c), c) for e in eta]) / 1e5
rho_sv = c.mass_density(pcsaft.eta_to_rho(eta_v, T, c))
rho_sl = c.mass_density(pcsaft.eta_to_rho(eta_l, T, c))
unstable = (eta > eta_v) & (eta < eta_l)

# what this package reports as "liquid"
P_grid = np.linspace(1.0, 42.0, 400)
rhoL_ours = [c.mass_density(pcsaft.density(T, p * 1e5, c, "liquid")) for p in P_grid]

# check: PC-SAFT pressure at Aspen's liquid densities
print("P asked [bar]   Aspen rhoL   P_EoS(rhoL) [bar]")
for p, r in zip(P_a, rhoL_a):
    print(f"{p:8.1f}      {r:9.3f}      {pcsaft.pressure(T, c.number_density(r), c)/1e5:8.3f}")
print(f"spinodal: rho = {rho_sv:.1f} and {rho_sl:.1f} kg/m3,  "
      f"P_max = {P_max/1e5:.2f} bar, P_min = {P_min/1e5:.2f} bar")

fig, ax = figure(1, 2, width=4.8, height=4.0)

# (a) where the points sit on the isotherm
ax[0].plot(np.where(unstable, np.nan, rho), P, color="black", label="PC-SAFT, 300 K")
ax[0].plot(np.where(unstable, rho, np.nan), P, color="black", linestyle="--", linewidth=0.9)
ax[0].plot([rho_sv, rho_sl], [P_max / 1e5, P_min / 1e5], "s", color="black",
           markersize=4, label="spinodal")
ax[0].plot(rhoL_a, P_a, "o", color="crimson", markerfacecolor="none", label="Aspen liquid")
ax[0].plot(rhoV_a, P_a, "^", color="tab:blue", markerfacecolor="none", label="Aspen vapour")
ax[0].set_xlim(0, 320)
ax[0].set_ylim(0, 50)
finish(ax[0], r"$\rho$ / kg m$^{-3}$", r"$P$ / bar", legend=True)

# (b) the reported liquid density against P
ax[1].plot(P_grid, rhoL_ours, color="black", label="this package")
ax[1].plot(P_a, rhoL_a, "o", color="crimson", markerfacecolor="none", label="Aspen")
ax[1].axvline(P_min / 1e5, color="gray", linestyle=":", linewidth=1)
ax[1].text(P_min / 1e5 - 0.5, 20, r"$P_{\min}$", ha="right", color="gray")
finish(ax[1], r"$P$ / bar", r"reported $\rho_L$ / kg m$^{-3}$", legend=True)

print(save(fig, "aspen_pseudoroot_ethane_300K"))
