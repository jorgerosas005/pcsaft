import warnings
warnings.filterwarnings("ignore")

import os
import numpy as np
import matplotlib.pyplot as plt
import pcsaft

COMPOUNDS = ["methane", "ethane", "propane", "butane", "pentane", "hexane",
             "heptane", "octane", "nonane", "decane", "dodecane"]
CARBONS = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12]
TR_MIN, TR_MAX = 0.45, 0.99     # T / Tc (Wagner)
N = 40                          # temperatures per compound

R = 8.314462618                 # J/(mol K)

fig, ax = plt.subplots(2, 2, figsize=(11, 8))
colors = plt.cm.viridis(np.linspace(0, 0.9, len(COMPOUNDS)))

for name, nc, col in zip(COMPOUNDS, CARBONS, colors):
    c = pcsaft.get(name)
    Tc = c.wagner25[5]
    T, P, rho_l, rho_v, dh = [], [], [], [], []
    for t in np.linspace(TR_MIN * Tc, TR_MAX * Tc, N):
        sat = pcsaft.saturation(t, c)
        up = pcsaft.saturation(1.0001 * t, c)
        down = pcsaft.saturation(0.9999 * t, c)
        if sat is None or up is None or down is None:
            continue
        T.append(t)
        P.append(sat[0])
        rho_l.append(c.mass_density(sat[1]))
        rho_v.append(c.mass_density(sat[2]))

        # Clapeyron: dh_vap = R T^2 (Z_V - Z_L) d ln P / dT, derivative by central difference
        Z_l = pcsaft.compressibility(t, sat[1], c)
        Z_v = pcsaft.compressibility(t, sat[2], c)
        dlnP_dT = (np.log(up[0]) - np.log(down[0])) / (0.0002 * t)
        dh.append(R * t**2 * (Z_v - Z_l) * dlnP_dT / 1000)              # kJ/mol
    T, P = np.array(T), np.array(P)

    label = f"C{nc}"
    ax[0, 0].plot(T, P / 1e5, color=col, label=label)
    ax[0, 1].plot(1000 / T, np.log(P / 1e5), color=col)
    ax[1, 0].plot(rho_l, T, color=col)
    ax[1, 0].plot(rho_v, T, color=col)
    ax[1, 1].plot(T, dh, color=col)
    print(f"{name:9s} {len(T):3d} of {N} temperatures solved")

ax[0, 0].set_yscale("log")
ax[0, 0].set_xlabel("T [K]")
ax[0, 0].set_ylabel(r"$P^{sat}$ [bar]")
ax[0, 0].set_title("vapour pressure")
ax[0, 0].legend(frameon=False, fontsize=8, ncol=2)

ax[0, 1].set_xlabel("1000 / T [1/K]")
ax[0, 1].set_ylabel(r"ln($P^{sat}$ / bar)")
ax[0, 1].set_title("Clausius-Clapeyron")

ax[1, 0].set_xlabel(r"$\rho$ [kg/m$^3$]")
ax[1, 0].set_ylabel("T [K]")
ax[1, 0].set_title("coexisting densities")

ax[1, 1].set_xlabel("T [K]")
ax[1, 1].set_ylabel(r"$\Delta h_{vap}$ [kJ/mol]")
ax[1, 1].set_title("enthalpy of vaporisation (Clapeyron)")

fig.suptitle(f"PC-SAFT, n-alkanes, T/Tc = {TR_MIN} to {TR_MAX}")
fig.tight_layout()
fig_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "figures")
os.makedirs(fig_dir, exist_ok=True)
file = os.path.join(fig_dir, "alkanes.png")
fig.savefig(file, dpi=200)
print("saved", os.path.normpath(file))
plt.show()
