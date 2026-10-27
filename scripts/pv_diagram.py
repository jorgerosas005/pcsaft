import warnings
warnings.filterwarnings("ignore")

import os
import numpy as np
import matplotlib.pyplot as plt
import pcsaft

COMPOUND = "hexane"
T_START = 360.0                         # K, bottom of the dome
ISOTHERMS = [420.0, 470.0, 500.0, 560.0]  # K; the critical isotherm is added automatically

c = pcsaft.get(COMPOUND)


def molar_volume(rho):
    return 6.02214076e23 / (rho * 1e24)       # molecules/A^3 -> cm3/mol


# ---- the dome: march up in T, each Psat starts from the previous one.
# Wagner stops at its own Tc; marching gets to the top of the PC-SAFT dome.
T_dome, P_dome, vL, vV = [], [], [], []
T, dT, P0 = T_START, 2.0, None
while dT > 0.005:
    sat = pcsaft.saturation(T, c, P0=P0)
    if sat is None:                     # past the top: go back, take a smaller step
        T = T - dT
        dT = dT / 2
        T = T + dT
        continue
    T_dome.append(T)
    P_dome.append(sat[0] / 1e5)
    vL.append(molar_volume(sat[1]))
    vV.append(molar_volume(sat[2]))
    P0 = sat[0]
    T = T + dT

Tc, Pc = T_dome[-1], P_dome[-1]
print(f"{COMPOUND}: top of the PC-SAFT dome  T = {Tc:.2f} K   P = {Pc:.2f} bar   "
      f"v_L = {vL[-1]:.0f}, v_V = {vV[-1]:.0f} cm3/mol")

# liquid branch going up, then vapour branch coming down: one closed curve
v_curve = vL + vV[::-1]
P_curve = P_dome + P_dome[::-1]
plt.plot(v_curve, P_curve, color="black", lw=2, label="saturation dome")

# ---- isotherms
etas = np.linspace(0.0005, 0.6, 4000)
all_T = sorted(ISOTHERMS + [Tc])
colors = plt.cm.plasma(np.linspace(0, 0.85, len(all_T)))
for T, col in zip(all_T, colors):
    rho = pcsaft.eta_to_rho(etas, T, c)
    P = np.array([pcsaft.pressure(T, r, c) for r in rho]) / 1e5
    label = f"{T:.1f} K (≈ Tc)" if T == Tc else f"{T:.0f} K"
    plt.plot(molar_volume(rho), P, color=col, lw=1, label=label)
    if T < Tc:                          # tie line at Psat
        sat = pcsaft.saturation(T, c, P0=np.interp(T, T_dome, P_dome) * 1e5)
        plt.plot([molar_volume(sat[1]), molar_volume(sat[2])], [sat[0] / 1e5] * 2,
                 "--", color=col, lw=1)

plt.plot(vL[-1], Pc, "o", color="black")
plt.axhline(0, color="gray", lw=0.5)
plt.xscale("log")
plt.xlim(100, 3e4)
plt.ylim(-20, 60)
plt.xlabel("v [cm$^3$/mol]")
plt.ylabel("P [bar]")
plt.title(f"{COMPOUND}, PC-SAFT")
plt.legend(frameon=False, fontsize=8)
fig_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "figures")
os.makedirs(fig_dir, exist_ok=True)
file = os.path.join(fig_dir, f"pv_diagram_{COMPOUND}.png")
plt.savefig(file, dpi=200)
print("saved", os.path.normpath(file))
plt.show()
