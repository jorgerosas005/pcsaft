import warnings
warnings.filterwarnings("ignore")

import os
import numpy as np
import matplotlib.pyplot as plt
import pcsaft

COMPOUND = "ethane"
T_MIN = 150.0   # K
N = 30          # number of temperatures, from T_MIN to 0.99 Tc (Wagner)

c = pcsaft.get(COMPOUND)
Tc = c.wagner25[5]
Ts = np.linspace(T_MIN, 0.99 * Tc, N)

print(f"{COMPOUND}   Tc (Wagner) = {Tc} K\n")
print("  T [K]     Psat [bar]    rho_L [kg/m3]   rho_V [kg/m3]")
T_ok, P_ok, rl_ok, rv_ok = [], [], [], []
for T in Ts:
    sat = pcsaft.saturation(T, c)
    if sat is None:
        print(f"{T:8.2f}     no solution")
        continue
    P, rho_l, rho_v = sat
    T_ok.append(T)
    P_ok.append(P / 1e5)
    rl_ok.append(c.mass_density(rho_l))
    rv_ok.append(c.mass_density(rho_v))
    print(f"{T:8.2f}   {P/1e5:12.6f}   {rl_ok[-1]:13.4f}   {rv_ok[-1]:13.4f}")

fig, ax = plt.subplots(1, 2, figsize=(9, 4))
ax[0].plot(T_ok, P_ok, "o-", color="black", ms=3)
ax[0].set_yscale("log")
ax[0].set_xlabel("T [K]")
ax[0].set_ylabel("Psat [bar]")
ax[0].set_title("vapour pressure")

ax[1].plot(rl_ok, T_ok, "o-", color="tab:blue", ms=3, label="liquid")
ax[1].plot(rv_ok, T_ok, "o-", color="tab:green", ms=3, label="vapour")
ax[1].set_xlabel(r"$\rho$ [kg/m$^3$]")
ax[1].set_ylabel("T [K]")
ax[1].set_title("coexisting densities")
ax[1].legend(frameon=False)

fig.suptitle(COMPOUND)
fig.tight_layout()
fig_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "figures")
os.makedirs(fig_dir, exist_ok=True)
file = os.path.join(fig_dir, f"saturation_curve_{COMPOUND}.png")
fig.savefig(file, dpi=200)
print("saved", os.path.normpath(file))
plt.show()
