import warnings
warnings.filterwarnings("ignore")

import os
import numpy as np
import matplotlib.pyplot as plt
import pcsaft

COMPOUND = "ethane"
T = 250.0       # K

c = pcsaft.get(COMPOUND)
etas = np.linspace(0.001, 0.36, 2000)
P = np.array([pcsaft.pressure(T, pcsaft.eta_to_rho(e, T, c), c) for e in etas]) / 1e5   # bar

# turning points of the isotherm (the spinodal), where dP/deta changes sign
turns = np.where(np.diff(np.sign(np.diff(P))) != 0)[0] + 1
print(f"{COMPOUND}   T = {T} K\n")
for i in turns:
    print(f"turning point   eta = {etas[i]:.4f}   P = {P[i]:10.3f} bar")

sat = pcsaft.saturation(T, c)
if sat is not None:
    Ps, rho_l, rho_v = sat
    eL = pcsaft.packing_fraction(T, rho_l, c)
    eV = pcsaft.packing_fraction(T, rho_v, c)
    print(f"\nPsat = {Ps/1e5:.4f} bar   eta_V = {eV:.5f}   eta_L = {eL:.5f}")

plt.plot(etas, P, color="black")
plt.axhline(0, color="gray", lw=0.5)
if sat is not None:
    plt.axhline(Ps / 1e5, color="tab:red", ls="--", lw=1, label=f"Psat = {Ps/1e5:.2f} bar")
    plt.plot([eV, eL], [Ps / 1e5, Ps / 1e5], "o", color="tab:red")
plt.plot(etas[turns], P[turns], "s", color="tab:blue", ms=5, label="turning points (spinodal)")
plt.ylim(min(P[turns].min() * 1.3, -5), max(3 * Ps / 1e5 if sat else 50, 20))
plt.xlabel(r"$\eta$")
plt.ylabel("P [bar]")
plt.title(f"{COMPOUND}, T = {T} K")
plt.legend(frameon=False)
os.makedirs("figures", exist_ok=True)
plt.savefig("figures/isotherm.png", dpi=200)
plt.show()
