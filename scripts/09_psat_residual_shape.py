"""What the Psat residual looks like.

    py scripts\\09_psat_residual_shape.py

    g(P; T) = ln phi(eta_L(P), T) - ln phi(eta_V(P), T)

g is a function of P. The two eta's are not unknowns of g: for each P they
are solved first (the density problem), then plugged in. So every P gives one
pair (eta_V, eta_L), and g can be drawn against P, or against either eta.
"""

import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pcsaft
from pcsaft.equilibrium import _both_roots
from pcsaft.plotting import figure, finish, save

COMPOUND = "ethane"
T = 250.0      # K

c = pcsaft.get(COMPOUND)
eta_sv, eta_sl, P_max, P_min = pcsaft.spinodal(T, c, n=3000)
P_sat, rho_l, rho_v = pcsaft.saturation(T, c)
eta_L_sat = pcsaft.packing_fraction(T, rho_l, c)
eta_V_sat = pcsaft.packing_fraction(T, rho_v, c)

# the isotherm
eta = np.logspace(-4, np.log10(0.45), 1500)
P_iso = np.array([pcsaft.pressure(T, pcsaft.eta_to_rho(e, T, c), c) for e in eta]) / 1e5

# g over the window where both roots exist
P = np.linspace(max(P_min, 0) + 0.2e5, 0.999 * P_max, 300)
eV, eL, g = [], [], []
for p in P:
    rv, rl = _both_roots(T, p, c, eta_sv, eta_sl)
    eV.append(pcsaft.packing_fraction(T, rv, c))
    eL.append(pcsaft.packing_fraction(T, rl, c))
    g.append(pcsaft.ln_fugacity_coefficient(T, rl, c) - pcsaft.ln_fugacity_coefficient(T, rv, c))

print(f"{COMPOUND}, T = {T} K")
print(f"  window where both roots exist: {max(P_min,0)/1e5:.2f} to {P_max/1e5:.2f} bar")
print(f"  Psat = {P_sat/1e5:.4f} bar,  eta_V = {eta_V_sat:.5f},  eta_L = {eta_L_sat:.5f}")

fig, ax = figure(1, 3, width=4.2, height=3.8)

# (a) the isotherm, the window, and Psat
ax[0].semilogx(eta, P_iso, color="black")
ax[0].axhspan(max(P_min, 0) / 1e5, P_max / 1e5, color="gray", alpha=0.15, label="both roots exist")
ax[0].axhline(P_sat / 1e5, color="crimson", linestyle="--", linewidth=1, label=r"$P^{sat}$")
ax[0].set_ylim(-20, 40)
finish(ax[0], r"$\eta$", r"$P$ / bar", legend=True)

# (b) g against P
ax[1].plot(np.array(P) / 1e5, g, color="black")
ax[1].axhline(0, color="gray", linewidth=0.6)
ax[1].axvline(P_sat / 1e5, color="crimson", linestyle="--", linewidth=1)
finish(ax[1], r"$P$ / bar", r"$g = \ln\varphi_L - \ln\varphi_V$")

# (c) the same g, drawn against each eta
ax[2].semilogx(eV, g, color="tab:blue", label=r"against $\eta_V(P)$")
ax[2].semilogx(eL, g, color="tab:orange", label=r"against $\eta_L(P)$")
ax[2].axhline(0, color="gray", linewidth=0.6)
ax[2].plot([eta_V_sat, eta_L_sat], [0, 0], "o", color="crimson", markersize=4)
finish(ax[2], r"$\eta$", r"$g$", legend=True)

print(save(fig, f"psat_residual_{COMPOUND}_{T:.0f}K"))
