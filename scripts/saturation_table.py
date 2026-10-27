import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pcsaft

COMPOUND = "hexane"
TEMPERATURES = np.arange(300.0, 510.0, 20.0)   # K

R = 8.314462618                                # J/(mol K)
c = pcsaft.get(COMPOUND)

print(f"{COMPOUND}, saturation properties (PC-SAFT)\n")
print("  T [K]   Psat [bar]   rho_L [kg/m3]  rho_V [kg/m3]     Z_L        Z_V       phi      dh_vap [kJ/mol]  ds_vap [J/mol K]")

for T in TEMPERATURES:
    sat = pcsaft.saturation(T, c)
    up = pcsaft.saturation(1.0001 * T, c)
    down = pcsaft.saturation(0.9999 * T, c)
    if sat is None or up is None or down is None:
        print(f"{T:7.2f}   no solution")
        continue
    P, rho_l, rho_v = sat
    Z_l = pcsaft.compressibility(T, rho_l, c)
    Z_v = pcsaft.compressibility(T, rho_v, c)
    phi = np.exp(pcsaft.ln_fugacity_coefficient(T, rho_v, c))     # same in both phases

    # Clapeyron: dh_vap = R T^2 (Z_V - Z_L) d ln Psat / dT
    dlnP_dT = (np.log(up[0]) - np.log(down[0])) / (0.0002 * T)
    dh = R * T**2 * (Z_v - Z_l) * dlnP_dT
    ds = dh / T

    print(f"{T:7.2f}   {P/1e5:10.5f}   {c.mass_density(rho_l):13.3f}  {c.mass_density(rho_v):13.4f}"
          f"   {Z_l:.6f}   {Z_v:.6f}   {phi:.6f}   {dh/1000:14.4f}   {ds:16.3f}")
