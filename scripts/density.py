import warnings
warnings.filterwarnings("ignore")

import pcsaft

COMPOUND = "ethane"
T = 273.15      # K
P = 1.0e5       # Pa

c = pcsaft.get(COMPOUND)
print(f"{COMPOUND}   T = {T} K   P = {P/1e5:g} bar\n")
print("phase      eta          rho [kg/m3]     Z")

for phase in ("liquid", "vapour"):
    rho = pcsaft.density(T, P, c, phase)                 # molecules / A^3
    if rho is None:
        print(f"{phase:8s} no {phase} root at this (T, P)")
        continue
    eta = pcsaft.packing_fraction(T, rho, c)
    Z = pcsaft.compressibility(T, rho, c)
    print(f"{phase:8s} {eta:.8f}   {c.mass_density(rho):11.4f}    {Z:.6f}")
