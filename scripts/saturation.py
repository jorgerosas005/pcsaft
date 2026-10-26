import warnings
warnings.filterwarnings("ignore")

import pcsaft

COMPOUND = "ethane"
T = 250.0       # K

c = pcsaft.get(COMPOUND)
P0 = pcsaft.wagner25(T, c)
sat = pcsaft.saturation(T, c)

print(f"{COMPOUND}   T = {T} K\n")
if sat is None:
    print("No saturation at this T (above the Wagner Tc, or Newton failed).")
    raise SystemExit

P, rho_l, rho_v = sat
print(f"P0 Wagner = {P0/1e5:.6f} bar")
print(f"Psat      = {P/1e5:.6f} bar\n")
print("phase      eta          rho [kg/m3]     Z")

for phase, rho in (("liquid", rho_l), ("vapour", rho_v)):
    eta = pcsaft.packing_fraction(T, rho, c)
    Z = pcsaft.compressibility(T, rho, c)
    print(f"{phase:8s} {eta:.8f}   {c.mass_density(rho):11.4f}    {Z:.6f}")
