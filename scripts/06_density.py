"""Liquid and vapour density at a given T and P.

    py scripts\\06_density.py

Change only the three lines below.
"""

import warnings
warnings.filterwarnings("ignore")

import pcsaft

COMPOUND = "ethane"
T = 273.15      # K
P = 1.0e5       # Pa

c = pcsaft.get(COMPOUND)
print(f"{COMPOUND}   T = {T} K   P = {P/1e5:g} bar\n")
print("phase      eta          rho [kg/m3]   rho [mol/L]     Z")

for phase in ("liquid", "vapour"):
    rho = pcsaft.density(T, P, c, phase)                 # molecules / A^3
    eta = pcsaft.packing_fraction(T, rho, c)
    Z = pcsaft.compressibility(T, rho, c)
    mol_L = rho * 1e27 / 6.02214076e23                   # A^-3 -> mol/L
    print(f"{phase:8s} {eta:.8f}   {c.mass_density(rho):11.4f}   {mol_L:11.5f}   {Z:.6f}")

rho_l = pcsaft.density(T, P, c, "liquid")
rho_v = pcsaft.density(T, P, c, "vapour")
if rho_l == rho_v:
    print("\nOnly one root at this (T, P): there is no separate liquid and vapour.")
