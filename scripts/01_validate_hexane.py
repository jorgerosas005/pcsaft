"""Validate the pure-component implementation against n-hexane.

    python scripts/01_validate_hexane.py

Writes figures/hexane_saturation.{pdf,png} and results/hexane_saturation.csv.
"""

import pcsaft

comp = pcsaft.get("hexane")

# --- 1. is the code correct?  no experimental data involved -----------------

print("Internal consistency: analytical Z vs numerical derivative of a_res")
cons = pcsaft.check_consistency(comp)
print(cons.to_string(index=False, float_format=lambda v: f"{v:.8g}"))
print(f"\n  worst relative difference: {cons['rel_diff'].max():.2e}")

print("\nIdeal-gas limit:", pcsaft.check_ideal_gas_limit(comp))
print("Diameter limits:", pcsaft.check_diameter_limits(comp))

Tc, Pc = pcsaft.critical_point(comp)
print(f"\nCritical point predicted by the model: {Tc:.1f} K, {Pc/1e5:.1f} bar")
print("  (experimental hexane: 507.8 K, 30.3 bar. PC-SAFT is mean-field and")
print("   overpredicts Tc because it ignores critical fluctuations.)")

# --- 2. is the model any good?  reference data required ---------------------

ref = pcsaft.load_nist_saturation("hexane")
print("\nReference data checks:")
for k, v in pcsaft.check_saturation_data(ref).items():
    print(f"  {k:28s} {v}")

cmp = pcsaft.compare_saturation(comp, ref)
print("\n" + pcsaft.report(comp, cmp, T_max=450))

csv = pcsaft.save_results(cmp, "hexane_saturation")
path = pcsaft.plotting.plot_saturation(cmp, comp)
print(f"\nwrote {csv}\n      {path}")
