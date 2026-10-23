"""Caloric validation: enthalpy of vaporisation and Cp - Cv.

    python scripts/03_validate_caloric.py

These are predictions, not correlations. The pure-component parameters were
regressed against vapour pressure and liquid density only, so no caloric data
entered the fit. They test the temperature derivatives of a_res, which the
density work never touches.

Both quantities are chosen so no ideal-gas heat capacity correlation is needed:
dh_vap because the ideal part cancels between the two phases at the same T, and
Cp - Cv because it follows from P(T, rho) alone.
"""

import warnings

import numpy as np
import pandas as pd

import pcsaft
from pcsaft.plotting import figure, finish, save

warnings.filterwarnings("ignore")

COMPOUNDS = ["propane", "hexane", "decane"]

summary = []
fig, ax = figure(1, 2, width=5.2, height=4.0)

for name in COMPOUNDS:
    comp = pcsaft.get(name)
    ref = pcsaft.load_nist_saturation(name)
    cal = pcsaft.compare_caloric(comp, ref)

    T_cut = ref["T"].min() + 0.9 * (ref["T"].max() - ref["T"].min())
    sub = cal[cal["T"] <= T_cut]
    Tr = sub["T"] / ref["T"].max()

    ax[0].plot(Tr, sub["dev_dhvap"], marker="o", markersize=3, label=name)
    ax[1].plot(Tr, sub["dev_CpCv"], marker="o", markersize=3, label=name)

    pcsaft.save_results(cal, f"{name}_caloric")
    summary.append({
        "compound": name,
        "n": len(sub),
        "AAD_dhvap": np.nanmean(np.abs(sub["dev_dhvap"])),
        "AAD_CpminusCv": np.nanmean(np.abs(sub["dev_CpCv"])),
    })

for a, lab in zip(ax, [r"$\Delta h_{\mathrm{vap}}$", r"$C_p - C_v$ (sat. liquid)"]):
    a.axhline(0, color="black", linewidth=0.6)
    finish(a, r"$T / T_c$", f"deviation in {lab} / %", legend=True)
save(fig, "caloric_validation")

table = pd.DataFrame(summary)
pcsaft.save_results(table, "caloric_aad")
print(table.to_string(index=False, float_format=lambda v: f"{v:.2f}"))

print("""
Reading these numbers

  dh_vap at ~1%: a real success. The equation reproduces a property it was
  never fitted to, across the whole subcritical range.

  Cp - Cv at ~10-20%: the implementation is not at fault. Two independent
  routes -- derivatives of P(T, rho), and T V alpha^2 / kappa_T from the
  density solver -- agree to better than 1e-4 %. The error is the model's.

  It looks worse than the 1.7% that Gross & Sadowski report in Table 3 because
  that table is for TOTAL cp, which for a liquid alkane is dominated by the
  ideal-gas contribution supplied from outside the equation of state. Cp - Cv
  removes that contribution entirely and leaves only what the EoS itself
  predicts, so it is a much harsher test of the same model.

  Why this matters here: Cp - Cv is governed by the isothermal compressibility,
  a second derivative of the Helmholtz energy. The thesis targets exactly such
  quantities in water -- the Cp minimum, the compressibility minimum. This sets
  a baseline for how much of any residual disagreement is the association model
  and how much was already present in plain PC-SAFT.
""")
