"""Validate the whole n-alkane series, which is what the protocol asks for.

    python scripts/02_validate_series.py

Runs on whatever NIST files are present in data/nist/ and says which are
missing. To add one: NIST Webbook -> Fluid Properties -> the compound ->
units K / bar / kg per cubic metre -> "Saturation properties - temperature
increments" -> download as tab-delimited -> data/nist/<name>_sat.txt
"""

import pandas as pd

import pcsaft

rows = []
missing = []

for name in pcsaft.N_ALKANES:
    comp = pcsaft.get(name)
    try:
        ref = pcsaft.load_nist_saturation(name)
    except FileNotFoundError:
        missing.append(name)
        continue

    cmp = pcsaft.compare_saturation(comp, ref)
    pcsaft.save_results(cmp, f"{name}_saturation")
    pcsaft.plotting.plot_saturation(cmp, comp)

    # exclude the top 10% of the temperature range, where every mean-field
    # equation of state fails for reasons that are not the implementation's
    T_cut = ref["T"].min() + 0.9 * (ref["T"].max() - ref["T"].min())
    stats = pcsaft.aad(cmp, T_max=T_cut)

    rows.append({
        "compound": name,
        "m": comp.m,
        "sigma": comp.sigma,
        "eps_k": comp.eps_k,
        "n": stats["n"],
        "AAD_Psat": stats["P_sat"],
        "AAD_rhoL": stats["rho_L"],
        "AAD_Psat_published": comp.aad_psat,
        "AAD_v_published": comp.aad_v,
    })

if rows:
    table = pd.DataFrame(rows)
    path = pcsaft.save_results(table, "alkane_series_aad")
    print(table.to_string(index=False, float_format=lambda v: f"{v:.3f}"))
    print(f"\nwrote {path}")

if missing:
    print("\nNo NIST file yet for: " + ", ".join(missing))
    print("See the docstring of this script for how to download them.")
