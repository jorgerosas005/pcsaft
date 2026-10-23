"""Checking the code and checking the model.

Two different questions live here, and keeping them apart matters:

  1. Is this a correct implementation of PC-SAFT?
     Answered without any experimental data, by internal consistency:
     the analytical Z must equal the numerical derivative of a_res.
     `check_consistency` does this.

  2. Is PC-SAFT a good model of reality?
     Answered against reference data. `compare_saturation` does this.

A wrong answer to (2) with a right answer to (1) is the model's fault.
Both wrong means the code is at fault. Never mix them up.
"""

import numpy as np
import pandas as pd

from .eos import (a_res, compressibility, compressibility_numeric,
                  eta_to_rho, segment_diameter)
from .equilibrium import density, saturation_curve
from .parameters import get


# ------------------------------------------------- (1) is the code correct?

def check_consistency(comp, T_list=(200.0, 300.0, 450.0), eta_list=(0.05, 0.2, 0.35, 0.5)):
    """Compare the analytical Z with a central difference of a_res.

    a_res and Z are coded from different equations in the appendix. If they
    agree to 6+ digits, a typo in either is very unlikely. Needs no data.
    """
    rows = []
    for T in T_list:
        for eta in eta_list:
            rho = eta_to_rho(eta, T, comp)
            Za = compressibility(T, rho, comp)
            Zn = compressibility_numeric(T, rho, comp)
            rows.append({"T": T, "eta": eta, "Z_analytic": Za, "Z_numeric": Zn,
                         "rel_diff": abs(Za - Zn) / max(abs(Za), 1e-30)})
    return pd.DataFrame(rows)


def check_ideal_gas_limit(comp, T=300.0, P=1.0):
    """At vanishing pressure Z must go to 1 and a_res to 0."""
    rho = density(T, P, comp, "vapour")
    return {"Z": compressibility(T, rho, comp), "a_res": a_res(T, rho, comp)}


def check_diameter_limits(comp):
    """d must stay below sigma, rise as T falls, and bracket [0.88 sigma, sigma]."""
    lo = segment_diameter(1.0, comp)
    hi = segment_diameter(1e8, comp)
    return {"d(T->0)/sigma": lo / comp.sigma, "d(T->inf)/sigma": hi / comp.sigma}


# ----------------------------------------------- (2) is the model any good?

def compare_saturation(comp, ref, T_max_subcritical=None):
    """PC-SAFT saturation against a reference table.

    `ref` is a DataFrame from data.load_nist_saturation. Returns a DataFrame
    with both sets of values and the percentage deviations.
    """
    T = ref["T"].values
    P_calc, rl_calc, rv_calc = saturation_curve(T, comp)

    out = pd.DataFrame({
        "T": T,
        "P_ref": ref["P"].values,
        "P_calc": P_calc,
        "rhoL_ref": ref["rho_L"].values,
        "rhoL_calc": comp.mass_density(rl_calc),
        "rhoV_ref": ref["rho_V"].values,
        "rhoV_calc": comp.mass_density(rv_calc),
    })
    for k in ["P", "rhoL", "rhoV"]:
        out[f"dev_{k}"] = 100 * (out[f"{k}_calc"] - out[f"{k}_ref"]) / out[f"{k}_ref"]
    return out


def aad(comparison, T_max=None):
    """Average absolute deviation, in percent.

    Pass T_max to exclude the critical region. PC-SAFT is mean-field and fails
    systematically near Tc, so an AAD over the whole range mixes a real model
    error with a known structural limitation. Report both.
    """
    df = comparison if T_max is None else comparison[comparison["T"] <= T_max]
    return {
        "n": int(df["dev_P"].notna().sum()),
        "P_sat": float(np.nanmean(np.abs(df["dev_P"]))),
        "rho_L": float(np.nanmean(np.abs(df["dev_rhoL"]))),
        "rho_V": float(np.nanmean(np.abs(df["dev_rhoV"]))),
    }


def report(comp, comparison, T_max=None):
    """A short text block with both AAD windows and the published values."""
    full = aad(comparison)
    lines = [f"{comp.name}   m={comp.m}  sigma={comp.sigma} A  eps/k={comp.eps_k} K",
             f"source: {comp.source}"]
    if T_max is not None:
        sub = aad(comparison, T_max)
        lines.append(f"\n  AAD% up to {T_max:.0f} K   ({sub['n']} points)")
        lines.append(f"    Psat  {sub['P_sat']:6.2f}   rho_L {sub['rho_L']:6.2f}   rho_V {sub['rho_V']:6.2f}")
    lines.append(f"\n  AAD% full range   ({full['n']} points)")
    lines.append(f"    Psat  {full['P_sat']:6.2f}   rho_L {full['rho_L']:6.2f}   rho_V {full['rho_V']:6.2f}")
    if comp.aad_psat is not None:
        lines.append(f"\n  published (Gross & Sadowski 2001)")
        lines.append(f"    Psat  {comp.aad_psat:6.2f}   volume {comp.aad_v:6.2f}")
    return "\n".join(lines)
