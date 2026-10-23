"""Caloric properties: residual enthalpy, entropy, and heat capacities.

Why this matters more than it looks. The pure-component parameters m, sigma and
eps/k were regressed against vapour pressure and liquid density ONLY. No caloric
data entered the fit. So every number in this module is a prediction, not a
correlation, and it tests something the density work cannot: whether the
TEMPERATURE dependence of the Helmholtz energy is right.

That is exactly the axis the thesis lives on. Water's anomalies -- the density
maximum, the Cp minimum, the compressibility minimum -- are temperature-derivative
phenomena. The known failure of the PC-SAFT-TS model of Novak et al. (2024) is
that its Cp minimum does not disappear at high pressure. If the second
temperature derivative of a_res is not trustworthy for a simple alkane, no
conclusion about a Cp minimum in water will mean anything.

Everything here uses numerical temperature derivatives of a_res. That is
deliberate: it works unchanged when the association term is added, with no new
derivations.

Two quantities need no ideal-gas correlation at all and are therefore the
cleanest validation targets:

  * dh_vap = h_V - h_L at one T. Both phases share the same ideal-gas enthalpy,
    so it cancels exactly and the difference is purely residual.
  * Cp - Cv. It follows from P(T, rho) alone. The ideal gas gives exactly R,
    and everything beyond that is the equation of state.
"""

import numpy as np

from .constants import R, kB
from .eos import a_res, compressibility, pressure


# ------------------------------------------------------- numerical derivatives

def _dadT(T, rho, comp, h=None):
    """(d a_res / dT) at constant rho, by central difference."""
    h = h or T * 1e-5
    return (a_res(T + h, rho, comp) - a_res(T - h, rho, comp)) / (2 * h)


def _d2adT2(T, rho, comp, h=None):
    """(d2 a_res / dT2) at constant rho, by central second difference."""
    h = h or T * 1e-4
    return (a_res(T + h, rho, comp) - 2 * a_res(T, rho, comp)
            + a_res(T - h, rho, comp)) / h**2


def _dPdT(T, rho, comp, h=None):
    """(dP/dT) at constant rho, in Pa/K."""
    h = h or T * 1e-5
    return (pressure(T + h, rho, comp) - pressure(T - h, rho, comp)) / (2 * h)


def _dPdrho(T, rho, comp, rel=1e-6):
    """(dP/drho) at constant T, in Pa per (molecule/A^3)."""
    h = rho * rel
    return (pressure(T, rho + h, comp) - pressure(T, rho - h, comp)) / (2 * h)


# -------------------------------------------------------- residual properties

def u_res(T, rho, comp):
    """Residual internal energy, J/mol.   U_res/(NkT) = -T (d a_res/dT)_rho"""
    return -R * T**2 * _dadT(T, rho, comp)


def h_res(T, rho, comp):
    """Residual enthalpy, J/mol.   H_res = U_res + NkT (Z - 1)"""
    Z = compressibility(T, rho, comp)
    return u_res(T, rho, comp) + R * T * (Z - 1)


def s_res(T, rho, comp):
    """Residual entropy at constant volume, J/(mol K)."""
    return -R * (a_res(T, rho, comp) + T * _dadT(T, rho, comp))


def cv_res(T, rho, comp):
    """Residual isochoric heat capacity, J/(mol K).

    Cv_res/(Nk) = -2 T (da/dT) - T^2 (d2a/dT2)
    """
    return -R * (2 * T * _dadT(T, rho, comp) + T**2 * _d2adT2(T, rho, comp))


def cp_minus_cv(T, rho, comp):
    """Cp - Cv, J/(mol K). Depends only on P(T, rho), never on ideal-gas data.

    Cp - Cv = -T (dP/dT)_V^2 / (dP/dV)_T, rewritten in density:
              = T (dP/dT)_rho^2 / (k rho^2 (dP/drho)_T)

    Densities go to molecules/m^3 here so the result comes out dimensionless
    before being multiplied by R. For an ideal gas this returns exactly R.
    """
    rho_si = rho * 1e30                       # molecules/A^3 -> molecules/m^3
    dPdT = _dPdT(T, rho, comp)                # Pa/K
    dPdrho = _dPdrho(T, rho, comp) / 1e30     # Pa per (molecule/m^3)
    return R * T * dPdT**2 / (kB * rho_si**2 * dPdrho)


def cp_res(T, rho, comp):
    """Residual isobaric heat capacity, J/(mol K).

    Cp_res = Cv_res + (Cp - Cv) - R, since the ideal gas already supplies R.
    """
    return cv_res(T, rho, comp) + cp_minus_cv(T, rho, comp) - R


def enthalpy_of_vaporisation(T, comp):
    """dh_vap in J/mol, from the two saturated densities.

    Both phases sit at the same temperature, so the ideal-gas enthalpy is
    identical in each and cancels: the difference of the residual enthalpies IS
    the enthalpy of vaporisation. No ideal-gas heat capacity is needed.
    """
    from .equilibrium import saturation
    r = saturation(T, comp)
    if r is None:
        return np.nan
    _, rho_l, rho_v = r
    return h_res(T, rho_v, comp) - h_res(T, rho_l, comp)


# ------------------------------------------------------------- validation

def compare_caloric(comp, ref):
    """PC-SAFT caloric predictions against a reference saturation table.

    Compares only the two quantities that need no ideal-gas correlation.
    """
    import pandas as pd
    from .equilibrium import saturation

    rows = []
    for _, row in ref.iterrows():
        T = row["T"]
        r = saturation(T, comp)
        if r is None:
            continue
        _, rho_l, rho_v = r
        rec = {"T": T}

        if "dh_vap" in ref.columns:
            rec["dhvap_ref"] = row["dh_vap"] * 1000.0          # kJ/mol -> J/mol
            rec["dhvap_calc"] = h_res(T, rho_v, comp) - h_res(T, rho_l, comp)

        if "Cp_minus_Cv_L" in ref.columns:
            rec["CpCv_ref"] = row["Cp_minus_Cv_L"]
            rec["CpCv_calc"] = cp_minus_cv(T, rho_l, comp)

        rows.append(rec)

    out = pd.DataFrame(rows)
    for k in ("dhvap", "CpCv"):
        if f"{k}_ref" in out:
            out[f"dev_{k}"] = 100 * (out[f"{k}_calc"] - out[f"{k}_ref"]) / out[f"{k}_ref"]
    return out
