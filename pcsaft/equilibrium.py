"""Density roots, the spinodal, and saturation.

The pressure of an analytical equation of state is not monotonic in density
below the critical temperature: the isotherm rises, turns over at a local
maximum, falls through a mechanically unstable stretch, turns again at a local
minimum, and then climbs steeply up the liquid branch. Everything here is about
handling those three roots without guessing.
"""

import numpy as np
from scipy.optimize import brentq

from .constants import ETA_MAX, kB
from .eos import eta_to_rho, pressure, ln_fugacity_coefficient, segment_diameter


def _P_of_eta(eta, T, comp):
    return pressure(T, eta_to_rho(eta, T, comp), comp)


def spinodal(T, comp, n=600):
    """Locate the two turning points of the isotherm.

    Returns (eta_vapour_side, eta_liquid_side, P_max, P_min), or None above the
    critical temperature where the isotherm is monotonic and no phase split
    exists. Psat must lie between P_min and P_max, which is what makes a
    bracketing solver possible without any initial guess from experiment.
    """
    etas = np.logspace(-10, np.log10(ETA_MAX), n)
    Ps = np.array([_P_of_eta(e, T, comp) for e in etas])
    dP = np.diff(Ps)
    turns = np.where(dP[:-1] * dP[1:] < 0)[0] + 1
    if len(turns) < 2:
        return None
    return etas[turns[0]], etas[turns[1]], Ps[turns[0]], Ps[turns[1]]


def density(T, P, comp, phase="liquid"):
    """Solve P(rho) = P for the density, in molecules/A^3.

    phase is "liquid" (largest root) or "vapour" (smallest). The middle root is
    mechanically unstable, dP/drho < 0, and is never returned.
    """
    f = lambda e: _P_of_eta(e, T, comp) - P

    # A floor derived from the ideal gas always sits below the true vapour root,
    # because attraction makes Z < 1 there. A fixed small number does not: at a
    # triple-point pressure the real root can be at eta ~ 1e-18.
    rho_ideal = P / (kB * T) / 1e30
    eta_floor = 0.1 * np.pi / 6 * rho_ideal * comp.m * segment_diameter(T, comp)**3
    eta_floor = min(max(eta_floor, 1e-20), 1e-8)

    etas = np.concatenate([[eta_floor],
                           np.logspace(np.log10(max(eta_floor, 1e-20)) + 0.1,
                                       np.log10(ETA_MAX), 1500)])
    fs = np.array([f(e) for e in etas])

    roots = []
    for k in range(len(etas) - 1):
        if fs[k] * fs[k + 1] < 0:
            roots.append(brentq(f, etas[k], etas[k + 1], xtol=1e-18))
    if not roots:
        return None

    eta = roots[-1] if phase.startswith("l") else roots[0]
    return eta_to_rho(eta, T, comp)


def _both_roots(T, P, comp, eta_v, eta_l):
    """Vapour and liquid roots at P, bracketed by the spinodal."""
    f = lambda e: _P_of_eta(e, T, comp) - P
    rho_ideal = P / (kB * T) / 1e30
    eta_floor = 0.1 * np.pi / 6 * rho_ideal * comp.m * segment_diameter(T, comp)**3
    ev = brentq(f, eta_floor, eta_v, xtol=1e-18)
    el = brentq(f, eta_l, ETA_MAX, xtol=1e-16)
    return eta_to_rho(ev, T, comp), eta_to_rho(el, T, comp)


def _delta_ln_phi(P, T, comp, eta_v, eta_l):
    """ln phi_liquid - ln phi_vapour. Zero at saturation.

    Far below Psat the liquid branch is almost vertical, Z loses its
    significant figures and ln Z can return NaN. The bracket search below
    expects that and filters on isfinite, so the warning is silenced here
    rather than printed for every probe.
    """
    with np.errstate(divide="ignore", invalid="ignore"):
        rho_v, rho_l = _both_roots(T, P, comp, eta_v, eta_l)
        return (ln_fugacity_coefficient(T, rho_l, comp)
                - ln_fugacity_coefficient(T, rho_v, comp))


def saturation(T, comp):
    """Saturation pressure and the two coexisting densities at T.

    Returns (P_sat in Pa, rho_liquid, rho_vapour) in molecules/A^3, or None
    above the critical temperature.

    No experimental starting guess is used: the bracket comes from the
    spinodal, which the equation of state supplies on its own.
    """
    sp = spinodal(T, comp)
    if sp is None:
        return None
    eta_v, eta_l, P_max, P_min = sp

    hi = P_max * 0.999999
    floor = max(P_min, 0.0)

    # Walk up from the bottom of the allowed window to the first pressure where
    # the liquid is the metastable phase. Starting far below Psat would put the
    # solver on the near-vertical part of the liquid branch, where Z loses all
    # its significant figures and ln Z returns NaN.
    lo = None
    for frac in np.logspace(-12, np.log10(0.5), 60):
        trial = floor + (hi - floor) * frac
        val = _delta_ln_phi(trial, T, comp, eta_v, eta_l)
        if np.isfinite(val) and val > 0:
            lo = trial
            break
    if lo is None:
        return None

    P_sat = brentq(lambda p: _delta_ln_phi(p, T, comp, eta_v, eta_l),
                   lo, hi, xtol=1e-12, rtol=1e-14)
    rho_v, rho_l = _both_roots(T, P_sat, comp, eta_v, eta_l)
    return P_sat, rho_l, rho_v


def saturation_curve(T_array, comp):
    """Saturation at many temperatures.

    Returns three arrays (P_sat in Pa, rho_liquid, rho_vapour), with NaN where
    the solver found no phase split.
    """
    P = np.full(len(T_array), np.nan)
    rl = np.full(len(T_array), np.nan)
    rv = np.full(len(T_array), np.nan)
    for i, T in enumerate(T_array):
        r = saturation(T, comp)
        if r is not None:
            P[i], rl[i], rv[i] = r
    return P, rl, rv


def critical_point(comp, T_lo=100.0, T_hi=1500.0):
    """Critical temperature and pressure, found where the spinodal closes.

    PC-SAFT is a mean-field theory: it ignores the long-range fluctuations that
    dominate near the critical point, so this overpredicts the real Tc,
    typically by several kelvin. Report it, do not hide it.
    """
    f = lambda T: 0 if spinodal(T, comp) is None else 1
    lo, hi = T_lo, T_hi
    if f(lo) == 0:
        return None
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if f(mid):
            lo = mid
        else:
            hi = mid
    Tc = 0.5 * (lo + hi)
    sp = spinodal(lo, comp)
    Pc = 0.5 * (sp[2] + sp[3]) if sp else np.nan
    return Tc, Pc
