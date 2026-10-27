"""Density roots and saturation pressure.

density()     Newton on the residual P_EoS / P - 1, in the packing fraction eta,
              from eta = 0.5 (liquid) or eta = 1e-10 (vapour).

saturation()  Initial estimate P0 from Wagner 25, then Newton on P with the
              residual g = ln phi_L - ln phi_V. Since (d ln phi / dP)_T = (Z - 1) / P,

                  dg/dP = (Z_L - Z_V) / P   ->   P_new = P (1 + g / (Z_V - Z_L))

Wagner 25 needs the constants of each compound (parameters.WAGNER25). Without
them, or above the Wagner Tc, saturation() returns None.
"""

import numpy as np

from .constants import ETA_MAX
from .eos import eta_to_rho, pressure, compressibility, ln_fugacity_coefficient


def density(T, P, comp, phase="liquid", verbose=False):
    """Density at (T, P) in molecules/A^3. phase is "liquid" or "vapour".

    Returns None if Newton does not converge or lands on the unstable branch
    (dP/deta < 0). verbose=True prints every iteration.
    """
    eta = 0.5 if phase.startswith("l") else 1e-10
    if verbose:
        print(f"  density, {phase}: Newton on P_EoS/P - 1 from eta = {eta:g}")
        print("    k   eta               P_EoS/P - 1")

    for k in range(100):
        f = pressure(T, eta_to_rho(eta, T, comp), comp) / P - 1
        h = 1e-6 * eta
        f_plus = pressure(T, eta_to_rho(eta + h, T, comp), comp) / P - 1
        f_minus = pressure(T, eta_to_rho(eta - h, T, comp), comp) / P - 1
        slope = (f_plus - f_minus) / (2 * h)

        if verbose:
            print(f"   {k:2d}   {eta:.10e}   {f: .3e}")
        step = -f / slope
        while not 0 < eta + step < ETA_MAX:          # stay inside 0 < eta < ETA_MAX
            step = step / 2
        eta = eta + step

        if abs(step) < 1e-12 * eta:
            if slope <= 0:                            # unstable root
                return None
            return eta_to_rho(eta, T, comp)
    return None


def wagner25(T, comp):
    """Vapour pressure in Pa from Wagner 25 (Aspen's WAGNER25), or None if T >= Tc.

        ln P = ln Pc + (C1 tau + C2 tau^1.5 + C3 tau^2.5 + C4 tau^5) / Tr
        Tr = T / Tc,   tau = 1 - Tr
    """
    C1, C2, C3, C4, ln_Pc, Tc = comp.wagner25
    Tr = T / Tc
    if Tr >= 1:
        return None
    tau = 1 - Tr
    return np.exp(ln_Pc + (C1 * tau + C2 * tau**1.5 + C3 * tau**2.5 + C4 * tau**5) / Tr)


def saturation(T, comp, tol=1e-10, verbose=False, P0=None):
    """Saturation pressure and the two coexisting densities at T.

    Returns (P_sat in Pa, rho_liquid, rho_vapour), densities in molecules/A^3,
    or None when there is no Wagner estimate or Newton fails.
    verbose=True prints every iteration. P0 (Pa) replaces the Wagner 25 estimate,
    e.g. the Psat of a nearby temperature when marching along the curve.
    """
    if P0 is None:
        if comp.wagner25 is None:
            raise ValueError(f"No Wagner 25 constants for {comp.name} in parameters.WAGNER25.")
        P0 = wagner25(T, comp)
        if P0 is None:
            return None
    P = P0

    if verbose:
        print(f"  saturation: P0 = {P/1e5:.6f} bar, Newton on g = ln phi_L - ln phi_V")
        print("    k   P [bar]            g            dP/P")
    step_old = 1.0
    for k in range(50):
        rho_l = density(T, P, comp, "liquid")
        rho_v = density(T, P, comp, "vapour")
        if rho_l is None or rho_v is None or abs(rho_l / rho_v - 1) < 1e-6:
            return None                               # only one root at this P

        g = ln_fugacity_coefficient(T, rho_l, comp) - ln_fugacity_coefficient(T, rho_v, comp)
        Z_l = compressibility(T, rho_l, comp)
        Z_v = compressibility(T, rho_v, comp)
        step = g / (Z_v - Z_l)
        if verbose:
            print(f"   {k:2d}   {P/1e5:.10f}   {g: .3e}   {step: .3e}")
        if abs(step) < tol:
            return P, rho_l, rho_v

        # Newton should shrink the step every iteration. At low T, Z_L ~ 1e-7 and
        # g cannot get below the round-off of double precision; once the step
        # stops shrinking, P is as good as it can get.
        if abs(step) >= abs(step_old) and abs(step) < 1e-5:
            return P, rho_l, rho_v

        step_old = step
        P = P * (1 + step)
    return None
