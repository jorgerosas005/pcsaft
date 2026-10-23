"""The PC-SAFT equation of state for a pure component.

Equation numbers refer to Appendix A of Gross & Sadowski (2001).

The residual Helmholtz energy is written as a sum of independent contributions

    a_res = a_hc + a_disp + a_assoc

and each contribution is its own function. That split is the point of the file:
the thesis replaces a_assoc with the CTS / two-state association term, and that
is a change to one function, not to the equation of state.

Every function takes (T, rho, comp) with
    T    temperature, K
    rho  number density, molecules/A^3
    comp a parameters.Component
"""

import numpy as np

from .constants import A_COEF, B_COEF, kB


# ------------------------------------------------------------ common pieces

def segment_diameter(T, comp):
    """Temperature-dependent segment diameter d, eq A.9, in A.

    The real repulsion is steep but soft, so the effective hard sphere shrinks
    as T rises and molecules interpenetrate further before turning around.
    """
    return comp.sigma * (1 - 0.12 * np.exp(-3 * comp.eps_k / T))


def zetas(T, rho, comp):
    """The four packing moments zeta_0..zeta_3, eq A.8, for a pure component.

    zeta_3 is the packing fraction eta.
    """
    d = segment_diameter(T, comp)
    c = np.pi / 6 * rho * comp.m
    return np.array([c, c * d, c * d**2, c * d**3])


def packing_fraction(T, rho, comp):
    return zetas(T, rho, comp)[3]


def eta_to_rho(eta, T, comp):
    """Invert the packing fraction back to a number density, eq A.20."""
    return 6 / np.pi * eta / (comp.m * segment_diameter(T, comp) ** 3)


def _dispersion_coefficients(m):
    """a_i(m) and b_i(m), eqs A.18 and A.19."""
    f1 = (m - 1) / m
    f2 = f1 * (m - 2) / m
    a = A_COEF[:, 0] + f1 * A_COEF[:, 1] + f2 * A_COEF[:, 2]
    b = B_COEF[:, 0] + f1 * B_COEF[:, 1] + f2 * B_COEF[:, 2]
    return a, b


def _C1_C2(eta, m):
    """Compressibility correction C1 (eq A.11) and its derivative C2 (eq A.31)."""
    C1 = 1 / (1
              + m * (8 * eta - 2 * eta**2) / (1 - eta)**4
              + (1 - m) * (20 * eta - 27 * eta**2 + 12 * eta**3 - 2 * eta**4)
              / ((1 - eta) * (2 - eta))**2)
    C2 = -C1**2 * (m * (-4 * eta**2 + 20 * eta + 8) / (1 - eta)**5
                   + (1 - m) * (2 * eta**3 + 12 * eta**2 - 48 * eta + 40)
                   / ((1 - eta) * (2 - eta))**3)
    return C1, C2


def radial_distribution(T, rho, comp):
    """Hard-sphere radial distribution function at contact, eq A.7.

    For a pure component d_i = d_j, so d_i*d_j/(d_i+d_j) collapses to d/2.
    """
    d = segment_diameter(T, comp)
    z = zetas(T, rho, comp)
    return (1 / (1 - z[3])
            + (d / 2) * 3 * z[2] / (1 - z[3])**2
            + (d / 2)**2 * 2 * z[2]**2 / (1 - z[3])**3)


# ------------------------------------------------- Helmholtz energy, by term

def a_hs(T, rho, comp):
    """Hard-sphere reference, eq A.6. Per segment, dimensionless (a/NkT)."""
    z = zetas(T, rho, comp)
    return 1 / z[0] * (3 * z[1] * z[2] / (1 - z[3])
                       + z[2]**3 / (z[3] * (1 - z[3])**2)
                       + (z[2]**3 / z[3]**2 - z[0]) * np.log(1 - z[3]))


def a_hc(T, rho, comp):
    """Hard-chain reference, eq A.4. Segments minus the bonds that join them."""
    return comp.m * a_hs(T, rho, comp) - (comp.m - 1) * np.log(radial_distribution(T, rho, comp))


def a_disp(T, rho, comp):
    """Dispersion, eqs A.10, A.12, A.13, A.16, A.17. Second-order perturbation."""
    eta = packing_fraction(T, rho, comp)
    a, b = _dispersion_coefficients(comp.m)
    i = np.arange(7)
    I1 = np.sum(a * eta**i)
    I2 = np.sum(b * eta**i)
    C1, _ = _C1_C2(eta, comp.m)

    m2es3 = comp.m**2 * (comp.eps_k / T) * comp.sigma**3
    m2e2s3 = comp.m**2 * (comp.eps_k / T)**2 * comp.sigma**3
    return -2 * np.pi * rho * I1 * m2es3 - np.pi * rho * comp.m * C1 * I2 * m2e2s3


def a_assoc(T, rho, comp):
    """Association. Returns 0 for a non-associating compound.

    THIS IS THE SLOT THE THESIS REPLACES.

    Plain PC-SAFT uses Wertheim's first-order thermodynamic perturbation theory
    (TPT1), the same term CPA carries:

        a_assoc = sum_A [ ln X_A - X_A/2 ] + M/2

    with the fraction of sites A not bonded, X_A, solved self-consistently from
    the association strength Delta_AB, which itself needs the radial
    distribution function above.

    The project replaces this with the two-state association model of the CTS
    equation of state (Medeiros 2008), so that water becomes a mixture of two
    species in chemical equilibrium instead of one species with four sites.

    When that term is added, do NOT immediately derive its analytical density
    derivative. Add the energy here, get Z from `compressibility_numeric`, and
    only derive the analytical Z once the physics is confirmed. The two must
    then agree, which is exactly the check in validation.check_consistency.
    """
    if not comp.associating:
        return 0.0
    raise NotImplementedError(
        "Association term not implemented yet. See the docstring for where the "
        "CTS / two-state term goes."
    )


def a_res(T, rho, comp):
    """Total residual Helmholtz energy, a_res/(N k T), dimensionless."""
    return a_hc(T, rho, comp) + a_disp(T, rho, comp) + a_assoc(T, rho, comp)


# -------------------------------------------- compressibility factor, by term

def Z_hc(T, rho, comp):
    """Hard-chain contribution to Z, eqs A.25 to A.27."""
    d = segment_diameter(T, comp)
    z = zetas(T, rho, comp)

    Z_hs = (z[3] / (1 - z[3])
            + 3 * z[1] * z[2] / (z[0] * (1 - z[3])**2)
            + (3 * z[2]**3 - z[3] * z[2]**3) / (z[0] * (1 - z[3])**3))

    g = radial_distribution(T, rho, comp)
    rho_dg = (z[3] / (1 - z[3])**2
              + (d / 2) * (3 * z[2] / (1 - z[3])**2 + 6 * z[2] * z[3] / (1 - z[3])**3)
              + (d / 2)**2 * (4 * z[2]**2 / (1 - z[3])**3 + 6 * z[2]**2 * z[3] / (1 - z[3])**4))

    return comp.m * Z_hs - (comp.m - 1) / g * rho_dg


def Z_disp(T, rho, comp):
    """Dispersion contribution to Z, eqs A.28 to A.31."""
    eta = packing_fraction(T, rho, comp)
    a, b = _dispersion_coefficients(comp.m)
    i = np.arange(7)
    I2 = np.sum(b * eta**i)
    detaI1 = np.sum(a * (i + 1) * eta**i)
    detaI2 = np.sum(b * (i + 1) * eta**i)
    C1, C2 = _C1_C2(eta, comp.m)

    m2es3 = comp.m**2 * (comp.eps_k / T) * comp.sigma**3
    m2e2s3 = comp.m**2 * (comp.eps_k / T)**2 * comp.sigma**3
    return (-2 * np.pi * rho * detaI1 * m2es3
            - np.pi * rho * comp.m * (C1 * detaI2 + C2 * eta * I2) * m2e2s3)


def compressibility(T, rho, comp):
    """Z = P v / (k T), eq A.24, from the analytical derivatives."""
    Z = 1 + Z_hc(T, rho, comp) + Z_disp(T, rho, comp)
    if comp.associating:
        raise NotImplementedError("Analytical Z for the association term is not implemented.")
    return Z


def compressibility_numeric(T, rho, comp, rel_step=1e-6):
    """Z from a central difference of a_res, Z = 1 + rho (d a_res / d rho).

    Slower and less precise than `compressibility`, but it needs no derivation.
    Use it to check a new term, and as the first working Z when a term is added.
    """
    h = rho * rel_step
    da = (a_res(T, rho + h, comp) - a_res(T, rho - h, comp)) / (2 * h)
    return 1 + rho * da


# ------------------------------------------------------ pressure and fugacity

def pressure(T, rho, comp):
    """Pressure in Pa, eq A.23. The 1e30 converts A^-3 to m^-3."""
    return compressibility(T, rho, comp) * kB * T * rho * 1e30


def ln_fugacity_coefficient(T, rho, comp):
    """ln phi, eqs A.32 and A.33.

    For a pure component every derivative with respect to mole fraction
    vanishes, so mu_res/kT reduces to a_res + (Z - 1).
    """
    Z = compressibility(T, rho, comp)
    return a_res(T, rho, comp) + (Z - 1) - np.log(Z)
