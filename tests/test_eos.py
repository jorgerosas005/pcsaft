"""Tests that must pass before any result is believed.

Run them with:  pytest -q

These are correctness tests for the code, not accuracy tests for the model.
Only test_hexane_matches_reference touches experimental data; everything else
is internal and needs no external source.
"""

import numpy as np
import pytest

import pcsaft


HEXANE = pcsaft.get("hexane")
METHANE = pcsaft.get("methane")
DECANE = pcsaft.get("decane")


# ------------------------------------------------------- segment diameter

def test_diameter_below_sigma():
    for comp in (METHANE, HEXANE, DECANE):
        for T in (100.0, 300.0, 700.0):
            assert pcsaft.segment_diameter(T, comp) < comp.sigma


def test_diameter_limits():
    # exp(-3 eps/kT) -> 0 as T -> 0, so d -> sigma: a cold collision stops at
    # the nominal diameter. As T -> infinity the exponential -> 1 and
    # d -> 0.88 sigma: a fast collision pushes 12% deeper into the soft wall.
    assert pcsaft.segment_diameter(1e-3, HEXANE) == pytest.approx(HEXANE.sigma, rel=1e-9)
    assert pcsaft.segment_diameter(1e9, HEXANE) == pytest.approx(0.88 * HEXANE.sigma, rel=1e-6)


def test_diameter_decreases_with_temperature():
    Ts = np.linspace(100, 800, 50)
    d = [pcsaft.segment_diameter(T, HEXANE) for T in Ts]
    assert np.all(np.diff(d) < 0)


# ------------------------------------------------- internal consistency

@pytest.mark.parametrize("comp", [METHANE, HEXANE, DECANE])
def test_analytical_Z_matches_numerical(comp):
    """a_res and Z come from different equations. They must still agree."""
    df = pcsaft.check_consistency(comp)
    assert df["rel_diff"].max() < 1e-6


def test_ideal_gas_limit():
    """Z -> 1 and a_res -> 0 at vanishing density."""
    rho = 1e-12
    assert pcsaft.compressibility(300.0, rho, HEXANE) == pytest.approx(1.0, abs=1e-8)
    assert pcsaft.a_res(300.0, rho, HEXANE) == pytest.approx(0.0, abs=1e-8)


def test_chain_term_vanishes_for_methane():
    """m = 1 means no bonds, so a_hc must reduce to m * a_hs."""
    rho = pcsaft.eta_to_rho(0.3, 200.0, METHANE)
    assert pcsaft.a_hc(200.0, rho, METHANE) == pytest.approx(
        METHANE.m * pcsaft.a_hs(200.0, rho, METHANE), rel=1e-12)


def test_association_is_zero_for_alkanes():
    assert pcsaft.a_assoc(300.0, 0.005, HEXANE) == 0.0


# ---------------------------------------------------- solvers behave

def test_density_roundtrip():
    """Solving P(rho)=P and putting rho back must return the same P."""
    P = 5e5
    rho = pcsaft.density(350.0, P, HEXANE, "liquid")
    assert pcsaft.pressure(350.0, rho, HEXANE) == pytest.approx(P, rel=1e-8)


def test_liquid_root_is_denser_than_vapour():
    P = 1e5
    rl = pcsaft.density(300.0, P, HEXANE, "liquid")
    rv = pcsaft.density(300.0, P, HEXANE, "vapour")
    assert rl > 100 * rv


def test_saturation_equalises_fugacity():
    """The definition of the answer: ln phi must match in both phases."""
    P, rl, rv = pcsaft.saturation(320.0, HEXANE)
    lnphi_l = pcsaft.ln_fugacity_coefficient(320.0, rl, HEXANE)
    lnphi_v = pcsaft.ln_fugacity_coefficient(320.0, rv, HEXANE)
    assert lnphi_l == pytest.approx(lnphi_v, abs=1e-9)


def test_saturation_pressure_increases_with_temperature():
    Ts = np.arange(200.0, 460.0, 20.0)
    P, _, _ = pcsaft.saturation_curve(Ts, HEXANE)
    assert np.all(np.diff(P) > 0)


def test_no_phase_split_above_critical():
    Tc, _ = pcsaft.critical_point(HEXANE)
    assert pcsaft.spinodal(Tc + 50.0, HEXANE) is None


# --------------------------------------------------- against the model

def test_hexane_matches_reference():
    """Model accuracy, not code correctness. Loose on purpose."""
    ref = pcsaft.load_nist_saturation("hexane")
    cmp = pcsaft.compare_saturation(HEXANE, ref)
    stats = pcsaft.aad(cmp, T_max=450)
    assert stats["P_sat"] < 1.5
    assert stats["rho_L"] < 1.5


# ------------------------------------------------------- caloric properties

def test_caloric_ideal_gas_limits():
    """At vanishing density: cv_res -> 0, cp_res -> 0, Cp - Cv -> R."""
    rho = pcsaft.density(400.0, 1.0, HEXANE, "vapour")
    assert pcsaft.cv_res(400.0, rho, HEXANE) == pytest.approx(0.0, abs=1e-4)
    assert pcsaft.cp_res(400.0, rho, HEXANE) == pytest.approx(0.0, abs=1e-4)
    assert pcsaft.cp_minus_cv(400.0, rho, HEXANE) == pytest.approx(pcsaft.R, rel=1e-4)


def test_cp_minus_cv_two_independent_routes():
    """Derivatives of P(T,rho) must agree with T V alpha^2 / kappa_T.

    Nothing is shared between the two calculations except the equation of
    state itself, so agreement rules out an algebra or unit error.
    """
    T = 300.0
    P = pcsaft.saturation(T, HEXANE)[0]
    rho = pcsaft.density(T, P, HEXANE, "liquid")
    route_a = pcsaft.cp_minus_cv(T, rho, HEXANE)

    dT, dP = 0.05, P * 1e-4
    V = 1.0 / rho
    alpha = (1 / pcsaft.density(T + dT, P, HEXANE, "liquid")
             - 1 / pcsaft.density(T - dT, P, HEXANE, "liquid")) / (2 * dT) / V
    kappa = -(1 / pcsaft.density(T, P + dP, HEXANE, "liquid")
              - 1 / pcsaft.density(T, P - dP, HEXANE, "liquid")) / (2 * dP) / V
    V_molar = V * 1e-30 * pcsaft.NA
    route_b = T * V_molar * alpha**2 / kappa

    assert route_a == pytest.approx(route_b, rel=1e-5)


def test_second_derivative_is_step_independent():
    """d2a/dT2 must not depend on the finite-difference step."""
    from pcsaft.caloric import _d2adT2
    rho = pcsaft.saturation(300.0, HEXANE)[1]
    vals = [_d2adT2(300.0, rho, HEXANE, h) for h in (0.3, 0.03, 0.003)]
    assert abs(vals[0] - vals[2]) / abs(vals[2]) < 1e-5


def test_enthalpy_of_vaporisation_predicts_reference():
    """dh_vap is a prediction: no caloric data entered the parameter fit."""
    ref = pcsaft.load_nist_saturation("hexane")
    cal = pcsaft.compare_caloric(HEXANE, ref)
    sub = cal[cal["T"] <= 450]
    assert np.nanmean(np.abs(sub["dev_dhvap"])) < 3.0


def test_enthalpy_of_vaporisation_vanishes_at_critical():
    """The two phases merge at Tc, so dh_vap must collapse toward zero."""
    Tc, _ = pcsaft.critical_point(HEXANE)
    near = pcsaft.enthalpy_of_vaporisation(Tc - 2.0, HEXANE)
    mid = pcsaft.enthalpy_of_vaporisation(400.0, HEXANE)
    assert 0 < near < 0.25 * mid


def test_associating_component_is_refused_until_implemented():
    """Water must not silently run without its association term."""
    water = pcsaft.get("water")
    assert water.associating
    with pytest.raises(NotImplementedError):
        pcsaft.a_res(373.15, 0.01, water)
    with pytest.raises(NotImplementedError):
        pcsaft.compressibility(373.15, 0.01, water)


def test_without_association_is_explicit_and_wrong():
    """The escape hatch works, and the error it leaves behind is enormous.

    Water with no association term overpredicts the boiling-point vapour
    pressure by more than an order of magnitude. That gap is the project.
    """
    w = pcsaft.get("water").without_association()
    assert not w.associating
    P, _, _ = pcsaft.saturation(373.15, w)
    assert P / 1e5 > 10 * 1.013
