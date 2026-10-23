# pcsaft

PC-SAFT equation of state for pure components, written to be extended with a
two-state association term.

Master's project, Posgrado en Ingeniería Química (Procesos), UNAM.
Advisor: Dr. Milton Medeiros, Departamento de Fisicoquímica, Facultad de Química.

## Where this is going

The end goal is **PC-SAFT-TSAM**: replace PC-SAFT's Wertheim association term
with the two-state association model of the CTS equation of state, and treat
water as a mixture of a high-density and a low-density species in chemical
equilibrium,

    LDW  ⇌  HDW

so that water's anomalies — the density maximum, the sound-speed maximum, the
compressibility minimum, the solubility minimum of hydrocarbons in water — come
out of the equation of state rather than being fitted around.

The order of work is deliberate. Non-associating n-alkanes first, because they
exercise every term except association and have excellent reference data. Water
only once that is verified.

| Stage                                             | Status             |
| ------------------------------------------------- | ------------------ |
| PC-SAFT, pure non-associating compounds           | done               |
| Density solver, saturation, fugacity coefficients | done               |
| Validation against NIST, n-alkane series          | done, 11 compounds |
| Caloric properties and their validation           | done               |
| Association term (Wertheim / CTS two-state)       | not started        |
| Binary mixtures                                   | not started        |
| LDW ⇌ HDW chemical equilibrium                    | not started        |
| Parameter regression                              | not started        |
|                                                   |                    |

## Install

```bash
pip install -e .
```

The `-e` matters: the package stays editable in place, so an edit to a source
file takes effect immediately with no reinstall.

## Use

```python
import pcsaft

hexane = pcsaft.get("hexane")

# density at a given state
rho = pcsaft.density(300.0, 1e5, hexane, "liquid")   # molecules/A^3
print(hexane.mass_density(rho), "kg/m3")

# saturation, with no experimental starting guess
P_sat, rho_L, rho_V = pcsaft.saturation(340.0, hexane)
print(P_sat / 1e5, "bar")

# validate against reference data
ref = pcsaft.load_nist_saturation("hexane")
cmp = pcsaft.compare_saturation(hexane, ref)
print(pcsaft.report(hexane, cmp, T_max=450))
pcsaft.plotting.plot_saturation(cmp, hexane)
```

## Layout

```
pcsaft/
  constants.py      physical constants, universal constants of Table 1
  parameters.py     Component dataclass, the library from Table 2
  eos.py            a_hc, a_disp, a_assoc, a_res, Z, P, ln phi
  equilibrium.py    density roots, spinodal, saturation, critical point
  data.py           reading NIST exports, sanity checks
  validation.py     consistency checks and AAD reports
  plotting.py       the figures, styled once

data/nist/          reference data, one file per compound
scripts/            runnable, numbered, each writes to results/ and figures/
tests/              pytest, run before believing any number
results/, figures/  generated output, not versioned
```

## The design decision that matters

`a_res` is a sum of independent terms:

```python
def a_res(T, rho, comp):
    return a_hc(T, rho, comp) + a_disp(T, rho, comp) + a_assoc(T, rho, comp)
```

The thesis replaces `a_assoc`. That is one function, not an equation of state.
When the CTS term goes in, do not derive its analytical density derivative
first — add the energy, take Z from `compressibility_numeric`, confirm the
physics, and only then derive the analytical Z. The two must agree afterwards,
which is exactly what `check_consistency` tests.

## Two different questions

Keeping these apart is the whole discipline of the validation module.

**Is the code a correct implementation of PC-SAFT?** Answered with no
experimental data at all. `a_res` and `Z` are coded from different equations in
Appendix A, so if the analytical `Z` matches a numerical derivative of `a_res`
to eight digits, a typo in either is very unlikely.

**Is PC-SAFT a good model of reality?** Answered against reference data, and
only meaningful once the first question is settled. A 5% error means nothing if
you cannot rule out a sign error.

## Current numbers

Saturation against the NIST reference equation of state, excluding the top 10%
of each compound's temperature range. `scripts/02_validate_series.py`.

| compound | m | AAD% P_sat | AAD% ρ_L | published P_sat | published v |
|---|---|---|---|---|---|
| methane | 1.000 | 0.33 | 0.52 | 0.36 | 0.67 |
| ethane | 1.607 | 0.31 | 0.71 | 0.30 | 0.57 |
| propane | 2.002 | 0.09 | 0.24 | 1.29 | 0.77 |
| butane | 2.332 | 0.98 | 0.24 | 0.75 | 1.59 |
| pentane | 2.690 | 1.19 | 1.04 | 1.45 | 0.78 |
| hexane | 3.058 | 0.36 | 0.52 | 0.31 | 0.76 |
| heptane | 3.483 | 0.66 | 0.54 | 0.34 | 2.10 |
| octane | 3.818 | 0.27 | 0.28 | 0.77 | 1.59 |
| nonane | 4.208 | 0.65 | 0.20 | 0.89 | 0.32 |
| decane | 4.663 | 1.68 | 0.47 | 0.24 | 1.18 |
| dodecane | 5.306 | 1.81 | 0.43 | 2.10 | 0.93 |

Internal consistency: analytical vs numerical Z agree to better than 1e-8
relative across T = 200–450 K and η = 0.05–0.5, for m from 1 to 4.7.

Predicted critical point for hexane, 519.3 K vs 507.8 K experimental. PC-SAFT
is a mean-field theory and ignores critical fluctuations, so it overpredicts
T_c and degrades above roughly 0.93 T_c. AAD is therefore reported over a
subcritical window as well as the full range.

### Caloric properties

Predictions, not correlations: no caloric data entered the parameter
regression. `scripts/03_validate_caloric.py`.

| compound | AAD% Δh_vap | AAD% C_p − C_v |
|---|---|---|
| propane | 2.39 | 8.26 |
| hexane | 1.25 | 13.59 |
| decane | 3.79 | 9.86 |

Δh_vap at ~1–4% is a real success. C_p − C_v at 8–14% is not an implementation
problem: two independent routes — derivatives of P(T, ρ), and T V α²/κ_T from
the density solver — agree to better than 1e-5 relative, and the deviation
curves are smooth and collapse on T/T_c across compounds, which is what model
error looks like rather than a bug.

It reads worse than the 1.7% Gross & Sadowski report in their Table 3 because
that table is for *total* c_p, which for a liquid alkane is dominated by the
ideal-gas contribution supplied from outside the equation of state. C_p − C_v
removes that contribution entirely, so it is a considerably harsher test of the
same model.

Both quantities were chosen because they need no ideal-gas heat-capacity
correlation: Δh_vap because the ideal part cancels between two phases at the
same temperature, and C_p − C_v because it follows from P(T, ρ) alone.

## On the reference data

The NIST Webbook does not serve raw experimental points. It serves a
multiparameter reference equation of state fitted to the experimental database,
with uncertainty around 0.05% in liquid density — far below PC-SAFT's ~1%, so
it is a fair stand-in for experiment when the question is whether the code is
right.

It is **not** the right source for regressing parameters. Fitting to a
correlation inherits its smoothing and its bias, with no honest way to report
uncertainty. Parameter work uses NIST-TDE or the primary literature.

Diego Wagner 25

## Reference

Gross, J.; Sadowski, G. Perturbed-Chain SAFT: An Equation of State Based on a
Perturbation Theory for Chain Molecules. *Ind. Eng. Chem. Res.* **2001**, 40,
1244–1260.

Medeiros, M.; Armas-Alemán, C. O. Cubic Two-State Equation of State for
Associating Fluids. *Ind. Eng. Chem. Res.* **2008**.

Novak, N.; Liang, X.; Kontogeorgis, G. M. Prediction of water anomalous
properties by introducing the two-state theory in SAFT. *J. Chem. Phys.*
**2024**, 160, 104505.
