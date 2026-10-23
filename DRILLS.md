# Drills

Practice using the package, not the theory. Short reps, one new thing at a time,
with the answer next to each so you know immediately whether you got it.

Work in a terminal: `py` then paste. Or a scratch file you keep overwriting.
Every drill starts from:

```python
import warnings; warnings.filterwarnings("ignore")
import numpy as np
import pcsaft
```

Levels 1-5 are fast and deliberately repetitive — change one word, run again.
Levels 6-10 add a piece at a time. Do a level, then move on; don't perfect it.

---

## Level 1 — the Component object  (5 min)

**The pattern.** One line, run it eight times, change only the name.

```python
c = pcsaft.get("ethane")
print(c.m, c.sigma, c.eps_k, c.M)
```

**Reps:** methane, ethane, propane, hexane, decane, nitrogen, carbon dioxide, water

**Check**

| compound | m | sigma | eps/k | M |
|---|---|---|---|---|
| methane | 1.0000 | 3.7039 | 150.03 | 16.043 |
| ethane | 1.6069 | 3.5206 | 191.42 | 30.07 |
| propane | 2.0020 | 3.6184 | 208.11 | 44.096 |
| hexane | 3.0576 | 3.7983 | 236.77 | 86.177 |
| decane | 4.6627 | 3.8384 | 243.87 | 142.285 |
| nitrogen | 1.2053 | 3.3130 | 90.96 | 28.01 |
| carbon dioxide | 2.0729 | 2.7852 | 169.21 | 44.01 |

**Then, one line each:**

```python
sorted(pcsaft.LIBRARY)              # 17 compounds
pcsaft.N_ALKANES                    # 12 in order
pcsaft.get("hexane").T_range        # (177, 503)
pcsaft.get("hexane").aad_psat       # 0.31  <- the paper's published value
pcsaft.get("water").associating     # True
```

**One trap, on purpose:**

```python
pcsaft.get("Ethane")      # works, case-insensitive
pcsaft.get("etano")       # KeyError, and it lists what IS available
```

---

## Level 2 — one function, many inputs  (5 min)

**The pattern.**

```python
c = pcsaft.get("ethane")
print(pcsaft.segment_diameter(300.0, c))
```

**Reps:** T = 100, 200, 300, 400, 500

**Check**

| T (K) | d (Å) | d/σ |
|---|---|---|
| 100 | 3.519245 | 0.999615 |
| 200 | 3.496677 | 0.993205 |
| 300 | 3.458302 | 0.982305 |
| 400 | 3.420068 | 0.971445 |
| 500 | 3.386632 | 0.961947 |

**New: converting between eta and rho.** Prove the roundtrip returns what you put in.

```python
rho = pcsaft.eta_to_rho(0.3, 250.0, c)
pcsaft.packing_fraction(250.0, rho, c)     # 0.3 exactly
c.mass_density(rho)                        # kg/m3
c.number_density(c.mass_density(rho))      # back to rho
```

---

## Level 3 — density  (10 min)

**The pattern.**

```python
c = pcsaft.get("ethane")
rho = pcsaft.density(250.0, 1.0e6, c, "liquid")
print(c.mass_density(rho))
```

**Reps:** (T, P) = (250, 1e5), (250, 1e6), (200, 1e5), (150, 1e5) — each for `"liquid"` and `"vapour"`

**Check** (kg/m³)

| T (K) | P (Pa) | liquid | vapour |
|---|---|---|---|
| 250 | 1e5 | 447.9032 | 1.4649 |
| 250 | 1e6 | 451.0167 | 16.7864 |
| 200 | 1e5 | 527.3093 | 1.8511 |
| 150 | 1e5 | 587.4407 | 2.5403 |

**Now a deliberate surprise.** Run this one:

```python
pcsaft.density(300.0, 1.0e6, c, "liquid")   # -> 13.0753 kg/m3
pcsaft.density(300.0, 1.0e6, c, "vapour")   # -> 13.0753 kg/m3, the SAME
```

Same answer for both phases. Find out why before reading on:

```python
pcsaft.spinodal(300.0, c)    # Pmax = 45.06 bar, Pmin = 39.18 bar
```

At 300 K ethane is at T/T_c ≈ 0.97, so the loop has shrunk to between 39 and
45 bar. A line at 10 bar passes entirely below it and crosses the isotherm only
once. One root, no two phases. Compare with 250 K, where the loop runs from
−92 to +24 bar and 10 bar sits inside it.

**New: the roundtrip check.** Any density you solve for must reproduce the pressure.

```python
rho = pcsaft.density(250.0, 1.0e6, c, "liquid")
pcsaft.pressure(250.0, rho, c)        # 1000000.0000...
```

---

## Level 4 — saturation  (10 min)

**The pattern.**

```python
Psat, rho_l, rho_v = pcsaft.saturation(250.0, c)
print(Psat / 1e5, c.mass_density(rho_l), c.mass_density(rho_v))
```

**Reps:** T = 150, 200, 250, 280, 300

**Check**

| T (K) | P_sat (bar) | ρ_L | ρ_V |
|---|---|---|---|
| 150 | 0.096484 | 587.371 | 0.2337 |
| 200 | 2.178007 | 527.465 | 4.1527 |
| 250 | 13.038828 | 452.026 | 23.2462 |
| 280 | 28.071033 | 382.949 | 54.6449 |
| 300 | 43.249322 | 301.098 | 103.9212 |

**New: the critical point and what's above it.**

```python
pcsaft.critical_point(c)        # (308.93 K, 51.60 bar)
pcsaft.saturation(320.0, c)     # None  <- above Tc, no phase split
pcsaft.spinodal(320.0, c)       # None, same reason
```

**New: the equilibrium condition.** Verify it rather than trusting it.

```python
Psat, rl, rv = pcsaft.saturation(250.0, c)
a = pcsaft.ln_fugacity_coefficient(250.0, rl, c)
b = pcsaft.ln_fugacity_coefficient(250.0, rv, c)
print(a, b, a - b)              # identical to ~1e-14
```

---

## Level 5 — loops over the library  (10 min)

Same calls as before, now inside a `for`. This is where it starts saving you time.

```python
for name in pcsaft.N_ALKANES:
    c = pcsaft.get(name)
    print(f"{name:10s} m = {c.m:.4f}")
```

**Rep 2 — P_sat at a fixed temperature across the series.** Note that some
return `None`: methane's critical point is 190 K, so at 250 K it has no
saturation pressure at all. Handle that.

```python
for name in pcsaft.N_ALKANES:
    c = pcsaft.get(name)
    r = pcsaft.saturation(250.0, c)
    print(f"{name:10s} {'-' if r is None else f'{r[0]/1e5:10.5f} bar'}")
```

**Check**

| compound | P_sat at 250 K |
|---|---|
| methane | none (above T_c) |
| ethane | 13.03883 bar |
| propane | 2.18184 bar |
| butane | 0.39332 bar |
| pentane | 0.07607 bar |
| hexane | 0.01528 bar |

Four orders of magnitude from ethane to hexane. That's the chain length.

**Rep 3 — critical point of every alkane.** Write it yourself; same shape.

**Rep 4 — an array of temperatures at once:**

```python
T = np.arange(150.0, 305.0, 10.0)
P, rl, rv = pcsaft.saturation_curve(T, c)
```

---

## Level 6 — reference data and AAD  (15 min)

```python
ref = pcsaft.load_nist_saturation("ethane")
print(len(ref), list(ref.columns))
```

18 rows, columns: `T, P, rho_L, rho_V, h_L, h_V, Cv_L, Cp_L, dh_vap, Cp_minus_Cv_L`

**Always check data before trusting it:**

```python
pcsaft.check_saturation_data(ref)
```

| check | value |
|---|---|
| no missing values | True |
| P strictly increasing | True |
| ρ_L decreasing | True |
| ρ_V increasing | True |
| ρ_L > ρ_V | True |
| normal boiling point | 133.62 K |
| Clausius–Clapeyron R² | 0.9996 |

**Then compare and score:**

```python
cmp = pcsaft.compare_saturation(c, ref)
pcsaft.aad(cmp)                  # full range
pcsaft.aad(cmp, T_max=290)       # subcritical only
print(pcsaft.report(c, cmp, T_max=290))
```

| window | n | P_sat | ρ_L | ρ_V |
|---|---|---|---|---|
| full | 18 | 0.3455 | 0.6495 | 3.4518 |
| ≤ 290 K | 14 | 0.2683 | 0.6263 | 2.1086 |

Note how ρ_V improves a lot when you cut the top. That's the critical region.

**Rep:** do the same for propane, hexane and decane. Four lines each, same shape.

---

## Level 7 — figures  (15 min)

Three functions, always in this order.

```python
from pcsaft.plotting import figure, finish, save

T = np.linspace(100, 305, 200)
d = [pcsaft.segment_diameter(t, c) for t in T]

fig, ax = figure()
ax.plot(T, d, label="ethane")
finish(ax, "T / K", "d / A", legend=True)
print(save(fig, "drill_diameter"))
```

**Reps:**
1. d(T) for three compounds on one axis
2. P_sat vs 1000/T on a log y axis (`ax.semilogy`)
3. The coexistence dome: ρ_L and ρ_V against T
4. Two panels: `figure(1, 2)`, then `ax[0]` and `ax[1]`

**Then the ready-made ones:**

```python
pcsaft.plotting.plot_saturation(cmp, c)        # the three-panel figure
pcsaft.save_results(cmp, "drill_ethane")       # -> results/drill_ethane.csv
```

Everything lands in `figures/` and `results/` regardless of where you launched from.

---

## Level 8 — caloric properties  (15 min)

```python
Psat, rl, rv = pcsaft.saturation(250.0, c)
pcsaft.h_res(250.0, rl, c)
pcsaft.s_res(250.0, rl, c)
pcsaft.cv_res(250.0, rl, c)
pcsaft.cp_minus_cv(250.0, rl, c)
```

**Check** (ethane, 250 K, J/mol or J/mol·K)

| | h_res | s_res | cv_res | C_p − C_v |
|---|---|---|---|---|
| liquid | −12496.49 | −22.134 | 10.287 | 43.542 |
| vapour | −1099.62 | −1.220 | 1.310 | 17.138 |

All the residual values are negative or small — they're differences from the
ideal gas at the same state, not absolute quantities.

```python
pcsaft.enthalpy_of_vaporisation(250.0, c) / 1000    # 11.3969 kJ/mol
```

**Rep:** Δh_vap across T = 150 to 300. It must fall toward zero at T_c. Plot it.

**Then the comparison layer:**

```python
cal = pcsaft.compare_caloric(c, ref)
cal[["T", "dhvap_ref", "dhvap_calc", "dev_dhvap"]]
```

---

## Level 9 — diagnostics  (10 min)

The functions you run when a number looks wrong.

```python
pcsaft.check_consistency(c)              # analytical Z vs numerical derivative
pcsaft.check_ideal_gas_limit(c)          # Z -> 1, a_res -> 0
pcsaft.check_diameter_limits(c)          # d/sigma -> 1 and 0.88
```

**Check:** worst `rel_diff` ≈ 4.4e−09, ideal-gas Z = 1.000, a_res ≈ −7.4e−08,
diameter limits 1.000000 and 0.880001.

**The one to internalise:**

```python
rho = pcsaft.eta_to_rho(0.35, 250.0, c)
pcsaft.compressibility(250.0, rho, c)            # analytical
pcsaft.compressibility_numeric(250.0, rho, c)    # from d(a_res)/d(rho)
```

They agree to 8+ digits. When you add the association term, this is the test
that tells you whether your new analytical Z is right.

**And the water guard:**

```python
w = pcsaft.get("water")
pcsaft.saturation(373.15, w)                      # NotImplementedError, on purpose
pcsaft.saturation(373.15, w.without_association()) # 41.703 bar vs 1.013 real
```

---

## Level 10 — build one yourself  (30 min)

No pattern given. Write a script from scratch that:

1. Takes a compound name and a temperature range
2. Loads the NIST file and checks it
3. Computes the saturation curve
4. Prints an AAD table against the reference, with a subcritical cut
5. Saves a CSV and a two-panel figure
6. Prints Δh_vap at three temperatures

You have used every piece. Now assemble them without looking at
`scripts/01_validate_hexane.py` — then open it and compare.

If your version is shorter, use yours.

---

## Quick reference

| What you want | Call |
|---|---|
| a compound | `pcsaft.get("ethane")` |
| its parameters | `c.m  c.sigma  c.eps_k  c.M  c.T_range` |
| segment diameter | `pcsaft.segment_diameter(T, c)` |
| eta from rho, rho from eta | `packing_fraction(T, rho, c)` / `eta_to_rho(eta, T, c)` |
| kg/m3 from rho | `c.mass_density(rho)` |
| the terms | `a_hs  a_hc  a_disp  a_assoc  a_res` |
| Z and P | `compressibility(T, rho, c)` / `pressure(T, rho, c)` |
| ln phi | `ln_fugacity_coefficient(T, rho, c)` |
| density at T, P | `density(T, P, c, "liquid")` |
| turning points | `spinodal(T, c)` |
| saturation | `saturation(T, c)` -> (Psat, rho_l, rho_v) |
| a whole curve | `saturation_curve(T_array, c)` |
| critical point | `critical_point(c)` |
| reference data | `load_nist_saturation("ethane")` |
| sanity-check it | `check_saturation_data(ref)` |
| model vs data | `compare_saturation(c, ref)` -> `aad(cmp, T_max=...)` -> `report(...)` |
| caloric | `h_res  s_res  cv_res  cp_res  cp_minus_cv  enthalpy_of_vaporisation` |
| is my code right | `check_consistency(c)` |
| a figure | `figure()` -> `finish(ax, ...)` -> `save(fig, "name")` |
| a csv | `save_results(df, "name")` |
