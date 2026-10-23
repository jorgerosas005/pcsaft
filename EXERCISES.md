
# Exercises

Ten exercises that walk through Gross & Sadowski (2001) with the code in your
hands. Each one names the part of the paper it belongs to, gives you a snippet
to run, and asks you to **predict the answer before running it**.

The prediction is not a formality. If you predict right, you understood the
paper. If you predict wrong, you have just found the exact paragraph you need
to reread, which is worth more than any amount of passive highlighting.

Work through them in order; each builds on the last. Two or three per sitting
is plenty. Put your answers in a file next to this one.

Every snippet assumes you are in the project folder with the package installed:

```python
import warnings; warnings.filterwarnings("ignore")
import numpy as np
import pcsaft
```

---

## Block A — the molecular model

### 1. The temperature-dependent diameter
*Paper: §2.1, eq 3*

**Predict first.** You have methane (ε/k = 150 K) and decane (ε/k = 244 K). At
300 K, which one has shrunk further below its own σ, in relative terms? Why?

```python
for name in ["methane", "decane"]:
    c = pcsaft.get(name)
    d = pcsaft.segment_diameter(300.0, c)
    print(f"{name:8s} d/sigma = {d/c.sigma:.4f}")

# now plot the whole curve against reduced temperature
import matplotlib.pyplot as plt
for name in ["methane", "decane"]:
    c = pcsaft.get(name)
    T = np.linspace(20, 2000, 300)
    d = np.array([pcsaft.segment_diameter(t, c) for t in T])
    plt.plot(T / c.eps_k, d / c.sigma, label=name)
plt.xlabel("T / (eps/k)"); plt.ylabel("d / sigma"); plt.legend(); plt.show()
```

**Then answer.** The two curves land exactly on top of each other when plotted
this way. What does that tell you about how many independent parameters the
diameter really has? What is ε/k doing in a formula that is supposed to be
about size?

---

### 2. Which term is doing the work
*Paper: §2.2, eq 4 and 10*

**Predict first.** For saturated liquid hexane at 300 K: will a_hc be positive
or negative? a_disp? Which is larger in magnitude? What would it mean physically
if they had the same sign?

```python
c = pcsaft.get("hexane")
rho = pcsaft.saturation(300.0, c)[1]
print("a_hc   =", round(pcsaft.a_hc(300.0, rho, c), 4))
print("a_disp =", round(pcsaft.a_disp(300.0, rho, c), 4))
print("a_res  =", round(pcsaft.a_res(300.0, rho, c), 4))
```

**Then answer.** The rv  result is a large positive number and a larger negative
one, nearly cancelling. Why is a liquid a near-cancellation of two big effects,
and what does that imply about how accurate each term has to be?

---

## Block B — the hard-chain reference

### 3. What the chain term actually buys you
*Paper: §2.2, eq 4*

**Predict first.** Take hexane's σ and ε/k but set m = 1, so it becomes one
fat sphere instead of a three-segment chain. Will its predicted saturation
pressure go up or down? By roughly how much — 10%, or orders of magnitude?

```python
c = pcsaft.get("hexane")
sphere = pcsaft.Component("fake hexane", c.M, 1.0, c.sigma, c.eps_k)
for comp in (c, sphere):
    P, rl, rv = pcsaft.saturation(300.0, comp)
    print(f"{comp.name:12s} Psat = {P/1e5:10.4f} bar   rhoL = {comp.mass_density(rl):7.1f} kg/m3")
```

**Then answer.** Connect the size of that change to Figure 1 in the paper,
where PC-SAFT beats SAFT mainly for chain molecules. What physical effect are
you removing when you set m = 1?

---

### 4. The radial distribution function
*Paper: §2.2, eq 8 and 9*

**Predict first.** g_hs at contact measures how much more likely two segments
are to touch than in an ideal gas. As packing increases toward η → 0.74, should
it go to 0, to 1, or to infinity?

```python
c = pcsaft.get("hexane")
eta = np.linspace(0.01, 0.6, 200)
g = [pcsaft.radial_distribution(300.0, pcsaft.eta_to_rho(e, 300.0, c), c) for e in eta]
plt.plot(eta, g); plt.xlabel("eta"); plt.ylabel("g_hs at contact"); plt.show()
```

**Then answer.** The chain term is −(m−1)·ln(g_hs). As η grows, g_hs grows, so
ln(g_hs) grows, so the chain term becomes *more negative*. Explain in words why
forming a chain lowers the free energy more in a dense fluid than in a dilute
one.

---

## Block C — dispersion

### 5. The integrals I₁ and I₂
*Paper: §2.2, eq 16 and 17*

**Predict first.** I₁ and I₂ are power series in η with coefficients that depend
on m. At η = 0 (infinite dilution) what must they equal, and why?

```python
from pcsaft.eos import _dispersion_coefficients
a, b = _dispersion_coefficients(pcsaft.get("hexane").m)
eta = np.linspace(0, 0.6, 200)
i = np.arange(7)
I1 = [np.sum(a * e**i) for e in eta]
I2 = [np.sum(b * e**i) for e in eta]
plt.plot(eta, I1, label="I1"); plt.plot(eta, I2, label="I2")
plt.xlabel("eta"); plt.legend(); plt.show()
print("a_i(m) =", np.round(a, 4))
```

**Then answer.** The paper fitted those 42 universal constants to n-alkane data
once, and they are never refitted. What does that mean for a compound that is
*not* an n-alkane — say, water? Is that a strength or a hidden assumption?

---

### 6. The compressibility correction C₁
*Paper: §2.2, eq 13*

**Predict first.** C₁ appears only in the second-order dispersion term. At
η → 0 what should it be? Does it grow or shrink as the fluid gets dense?

```python
from pcsaft.eos import _C1_C2
c = pcsaft.get("hexane")
eta = np.linspace(0.001, 0.6, 200)
C1 = [_C1_C2(e, c.m)[0] for e in eta]
plt.plot(eta, C1); plt.xlabel("eta"); plt.ylabel("C1"); plt.show()
```

**Then answer.** C₁ is defined as (1 + Z_hc + ρ dZ_hc/dρ)⁻¹. Why does a
*second-order perturbation* term need to know about the compressibility of the
*reference* fluid? (Hint: Barker–Henderson, and what fluctuates in a
compressible fluid.)

---

### 7. Kill the attraction
*Paper: §2.2 — the point of the whole dispersion term*

**Predict first.** If you delete a_disp entirely, leaving only hard chains,
does hexane still have a liquid–vapour transition? Will it just become a worse
prediction, or will something qualitative break?

```python
from pcsaft import eos
c = pcsaft.get("hexane")

keep_a, keep_Z = eos.a_disp, eos.Z_disp
eos.a_disp = lambda T, rho, comp: 0.0
eos.Z_disp = lambda T, rho, comp: 0.0

print("spinodal without dispersion:", pcsaft.spinodal(300.0, c))
print("saturation without dispersion:", pcsaft.saturation(300.0, c))

eos.a_disp, eos.Z_disp = keep_a, keep_Z          # always restore
print("saturation restored:", round(pcsaft.saturation(300.0, c)[0], 1), "Pa")
```

**Then answer.** Both come back `None`. Explain, in one sentence, what `None`
means physically here. This is the single most important exercise in the set:
it shows which term is responsible for the existence of a liquid at all.

*(Note why both `a_disp` and `Z_disp` had to be switched off. What would have
gone wrong if you had only switched off one?)*

---

## Block D — solvers and phase equilibrium

### 8. Look at a real isotherm
*Paper: Appendix A, the density algorithm*

**Predict first.** Below T_c, sketch on paper what P versus η looks like. How
many times does a horizontal line at P = P_sat cross it? Which crossings are
physical?

```python
c = pcsaft.get("hexane")
T = 400.0
eta = np.logspace(-4, np.log10(0.6), 500)
P = [pcsaft.pressure(T, pcsaft.eta_to_rho(e, T, c), c) for e in eta]
Psat = pcsaft.saturation(T, c)[0]

plt.plot(eta, np.array(P)/1e5)
plt.axhline(Psat/1e5, color="k", ls="--", label="Psat")
plt.ylim(-30, 60); plt.xlabel("eta"); plt.ylabel("P / bar"); plt.legend(); plt.show()
print("spinodal (eta_v, eta_l, Pmax, Pmin):", pcsaft.spinodal(T, c))
```

**Then answer.** Part of that curve has negative pressure. Is that a bug, a
numerical artefact, or real physics? What is the middle root, and what is
`dP/dρ` doing there? Why can the code use the two turning points as a rigorous
bracket for P_sat without ever consulting experimental data?

---

### 9. Parameter sensitivity
*Paper: §3.1 and Figure 5*

**Predict first.** Change ε/k by +1% and by −1%. Which changes more, P_sat or
ρ_L? By how much, roughly?

```python
c = pcsaft.get("hexane")
base = pcsaft.saturation(300.0, c)
for f in (0.99, 1.00, 1.01):
    t = pcsaft.Component("t", c.M, c.m, c.sigma, c.eps_k * f)
    P, rl, _ = pcsaft.saturation(300.0, t)
    print(f"eps/k x{f}:  Psat {P/1e5:8.5f} bar ({100*(P/base[0]-1):+6.2f}%)   "
          f"rhoL {t.mass_density(rl):6.1f} ({100*(rl/base[1]-1):+5.2f}%)")
```

Repeat for σ and for m.

**Then answer.** One parameter moves P_sat enormously while barely moving ρ_L.
That is *why* the regression works: the two properties constrain different
parameters. Now the harder question — when you fit water's two-state parameters
later, which properties will you need in the objective function so that the
parameters are actually identifiable?

---

## Block E — toward the thesis

### 10. Why water needs association
*Paper: nothing. This is the gap the project fills.*

**Predict first.** Water's PC-SAFT parameters without association are
m = 1.0656, σ = 3.0007 Å, ε/k = 366.51 K. At the normal boiling point, 373.15 K,
the real saturation pressure is 1.013 bar and the liquid density is 958 kg/m³.
How far off will a model with no association term be — 10%? A factor of 2? More?

```python
# Asking for water directly raises NotImplementedError, on purpose: the
# association term is still a stub, and silently ignoring hydrogen bonding
# would give numbers that look plausible and are badly wrong. Dropping it has
# to be deliberate.
c = pcsaft.get("water").without_association()
P, rl, rv = pcsaft.saturation(373.15, c)
print(f"Psat = {P/1e5:7.3f} bar   (real 1.013)")
print(f"rhoL = {c.mass_density(rl):7.1f} kg/m3   (real 958)")

# try it without the escape hatch and read the error
pcsaft.saturation(373.15, pcsaft.get("water"))
```

**Then answer.** The prediction comes out around 42 bar instead of 1 bar — the
model says water at 100 °C should be a gas at forty atmospheres. That factor of
forty is hydrogen bonding, and it is the entire reason your thesis exists.

Now write down, in your own words and in three sentences:

1. What the association term adds that dispersion cannot capture.
2. Why the two-state (CTS / TSAM) treatment is different from Wertheim's TPT1.
3. Which of water's anomalies you expect the two-state picture to reproduce
   that a four-site Wertheim model cannot.

If you can answer those three without looking anything up, you are ready to read
the Medeiros CTS paper properly.

---

## Where to go after

- Reread §2.2 with exercise 7's result in mind. The paper's structure —
  reference fluid plus perturbation — stops being an abstraction once you have
  seen the reference fluid fail to condense on its own.
- Open `pcsaft/eos.py` and match every function to its equation number. You have
  now used each one.
- Then, and only then, open the CTS paper. The association slot in `a_assoc` is
  waiting.
