# PC-SAFT thesis package — context for Claude

Jorge Rosas, MSc (Ing. Química, Procesos), UNAM. Advisor: Dr. Milton Medeiros.
Previous student on the same line: Diego de la Cruz (PC-SAFT up to binaries).

## The project (from the protocol)

- Goal: describe water's anomalies (density maximum at 3.98 °C, minimum in alkane
  solubility near 25 °C) with **PC-SAFT** and **PC-SAFT-TST**, treating liquid water
  as two structures, LAD (high density) ⇌ LBD (low density), in chemical equilibrium.
- PC-SAFT = Gross & Sadowski 2001, with Wertheim (TPT1) association.
- PC-SAFT-TST = Wertheim replaced by the two-state theory (Medeiros & Téllez-Arredondo,
  cubic two-state). The PC-SAFT version is still in development in the group.
- Original contribution: the LAD ⇌ LBD equilibrium written with fugacity coefficients
  from the EoS, K = x2 Φ2 / (x1 Φ1), instead of Novak, Liang & Kontogeorgis (JCP 2024),
  who assume an ideal solution and fit ΔH, ΔS, ΔP.
- Semester I: pure compounds with PC-SAFT **and** PC-SAFT-TST; NIST data for water.

## How we work

- Replies to Jorge in English, casual and direct, even when he writes in Spanish.
  Math in LaTeX. Short texts; more explicit steps when the topic is new to him.
- Code: Claude writes it, Jorge reviews it. Simple style: no rcParams, no config
  dictionaries, options written where they are used. Scripts have a small CONFIG block.
- Git: GitHub `jorgerosas005/pcsaft` is the source of truth. Claude pushes from the
  cloud clone; Jorge runs `git pull` in his OneDrive folder. Claude never runs git in
  OneDrive. Always pull before editing. Say what is in each commit. Run `pytest` first.
- Commits: author "Jorge Rosas <jorgerosas005@gmail.com>" plus the Co-Authored-By line.
- Never report a number that was not computed or read from a source.

## Layout

    pcsaft/       the engine: eos.py, equilibrium.py, caloric.py, data.py, ...
    tests/        pytest, 22 tests
    scripts/      01-03 validation (NIST), 04 Aspen pseudo-root, 05 two methods,
                  06 density, 07 residual analysis, 08 Wagner P0, 09 Psat residual shape
    exercises/    learning scripts (compare_diego.py, ex.py, plot_residuals.py)
    docs/         saturation.tex / .pdf (how saturation() works, with flowchart)
    data/nist/    saturation data, 11 n-alkanes
    EXERCISES.md  theory exercises      DRILLS.md  usage drills

## What is established

- PC-SAFT pure, non-associating: validated on 11 n-alkanes (AAD Psat 0.09-1.8 %),
  caloric (Δh_vap 1-4 %, Cp-Cv 8-14 %, model error not a bug).
- Aspen Plus, same parameters (ethane): Psat identical to every digit shown, densities
  to ~4 ppm. Where no liquid root exists Aspen returns a pseudo-root inside the unstable
  region, continuous with the real branch at P_min (see scripts/04).
- saturation(): brentq on ln φ_L − ln φ_V inside the spinodal window
  max(P_min, 0) < P < P_max; inner densities bracketed by the spinodal too.
- Density residuals from the paper's starts (η = 0.5 / 1e-10), Newton, 7 alkanes:
  paper's P_EoS/P − 1 is best (100 % success, cheapest). Milton's 1 − Z_PV/Z_EoS fails
  the liquid ~94 % (flat at 0.5, pole where Z_EoS = 0 just below the liquid root).
  Any Z form is ~6× slower for the vapour. See scripts/07.
- Psat with Wagner as P0 + Newton on P, P ← P(1 + g/(Z_V − Z_L)): 3-5 iterations,
  3-10× cheaper than the spinodal route. Spinodal still needed as fallback near Tc and
  between real Tc and EoS Tc (305.3-308.9 K for ethane).
- Milton's textbook (Equilibrio de Fases y Químico, 2nd ed., ch. 14, pp. 190-194):
  cubic EoS, analytic roots, "valor inicial" + Solver for Psat, Wagner as reference.
  Water Wagner (1.5-3-6 form) in table A.3.

## Open items

1. Fix `density()`: near P_max it can return the liquid when asked for vapour
   (ethane 200 K, P = (1−1e-6) P_max). Use spinodal brackets when P is in the window;
   add a test for that case.
2. Push the fixed `pcsaft/plotting.py` (repo still forces Agg, so plt.show() is dead;
   save() lacks show=True). Restore the "Making a figure" section of EXERCISES.md.
3. Reorganize folders (agreed): src/pcsaft/, tests/, scripts/, data/, docs/,
   learning/ (EXERCISES, DRILLS, exercises/), sandbox/ and outputs/ gitignored.
   Rename scripts/test_playgroud.py (pytest would collect it).
4. scripts/08: needs Aspen WAGNER25 constants for ethane (8 elements, unit of ln Pc).
   Then add the Newton-on-P column.
5. Ask Diego: which solver, and how he gets the liquid root with the ratio residual.
6. Ask Milton: TPT1 first then TST? Equations/notes for PC-SAFT-TST? Where does the
   Z-ratio residual come from?
7. TPT1 association in `a_assoc` (2B and 4C), after reading Gross & Sadowski 2002.
   Validate with water; exercise 10 should go from ~42 bar to ~1 bar.
8. To read: Segura & Wisniak 1997 (spinodal initialisation of Psat);
   Mathias, Boston & Watanasiri 1984 (pseudo-roots, verify reference).
9. Question for Pizio (Stat. Thermo): can theory beyond TPT1 (rings, cooperativity,
   bond angles) make the tetrahedral species less dense, unlike TPT1?
10. Unanswered from earlier: "problema 4 — entropía negativa del nitrógeno"
    (ask which course: Termo Química or Termo Estadística).
