# pcsaft

PC-SAFT (Gross & Sadowski, 2001) for pure components: density and saturation pressure.

Master's project, Posgrado en Ingeniería Química (Procesos), UNAM.
Advisor: Dr. Milton Medeiros, Facultad de Química.

## Install

    py -m pip install -e .

or double-click `install.bat`.

## Use

    py scripts\density.py       liquid and vapour density at (T, P)
    py scripts\saturation.py    Psat and coexisting densities at T

## Method

- **Density:** Newton on $P_{EoS}/P - 1$ in the packing fraction $\eta$, from
  $\eta = 0.5$ (liquid) or $\eta = 10^{-10}$ (vapour).
- **Saturation:** $P_0$ from Wagner 25, then Newton on $P$ with
  $g = \ln\varphi_L - \ln\varphi_V$:

  $$P_{k+1} = P_k\left(1 + \frac{g}{Z_V - Z_L}\right)$$

Wagner 25 constants are in `pcsaft/parameters.py` (provisional, to be replaced by Aspen's WAGNER25).

## Layout

    pcsaft/constants.py     physical and universal PC-SAFT constants
    pcsaft/parameters.py    compound parameters and Wagner 25 constants
    pcsaft/eos.py           the equation of state: a_res, Z, P, ln phi
    pcsaft/equilibrium.py   density() and saturation()
    scripts/                density.py, saturation.py
