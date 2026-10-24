"""Physical constants and the PC-SAFT universal constants.

Universal constants are Table 1 of
Gross, J.; Sadowski, G. Ind. Eng. Chem. Res. 2001, 40, 1244-1260.
They were regressed against n-alkane data and are the same for every compound.
"""

import numpy as np

kB = 1.380649e-23      # Boltzmann constant, J/K
NA = 6.02214076e23     # Avogadro number, 1/mol
R = kB * NA            # gas constant, J/(mol K)

# Working units inside the equation of state:
#   T     [K]
#   sigma [A], d [A]
#   eps_k [K]           (this is epsilon/k, an energy expressed as a temperature)
#   rho   [molecules/A^3]
# Conversions to SI happen only at the edges, in eos.pressure and
# parameters.Component.mass_density.

# Table 1, dispersion term. Row i of a_coef is (a0i, a1i, a2i).
A_COEF = np.array([
    [0.9105631445, -0.3084016918, -0.0906148351],
    [0.6361281449,  0.1860531159,  0.4527842806],
    [2.6861347891, -2.5030047259,  0.5962700728],
    [-26.547362491, 21.419793629, -1.7241829131],
    [97.759208784, -65.255885330, -4.1302112531],
    [-159.59154087, 83.318680481, 13.776631870],
    [91.297774084, -33.746922930, -8.6728470368],
])

B_COEF = np.array([
    [0.7240946941, -0.5755498075,  0.0976883116],
    [2.2382791861,  0.6995095521, -0.2557574982],
    [-4.0025849485, 3.8925673390, -9.1558561530],
    [-21.003576815, -17.215471648, 20.642075974],
    [26.855641363, 192.67226447, -38.804430052],
    [206.55133841, -161.82646165, 93.626774077],
    [-355.60235612, -165.20769346, -29.666905585],
])

# Random close packing of hard spheres. eta cannot physically exceed this,
# and the density solver uses it as the upper bound of the liquid branch.
ETA_MAX = 0.7405
