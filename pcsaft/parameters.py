"""Pure-component parameters.

Values are Table 2 of Gross & Sadowski (2001), transcribed as printed.
Every record also carries the temperature range of the original regression and
the AAD% the paper reports, so a validation run can be compared with the
published numbers instead of with a guess.

Association parameters (eps_assoc_k, kappa_assoc, scheme) are unused by the
2001 equation and stay None for the non-associating compounds. They are here so
that water and the CTS/two-state work of the thesis slot in without changing
any function signature.
"""

from dataclasses import dataclass, replace

from .constants import NA


@dataclass
class Component:
    name: str
    M: float            # molar mass, g/mol
    m: float            # number of segments per chain, -
    sigma: float        # segment diameter, A
    eps_k: float        # segment energy / k, K

    # association, not used by plain PC-SAFT
    eps_assoc_k: float = None    # association energy / k, K
    kappa_assoc: float = None    # association volume, -
    scheme: str = None           # "2B", "4C", ...

    # provenance, for validation reports
    T_range: tuple = None        # regression range of the source, K
    aad_psat: float = None       # AAD% in Psat reported by the source
    aad_v: float = None          # AAD% in liquid volume reported by the source
    source: str = "Gross & Sadowski 2001, Table 2"

    # Wagner 25 vapour pressure: (C1, C2, C3, C4, ln(Pc/Pa), Tc/K), filled from WAGNER25 below
    wagner25: tuple = None

    @property
    def associating(self):
        return self.eps_assoc_k is not None

    def without_association(self):
        """A copy with the association parameters stripped off.

        For an associating compound the equation of state refuses to run while
        a_assoc is still a stub, which is the correct behaviour: silently
        dropping hydrogen bonding would produce numbers that look fine and are
        badly wrong. This makes dropping it deliberate and visible, so the
        error can be measured on purpose -- see exercise 10.
        """
        return replace(self, name=f"{self.name} (no association)",
                       eps_assoc_k=None, kappa_assoc=None, scheme=None,
                       aad_psat=None, aad_v=None,
                       source=f"{self.source}, association term removed")

    def mass_density(self, rho):
        """molecules/A^3 -> kg/m^3."""
        return rho * 1e30 / NA * self.M / 1000

    def number_density(self, rho_mass):
        """kg/m^3 -> molecules/A^3."""
        return rho_mass * 1000 / self.M * NA / 1e30


# ---------------------------------------------------------------- n-alkanes

LIBRARY = {
    "methane":  Component("methane",  16.043, 1.0000, 3.7039, 150.03, T_range=(97, 300),  aad_psat=0.36, aad_v=0.67),
    "ethane":   Component("ethane",   30.070, 1.6069, 3.5206, 191.42, T_range=(90, 305),  aad_psat=0.30, aad_v=0.57),
    "propane":  Component("propane",  44.096, 2.0020, 3.6184, 208.11, T_range=(85, 523),  aad_psat=1.29, aad_v=0.77),
    "butane":   Component("butane",   58.123, 2.3316, 3.7086, 222.88, T_range=(135, 573), aad_psat=0.75, aad_v=1.59),
    "pentane":  Component("pentane",  72.146, 2.6896, 3.7729, 231.20, T_range=(143, 469), aad_psat=1.45, aad_v=0.78),
    "hexane":   Component("hexane",   86.177, 3.0576, 3.7983, 236.77, T_range=(177, 503), aad_psat=0.31, aad_v=0.76),
    "heptane":  Component("heptane", 100.203, 3.4831, 3.8049, 238.40, T_range=(182, 623), aad_psat=0.34, aad_v=2.10),
    "octane":   Component("octane",  114.231, 3.8176, 3.8373, 242.78, T_range=(216, 569), aad_psat=0.77, aad_v=1.59),
    "nonane":   Component("nonane",  128.250, 4.2079, 3.8448, 244.51, T_range=(219, 595), aad_psat=0.89, aad_v=0.32),
    "decane":   Component("decane",  142.285, 4.6627, 3.8384, 243.87, T_range=(243, 617), aad_psat=0.24, aad_v=1.18),
    "undecane": Component("undecane", 156.312, 4.9082, 3.8893, 248.82, T_range=(247, 639), aad_psat=2.02, aad_v=0.69),
    "dodecane": Component("dodecane", 170.338, 5.3060, 3.8959, 249.21, T_range=(263, 658), aad_psat=2.10, aad_v=0.93),

    # Water. The association parameters are carried but NOT used yet, because
    # a_assoc is still a stub. Treating water with these three parameters alone
    # is physically wrong on purpose -- see exercise 10 in EXERCISES.md, where
    # the size of that error is the motivation for the whole project.
    # CHECK these against Gross & Sadowski, Ind. Eng. Chem. Res. 2002, 41, 5510
    # (the associating-systems paper) before using them for anything real.
    "water": Component("water", 18.015, 1.0656, 3.0007, 366.51,
                       eps_assoc_k=2500.7, kappa_assoc=0.034868, scheme="2B",
                       source="Gross & Sadowski 2002 (VERIFY before use)"),

    # useful non-alkanes
    "carbon dioxide": Component("carbon dioxide", 44.010, 2.0729, 2.7852, 169.21, T_range=(216, 304), aad_psat=2.78, aad_v=2.73),
    "nitrogen":       Component("nitrogen",       28.010, 1.2053, 3.3130,  90.96, T_range=(63, 126),  aad_psat=0.34, aad_v=1.50),
    "benzene":        Component("benzene",        78.114, 2.4653, 3.6478, 287.35, T_range=(278, 562), aad_psat=0.64, aad_v=1.42),
    "cyclohexane":    Component("cyclohexane",    84.147, 2.5303, 3.8499, 278.11, T_range=(279, 553), aad_psat=0.53, aad_v=3.12),
}

# ---------------------------------------------------------------- Wagner 25

# Initial estimate of Psat in equilibrium.saturation(). Same order as Aspen's
# WAGNER25/1..6: C1, C2, C3, C4, ln(Pc / Pa), Tc / K.
# PROVISIONAL: Tc and Pc are NIST values to be checked; C1..C4 are a least-squares
# fit to NIST saturation data (AAD < 0.1 %). Replace with the Aspen values.
WAGNER25 = {
    "methane":  ( -6.02563,   1.27553,  -0.58561,  -1.33956, 15.341393, 190.564),
    "ethane":   ( -6.47597,   1.40506,  -1.13514,  -1.83357, 15.399056, 305.322),
    "propane":  ( -6.75432,   1.51780,  -1.51501,  -2.26582, 15.262712, 369.890),
    "butane":   ( -7.08718,   1.80518,  -2.04055,  -2.22511, 15.149458, 425.125),
    "pentane":  ( -7.33034,   1.86168,  -2.39243,  -2.44170, 15.030423, 469.700),
    "hexane":   ( -7.51273,   1.78936,  -2.56396,  -3.02344, 14.925392, 507.820),
    "heptane":  ( -7.83147,   2.01391,  -3.02631,  -3.31015, 14.822008, 540.130),
    "octane":   ( -7.91418,   1.62396,  -2.61957,  -4.52734, 14.730601, 569.320),
    "nonane":   ( -8.44693,   2.58323,  -4.22813,  -3.43963, 14.640125, 594.550),
    "decane":   ( -8.63638,   2.56779,  -4.48981,  -3.69907, 14.558875, 617.700),
    "dodecane": ( -8.98105,   2.43074,  -4.65153,  -5.30642, 14.412697, 658.100),
}

for name in WAGNER25:
    LIBRARY[name].wagner25 = WAGNER25[name]


# The n-alkane series in order, which is the homologous set the protocol asks
# for. Used by scripts/02_validate_series.py.
N_ALKANES = ["methane", "ethane", "propane", "butane", "pentane", "hexane",
             "heptane", "octane", "nonane", "decane", "undecane", "dodecane"]


def get(name):
    """Look a component up by name, case-insensitively."""
    key = name.strip().lower()
    if key not in LIBRARY:
        raise KeyError(f"{name!r} is not in the library. Available: {sorted(LIBRARY)}")
    return LIBRARY[key]
