"""Loading reference data.

What the NIST Webbook serves is not raw experimental points: it is a
multiparameter reference equation of state fitted to the experimental database,
accurate to roughly experimental uncertainty (~0.05% in liquid density). That is
far below the ~1% error of PC-SAFT, so it is a legitimate stand-in for
experiment when the question is "is my implementation correct".

It is NOT the right source for regressing parameters. Fitting to a correlation
means inheriting its smoothing and its bias with no honest way to report
uncertainty. Parameter work uses NIST-TDE or the primary literature, with
references and stated uncertainties.
"""

from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = ROOT / "data"
RESULTS_DIR = ROOT / "results"


def save_results(df, name):
    """Write a DataFrame to results/<name>.csv, creating the folder if needed.

    Paths are resolved from the package location, not the working directory,
    so a script lands its output in the same place whether it was launched
    from the project root, from an editor, or from anywhere else.
    """
    RESULTS_DIR.mkdir(exist_ok=True)
    path = RESULTS_DIR / f"{name}.csv"
    df.to_csv(path, index=False)
    return path


def _find(name):
    """Locate a data file by component name, or accept a path as given."""
    path = Path(name)
    if path.exists():
        return path
    stem = str(name).strip().lower()
    for ext in (".tsv", ".txt"):
        cand = DATA_DIR / "nist" / f"{stem}_sat{ext}"
        if cand.exists():
            return cand
    raise FileNotFoundError(
        f"No NIST file for {name!r}. Expected data/nist/{stem}_sat.tsv or .txt")


# The full Webbook export has 25 columns; the trimmed files keep 8. Both are
# read into the same tidy frame, so nothing downstream cares which it got.
_ALIASES = {
    "T": ["Temperature (K)", "T (K)"],
    "P": ["Pressure (bar)", "P (bar)"],
    "rho_L": ["Density (l, kg/m3)", "rho_L (kg/m3)"],
    "rho_V": ["Density (v, kg/m3)", "rho_V (kg/m3)"],
    "h_L": ["Enthalpy (l, kJ/mol)", "h_L (kJ/mol)"],
    "h_V": ["Enthalpy (v, kJ/mol)", "h_V (kJ/mol)"],
    "Cv_L": ["Cv (l, J/mol*K)", "Cv_L (J/mol*K)"],
    "Cp_L": ["Cp (l, J/mol*K)", "Cp_L (J/mol*K)"],
}


def load_nist_saturation(name):
    """Read a NIST Webbook 'saturation properties' export.

    Takes a component name (looked up under data/nist/) or a path. Returns a
    DataFrame with columns T (K), P (Pa), rho_L, rho_V (kg/m3), and where the
    file has them h_L, h_V (kJ/mol), Cv_L, Cp_L (J/mol/K), dh_vap, and
    Cp_minus_Cv_L.
    """
    raw = pd.read_csv(_find(name), sep="\t")

    out = {}
    for key, options in _ALIASES.items():
        for col in options:
            if col in raw.columns:
                out[key] = raw[col].values
                break
    if "T" not in out or "rho_L" not in out:
        raise ValueError(f"Unrecognised column layout in {name!r}: {list(raw.columns)}")

    df = pd.DataFrame(out)
    df["P"] = df["P"] * 1e5                       # bar -> Pa
    if "h_L" in df:
        df["dh_vap"] = df["h_V"] - df["h_L"]      # ideal part cancels exactly
    if "Cp_L" in df and "Cv_L" in df:
        df["Cp_minus_Cv_L"] = df["Cp_L"] - df["Cv_L"]
    return df.dropna(subset=["T", "P", "rho_L", "rho_V"]).reset_index(drop=True)


def check_saturation_data(df):
    """Sanity-check a saturation table before trusting it.

    Returns a dict of check name -> result. Everything here is cheap and
    catches the usual import mistakes: wrong unit column, a truncated download,
    columns swapped between liquid and vapour.
    """
    import numpy as np

    T, P = df["T"].values, df["P"].values
    rl, rv = df["rho_L"].values, df["rho_V"].values

    checks = {
        "no missing values": bool(df.notna().all().all()),
        "P strictly increasing": bool(np.all(np.diff(P) > 0)),
        "rho_L decreasing": bool(np.all(np.diff(rl) < 0)),
        "rho_V increasing": bool(np.all(np.diff(rv) > 0)),
        "rho_L > rho_V": bool(np.all(rl > rv)),
    }
    if P.min() < 1.01325e5 < P.max():
        checks["normal boiling point, K"] = float(np.interp(1.01325e5, P, T))

    # ln P against 1/T should be close to a straight line (Clausius-Clapeyron)
    c = np.polyfit(1 / T, np.log(P), 1)
    resid = np.log(P) - np.polyval(c, 1 / T)
    checks["Clausius-Clapeyron R2"] = float(
        1 - np.sum(resid**2) / np.sum((np.log(P) - np.log(P).mean()) ** 2))
    return checks
