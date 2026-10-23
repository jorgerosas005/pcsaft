"""Figures, with one style in one place.

Deliberately no rcParams and no style dictionary: every option is written
literally where it is used, so a figure can be read top to bottom without
chasing a global setting. The repetition is avoided by reusing these three
functions, not by hiding the settings.

save() writes both a .pdf and a .png: the PDF goes into LaTeX as vector art,
the PNG is for a quick look or a slide.
"""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

FIG_DIR = Path(__file__).resolve().parent.parent / "figures"

REF_STYLE = dict(color="black", linewidth=1.2, zorder=1)
CALC_STYLE = dict(linestyle="none", marker="o", markersize=4,
                  markerfacecolor="none", color="crimson", zorder=2)


def figure(nrows=1, ncols=1, width=6.0, height=4.0, **kw):
    """A figure sized for a one-column LaTeX layout."""
    fig, ax = plt.subplots(nrows, ncols, figsize=(width * ncols, height * nrows), **kw)
    return fig, ax


def finish(ax, xlabel=None, ylabel=None, legend=False):
    """Axis cosmetics applied the same way everywhere."""
    if xlabel:
        ax.set_xlabel(xlabel)
    if ylabel:
        ax.set_ylabel(ylabel)
    ax.tick_params(direction="in", top=True, right=True)
    if legend:
        ax.legend(frameon=False)


def save(fig, name, tight=True):
    """Write figures/<name>.pdf and figures/<name>.png."""
    if tight:
        fig.tight_layout()
    FIG_DIR.mkdir(exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(FIG_DIR / f"{name}.{ext}", dpi=200)
    plt.close(fig)
    return FIG_DIR / f"{name}.pdf"


# --------------------------------------------------- the recurring figures

def plot_saturation(comparison, comp, name=None, every=2):
    """Psat, the coexistence dome, and the deviations. Three panels."""
    df = comparison
    fig, ax = figure(1, 3, width=4.3, height=3.8)

    ax[0].semilogy(1000 / df["T"], df["P_ref"] / 1e5, label="reference", **REF_STYLE)
    ax[0].semilogy(1000 / df["T"][::every], df["P_calc"][::every] / 1e5,
                   label="PC-SAFT", **CALC_STYLE)
    finish(ax[0], r"1000/$T$  /  K$^{-1}$", r"$P^{\mathrm{sat}}$ / bar", legend=True)

    ax[1].plot(df["rhoL_ref"], df["T"], **REF_STYLE)
    ax[1].plot(df["rhoV_ref"], df["T"], label="reference", **REF_STYLE)
    ax[1].plot(df["rhoL_calc"][::every], df["T"][::every], **CALC_STYLE)
    ax[1].plot(df["rhoV_calc"][::every], df["T"][::every], label="PC-SAFT", **CALC_STYLE)
    finish(ax[1], r"$\rho$ / kg m$^{-3}$", r"$T$ / K", legend=True)

    ax[2].axhline(0, color="black", linewidth=0.6)
    ax[2].plot(df["T"], df["dev_P"], marker="o", markersize=3, label=r"$P^{\mathrm{sat}}$")
    ax[2].plot(df["T"], df["dev_rhoL"], marker="s", markersize=3, label=r"$\rho_L^{\mathrm{sat}}$")
    ax[2].set_ylim(-10, 10)
    finish(ax[2], r"$T$ / K", "deviation / %", legend=True)

    return save(fig, name or f"{comp.name}_saturation")


def plot_property(x, y_ref, y_calc, xlabel, ylabel, name, every=2, logy=False):
    """One property against a reference, with the deviation underneath."""
    fig, ax = figure(2, 1, width=6.0, height=3.4, sharex=True)

    draw = ax[0].semilogy if logy else ax[0].plot
    draw(x, y_ref, label="reference", **REF_STYLE)
    draw(x[::every], y_calc[::every], label="PC-SAFT", **CALC_STYLE)
    finish(ax[0], None, ylabel, legend=True)

    dev = 100 * (y_calc - y_ref) / y_ref
    ax[1].axhline(0, color="black", linewidth=0.6)
    ax[1].plot(x, dev, marker="o", markersize=3, color="crimson")
    finish(ax[1], xlabel, "deviation / %")

    return save(fig, name)
