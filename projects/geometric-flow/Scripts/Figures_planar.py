#!/usr/bin/env python3
"""The phase space of the planar membrane in the style of the manuscript (figstyle.py), written as png and eps to Results/Paper_figures/<name>_vNN.*.

usage:  python3 Figures_planar.py                    the phase space, one column wide
        python3 Figures_planar.py --cols 1.5         another width (columns); --height_cm 6 another height
        python3 Figures_planar.py --no_theory        without the lines of Deserno's Fig. 2

x = reduced adhesion  wbar = 4 K_I r^2 / K_B,  y = reduced surface tension  sigmabar = 2 K_A r^2 / K_B  (r the bead radius; kappa = K_B / 2, w = K_I, sigma = K_A: Plots.py, planar_phase_deserno),
the colour the coverage (CoverageUnion of PostProcessing). Every point takes the run of lower Total_E of the flat start (Coverage_data.txt) and of the full, wrapped start (Coverage_data_full.txt).
The lines: W (wbar = 4, dashed), the zero-energy line of the enveloped state (wbar = 4 + 2 sigmabar, dotted), E (solid) and the spinodals S1, S2 (dashed) of Deserno (deserno_theory.py).
"""
import argparse
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import figstyle                                    # noqa: E402  the style of the manuscript

RESULTS = os.path.join(os.path.dirname(HERE), "Results")
COLS = 1.0                                         # width in columns
FLAT = os.path.join(RESULTS, "Wrapping_planar_excess", "Coverage_data.txt")
FULL = os.path.join(RESULTS, "Wrapping_planar_excess", "Coverage_data_full.txt")


def load_phase(path):
    """x = wbar = 4 KI r^2/KB, y = sigmabar = 2 KA r^2/KB, the coverage and Total_E of the runs of one Coverage_data.txt of PostProcessing (columns by the header names)."""
    with open(path) as f:
        names = f.readline().lstrip("#").split()
    col = {n: i for i, n in enumerate(names[1:])}                  # the DIR column is a string: usecols skips it
    data = np.loadtxt(path, comments="#", usecols=range(1, len(names)), ndmin=2)
    KA, KB, KI, r = (data[:, col[n]] for n in ("KA", "KB", "KI", "BeadRadius"))
    return np.round(4 * KI * r ** 2 / KB, 6), np.round(2 * KA * r ** 2 / KB, 6), data[:, col["CoverageUnion"]], data[:, col["Total_E"]]


def edges(c):
    """Cell edges halfway between the centres of a sorted grid, so that every cell is centred on its run."""
    if len(c) == 1:
        return np.array([c[0] - 0.5, c[0] + 0.5])
    mid = 0.5 * (c[1:] + c[:-1])
    return np.concatenate(([c[0] - (mid[0] - c[0])], mid, [c[-1] + (c[-1] - mid[-1])]))


def fig_phase_space(cols=COLS, height_cm=None, theory=True, xlim=(2.0, 10.0), ylim=(0.0, 2.0), flat=FLAT, full=FULL):
    """Coverage over the (wbar, sigmabar) plane, the lower-energy run of the flat and the full start at every point, with the lines of Deserno's Fig. 2."""
    wa, sa, ca, Ea = load_phase(flat)
    wb, sb, cb, Eb = load_phase(full)
    ws, ss = np.unique(np.concatenate((wa, wb))), np.unique(np.concatenate((sa, sb)))

    def grid(w, s, v):
        G = np.full((len(ss), len(ws)), np.nan)
        G[np.searchsorted(ss, s), np.searchsorted(ws, w)] = v
        return G

    Ea_g, Eb_g = grid(wa, sa, Ea), grid(wb, sb, Eb)
    full_wins = np.where(np.isnan(Eb_g), np.inf, Eb_g) < np.where(np.isnan(Ea_g), np.inf, Ea_g)
    Z = np.where(full_wins, grid(wb, sb, cb), grid(wa, sa, ca))

    figstyle.use()
    fig, ax = figstyle.figure(cols=cols, aspect=0.8, height_cm=height_cm)
    mesh = ax.pcolormesh(edges(ws), edges(ss), Z, cmap="viridis", shading="flat", rasterized=True)
    cbar = fig.colorbar(mesh, ax=ax, pad=0.02, aspect=22)
    cbar.set_label("Coverage")
    cbar.outline.set_linewidth(figstyle.AXES_LW)
    cbar.ax.tick_params(width=figstyle.AXES_LW, direction="out")
    if theory:
        ax.axvline(4.0, ls="--", color="black", lw=1.0)                                       # W
        sig = np.linspace(0.0, ylim[1], 50)
        ax.plot(4.0 + 2.0 * sig, sig, ls=":", color="black", lw=1.0)                          # the enveloped state has zero energy
        import deserno_theory
        deserno_theory.plot_fig2_lines(ax, ylim=ylim, lw_E=1.6, lw_S=1.0)                    # E, S1, S2
    ax.set_xlabel(figstyle.label("Reduced adhesion", r"\bar{w} = 4K_Ir^2/K_B"))
    ax.set_ylabel(figstyle.label("Reduced surface tension", r"\bar{\sigma}"))
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.xaxis.set_major_locator(figstyle.mpl.ticker.MultipleLocator(2))
    ax.yaxis.set_major_locator(figstyle.mpl.ticker.MultipleLocator(0.5))
    ax.set_aspect("auto")
    return fig


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cols", type=float, default=COLS, help="width in columns (default 1)")
    ap.add_argument("--height_cm", type=float, default=None, help="height in cm (default: 0.8 of the width)")
    ap.add_argument("--no_theory", action="store_true", help="without the lines of Deserno's Fig. 2")
    args = ap.parse_args()
    import matplotlib
    matplotlib.use("Agg")
    fig = fig_phase_space(cols=args.cols, height_cm=args.height_cm, theory=not args.no_theory)
    figstyle.save_figure(fig, "planar_phase_space")


if __name__ == "__main__":
    main()
