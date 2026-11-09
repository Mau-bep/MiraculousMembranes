#!/usr/bin/env python3
"""The interaction potential of two beads on a vesicle against the distance between the bead centres, at coverage 0.9 and the three surface tensions of WrappingTwoVesicle, in the style of the
manuscript (figstyle.py); png and eps to Results/Paper_figures/two_vesicle_<energy>_<leg>_vNN.*.

usage:  python3 Figures_two_vesicle.py                       the six figures: bending, tension and total energy, each for the way in and the way out
        python3 Figures_two_vesicle.py bending total         only these energies (bending, tension, total)
        python3 Figures_two_vesicle.py --legs in             only the way in (the beads approach), out (they move apart)
        python3 Figures_two_vesicle.py --cols 1.5            another width (columns); --height_cm another height; --coverage 0.8 another coverage

The data are the chains TwoInteraction_r_0.5_coverage_<c>_tension_<t>_replica_0_SurfaceTension of Results/WrappingTwoVesicle, as analysed by Analyze_two_vesicle.py
(two_beads_analysis_SurfaceTension/csv/*.csv, one row per relaxed link: leg in / out, the distance d between the bead centres, the energies in kappa = K_B / 2). u is the energy of a link minus the
energy far away, in the same way as the csv's u_kappa (the mean over the links of the way in beyond d = 4 D_p; the way out has the same zero), for
    bending   the bending energy                       (u_kappa of the csv)
    tension   the surface-tension energy sigma A       (the Surface_tension column of the run's Output_data.txt)
    total     the whole energy: bending + tension + constraints + adhesion  (u_tot_kappa of the csv)
The axis is written in K_B as Mau asked, the numbers are those of the old plots (units of kappa = K_B / 2, not converted). D_p = 2 a is the bead diameter (a = 0.5).
"""
import argparse
import csv
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import figstyle                                    # noqa: E402  the style of the manuscript

RESULTS = os.path.join(os.path.dirname(HERE), "Results")
ANALYSIS = os.path.join(RESULTS, "WrappingTwoVesicle", "two_beads_analysis_SurfaceTension", "csv")
COLS = 1.0                                         # width in columns
TENSIONS = ["0.0", "0.02", "1.0"]                  # the surface tensions (KA of Surface_tension) of the three chains
COLORS = ["#1b4fb8", "#9a66e0", "#e0559c"]         # blue, purple, pink (as the tube figures)
MARKERS = ["o", "s", "^"]
D_FAR = 4.0                                        # D_p: the links of the way in beyond this distance form the zero of u (the zero of the csv: --ref mean_far --s_far 4 --xaxis chord)
ENERGIES = {"bending": "Ebend_kappa", "tension": "Esigma_kappa", "total": "Etot_kappa"}


def chain_csv(tension, coverage):
    return os.path.join(ANALYSIS, "TwoInteraction_r_0.5_coverage_%s_tension%s_replica_0_SurfaceTension.csv" % (coverage, tension))


def sigma_energy_kappa(run_rel):
    """The surface-tension energy of a run in kappa: the Surface_tension column of the last line of its Output_data.txt over kappa = (the constant of Bending) / 2."""
    run = os.path.normpath(os.path.join(HERE, run_rel))
    cfg = json.load(open(os.path.join(run, "Input_file.json")))
    constants = {e["Name"]: e["constants"] for e in cfg["Energies"]}
    kappa = constants["Bending" if "Bending" in constants else "Bending_tan"][0] / 2.0
    with open(os.path.join(run, "Output_data.txt")) as f:
        header = f.readline().split()
        last = None
        for line in f:
            if line.strip():
                last = line.split()
    return float(last[header.index("Surface_tension")]) / kappa if "Surface_tension" in header else 0.0


def load_chain(tension, coverage):
    """The links of one chain: dict of arrays d, leg, and the energies of ENERGIES (kappa) already referred to the far energy (u)."""
    with open(chain_csv(tension, coverage), newline="") as f:
        rows = list(csv.DictReader(f))
    d = np.array([float(r["d"]) for r in rows]) / 1.0                       # D_p = 2 a = 1 (a = 0.5): d is already d / D_p
    leg = np.array([r["leg"] for r in rows])
    E = {"Ebend_kappa": np.array([float(r["Ebend_kappa"]) for r in rows]), "Etot_kappa": np.array([float(r["Etot_kappa"]) for r in rows]),
         "Esigma_kappa": np.array([sigma_energy_kappa(r["run"]) for r in rows])}
    far = (leg == "in") & (d > D_FAR)
    u = {name: E[col] - E[col][far].mean() for name, col in ENERGIES.items()}
    # the bending and the total u of the csv were made in the same way: check
    for name, col in (("bending", "u_kappa"), ("total", "u_tot_kappa")):
        ref = np.array([float(r[col]) for r in rows])
        if np.max(np.abs(u[name] - ref)) > 2e-3:
            print("note: tension %s, %s: my u differs from the csv's %s by up to %.3g kappa (the csv used another zero)" % (tension, name, col, np.max(np.abs(u[name] - ref))))
    return {"d": d, "leg": leg, "u": u, "coverage": [float(r["frac0"]) for r in rows]}


def fig_potential(energy, leg, cols=COLS, height_cm=None, coverage="0.9", tensions=TENSIONS, xlim=(0.8, 6.2)):
    """u(d / D_p) of the three tensions on the way in or the way out: bending, tension or total energy."""
    figstyle.use()
    fig, ax = figstyle.figure(cols=cols, aspect=0.75, height_cm=height_cm)
    ax.axhline(0.0, color="0.55", lw=0.5, zorder=0)
    for t, c, mk in zip(tensions, COLORS, MARKERS):
        ch = load_chain(t, coverage)
        sel = ch["leg"] == leg
        order = np.argsort(ch["d"][sel])
        ax.plot(ch["d"][sel][order], ch["u"][energy][sel][order], ls="-", lw=1.0, color=c, marker=mk, ms=figstyle.MARKER, mew=0,
                label=r"Surface tension $\sigma = %s$" % ("%g" % float(t)), zorder=2)
    ax.set_xlim(*xlim)
    ax.set_xlabel(r"Distance between bead centers $d/D_p$")
    ax.set_ylabel(r"Interaction potential $u\,[K_B]$")
    ax.legend(loc="lower right", handlelength=1.6, handletextpad=0.4)
    return fig


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("energies", nargs="*", help="bending, tension, total (default: all three)")
    ap.add_argument("--legs", default="in,out", help="in (the beads approach), out (they move apart), or both (default)")
    ap.add_argument("--coverage", default="0.9", help="the coverage of the chains (default 0.9)")
    ap.add_argument("--cols", type=float, default=COLS, help="width in columns (default 1)")
    ap.add_argument("--height_cm", type=float, default=None, help="height in cm (default 0.75 of the width)")
    args = ap.parse_args()
    bad = [e for e in args.energies if e not in ENERGIES]
    if bad:
        raise SystemExit("energies are bending, tension, total, not %s" % bad)
    import matplotlib
    matplotlib.use("Agg")
    for energy in args.energies or list(ENERGIES):
        for leg in args.legs.split(","):
            fig = fig_potential(energy, leg, cols=args.cols, height_cm=args.height_cm, coverage=args.coverage)
            figstyle.save_figure(fig, "two_vesicle_%s_%s" % (energy, "way_" + leg))


if __name__ == "__main__":
    main()
