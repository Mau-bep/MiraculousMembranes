#!/usr/bin/env python3
"""Figures of the pilot (README section 5) from the small csv files in data/ (written by Scripts/Analyze_two_vesicle.py).

usage: micromamba run -n mir_mem python make_figures.py       (writes figures/u_of_s_pilot.png and figures/u_of_s_round2.png)
u is the change of the TOTAL energy (bending + tension + constraints + adhesion) in kappa, zero = mean of the runs beyond 4.3 D_p.
"""
import csv
import os

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))


def read(name):
    path = os.path.join(HERE, "data", name + ".csv")
    if not os.path.exists(path):
        return []
    rows = list(csv.DictReader(open(path)))
    out = []
    for r in rows:
        out.append({"s": float(r["s_over_Dp"]), "u_tot": float(r["u_tot_kappa"]), "u_bend": float(r["u_kappa"]),
                    "overlap": r["overlap"] == "True", "cov": 0.5 * (float(r["cov0"]) + float(r["cov1"])), "d": float(r["d"]),
                    "chord": float(r["s_chord"]) if r.get("s_chord") else float(r["d"]) / 2.0})
    return sorted(out, key=lambda p: p["s"])


def paper():
    rows = {}
    for line in open(os.path.join(HERE, "paper_fig4a_digitized.csv")):
        if line.startswith("#") or not line.strip():
            continue
        name, s, u = line.strip().split(",")
        rows.setdefault(name, []).append((float(s), float(u)))
    return {k: np.array(sorted(v)) for k, v in rows.items()}


def draw(ax, names, labels, colors, key="u_tot"):
    for name, label, color in zip(names, labels, colors):
        pts = read(name)
        good = [p for p in pts if not p["overlap"]]
        bad = [p for p in pts if p["overlap"]]
        ax.plot([p["s"] for p in good], [p[key] for p in good], "o-", ms=3.5, lw=1.0, color=color, label=label)
        if bad:
            ax.plot([p["s"] for p in bad], [p[key] for p in bad], "x", ms=6, mew=1.3, color=color)


def main():
    pap = paper()
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.8))
    for ax in axes:
        if "experiment_floppy" in pap:
            ax.plot(pap["experiment_floppy"][:, 0], pap["experiment_floppy"][:, 1], "o", color="tab:red", ms=3, alpha=0.4, label="experiment, floppy (paper Fig. 4a)")
        if "simulation" in pap:
            ax.plot(pap["simulation"][:, 0], pap["simulation"][:, 1], "D", color="k", ms=3.5, label="paper simulation (Fig. 4a)")
        ax.axhspan(-0.03, 0.03, color="0.85", zorder=0, label="restart noise $\\pm$0.03 $\\kappa$")
        ax.axhline(0, color="0.5", lw=0.6)
        ax.set_xlabel(r"$s\ [D_p]$  (arc of a sphere of radius $R_v = 6.25\,a$)")
        ax.set_ylabel(r"$u_{tot}\ [\kappa]$")
    draw(axes[0], ["floppy_90", "area5_90", "pipeline_check_floppy_90"],
         ["floppy ($\\tilde\\sigma$=0.02), 90 % wrapped", "area 5 %, 90 % (taut, $\\tilde\\sigma\\approx$1.1)", "default pipeline run, floppy 90 %"],
         ["tab:blue", "tab:green", "tab:cyan"])
    axes[0].set_xlim(1.2, 4.4)
    axes[0].set_ylim(-0.3, 1.1)
    axes[0].set_title("90 % wrapped (x: beads overlap)", fontsize=10)
    axes[0].legend(fontsize=7, loc="upper right")
    draw(axes[1], ["floppy_50", "floppy_40", "floppy_30", "floppy_60", "area5_60"],
         ["floppy 50 %", "floppy 40 %", "floppy 30 %", "floppy 60 %", "area 5 %, 60 %"],
         ["tab:orange", "tab:purple", "tab:brown", "tab:red", "tab:olive"])
    axes[1].set_xlim(1.2, 4.4)
    axes[1].set_ylim(-1.5, 8.0)
    axes[1].set_title("shallower wrapping (noise 0.1-0.5 $\\kappa$)", fontsize=10)
    axes[1].legend(fontsize=7, loc="upper right")
    fig.tight_layout()
    os.makedirs(os.path.join(HERE, "figures"), exist_ok=True)
    out = os.path.join(HERE, "figures", "u_of_s_pilot.png")
    fig.savefig(out, dpi=150)
    print("wrote", out)


def round2():
    """Second round (README section 6): beads inside, theta_Wr = 98 %, planar disk. u_tot against s (arc of the vesicle, or the distance d/D_p on the disk)."""
    pap = paper()
    fig, axes = plt.subplots(1, 3, figsize=(15.5, 4.6))
    panels = [
        ("beads outside, 98 % wrapped (budded)", "s [$D_p$] arc of the vesicle",
         [("out_floppy98", "floppy ($\\tilde\\sigma$=0.02)", "tab:blue"), ("out_area5_98", "area 5 % (taut, $\\tilde\\sigma\\approx$1.9)", "tab:green"),
          ("floppy_90", "floppy, 90 % (pilot)", "tab:cyan")]),
        ("beads INSIDE, 98 % wrapped", "s [$D_p$] arc of the vesicle",
         [("in_floppy98", "floppy ($\\tilde\\sigma$=0.02)", "tab:blue"), ("in_area5_98", "paper's ensemble: area 5 % ($\\tilde\\sigma\\approx$0.05)", "tab:red"),
          ("in_floppy98_rep1", "floppy, replica 1", "tab:purple"), ("in_floppy90", "floppy, 90 %", "tab:orange"),
          ("tense_in_98", "tense ($\\tilde\\sigma$=3, W=20)", "tab:gray")]),
        ("flat disk (R = 24 a), floppy", "d / $D_p$ distance of the bead centres",
         [("planar_floppy60", "60 %", "tab:brown"), ("planar_floppy75", "75 %", "tab:olive"), ("planar_floppy90", "90 %, inward", "tab:blue"),
          ("planar_floppy90_out", "90 %, outward", "tab:cyan"), ("planar_floppy98", "98 %", "tab:green")]),
    ]
    for ax, (title, xlabel, series) in zip(axes, panels):
        if "experiment_floppy" in pap:
            ax.plot(pap["experiment_floppy"][:, 0], pap["experiment_floppy"][:, 1], "o", color="tab:red", ms=3, alpha=0.3, label="experiment, floppy (Fig. 4a)")
        if "simulation" in pap:
            ax.plot(pap["simulation"][:, 0], pap["simulation"][:, 1], "D", color="k", ms=3.5, label="paper simulation (Fig. 4a)")
        ax.axhspan(-0.03, 0.03, color="0.85", zorder=0)
        ax.axhline(0, color="0.5", lw=0.6)
        for name, label, color in series:
            pts = read(name)
            if not pts:
                continue
            good = [p for p in pts if not p["overlap"]]
            ax.plot([p["s"] for p in good], [p["u_tot"] for p in good], "o-", ms=3, lw=0.9, color=color, label=label)
        ax.set_xlim(0.8, 4.4)
        ax.set_ylim(-0.3, 0.9)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(r"$u_{tot}\ [\kappa]$")
        ax.set_title(title, fontsize=10)
        ax.legend(fontsize=6.5, loc="upper right")
    fig.tight_layout()
    out = os.path.join(HERE, "figures", "u_of_s_round2.png")
    fig.savefig(out, dpi=150)
    print("wrote", out)


if __name__ == "__main__":
    main()
    round2()
