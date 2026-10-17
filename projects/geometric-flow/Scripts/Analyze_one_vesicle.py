#!/usr/bin/env python3
"""Curves of a one-bead chain (Create_subjob_one_vesicle.py, Run_one_vesicle_chain.py): energies, wrapped fraction and shape against the scanned parameter.

usage:  python3 Analyze_one_vesicle.py [options] [label=]RESULTS_DIR [[label=]RESULTS_DIR ...]
        (label=DIR1,DIR2 puts the runs of several folders in one series; give each replica its own label to compare them)

Every folder below a RESULTS_DIR that holds an Input_file.json, an Output_data.txt and the bead file is one point of the curve. The scanned parameter and its
value come from chain_summary.csv of the RESULTS_DIR (columns scan_var, scan_value, leg: in = the way along the list, out = the way back), or, without it, from
the config of the run (--var says which parameter: f, theta, wt, W, KB, st, sigma, excess, volume or cov). The Coverage pre-wrap and the cold start of a chain are skipped
(--keep_cold keeps the cold start): they are only the starting mesh of the chain.

Per run (from the last saved state, Final_state.obj if the run finished), energies in units of kappa = KB/2:
    E_bend, E_adh (the adhesion energy of the bead), E_total, E_el = E_total - E_adh (what the membrane pays: bending, tension, constraints)
    cap_cfg    the solid-angle fraction of the bead covered by membrane faces facing it, in the shell of the interaction of the config: the wrapped fraction the
               bead sees (Adhesion: r within shell_width of the radius; Frenkel_Normal_nopush: r within rc - a of the radius, the well of the potential)
    cap25/cov25 the same in a shell of +-25 percent (plain / cosine weighted: the quantity the Coverage energy constrains; includes the neck)
    f_patch    the active area of the patch, (1 - cos theta_Wr)/2
    R_c        distance of the bead centre from the volume centre of the vesicle (the height above the disk for --planar): the depth of the pocket
    nu         reduced volume 6 sqrt(pi) V / A^(3/2) of the membrane; sigma_eff the tension (the multiplier of the area constraint, or the Surface_tension), kappa/a^2
    W, wt      the adhesion strength of the bead in the energy unit of the config per a^2, and the dimensionless w~ = 4 W a^2/KB of the Deserno mapping
    G, f_ad    the coverage the bead energy measures: G = -E_bead/(2 pi W r0^2) runs from 0 to 2 (a closed sphere), f_ad = G/2 is the wrapped fraction (the z_ad of the
               Deserno comparisons, and the G of Studies/deserno_shell/theory/softshell_direct.py). r0 = a for Adhesion, where G = (1/2 pi) sum_f w(r_f) |omega_f|; for the
               Frenkel potential r0 = r_min, the distance of its minimum (1.055 a for rc = 1.2): a membrane lying in the minimum has f_ad = its solid angle fraction, one
               off the minimum has less, so f_ad is a lower bound there and the solid angle in the well (cap_cfg) is the wrapped fraction of the plot
    frac       the wrapped fraction of the coverage plot: cap_cfg for the Frenkel potential, f_ad for Adhesion
Writes <out>.csv, <out>.png (four panels: the energies, the wrapped fraction, the elastic energy against the wrapped fraction, tension and depth) and, for a scan of the
adhesion (wt or W), <out>_coverage.png: the wrapped fraction against the interaction strength with the sharp-contact transition of Deserno for comparison (w~ = 4, or
w~ = 4 (a/r_min)^2 = 3.59 for the Frenkel potential of rc = 1.2, whose membrane sits at r_min) and, for the Adhesion cosine shell of half width 0.25, the zero-tension
soft-shell theory of the planar disk (Studies/deserno_shell/data/softshell_sigma0.csv). There is no continuum theory for the Frenkel potential.
"""
import argparse
import csv
import json
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import Analyze_two_vesicle as A2  # noqa: E402  the mesh helpers (solid angle, coverage in a shell, run folders) are those of the two-bead analysis


def frenkel_rmin(sigma, rc):
    """Distance of the minimum of the Frenkel potential V = alpha (sigma^2/r^2 - 1)(rc^2/r^2 - 1)^2 (zero at sigma, zero value and slope at rc, depth 1)."""
    return math.sqrt(3.0 * sigma ** 2 * rc ** 2 / (sigma ** 2 + 2.0 * rc ** 2))


def read_summary(root):
    path = os.path.join(root, "chain_summary.csv")
    out = {}
    if os.path.exists(path):
        for r in csv.DictReader(open(path)):
            out[r["run"]] = r
    return out


def value_from_cfg(cfg, var):
    b = cfg["Beads"][0]
    energies = {e["Name"]: e["constants"] for e in cfg["Energies"]}
    Rv = cfg["Rv"]
    if var in ("theta", "f"):
        theta = b.get("patch_angle")
        if theta is None:
            return float("nan")
        return theta if var == "theta" else (1.0 - math.cos(math.radians(theta))) / 2.0
    KB = [energies[n][0] for n in ("Bending", "Bending_tan") if n in energies][0]
    if var == "KB":
        return KB
    if var == "W":
        return b["inter_str"]
    if var == "wt":
        return 4.0 * b["inter_str"] * b["radius"] ** 2 / KB
    if var in ("sigma", "st"):
        for name in ("Surface_tension", "Excess_tension"):
            if name in energies:
                return energies[name][0] if var == "sigma" else 2.0 * energies[name][0] * b["radius"] ** 2 / KB
        return float("nan")
    if var == "cov":
        return energies["Coverage"][1] if "Coverage" in energies else float("nan")
    if var == "excess":
        if "Area_constraint" in energies:
            return energies["Area_constraint"][1] / (4.0 * math.pi * Rv ** 2) - 1.0
        return float("nan")
    if var == "volume":
        if "Volume_constraint" in energies:
            return energies["Volume_constraint"][1] / (4.0 / 3.0 * math.pi * Rv ** 3)
        return float("nan")
    raise SystemExit("unknown parameter " + var)


def analyze_run(run, summary_row, var_hint):
    cfg = json.load(open(os.path.join(run, "Input_file.json")))
    energies = {e["Name"]: e["constants"] for e in cfg["Energies"]}
    bend_name = "Bending" if "Bending" in energies else "Bending_tan"
    kappa = energies[bend_name][0] / 2.0
    bead = cfg["Beads"][0]
    a = bead["radius"]
    side = -1.0 if bead.get("outside", 1) < 0 else 1.0
    planar = bool(cfg.get("boundary"))
    col, data = A2.read_output(run)
    last = data[-1]
    mesh, finished = A2.last_mesh(run)
    V, F = A2.read_obj(mesh)
    p = A2.read_bead_positions(run, 1)[0]
    area, vol, vc = A2.mesh_geometry(V, F)

    def energy(name):
        return last[col[name]] / kappa if name in col else 0.0

    E_adh_code = sum(last[i] for i in col["_bead_columns"])
    E_adh = E_adh_code / kappa
    E_total = energy("Total_E")
    KB = energies[bend_name][0]
    W = bead.get("inter_str", 0.0)
    out = {"run": run, "finished": finished, "nV": len(V), "step": last[col["step"]], "kappa": kappa,
           "E_bend": energy(bend_name), "E_adh": E_adh, "E_total": E_total, "E_el": E_total - E_adh,
           "E_tension": energy("Surface_tension") + energy("Excess_tension"), "E_constraints": energy("Volume_constraint") + energy("Area_constraint"),
           "Area": area, "Volume": float("nan") if planar else vol,
           "nu": float("nan") if planar else 6.0 * math.sqrt(math.pi) * vol / area ** 1.5}
    mem_inter = bead.get("mem_inter", "Adhesion")
    frenkel = mem_inter in ("Frenkel_Normal_nopush", "Frenkel")
    if frenkel:
        rc_cfg = bead.get("rc", 2.0 * a)            # the cutoff, an absolute length; sigma of the potential is the bead radius
        shell, power, r0 = rc_cfg / a - 1.0, 1.0, frenkel_rmin(a, rc_cfg)
    else:
        rc_cfg = float("nan")
        shell, power, r0 = bead.get("shell_width", 0.25), bead.get("shell_power", 1.0), a
    w25, plain25, _ = A2.coverage(V, F, p, a, side=side)
    wcfg, plaincfg, rmean = A2.coverage(V, F, p, a, width=shell, power=power, side=side)
    out.update({"cov25": w25, "cap25": plain25, "cov_cfg": wcfg, "cap_cfg": plaincfg, "rmean": rmean})
    out["W"] = W
    out["wt"] = 4.0 * W * a * a / KB if W > 0 else float("nan")
    out["wt_eff"] = 4.0 * W * r0 * r0 / KB if W > 0 else float("nan")       # with the radius of contact r0: the w~ of Deserno's geometry
    out["G"] = -E_adh_code / (2.0 * math.pi * W * r0 * r0) if W > 0 else float("nan")
    out["f_ad"] = out["G"] / 2.0
    out["potential"] = "frenkel" if frenkel else ("adhesion" if mem_inter == "Adhesion" else mem_inter)
    out["rmin"], out["rc"] = r0, rc_cfg
    out["frac"] = out["cap_cfg"] if frenkel else out["f_ad"]
    out["KB"] = KB
    out["shell_width"], out["shell_power"] = shell, power
    theta = bead.get("patch_angle")
    out["f_patch"] = (1.0 - math.cos(math.radians(theta))) / 2.0 if theta is not None else float("nan")
    out["Rc"] = float(p[0]) if planar else float(np.linalg.norm(p - vc))
    out["sigma_eff"] = 0.0
    if "Area_constraint" in energies:
        # E_A = KA (A - A_bar)^2 / (2 A_bar^2) and the tension is dE_A/dA = KA (A - A_bar)/A_bar^2, so |sigma| = KA sqrt(2 E_A/KA)/A_bar (the sign is lost)
        KA = energies["Area_constraint"][0]
        out["sigma_eff"] = KA * math.sqrt(max(2.0 * last[col["Area_constraint"]] / KA, 0.0)) / area / kappa * a * a
    else:
        for name in ("Surface_tension", "Excess_tension"):
            if name in energies:
                out["sigma_eff"] = energies[name][0] / kappa * a * a
    if summary_row is not None:
        out["scan_var"], out["scan_value"], out["leg"], out["link"] = summary_row["scan_var"], float(summary_row["scan_value"]), summary_row["leg"], int(summary_row["link"])
    else:
        out["scan_var"], out["scan_value"], out["leg"], out["link"] = var_hint, value_from_cfg(cfg, var_hint), "?", -1
    return out


def read_theory():
    """The zero-tension soft-shell theory of the planar disk (cosine shell s = 0.25, p = 1, R = 10): coverage G(w~), Studies/deserno_shell/data/softshell_sigma0.csv."""
    path = os.path.join(HERE, "..", "Studies", "deserno_shell", "data", "softshell_sigma0.csv")
    if not os.path.exists(path):
        return None
    rows = [l.split() for l in open(path) if l.strip() and not l.startswith("#")]
    return np.array([float(r[0]) for r in rows]), np.array([float(r[1]) for r in rows])


def coverage_figure(series_points, path, plt):
    """Wrapped fraction against the interaction strength, the way along the list and the way back, with the sharp-contact transition of Deserno (and the theory of the
    planar soft shell for the cosine shell of Adhesion)."""
    fig, ax = plt.subplots(figsize=(7.4, 5.4))
    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    allq = [q for pts in series_points.values() for q in pts]
    KB, a = allq[0]["KB"], 1.0
    frenkel = all(q["potential"] == "frenkel" for q in allq)
    if frenkel:
        # the membrane lies in the minimum of the well, r_min from the bead centre: Deserno's geometry has that radius, w~ = 4 W r_min^2/KB = 4 at the sharp transition
        rmin = allq[0]["rmin"]
        ax.axvline(4.0 / rmin ** 2, color="0.5", ls=":", lw=1.0,
                   label=r"Deserno, sharp contact at $r_{min}=%.3f\,a$, $\tilde\sigma=0$: $\tilde w=%.2f$" % (rmin, 4.0 / rmin ** 2))
        print("the Frenkel potential has no continuum theory: the coverage figure shows the sharp-contact line of Deserno only")
    else:
        theory = read_theory()
        soft = all(abs(q["shell_width"] - 0.25) < 1e-9 and abs(q["shell_power"] - 1.0) < 1e-9 for q in allq)
        if theory is not None and soft:
            ax.plot(theory[1], theory[0] / 2.0, "k--", lw=1.1, label=r"planar disk, soft-shell theory ($\tilde\sigma=0$, $s=0.25$)")
        elif theory is not None:
            print("the shell of the runs is not the cosine shell of half width 0.25: no theory curve in the coverage figure")
        ax.axvline(4.0, color="0.5", ls=":", lw=1.0, label=r"Deserno, sharp contact, $\tilde\sigma=0$: $\tilde w=4$")
    ax.axhspan(0.5, 0.6, color="0.9", zorder=0, label="50-60 % wrapped")
    main_name, alt_name = (("solid angle in the well", "bead energy"), ("adhesion energy", "solid angle in the shell"))[0 if frenkel else 1]
    alt_key = "f_ad" if frenkel else "cap_cfg"
    for si, (label, pts) in enumerate(series_points.items()):
        c = colors[si % len(colors)]
        for leg, marker, alpha in (("in", "o", 1.0), ("out", "s", 0.6), ("?", "o", 1.0)):
            sel = sorted([q for q in pts if q["leg"] == leg and not math.isnan(q["wt"])], key=lambda q: q["wt"])
            if not sel:
                continue
            name = label + (" (way back)" if leg == "out" else "")
            ax.plot([q["wt"] for q in sel], [q["frac"] for q in sel], marker + "-", color=c, alpha=alpha, ms=4, lw=1.0, label=name + ": " + main_name)
            ax.plot([q["wt"] for q in sel], [q[alt_key] for q in sel], marker, color=c, alpha=0.5 * alpha, ms=3, mfc="none", lw=0, label=name + ": " + alt_name)
    ax.set_xlabel(r"interaction strength $\tilde w = 4Wa^2/K_B = 2Wa^2/\kappa$ ($W$: depth of the potential per area)")
    if frenkel:
        ax.set_ylabel(r"wrapped fraction of the bead: solid angle covered within the well, $r<r_c$")
    else:
        ax.set_ylabel(r"wrapped fraction of the bead, $f=G/2$  ($G=-E_{adh}/(2\pi W a^2)$)")
    ax.set_ylim(-0.02, 1.02)
    ax.grid(alpha=0.25, lw=0.5)
    sec = ax.secondary_xaxis("top", functions=(lambda x: x * KB / (4.0 * a * a), lambda x: x * 4.0 * a * a / KB))
    sec.set_xlabel(r"$W$ [energy unit of the config per $a^2$; $k_BT/a^2$ for $K_B=42$]")
    if not frenkel:
        sec2 = ax.secondary_yaxis("right", functions=(lambda y: 2.0 * y, lambda y: y / 2.0))
        sec2.set_ylabel(r"$G$ (2 = a closed sphere)")
    ax.legend(fontsize=6.5, loc="lower right")
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    print("wrote %s" % path)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("series", nargs="+", help="[label=]results folder")
    ap.add_argument("--var", default=None, help="the scanned parameter when a folder has no chain_summary.csv: f, theta, wt, W, KB, st, sigma, excess, volume or cov")
    ap.add_argument("--out", default="one_vesicle", help="output prefix (<out>.csv, <out>.png)")
    ap.add_argument("--ref", default="none", help='energy zero of the energy panel: "none" (absolute), "first" or "last" (the point with the first / last value of the way along the list), or a number (kappa)')
    ap.add_argument("--only_finished", action="store_true", help="skip runs that did not finish")
    ap.add_argument("--keep_cold", action="store_true", help="keep the cold start (link00_cold) of a chain as a point")
    ap.add_argument("--waypoints", action="store_true", help="keep the waypoint runs (link03_f0.5_w1, ...): by default only the links of the list are points")
    args = ap.parse_args()
    A2.KEEP_COLD[0] = args.keep_cold

    series_points = {}
    for item in args.series:
        label, _, path = item.rpartition("=") if "=" in item else (os.path.basename(os.path.normpath(item)), "", item)
        pts = []
        for one_path in path.split(","):
            summary = read_summary(one_path)
            for run in A2.find_runs(one_path):
                name = os.path.basename(os.path.dirname(run.rstrip("/")))
                if not args.waypoints and ("_w" in name.split("_", 1)[-1] or name.endswith("_mid")):
                    continue
                try:
                    if name not in summary and not args.var:
                        raise RuntimeError("not in chain_summary.csv and no --var")
                    p = analyze_run(run, summary.get(name), args.var)
                except Exception as exc:  # unfinished or odd folder
                    print("skipped %s: %s" % (run, exc), file=sys.stderr)
                    continue
                if args.only_finished and not p["finished"]:
                    continue
                p["series"] = label
                pts.append(p)
        if not pts:
            print("no runs under %s" % path, file=sys.stderr)
            continue
        pts.sort(key=lambda q: (q["leg"] != "in", q["link"], q["scan_value"]))
        along = [q for q in pts if q["leg"] in ("in", "cold", "?")] or pts
        if args.ref == "first":
            ref = along[0]
        elif args.ref == "last":
            ref = along[-1]
        else:
            ref = None
        for q in pts:
            for key in ("E_bend", "E_adh", "E_total", "E_el"):
                q["d" + key] = q[key] - (ref[key] if ref is not None else (float(args.ref) if args.ref not in ("none", "first", "last") and key == "E_total" else 0.0))
        series_points[label] = pts
        var = pts[0]["scan_var"]
        print("\n%s  (scan %s%s)" % (label, var, "" if args.ref == "none" else ", energies relative to %s" % ("the %s point" % args.ref if ref is not None else args.ref)))
        print("  leg  %9s  E_total    E_bend      E_adh      E_el     cap_cfg  cap25  f_patch   R_c    nu    sigma~  finished  nV   run" % var)
        for q in pts:
            print("  %-4s %9.4g %9.3f %9.3f %10.3f %9.3f   %5.3f  %5.3f  %6.3f  %6.3f %5.3f %7.3f  %-8s %5d  %s" % (
                q["leg"], q["scan_value"], q["dE_total"], q["dE_bend"], q["dE_adh"], q["dE_el"], q["cap_cfg"], q["cap25"], q["f_patch"], q["Rc"], q["nu"],
                q["sigma_eff"], q["finished"], q["nV"], os.path.relpath(q["run"])))

    if not series_points:
        raise SystemExit("nothing to analyze")
    allrows = [q for pts in series_points.values() for q in pts]
    fields = ["series", "run", "link", "leg", "scan_var", "scan_value", "E_total", "E_bend", "E_adh", "E_el", "E_tension", "E_constraints", "f_patch", "cap_cfg", "cov_cfg", "cap25", "cov25", "rmean", "W", "wt", "wt_eff", "G", "f_ad", "frac", "potential", "rmin", "rc", "KB", "shell_width", "shell_power",
              "Rc", "Area", "Volume", "nu", "sigma_eff", "nV", "step", "finished"]
    with open(args.out + ".csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(fields)
        for q in allrows:
            w.writerow([q[k] for k in fields])
    print("\nwrote %s.csv" % args.out)

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib is not installed here: the table and the csv are written, no png (the mir_mem or mir_membranes environment has it)")
        return
    var = allrows[0]["scan_var"]
    shell_name = "well of the potential" if allrows[0]["potential"] == "frenkel" else "shell of the adhesion"
    if var in ("wt", "W"):
        coverage_figure(series_points, args.out + "_coverage.png", plt)
    fig, ax = plt.subplots(2, 2, figsize=(11.0, 8.2))
    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    energy_style = (("dE_total", "total", "-"), ("dE_bend", "bending", "--"), ("dE_adh", "adhesion", ":"), ("dE_el", "total - adhesion", "-."))
    for si, (label, pts) in enumerate(series_points.items()):
        c = colors[si % len(colors)]
        for leg, marker, alpha in (("in", "o", 1.0), ("out", "s", 0.55), ("?", "o", 1.0), ("cold", "D", 1.0)):
            sel = [q for q in pts if q["leg"] == leg]
            if not sel:
                continue
            x = [q["scan_value"] for q in sel]
            for key, name, ls in energy_style:
                ax[0, 0].plot(x, [q[key] for q in sel], marker=marker, ls=ls, color=c, alpha=alpha, ms=3.5, lw=1.0, label=("%s %s" % (label, name)) if leg in ("in", "?") else None)
            ax[0, 1].plot(x, [q["cap_cfg"] for q in sel], marker=marker, ls="-", color=c, alpha=alpha, ms=3.5, lw=1.0, label=("%s, %s" % (label, shell_name)) if leg in ("in", "?") else None)
            ax[0, 1].plot(x, [q["cap25"] for q in sel], marker=marker, ls="--", color=c, alpha=0.5 * alpha, ms=3, lw=0.8, label=("%s, shell +-25 percent" % label) if leg in ("in", "?") else None)
            if not all(math.isnan(q["f_patch"]) for q in sel):
                ax[0, 1].plot(x, [q["f_patch"] for q in sel], ls=":", color="k", lw=0.8, label="active area of the patch" if si == 0 and leg in ("in", "?") else None)
            ax[1, 0].plot([q["cap_cfg"] for q in sel], [q["dE_el"] for q in sel], marker=marker, ls="-", color=c, alpha=alpha, ms=3.5, lw=1.0, label=("%s, total - adhesion" % label) if leg in ("in", "?") else None)
            ax[1, 0].plot([q["cap_cfg"] for q in sel], [q["dE_bend"] for q in sel], marker=marker, ls="--", color=c, alpha=0.6 * alpha, ms=3, lw=0.8, label=("%s, bending" % label) if leg in ("in", "?") else None)
            ax[1, 1].plot(x, [q["sigma_eff"] for q in sel], marker=marker, ls="-", color=c, alpha=alpha, ms=3.5, lw=1.0, label=("%s, tension" % label) if leg in ("in", "?") else None)
    ax2 = ax[1, 1].twinx()
    for si, (label, pts) in enumerate(series_points.items()):
        sel = [q for q in pts if q["leg"] in ("in", "?")]
        ax2.plot([q["scan_value"] for q in sel], [q["Rc"] for q in sel], marker="^", ls=":", color=colors[si % len(colors)], ms=3.5, lw=0.9)
    ax2.set_ylabel(r"$R_c$: bead centre to volume centre (dotted) [$a$]")
    ax[0, 0].set_ylabel(r"energy [$\kappa$]" + ("" if args.ref == "none" else " relative to the reference"))
    ax[0, 1].set_ylabel("wrapped fraction of the bead")
    ax[1, 0].set_xlabel("wrapped fraction (solid angle in the %s)" % shell_name)
    ax[1, 0].set_ylabel(r"energy [$\kappa$]")
    ax[1, 1].set_ylabel(r"tension $\sigma a^2/\kappa$")
    for a_ in (ax[0, 0], ax[0, 1], ax[1, 1]):
        a_.set_xlabel(var)
    if args.var == "sigma" or var in ("sigma", "W", "KB"):
        for a_ in (ax[0, 0], ax[0, 1], ax[1, 1]):
            if all(q["scan_value"] > 0 for q in allrows) and max(q["scan_value"] for q in allrows) / min(q["scan_value"] for q in allrows) > 20 and var in ("sigma", "KB"):
                a_.set_xscale("log")
    for a_ in (ax[0, 0], ax[0, 1], ax[1, 0], ax[1, 1]):
        a_.grid(alpha=0.25, lw=0.5)
    ax[0, 0].legend(fontsize=6, ncol=2)
    ax[0, 1].legend(fontsize=6)
    ax[1, 0].legend(fontsize=6)
    fig.suptitle("one bead on a vesicle, scan of %s (circles: the way along the list, squares: the way back)" % var, fontsize=10)
    fig.tight_layout()
    fig.savefig(args.out + ".png", dpi=150)
    print("wrote %s.png" % args.out)


if __name__ == "__main__":
    main()
