#!/usr/bin/env python3
"""Force, energy and shape of the tube against the position of the bead, the way out and the way back (Create_subjob_tube_material.py, Run_tube_material_chain.py).

usage:  python3 Analyze_tube_material.py [--out tube_material] [label=]RESULTS_DIR [[label=]RESULTS_DIR ...]      (matplotlib for the figure: micromamba run -n mir_mem)

RESULTS_DIR is the batch folder of a chain (../Results/Tube_material/KA_16.0_KB_4.0_r_0.3): chain_summary.csv and the run folders link000_x0.2/1/, ... Give one
RESULTS_DIR per (KA, KB) to overlay them; label=DIR names the series.

Per run (the waypoints between two positions of the list are left out, --waypoints keeps them): x (the bead), E_total, E_bend, E_tension, E_bead (energy unit of the
config), Fx (the force on the bead along x, last row of Bead_0_data.txt: negative = the membrane pulls the bead back), x_tip (largest x of the mesh), d_bead (distance of the
bead centre to the mesh, about its radius when the bead rests on the membrane), the radius of the tube (median distance from the x axis of the vertices between 1 and x_tip - 1,
only for x_tip > 2.5) and the area of the membrane. capped = the relaxation hit its step cap (not converged): marked in the table and in the figure. The theory of the infinite
tube (kappa = KB/2, sigma = KA): force f0 = 2 pi sqrt(KB KA), radius R0 = sqrt(KB/KA)/2, energy slope f0; the figure draws them. The cluster scan of the old method
(Results/Pulling_and_relaxing_2/Force_data.csv, one bead pulled by a spring to every x: dir KA r x Fx Fy Fz) is drawn as grey crosses for the same KA and radius when the file
exists (--old_scan other.csv, --old_scan none). Writes <out>.csv and <out>.png (the force, the energy, the radius and the area against x; circles: the way out, squares: the way
back) and prints the hysteresis (the way back minus the way out at the same x).
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
from Run_tube_material_chain import last_mesh, read_obj  # noqa: E402

OLD_SCAN_NAME = os.path.join("Results", "Pulling_and_relaxing_2", "Force_data.csv")
OLD_SCAN = next((c for c in (os.path.join("..", OLD_SCAN_NAME), os.path.join(os.path.dirname(os.path.realpath(__file__)), "..", OLD_SCAN_NAME)) if os.path.exists(c)),
                os.path.join("..", OLD_SCAN_NAME))


def tube_radius(V, x_tip):
    if x_tip < 2.5:
        return float("nan")
    sel = (V[:, 0] > 1.0) & (V[:, 0] < x_tip - 1.0)
    if sel.sum() < 5:
        return float("nan")
    return float(np.median(np.hypot(V[sel, 1], V[sel, 2])))


def load(root, waypoints=False):
    path = os.path.join(root, "chain_summary.csv")
    if not os.path.exists(path):
        raise SystemExit("no chain_summary.csv in %s" % root)
    rows = [r for r in csv.DictReader(open(path)) if r.get("x") not in (None, "")]
    if not waypoints:
        # the waypoints (kind waypoint) and the half move of a retry (mid) are relaxations on the way, not positions of the list
        rows = [r for r in rows if r.get("kind", "point") == "point"]
    if not rows:
        raise SystemExit("no finished positions in %s" % path)
    cfg = None
    pts = []
    for r in rows:
        run = os.path.join(root, r["rundir"]) if r.get("rundir") else None
        if run is None or not os.path.isdir(run):
            print("skipped %s: no run folder" % r["run"], file=sys.stderr)
            continue
        if cfg is None and os.path.exists(os.path.join(run, "Input_file.json")):
            cfg = json.load(open(os.path.join(run, "Input_file.json")))
        V, _ = read_obj(last_mesh(run))
        p = {k: (v if k in ("run", "rundir", "leg", "kind") else (float(v) if v not in (None, "") else float("nan"))) for k, v in r.items()}
        p["R_tube"] = tube_radius(V, p["x_tip"])
        p["path"] = run
        pts.append(p)
    if not pts or cfg is None:
        raise SystemExit("no usable runs (with an Input_file.json) in %s" % root)
    en = {e["Name"]: e["constants"] for e in cfg["Energies"]}
    KA, KB = en["Surface_tension"][0], [en[n][0] for n in ("Bending", "Bending_tan") if n in en][0]
    return pts, KA, KB, cfg["Beads"][0]["radius"]


def old_scan(path, KA, radius):
    if path in (None, "none"):
        return None
    if not os.path.exists(path):
        print("note: the old scan %s does not exist, not drawn" % path)
        return None
    out = []
    for line in open(path):
        parts = line.split()
        if len(parts) < 5 or parts[0] in ("dir", "#"):
            continue
        try:
            if abs(float(parts[1]) - KA) < 1e-9 and abs(float(parts[2]) - radius) < 1e-9:
                out.append((float(parts[3]), float(parts[4])))
        except ValueError:
            continue
    return sorted(out) or None


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("series", nargs="+", help="[label=]results folder of a chain")
    ap.add_argument("--out", default="tube_material", help="output prefix (<out>.csv, <out>.png)")
    ap.add_argument("--waypoints", action="store_true", help="keep the waypoint runs: by default only the positions of the list are points")
    ap.add_argument("--old_scan", default=OLD_SCAN, help="Force_data.csv of the old method (dir KA r x Fx Fy Fz), 'none' to leave it out")
    args = ap.parse_args()

    data = {}
    for item in args.series:
        head, eq, tail = item.partition("=")
        label, path = (head, tail) if eq and "/" not in head else (os.path.basename(os.path.normpath(item)), item)
        if label in data:
            raise SystemExit("two series are called %r: give them labels, label=folder" % label)
        pts, KA, KB, radius = load(path, args.waypoints)
        data[label] = (pts, KA, KB, radius)
        skipped = [r for r in csv.DictReader(open(os.path.join(path, "chain_summary.csv"))) if r.get("kind") in ("waypoint", "mid") and str(r.get("capped")) == "1"]
        if skipped and not args.waypoints:
            print("note: %d waypoint runs of %s hit their step cap (not converged); --waypoints shows them" % (len(skipped), label))
        f0, R0 = 2.0 * math.pi * math.sqrt(KB * KA), 0.5 * math.sqrt(KB / KA)
        print("\n%s  (KA %g, KB %g, bead radius %g; theory of the tube: force %.3f, radius %.3f)" % (label, KA, KB, radius, f0, R0))
        print("  leg    x      E_total    E_bend  E_tension   E_bead     Fx      x_tip  d_bead  R_tube    Area     steps  seconds  run")
        for q in pts:
            print("  %-4s %6.2f %10.3f %9.3f %9.3f %9.3f %9.3f %7.3f %7.3f %7.3f %9.3f %7d %7d%s  %s" % (
                q["leg"], q["x"], q["E_total"], q["E_bend"], q["E_tension"], q["E_bead"], q["Fx"], q["x_tip"], q["d_bead"], q["R_tube"], q["Area"], q["step"], q["seconds"],
                " CAPPED" if q.get("capped") else "       ", q["run"]))
        # hysteresis: the way back against the way out at the same x
        out = {round(q["x"], 6): q for q in pts if q["leg"] == "in"}
        back = {round(q["x"], 6): q for q in pts if q["leg"] == "out"}
        common = sorted(set(out) & set(back))
        if common:
            dF = [back[x]["Fx"] - out[x]["Fx"] for x in common]
            dE = [back[x]["E_total"] - out[x]["E_total"] for x in common]
            i = int(np.argmax(np.abs(dF)))
            trap = getattr(np, "trapezoid", None) or np.trapz
            print("  hysteresis (way back - way out): largest |dFx| %.3f at x = %g, largest |dE_total| %.3f at x = %g; the loop area int dFx dx = %.3f" % (
                abs(dF[i]), common[i], max(map(abs, dE)), common[int(np.argmax(np.abs(dE)))], trap(dF, common)))
    allrows = [(l, q) for l, (pts, *_rest) in data.items() for q in pts]
    fields = ["series", "link", "kind", "leg", "x", "E_total", "E_bend", "E_tension", "E_bead", "Fx", "Fy", "Fz", "x_tip", "d_bead", "R_tube", "Area", "nV", "step", "seconds", "capped", "run"]
    if args.out.endswith("/") or os.path.isdir(args.out):
        args.out = os.path.join(args.out, "tube_material")
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    with open(args.out + ".csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(fields)
        for l, q in allrows:
            w.writerow([l] + [q.get(k, "") for k in fields[1:]])
    print("\nwrote %s.csv" % args.out)

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib is not installed here: the table and the csv are written, no png (the mir_mem or mir_membranes environment has it)")
        return
    fig, ax = plt.subplots(2, 2, figsize=(11.0, 8.0))
    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    any_capped = False
    for si, (label, (pts, KA, KB, radius)) in enumerate(data.items()):
        c = colors[si % len(colors)]
        f0, R0 = 2.0 * math.pi * math.sqrt(KB * KA), 0.5 * math.sqrt(KB / KA)
        E0 = pts[0]["E_total"]
        A0 = pts[0]["Area"]
        for leg, marker, alpha in (("in", "o", 1.0), ("out", "s", 0.6)):
            sel = [q for q in pts if q["leg"] == leg]
            if not sel:
                continue
            x = [q["x"] for q in sel]
            name = "%s (KA %g, KB %g)%s" % (label, KA, KB, "" if leg == "in" else ", way back")
            ax[0, 0].plot(x, [q["Fx"] for q in sel], marker + "-", color=c, alpha=alpha, ms=3, lw=0.9, label=name)
            ax[0, 1].plot(x, [q["E_total"] - E0 for q in sel], marker + "-", color=c, alpha=alpha, ms=3, lw=0.9, label=name)
            ax[1, 0].plot(x, [q["R_tube"] for q in sel], marker + "-", color=c, alpha=alpha, ms=3, lw=0.9)
            ax[1, 1].plot(x, [q["Area"] - A0 for q in sel], marker + "-", color=c, alpha=alpha, ms=3, lw=0.9)
        capped = [q for q in pts if q.get("capped")]
        if capped:
            ax[0, 0].plot([q["x"] for q in capped], [q["Fx"] for q in capped], "x", color="red", ms=8, label="hit the step cap (not converged)" if not any_capped else None)
            any_capped = True
        xs = np.array(sorted({q["x"] for q in pts if q["x"] >= pts[0]["x"] - 1e-9}))         # the theory lines from the first position of the way out on
        ax[0, 0].axhline(-f0, color=c, ls=":", lw=0.9)
        ax[0, 1].plot(xs, f0 * (xs - xs[0]), color=c, ls=":", lw=0.9)
        ax[1, 0].axhline(R0, color=c, ls=":", lw=0.9)
        ax[1, 1].plot(xs, 2.0 * math.pi * R0 * (xs - xs[0]), color=c, ls=":", lw=0.9)
        old = old_scan(args.old_scan, KA, radius)
        if old:
            ax[0, 0].plot([o[0] for o in old], [o[1] for o in old], "+", color="0.45", ms=6, label="old method, KA %g (spring pulled to each x)" % KA)
    ax[0, 0].set_ylabel("force on the bead, x (dotted: -f0)", fontsize=9)
    ax[0, 1].set_ylabel("E_total - E_total(x0) (dotted: f0 (x - x0))", fontsize=9)
    ax[1, 0].set_ylabel("radius of the tube (dotted: R0)", fontsize=9)
    ax[1, 1].set_ylabel("area - area(x0) (dotted: 2 pi R0 (x - x0))", fontsize=9)
    for a_ in ax.flat:
        a_.set_xlabel("x of the bead")
        a_.grid(alpha=0.25, lw=0.5)
    ax[0, 0].legend(fontsize=6.5)
    fig.suptitle("bead through the planar membrane, circles: way out, squares: way back", fontsize=10)
    fig.tight_layout()
    fig.savefig(args.out + ".png", dpi=150)
    print("wrote %s.png" % args.out)


if __name__ == "__main__":
    main()
