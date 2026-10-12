#!/usr/bin/env python3
"""Interaction energy u(s) between two wrapped beads on a closed vesicle (van der Wel et al., Sci. Rep. 6, 32825, Fig. 3e / 4a).

usage:  python3 Analyze_two_vesicle.py [options] [label=]RESULTS_DIR [[label=]RESULTS_DIR ...]
        (label=DIR1,DIR2 puts the runs of several folders, e.g. the inward and the outward leg of a scan, in one series)

Every folder below a RESULTS_DIR that holds an Input_file.json, an Output_data.txt and the bead files is one
point of the curve (one separation). Folders of the same RESULTS_DIR form one series; give each series a label
(label=path) to name it in the legend. Chains written by Run_two_vesicle_chain.py are read the same way.

Per run (from the last saved state, Final_state.obj if the run finished):
    d          distance between the bead centres
    s          arc of a sphere: s = R_v * phi, phi the angle between the two bead centres seen from the volume
               centre of the vesicle, R_v the radius of the vesicle of the template (key "Rv" of the config, else "rescale")
               (not a geodesic on the mesh: the paper's s is the distance along the membrane between the two contact points
               of a nearly spherical vesicle)
    E_bend     Bending_tan, in units of kappa = KB/2 (the code's bending energy is KB * H^2, H = (c1 + c2)/2)
    u(s)       E(s) - E(reference), reference = the run with the largest s of the series (--ref far) or the mean of the runs
               with s above --s_far in D_p (--ref mean_far)
    u in kBT   u * kappa_kbt, with kappa = 20.9 kBT (Sci. Rep. Fig. 4)

Also printed/written: the weighted coverage (the quantity the Coverage energy constrains), the plain solid-angle fraction
of the facing faces in the shell, the area, volume and the area-penalty tension sigma_eff = KA (A - A_bar)/A_bar^2 in
units of kappa/a^2 (sigma~ of the paper's regime: floppy < 0.03, tense > 2.7), and whether the run finished.
Writes <out>.csv and <out>.png (the paper's digitized Fig. 4a, Studies/two_beads_vesicle/paper_fig4a_digitized.csv, is overlaid).
"""
import argparse
import csv
import json
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
PAPER_CSV = os.path.join(HERE, "..", "Studies", "two_beads_vesicle", "paper_fig4a_digitized.csv")


# ----------------------------------------------------------------------------------------------------------------------
# mesh helpers (the shell and solid-angle formulas are those of src/BeadGeometry.cpp)
# ----------------------------------------------------------------------------------------------------------------------
def read_obj(path):
    V, F = [], []
    with open(path) as f:
        for line in f:
            p = line.split()
            if not p:
                continue
            if p[0] == "v":
                V.append([float(x) for x in p[1:4]])
            elif p[0] == "f":
                F.append([int(x.split("/")[0]) - 1 for x in p[1:4]])
    return np.array(V), np.array(F, dtype=int)


def mesh_geometry(V, F):
    a, b, c = V[F[:, 0]], V[F[:, 1]], V[F[:, 2]]
    cr = np.cross(b - a, c - a)
    area = 0.5 * np.linalg.norm(cr, axis=1).sum()
    trip = (a * np.cross(b, c)).sum(1)
    vol = trip.sum() / 6.0
    vol_centre = (((a + b + c) / 4.0) * (trip / 6.0)[:, None]).sum(0) / vol
    return area, vol, vol_centre


def solid_angle(a, b, c):
    ra, rb, rc = np.linalg.norm(a, axis=1), np.linalg.norm(b, axis=1), np.linalg.norm(c, axis=1)
    num = (a * np.cross(b, c)).sum(1)
    den = ra * rb * rc + (a * b).sum(1) * rc + (b * c).sum(1) * ra + (c * a).sum(1) * rb
    return 2.0 * np.arctan2(num, den)


def coverage(V, F, pos, sigma, width=0.25, power=1.0, side=1.0):
    """(weighted coverage, plain facing solid-angle fraction in the shell, solid-angle weighted mean radius).
    side = -1: a bead inside the closed membrane (the faces that count have their outward normal pointing away from it)."""
    a, b, c = V[F[:, 0]] - pos, V[F[:, 1]] - pos, V[F[:, 2]] - pos
    cen = (a + b + c) / 3.0
    r = np.linalg.norm(cen, axis=1)
    n = np.cross(b - a, c - a)
    n /= np.linalg.norm(n, axis=1)[:, None]
    facing = side * ((cen / r[:, None]) * n).sum(1) < 0
    x = (r - sigma) / (width * sigma)
    inside = np.abs(x) < 1.0
    w = np.where(inside, 0.5 * (1.0 + np.cos(np.pi * x)), 0.0) ** power
    om = np.abs(solid_angle(a, b, c))
    sel = facing & inside
    weighted = (w * facing * om).sum() / (4 * math.pi)
    plain = (om * sel).sum() / (4 * math.pi)
    rmean = (om[sel] * r[sel]).sum() / max(om[sel].sum(), 1e-30)
    return weighted, plain, rmean


# ----------------------------------------------------------------------------------------------------------------------
# run folders
# ----------------------------------------------------------------------------------------------------------------------
KEEP_COLD = [False]


def find_runs(root):
    runs = []
    for dirpath, dirnames, filenames in os.walk(root):
        parts = dirpath.rstrip("/").split("/")
        if any("prewrap" in x for x in parts[-2:]) or (not KEEP_COLD[0] and any(x == "link00_cold" for x in parts[-2:])):
            # the Coverage pre-wrap of a chain is not a point of the curve, and neither is its cold start: it stops in a state 0.3-0.5 kappa above
            # what the next link finds from it after a small move (it is only the starting mesh of the chain); --keep_cold keeps it
            dirnames[:] = []
            continue
        if "Input_file.json" in filenames and "Output_data.txt" in filenames:
            runs.append(dirpath)
            dirnames[:] = []
    return sorted(runs)


def last_mesh(run):
    final = os.path.join(run, "Final_state.obj")
    if os.path.exists(final):
        return final, True
    meshes = sorted((f for f in os.listdir(run) if f.startswith("membrane_") and f.endswith(".obj")),
                    key=lambda s: int(s[9:-4]))
    return (os.path.join(run, meshes[-1]), False) if meshes else (None, False)


def read_output(run):
    with open(os.path.join(run, "Output_data.txt")) as f:
        header = f.readline().split()
        rows = [[float(x) for x in l.split()] for l in f if l.strip()]
    col = {}
    for i, name in enumerate(header):
        col.setdefault(name, i)
    col["_bead_columns"] = [i for i, name in enumerate(header) if name == "Bead"]  # the adhesion energy of each bead
    return col, np.array(rows)


def read_bead_positions(run, n=2):
    out = []
    for i in range(n):
        rows = [l.split() for l in open(os.path.join(run, "Bead_%d_data.txt" % i)) if l.strip() and not l.startswith("#")]
        out.append(np.array([float(x) for x in rows[-1][1:4]]))
    return out


def analyze_run(run):
    cfg = json.load(open(os.path.join(run, "Input_file.json")))
    energies = {e["Name"]: e["constants"] for e in cfg["Energies"]}
    kappa = energies["Bending_tan"][0] / 2.0
    a = cfg["Beads"][0]["radius"]
    Rv = cfg.get("Rv", cfg.get("rescale", 1.0))
    col, data = read_output(run)
    last = data[-1]
    out = {"run": run, "kappa": kappa, "a": a, "Rv": Rv, "step": last[col["step"]]}
    for name in ("Volume", "Area", "Bending_tan", "Volume_constraint", "Area_constraint", "Surface_tension",
                 "Excess_tension", "Coverage", "Total_E"):
        out[name] = last[col[name]] if name in col else 0.0
    mesh, finished = last_mesh(run)
    out["finished"] = finished
    V, F = read_obj(mesh)
    pos = read_bead_positions(run)
    out["nV"] = len(V)
    out["d"] = float(np.linalg.norm(pos[0] - pos[1]))
    # the beads have no excluded volume: below 2 a they overlap (deep pockets: the arc s = R_v phi reaches 2 a at 1.8 D_p, not at 1 D_p)
    out["overlap"] = out["d"] < 2.0 * a * 0.999
    out["s_chord"] = out["d"] / (2.0 * a)  # the straight distance between the bead centres in D_p
    if cfg.get("boundary"):
        # planar disk with a clamped edge (--planar): no volume centre, s is the distance between the bead centres; Rc is the height of the bead centre above the plane
        area = 0.5 * np.linalg.norm(np.cross(V[F[:, 1]] - V[F[:, 0]], V[F[:, 2]] - V[F[:, 0]]), axis=1).sum()
        out["phi"] = float("nan")
        out["s_over_Dp"] = out["s_chord"]
        out["Rc0"], out["Rc1"] = (float(p[0]) for p in pos)
    else:
        area, vol, vc = mesh_geometry(V, F)
        u = [(p - vc) / np.linalg.norm(p - vc) for p in pos]
        out["phi"] = math.acos(max(-1.0, min(1.0, float(u[0] @ u[1]))))
        out["s_over_Dp"] = Rv * out["phi"] / (2.0 * a)
        out["Rc0"], out["Rc1"] = (float(np.linalg.norm(p - vc)) for p in pos)
    for i, p in enumerate(pos):
        w, plain, rmean = coverage(V, F, p, cfg["Beads"][i]["radius"], side=-1.0 if cfg["Beads"][i].get("outside", 1) < 0 else 1.0)
        out["cov%d" % i], out["cap%d" % i], out["rmean%d" % i] = w, plain, rmean
    out["Ebead_kappa"] = sum(last[i] for i in col["_bead_columns"]) / kappa
    out["Ebend_kappa"] = out["Bending_tan"] / kappa
    out["Esigma_kappa"] = (out["Surface_tension"] + out["Excess_tension"]) / kappa
    out["Etot_kappa"] = out["Total_E"] / kappa
    out["sigma_eff"] = 0.0
    if "Area_constraint" in energies:
        # E_A = KA (A - A_bar)^2 / (2 A_bar^2) and the tension is dE_A/dA = KA (A - A_bar)/A_bar^2, so
        # |sigma| = KA sqrt(2 E_A/KA)/A_bar (the sign is lost, A_bar is taken as the mesh area)
        KA = energies["Area_constraint"][0]
        out["sigma_eff"] = KA * math.sqrt(2.0 * out["Area_constraint"] / KA) / area / kappa * a * a
    elif "Surface_tension" in energies:
        out["sigma_eff"] = energies["Surface_tension"][0] / kappa * a * a
    elif "Excess_tension" in energies:
        out["sigma_eff"] = energies["Excess_tension"][0] / kappa * a * a
    return out


def reference(points, mode, s_far):
    if mode == "far":
        return max(points, key=lambda p: p["s_over_Dp"])["Ebend_kappa"]
    if mode == "mean_far":
        far = [p["Ebend_kappa"] for p in points if p["s_over_Dp"] > s_far]
        if not far:
            raise SystemExit("--ref mean_far: no run with s > %g D_p" % s_far)
        return float(np.mean(far))
    return float(mode)


def load_paper():
    rows = {}
    if os.path.exists(PAPER_CSV):
        for line in open(PAPER_CSV):
            if line.startswith("#") or not line.strip():
                continue
            name, s, u = line.strip().split(",")
            rows.setdefault(name, []).append((float(s), float(u)))
    return {k: np.array(sorted(v)) for k, v in rows.items()}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("series", nargs="+", help="[label=]results folder")
    ap.add_argument("--ref", default="far", help='energy zero: "far" (largest s), "mean_far", or a number (kappa)')
    ap.add_argument("--s_far", type=float, default=3.5, help="runs with s above this (D_p) form the reference of --ref mean_far")
    ap.add_argument("--kappa_kbt", type=float, default=20.9, help="kappa in kBT (the paper: 20.9)")
    ap.add_argument("--out", default="two_vesicle", help="output prefix (<out>.csv, <out>.png)")
    ap.add_argument("--only_finished", action="store_true", help="skip runs that did not finish")
    ap.add_argument("--keep_cold", action="store_true", help="keep the cold start (link00_cold) of a chain as a point (default: skipped, it is not relaxed like the links)")
    ap.add_argument("--ylim", type=float, nargs=2, default=None, metavar=("LO", "HI"), help="y range of the plot in kappa (default: from the data, at least -0.3 to 0.15)")
    ap.add_argument("--xlim", type=float, nargs=2, default=[0.8, 4.2], metavar=("LO", "HI"), help="x range of the plot in D_p")
    ap.add_argument("--xaxis", choices=["arc", "chord"], default="arc",
                    help="arc (default): s = R_v phi, the arc of the sphere; chord: the straight distance d/D_p between the bead centres. The two differ by the depth of the "
                         "beads: R_c/R_v (0.6 for deep pockets, 1.1 for beads inside a vesicle). The paper's s is a distance between particle centres")
    ap.add_argument("--curves", choices=["both", "total", "bending"], default="both", help="which energy differences the plot shows (default both: bending solid, total dashed)")
    ap.add_argument("--noise", type=float, default=0.03,
                    help="half width (kappa) of the grey band around u = 0 in the plot: the change of E_bend when a converged state is re-meshed "
                         "and relaxed again at the same separation (measured 0.01-0.04 kappa, Studies/two_beads_vesicle/README.md); 0 = no band")
    args = ap.parse_args()
    KEEP_COLD[0] = args.keep_cold

    allrows = []
    series_points = {}
    for item in args.series:
        label, _, path = item.rpartition("=") if "=" in item else (os.path.basename(os.path.normpath(item)), "", item)
        pts = []
        for one_path in path.split(","):  # several folders (the legs of one scan: inward, outward, extension) can form one series: label=a,b,c
            for run in find_runs(one_path):
                try:
                    p = analyze_run(run)
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
        pts.sort(key=lambda p: p["s_over_Dp"])
        e0 = reference(pts, args.ref, args.s_far)
        # u_tot: the change of the whole energy (bending + tension + constraints + adhesion), the u_T of the paper's Fig. S5b. The bending
        # energy alone trades against the adhesion energy (tilting the patch of a bead by 4 degrees moves 0.06 kappa between them and
        # leaves the total unchanged to 0.003), so u_tot is the robust curve; the paper finds u_Curv ~ u_T
        # the zero of u_tot: the total energy of the run with the largest s (--ref far, numbers) or the mean over s > --s_far (--ref mean_far)
        if args.ref == "mean_far":
            e0_tot = float(np.mean([p["Etot_kappa"] for p in pts if p["s_over_Dp"] > args.s_far]))
        else:
            e0_tot = max(pts, key=lambda p: p["s_over_Dp"])["Etot_kappa"]
        for p in pts:
            p["u_kappa"] = p["Ebend_kappa"] - e0
            p["u_kbt"] = p["u_kappa"] * args.kappa_kbt
            p["u_tot_kappa"] = p["Etot_kappa"] - e0_tot
        series_points[label] = pts
        allrows += pts
        print("\n%s  (reference E_bend = %.4f kappa)" % (label, e0))
        print("   d      s/Dp   Ebend/k     u/k   u_tot/k  u(kBT)  cov0   cov1   cap0   cap1  sigma~   finished  nV   run")
        for p in pts:
            print("  %5.2f  %5.2f  %9.4f  %7.4f  %7.4f  %7.3f  %5.3f  %5.3f  %5.3f  %5.3f  %6.3f  %-8s %5d  %s" % (
                p["d"], p["s_over_Dp"], p["Ebend_kappa"], p["u_kappa"], p["u_tot_kappa"], p["u_kbt"], p["cov0"], p["cov1"], p["cap0"], p["cap1"],
                p["sigma_eff"], p["finished"], p["nV"], os.path.relpath(p["run"])) + ("   BEADS OVERLAP (d < 2a)" if p["overlap"] else ""))

    if not allrows:
        raise SystemExit("nothing to analyze")
    fields = ["series", "run", "d", "s_over_Dp", "s_chord", "phi", "Rc0", "Rc1", "Ebend_kappa", "Esigma_kappa", "Ebead_kappa", "Etot_kappa", "u_kappa", "u_kbt", "u_tot_kappa",
              "cov0", "cov1", "cap0", "cap1", "rmean0", "rmean1", "Area", "Volume", "sigma_eff", "nV", "step", "finished", "overlap"]
    with open(args.out + ".csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(fields)
        for p in allrows:
            w.writerow([p[k] for k in fields])

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    paper = load_paper()
    fig, ax = plt.subplots(figsize=(6.4, 5.6))
    if "experiment_floppy" in paper:
        e = paper["experiment_floppy"]
        ax.plot(e[:, 0], e[:, 1], "o", color="tab:red", ms=3.5, alpha=0.45, label="experiment, floppy (paper Fig. 4a)")
    if "simulation" in paper:
        e = paper["simulation"]
        ax.plot(e[:, 0], e[:, 1], "D", color="k", ms=4, label="paper simulation (Fig. 4a)")
    xkey = "s_over_Dp" if args.xaxis == "arc" else "s_chord"
    for label, pts in series_points.items():
        color = None
        if args.curves in ("both", "bending"):
            line, = ax.plot([p[xkey] for p in pts], [p["u_kappa"] for p in pts], "s-", ms=4.5, lw=1.2,
                            label=label + (" (bending)" if args.curves == "both" else ""))
            color = line.get_color()
        if args.curves in ("both", "total"):
            line, = ax.plot([p[xkey] for p in pts], [p["u_tot_kappa"] for p in pts], "o--" if args.curves == "both" else "o-", ms=3.5, lw=0.9,
                            color=color, mfc="none" if args.curves == "both" else None, label=label + (" (total)" if args.curves == "both" else ""))
    for label, pts in series_points.items():
        bad = [p for p in pts if p["overlap"]]
        if bad:
            ax.plot([p[xkey] for p in bad], [p["u_tot_kappa"] for p in bad], "x", color="tab:red", ms=7, mew=1.5, label="beads overlap (d < 2a)" if label == list(series_points)[0] else None)
    ax.axhline(0, color="0.5", lw=0.6)
    if args.noise > 0:
        ax.axhspan(-args.noise, args.noise, color="0.85", zorder=0, label="restart noise $\\pm$%.2g $\\kappa$" % args.noise)
    ax.set_xlabel(r"$s\ [D_p]$  (arc of a sphere of radius $R_v$)" if args.xaxis == "arc" else r"$d/D_p$  (distance between the bead centres)")
    ax.set_ylabel(r"$u\ [\kappa]$")
    ax.set_xlim(*args.xlim)
    if args.ylim:
        lo, hi = args.ylim
    else:
        lo = min(-0.3, min(min(p["u_kappa"], p["u_tot_kappa"]) for p in allrows) * 1.1)
        hi = max(0.15, max(max(p["u_kappa"], p["u_tot_kappa"]) for p in allrows) * 1.1)
    ax.set_ylim(lo, hi)
    sec = ax.secondary_yaxis("right", functions=(lambda y: y * args.kappa_kbt, lambda y: y / args.kappa_kbt))
    sec.set_ylabel(r"$u\ [k_BT]$ for $\kappa=%.1f\,k_BT$" % args.kappa_kbt)
    ax.legend(fontsize=6.5, loc="upper center", bbox_to_anchor=(0.5, -0.14), ncol=2, frameon=False)
    fig.tight_layout()
    fig.savefig(args.out + ".png", dpi=170)
    print("\nwrote %s.csv and %s.png" % (args.out, args.out))


if __name__ == "__main__":
    main()
