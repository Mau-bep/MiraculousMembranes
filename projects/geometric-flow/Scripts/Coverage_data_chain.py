#!/usr/bin/env python3
"""The coverage of the points of the one-bead chains of a batch, one csv per chain: <batch>/Coverage_data_<Tag>.csv

chain_summary.csv has no coverage (it is measured on the meshes), so this reads the points of every chain folder <batch>/<Tag>/ (the links of the list, both legs;
not the waypoints _w1 nor the half step _mid of a retry; the same points as Collect_chain_results.sh), measures the coverage on the Final_state.obj of each and writes

    <batch>/Coverage_data_<Tag>.csv
        DIR            the run folder, relative to the batch folder
        link, leg      the link of the chain and the leg: cold (the cold start), in (the way along the list), out (the way back)
        scan_var, scan_value
        KA, KB, KI, BeadRadius, rc     the constants of the run (KA = Surface_tension / Excess_tension, KB = Bending, KI = inter_str of the bead): those of the
                                       Coverage_data.txt of PostProcessing, so that the plots of Plots.py (deserno_fig5) read them
        CoveredArea    the coverage of the Coverage energy: sum over the faces in the shell of weight * solid angle / 4 pi (PostProcessing: CoveredArea)
        CoverageUnion  the fraction of 20000 directions from the bead centre (Fibonacci lattice) hit by at least one face of the shell (PostProcessing: CoverageUnion)
        MultilayerFrac the fraction hit by two or more (membrane folded over itself in the shell)
        cap_cfg        the plain solid-angle fraction of the faces in the shell (Analyze_one_vesicle: the wrapped fraction of its coverage plot for the Frenkel potential)
        f_ad           the coverage the adhesion energy of the bead measures, -E_bead/(4 pi W r0^2) (Analyze_one_vesicle: the wrapped fraction for the Adhesion potential)
        shell_width    the shell of the three coverages: the half width of the cosine shell of the config (Adhesion), rc/a - 1 for the Frenkel potential
        nV, finished   vertices of the mesh, 1 if the run wrote Final_state.obj
The points of a chain that goes back appear twice, once for each leg. Tube chains (chain_summary.csv with a column "kind") have no bead to wrap and are skipped.

usage:  python3 Coverage_data_chain.py BATCH_DIR [TAG ...]          (BATCH_DIR holds the chain folders, e.g. ../Results/Phase_one_vesicle; needs numpy only)
        Collect_chain_results.sh calls it, so that one command brings the csv's, the coverage and the final states of a batch.
"""
import csv
import math
import os
import re
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import Analyze_one_vesicle as A1  # noqa: E402  analyze_run: the energies, the shell of the config and the coverages in it
import Analyze_two_vesicle as A2  # noqa: E402

N_DIRECTIONS = 20000   # as PostProcessing
EPS = 1e-12            # a face seen edge-on from the bead has no solid angle

FIELDS = ["DIR", "link", "leg", "scan_var", "scan_value", "KA", "KB", "KI", "BeadRadius", "rc", "CoveredArea", "CoverageUnion", "MultilayerFrac",
          "cap_cfg", "f_ad", "shell_width", "nV", "finished"]


def chain_points(chain):
    """The rows of chain_summary.csv that are points of the chain (None for a tube chain). Same rule as Collect_chain_results.sh."""
    with open(os.path.join(chain, "chain_summary.csv")) as f:
        rows = list(csv.DictReader(f))
    if rows and "kind" in rows[0]:
        return None
    return [r for r in rows if not re.search(r"_w[0-9]+$", r["run"]) and not r["run"].endswith("_mid")]


def run_folder(chain, name):
    """The run folder main_cluster made inside <chain>/<name>/: the highest number."""
    d = os.path.join(chain, name)
    nums = sorted(int(x) for x in os.listdir(d) if x.isdigit() and os.path.isdir(os.path.join(d, x))) if os.path.isdir(d) else []
    return os.path.join(d, str(nums[-1])) + "/" if nums else None


def ray_union(V, F, pos, sigma, width, power, side, n_dir=N_DIRECTIONS):
    """(fraction of directions from pos hit by at least one face of the shell, fraction hit by two or more): CalculateCoverageRays of Post_processing.cpp.
    The faces of the shell are those of the coverage energy: the centroid in the cosine shell (weight > 0) and the face turned to the bead (side = -1: a bead
    inside the closed membrane)."""
    a, b, c = V[F[:, 0]] - pos, V[F[:, 1]] - pos, V[F[:, 2]] - pos
    cen = (a + b + c) / 3.0
    r = np.linalg.norm(cen, axis=1)
    n = np.cross(b - a, c - a)
    n /= np.linalg.norm(n, axis=1)[:, None]
    facing = side * ((cen / r[:, None]) * n).sum(1) < 0
    sel = facing & (np.abs((r - sigma) / (width * sigma)) < 1.0) & (power > 0)
    i = np.arange(n_dir)
    z = 1.0 - 2.0 * (i + 0.5) / n_dir
    rxy = np.sqrt(1.0 - z * z)
    phi = math.pi * (3.0 - math.sqrt(5.0)) * i
    D = np.stack([rxy * np.cos(phi), rxy * np.sin(phi), z], axis=1)
    unit = lambda x: x / np.linalg.norm(x, axis=1)[:, None]  # noqa: E731
    A, B, C = unit(a[sel]), unit(b[sel]), unit(c[sel])
    trip = (A * np.cross(B, C)).sum(1)
    ok = np.abs(trip) >= EPS
    A, B, C, s = A[ok], B[ok], C[ok], np.sign(trip[ok])[:, None]
    nab, nbc, nca = s * np.cross(A, B), s * np.cross(B, C), s * np.cross(C, A)
    M = unit(A + B + C)
    cos_cap = np.minimum(np.minimum((M * A).sum(1), (M * B).sum(1)), (M * C).sum(1)) - 1e-9
    hits = np.zeros(n_dir, dtype=int)
    for k in range(len(A)):
        idx = np.arange(n_dir) if cos_cap[k] <= 0.0 else np.nonzero(D @ M[k] >= cos_cap[k])[0]   # bounding cap around the face, to skip the far directions
        d = D[idx]
        inside = (d @ nab[k] >= 0.0) & (d @ nbc[k] >= 0.0) & (d @ nca[k] >= 0.0)               # on the inner side of the three great circles of its edges
        hits[idx[inside]] += 1
    return float((hits >= 1).sum()) / n_dir, float((hits >= 2).sum()) / n_dir


def coverage_row(run, point, chain_name):
    q = A1.analyze_run(run, point, point["scan_var"])
    cfg = A1.json.load(open(os.path.join(run, "Input_file.json")))
    energies = {e["Name"]: e["constants"] for e in cfg["Energies"]}
    bead = cfg["Beads"][0]
    KA = 0.0
    for name in ("Surface_tension", "Membrane_tension", "Excess_tension"):
        if name in energies:
            KA = energies[name][0]
    V, F = A2.read_obj(A2.last_mesh(run)[0])
    p = A2.read_bead_positions(run, 1)[0]
    side = -1.0 if bead.get("outside", 1) < 0 else 1.0
    union, multi = ray_union(V, F, p, bead["radius"], q["shell_width"], q["shell_power"], side)
    return {"DIR": os.path.relpath(run, os.path.dirname(chain_name)) + "/", "link": q["link"], "leg": q["leg"], "scan_var": q["scan_var"], "scan_value": q["scan_value"],
            "KA": KA, "KB": q["KB"], "KI": bead.get("inter_str", 0.0), "BeadRadius": bead["radius"], "rc": 0.0 if q["rc"] != q["rc"] else q["rc"],
            "CoveredArea": q["cov_cfg"], "CoverageUnion": union, "MultilayerFrac": multi, "cap_cfg": q["cap_cfg"], "f_ad": q["f_ad"],
            "shell_width": q["shell_width"], "nV": q["nV"], "finished": int(q["finished"])}


def main(batch, tags):
    batch = batch.rstrip("/")
    if not tags:
        tags = sorted(d for d in os.listdir(batch) if os.path.exists(os.path.join(batch, d, "chain_summary.csv")))
    if not tags:
        raise SystemExit("no <tag>/chain_summary.csv in %s (BATCH_DIR is the folder that holds the chain folders)" % batch)
    for tag in tags:
        chain = os.path.join(batch, tag)
        if not os.path.exists(os.path.join(chain, "chain_summary.csv")):
            print("%s: no chain_summary.csv, skipped" % tag)
            continue
        points = chain_points(chain)
        if points is None:
            print("%s: a tube chain (no bead to wrap), no coverage" % tag)
            continue
        rows, bad = [], 0
        for pt in points:
            run = run_folder(chain, pt["run"])
            try:
                rows.append(coverage_row(run, pt, chain))
            except Exception as e:  # a run that is missing or unreadable: reported, the chain goes on
                print("%s: link %s (%s): %s: %s" % (tag, pt["link"], pt["run"], type(e).__name__, e), file=sys.stderr)
                bad += 1
        out = os.path.join(batch, "Coverage_data_%s.csv" % tag)
        with open(out, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(FIELDS)
            for r in rows:
                w.writerow([("%.8g" % r[k]) if isinstance(r[k], float) else r[k] for k in FIELDS])
        print("%s: Coverage_data_%s.csv, %d points%s" % (tag, tag, len(rows), (" (%d unreadable)" % bad) if bad else ""))


if __name__ == "__main__":
    if len(sys.argv) < 2:
        raise SystemExit(__doc__)
    main(sys.argv[1], sys.argv[2:])
