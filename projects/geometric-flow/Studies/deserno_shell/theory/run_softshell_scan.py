#!/usr/bin/env python
"""Scan of the continuum soft-shell model: lines S1, E, S2 against shell width s, shell power p and sigma~.

    OMP_NUM_THREADS=1 micromamba run -n mir_mem python run_softshell_scan.py [--part N/M] [--set main|s0|R|sim]

Results are appended (one row per case, resumable) to ../data/softshell_lines.csv .  Sets
    main : s in {0.25,0.15,0.10,0.05,0.025} x p in {1,2,4,8} x sigma~ in {0.13,0.27,0.53,1.0}, R = 10 (the disk of the simulation)
    s0   : s -> 0 series, p = 1, s in {0.2,0.1,0.05,0.025,0.0125,0.00625}, R = 10 (+ R = 30 spot checks), validation against Deserno
    R    : R_dom study at s = 0.25 (R = 6, 10, 20, 40)
    sim  : sigma~ = 2k/15 (k = 1..15), the simulation grid, s = 0.25, p = 1, R = 10 (not used for the deliverable: compare_sim.py does this and writes softshell_vs_sim_lines.csv)
"""
import os
import sys
import csv
import time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "Scripts"))
import softshell_direct as sd
import deserno_theory as dt

OUT = os.path.join(HERE, "..", "data", "softshell_lines.csv")
COLS = ["set", "sigma_t", "s", "p", "R", "S1", "E", "S2", "G_S1", "G_S2", "G_E_partial", "G_E_env", "E_energy", "G_flat",
        "S1_coarse", "S2_coarse", "deserno_S1", "deserno_E", "deserno_S2", "grid_scan", "grid_ref", "note", "time_s"]
SIG = [0.13, 0.27, 0.53, 1.0]


def cases(which):
    out = []
    if which == "main":
        for p in (1, 2, 4, 8):
            for s in (0.25, 0.15, 0.10, 0.05, 0.025):
                for sg in SIG:
                    out.append(("main", sg, s, p, 10.0))
        # p = 1 first, then s = 0.25 with p > 1, then the rest
        out.sort(key=lambda c: (0 if c[3] == 1 else (1 if c[2] == 0.25 else 2), c[3], -c[2], c[1]))
    elif which == "s0":
        for s in (0.2, 0.1, 0.05, 0.025, 0.0125, 0.00625):
            for sg in SIG:
                out.append(("s0", sg, s, 1, 10.0))
        for s in (0.1, 0.025):                       # large-R check where the clamp matters most (smallest sigma~)
            out.append(("s0", 0.13, s, 1, 30.0))
            out.append(("s0", 1.0, s, 1, 30.0))
    elif which == "R":
        for R in (6.0, 10.0, 20.0, 40.0):
            for sg in SIG:
                out.append(("R", sg, 0.25, 1, R))
    elif which == "sim":
        for k in range(1, 16):
            out.append(("sim", round(2.0 * k / 15.0, 6), 0.25, 1, 10.0))
    return out


def grids_for(s, p):
    if s >= 0.1 and p <= 2:
        return "mid", "mid"
    return "mid", "mid"


def done_keys():
    keys = set()
    if os.path.exists(OUT):
        with open(OUT) as f:
            for r in csv.DictReader(f):
                keys.add((r["set"], round(float(r["sigma_t"]), 6), float(r["s"]), float(r["p"]), float(r["R"])))
    return keys


def append_row(row):
    new = not os.path.exists(OUT)
    with open(OUT, "a", newline="") as f:
        w = csv.DictWriter(f, fieldnames=COLS)
        if new:
            w.writeheader()
        w.writerow({k: row.get(k, "") for k in COLS})


def run_case(c, verbose=False):
    which, sg, s, p, R = c
    gs, gr = grids_for(s, p)
    row = dict(set=which, sigma_t=sg, s=s, p=p, R=R, grid_scan=gs, grid_ref=gr)
    try:
        dS1, dS2, dE = dt.spinodal_S1(sg), dt.spinodal_S2(sg), dt.w_E(sg)
    except Exception:
        dS1 = dS2 = dE = float("nan")
    row.update(deserno_S1=dS1, deserno_S2=dS2, deserno_E=dE)
    t0 = time.time()
    try:
        wt0 = 2.0
        out = sd.soft_lines(sg, s=s, p=p, R=R, grid_scan=gs, grid_ref=gr, wt0=wt0, verbose=verbose)
        row.update({k: out.get(k, "") for k in ("S1", "S2", "E", "G_S1", "G_S2", "S1_coarse", "S2_coarse")})
        row.update(G_E_partial=out.get("E_G_partial", ""), G_E_env=out.get("E_G_env", ""), E_energy=out.get("E_energy", ""),
                   G_flat=out.get("G0", ""), note=(out.get("note") or out.get("E_note") or ""))
    except Exception as exc:                                  # keep the scan going, report the failure
        row["note"] = "FAILED: %s: %s" % (type(exc).__name__, str(exc)[:120])
    row["time_s"] = round(time.time() - t0, 1)
    return row


if __name__ == "__main__":
    which = "main"
    part = (0, 1)
    for i, a in enumerate(sys.argv):
        if a == "--set":
            which = sys.argv[i + 1]
        if a == "--part":
            k, n = sys.argv[i + 1].split("/")
            part = (int(k), int(n))
    todo = cases(which)
    done = done_keys()
    mine = [c for j, c in enumerate(todo) if j % part[1] == part[0]]
    for c in mine:
        key = (c[0], round(c[1], 6), float(c[2]), float(c[3]), float(c[4]))
        if key in done:
            continue
        row = run_case(c)
        append_row(row)
        print("%s sigma~=%.3f s=%.4g p=%g R=%g  S1=%s E=%s S2=%s  %s (%.0fs)" % (
            c[0], c[1], c[2], c[3], c[4], row.get("S1"), row.get("E"), row.get("S2"), row.get("note", ""), row["time_s"]), flush=True)
