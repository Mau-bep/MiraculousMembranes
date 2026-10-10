#!/usr/bin/env python
"""Continuum soft-shell (s = 0.25, p = 1, R = 10) against the OLD simulation (data/baseline_lines.csv, baseline_partial_branch.csv).

    OMP_NUM_THREADS=1 micromamba run -n mir_mem python compare_sim.py

writes ../data/softshell_vs_sim_lines.csv, ../data/softshell_vs_sim_partial.csv, ../data/softshell_branches.npz
and the figures ../figures/softshell_vs_sim_lines.png, softshell_vs_sim_partial.png, softshell_branch_w_of_G.png
"""
import os
import sys
import csv
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "Scripts"))
import softshell_direct as sd
import deserno_theory as dt
from scipy.interpolate import CubicHermiteSpline

DATA = os.path.join(HERE, "..", "data")
FIG = os.path.join(HERE, "..", "figures")


def read_csv(path):
    with open(path) as f:
        return list(csv.DictReader(f))


def fnum(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return float("nan")


SIMDIR = "/home/mrojasve/Documents/DDG/MiraculousMembranes/projects/geometric-flow/Results/Wrapping_planar_excess/"


def read_sim(name):
    """Coverage_data*.txt of the old campaign -> list of dicts (KA, KI, z_ad = -Bead/(2 pi KI), zu = 2 CoverageUnion, Total_E)"""
    out = []
    with open(SIMDIR + name) as f:
        hdr = f.readline().split()[1:]
        for line in f:
            t = line.split()
            if len(t) != len(hdr):
                continue
            d = dict(zip(hdr, t))
            ki, ka = float(d["KI"]), float(d["KA"])
            out.append(dict(KA=ka, KI=ki, sigma_t=2.0 * ka, w=4.0 * ki, zad=-float(d["Bead"]) / (2 * np.pi * ki),
                            zu=2.0 * float(d["CoverageUnion"]), E=2.0 * float(d["Total_E"]) / np.pi))
    return out


def sim_at(rows, sg, w):
    m = [r for r in rows if abs(r["sigma_t"] - sg) < 2e-3 and abs(r["w"] - w) < 0.05]
    return m[0] if m else None


def lower_state(br_arr, i1, w):
    """coverage and energy of the stationary state with w~(G) = w on the lower branch (G < G_S1), Hermite interpolation of e(G)"""
    G, e, wt = br_arr["G"][: i1 + 1], br_arr["e"][: i1 + 1], br_arr["wt"][: i1 + 1]
    if w > wt.max() or w < wt[0]:
        return np.nan, np.nan
    k = int(np.nonzero(wt >= w)[0][0])                    # wt rises monotonically on the lower branch
    if k == 0:
        return G[0], e[0] - w * G[0]
    Gq = G[k - 1] + (G[k] - G[k - 1]) * (w - wt[k - 1]) / (wt[k] - wt[k - 1])
    eq = float(CubicHermiteSpline(G[max(0, k - 2): k + 2], e[max(0, k - 2): k + 2], wt[max(0, k - 2): k + 2])(Gq))
    return Gq, eq - w * Gq


def upper_state(br_arr, i2, w):
    G, e, wt = br_arr["G"][i2:], br_arr["e"][i2:], br_arr["wt"][i2:]
    if w < wt.min() or w > wt.max():
        return np.nan, np.nan
    k = int(np.nonzero(wt >= w)[0][0])
    if k == 0:
        return G[0], e[0] - w * G[0]
    Gq = G[k - 1] + (G[k] - G[k - 1]) * (w - wt[k - 1]) / (wt[k] - wt[k - 1])
    eq = float(CubicHermiteSpline(G[max(0, k - 2): k + 2], e[max(0, k - 2): k + 2], wt[max(0, k - 2): k + 2])(Gq))
    return Gq, eq - w * Gq


if __name__ == "__main__":
    base = read_csv(os.path.join(DATA, "baseline_lines.csv"))
    part = read_csv(os.path.join(DATA, "baseline_partial_branch.csv"))
    sigmas = sorted({round(fnum(r["sigma_t"]), 4) for r in base if fnum(r["sigma_t"]) > 0})
    FLAT, FULL = read_sim("Coverage_data.txt"), read_sim("Coverage_data_full.txt")
    lines, partial, store = [], [], {}
    for sg in sigmas:
        out = sd.soft_lines(sg, s=0.25, p=1.0, R=10.0, grid_scan="mid", grid_ref="fine")
        br = out["branch"]
        a = br.arrays()
        store["G_%.4f" % sg], store["e_%.4f" % sg], store["wt_%.4f" % sg] = a["G"], a["e"], a["wt"]
        brow = [r for r in base if abs(fnum(r["sigma_t"]) - sg) < 1e-3][0]
        row = dict(sigma_t=sg, S1_cont=out["S1"], E_cont=out["E"], S2_cont=out["S2"],
                   G_S1=out.get("G_S1"), G_S2=out.get("G_S2"), G_E_partial=out.get("E_G_partial"), G_E_env=out.get("E_G_env"),
                   S1_sim=fnum(brow["S1_sim"]), E_sim=fnum(brow["E_sim"]), S2_sim=fnum(brow["S2_sim"]),
                   S1_sim_lo=fnum(brow["S1_lo"]), S1_sim_hi=fnum(brow["S1_hi"]), S2_sim_lo=fnum(brow["S2_lo"]), S2_sim_hi=fnum(brow["S2_hi"]),
                   E_sim_lo=fnum(brow["E_lo"]), E_sim_hi=fnum(brow["E_hi"]),
                   S1_deserno=fnum(brow["theory_S1"]), E_deserno=fnum(brow["theory_E"]), S2_deserno=fnum(brow["theory_S2"]), note=out.get("note", ""))
        lines.append(row)
        print("sigma~=%.3f cont S1 %.3f E %.3f S2 %.3f | sim %.2f %.2f %.2f | Deserno %.3f %.3f %.3f" % (
            sg, row["S1_cont"], row["E_cont"], row["S2_cont"], row["S1_sim"], row["E_sim"], row["S2_sim"],
            row["S1_deserno"], row["E_deserno"], row["S2_deserno"]), flush=True)
        # coverage before / after the jumps: continuum vs simulation (bracket ends of the simulated jumps)
        imax, imin = out.get("imax"), out.get("imin")
        if imax is not None and imin is not None:
            Gp_S1, _ = lower_state(a, imax, out["S1"] - 1e-3)
            Ge_S1, _ = upper_state(a, imin, out["S1"] + 1e-3)           # enveloped state just above S1 (flat start jumps to it)
            Gp_S2, _ = lower_state(a, imax, out["S2"] + 1e-3)           # partially wrapped state the unwrapping enveloped state falls into
            row.update(G_before_S1=out.get("G_S1"), G_after_S1=Ge_S1, G_at_S2_env=out.get("G_S2"), G_after_S2=Gp_S2)
            lo, hi = fnum(brow["S1_lo"]), fnum(brow["S1_hi"])
            fl_, fu_ = sim_at(FLAT, sg, lo), sim_at(FLAT, sg, hi)
            if fl_ and fu_:
                row.update(sim_zad_before_S1=fl_["zad"], sim_zad_after_S1=fu_["zad"], sim_zu_before_S1=fl_["zu"], sim_zu_after_S1=fu_["zu"])
            lo, hi = fnum(brow["S2_lo"]), fnum(brow["S2_hi"])
            gl_, gu_ = sim_at(FULL, sg, lo), sim_at(FULL, sg, hi)
            if gl_ and gu_:
                row.update(sim_zad_after_S2=gl_["zad"], sim_zad_before_S2=gu_["zad"], sim_zu_after_S2=gl_["zu"], sim_zu_before_S2=gu_["zu"])
        # partially wrapped branch: sim flat-start energies vs continuum lower-branch energies
        if imax is None:
            continue
        for r in part:
            if r["file"] != "flat" or abs(fnum(r["sigma"]) - sg) > 1e-3:
                continue
            w, Es, zs = fnum(r["w"]), fnum(r["E"]), fnum(r["zad"])
            Gc, Ec = lower_state(a, imax, w)
            partial.append(dict(sigma_t=sg, w=w, E_sim=Es, zad_sim=zs, G_cont=Gc, E_cont=Ec, dE=Es - Ec, dG=zs - Gc))
    keys = []
    for r in lines:
        keys += [k for k in r if k not in keys]
    with open(os.path.join(DATA, "softshell_vs_sim_lines.csv"), "w", newline="") as f:
        wr = csv.DictWriter(f, fieldnames=keys, restval='')
        wr.writeheader()
        wr.writerows(lines)
    keys = list(partial[0].keys())
    with open(os.path.join(DATA, "softshell_vs_sim_partial.csv"), "w", newline="") as f:
        wr = csv.DictWriter(f, fieldnames=keys)
        wr.writeheader()
        wr.writerows(partial)
    np.savez(os.path.join(DATA, "softshell_branches.npz"), **store)
    ok = [r for r in partial if np.isfinite(r["E_cont"])]
    print("partial-branch rows with a continuum partial state: %d of %d" % (len(ok), len(partial)))
    if ok:
        dE = np.array([r["dE"] for r in ok])
        dG = np.array([r["dG"] for r in ok])
        print("E_sim - E_cont: mean %.3f  rms %.3f  (pi kappa);  zad_sim - G_cont: mean %.3f rms %.3f" % (dE.mean(), np.sqrt((dE ** 2).mean()), dG.mean(), np.sqrt((dG ** 2).mean())))
