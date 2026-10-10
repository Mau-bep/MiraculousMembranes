#!/usr/bin/env python
"""Simulation campaigns (run_condition.py tags) against the continuum soft-shell model (theory/softshell_direct.py).

    OMP_NUM_THREADS=2 micromamba run -n mir_mem python compare_campaign.py --tags p1,p2,p4,p8 [--out ../data/campaign_pXX.csv]

For every finished run (state done in Results/_sim_ctl/<tag>/summary.csv, so run `run_condition.py collect` first):
  * w~ = 4 KI r^2/KB, sigma~ = 2 KA r^2/KB, coverage G_sim = -E_bead/(2 pi KI r^2) (the 'z_ad' of the continuum model),
    E~_sim = Total_E / (pi kappa) = 2 Total_E/(pi KB)  (energy in units pi kappa);
  * the continuum model's equilibrium curve (e(G), w~(G)) for the same (sigma~, s, p) gives, at this w~, every state
    (stable = local minimum at fixed w~): its coverage and energy E~ = e(G) - w~ G. The simulated state is matched with the
    stable state of nearest coverage; dG, dE are sim minus continuum.
  * per tag the simulated lines are bracketed from the runs (flat start jump = S1, full start unwrapping = S2, flat/full
    energy crossing = E) and compared with the continuum and Deserno's lines.
Theory curves are cached in ../data/campaign_theory_cache.json (keyed by sigma~, s, p).
"""
import argparse, csv, json, math, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
STUDY = os.path.dirname(HERE)
P = os.path.dirname(os.path.dirname(STUDY))
sys.path.insert(0, os.path.join(STUDY, "theory"))
sys.path.insert(0, os.path.join(P, "Scripts"))
CTL = os.path.join(P, "Results", "_sim_ctl")
CACHE = os.path.join(STUDY, "data", "campaign_theory_cache.json")
GJUMP = 1.7                                    # coverage above which a state counts as enveloped


def theory(sigma, s, p):
    key = "%.4f_%.4f_%.3f" % (sigma, s, p)
    cache = json.load(open(CACHE)) if os.path.exists(CACHE) else {}
    if key not in cache:
        import softshell_direct as sd
        out = sd.soft_lines(sigma, s=s, p=p)
        a = out["branch"].arrays()
        cache[key] = dict(G=a["G"].tolist(), e=a["e"].tolist(), wt=a["wt"].tolist(),
                          S1=out.get("S1"), S2=out.get("S2"), E=out.get("E"))
        json.dump(cache, open(CACHE, "w"))
    c = cache[key]
    return {k: (np.array(v) if isinstance(v, list) else v) for k, v in c.items()}


def states_at(th, w):
    G, e, wt = th["G"], th["e"], th["wt"]
    res = []
    for i in range(len(G) - 1):
        a, b = wt[i] - w, wt[i + 1] - w
        if a == 0 or a * b < 0:
            t = a / (a - b)
            Gx = G[i] + t * (G[i + 1] - G[i])
            ex = e[i] + t * (e[i + 1] - e[i])
            res.append((Gx, ex - w * Gx, wt[i + 1] > wt[i]))
    return res


def deserno_lines(sigma):
    import deserno_theory as dt
    return dt.spinodal_S1(sigma), dt.w_E(sigma), dt.spinodal_S2(sigma)


def bracket(pairs, pred):
    """pairs: sorted (w, value). last w where pred false and first later w where pred true (None if no change)."""
    pairs = sorted(pairs)
    for (w0, v0), (w1, v1) in zip(pairs[:-1], pairs[1:]):
        if (not pred(v0)) and pred(v1):
            return (w0, w1)
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tags", required=True)
    ap.add_argument("--out", default=None)
    args = ap.parse_args()
    rows = []
    for tag in args.tags.split(","):
        man = json.load(open(os.path.join(CTL, tag, "manifest.json")))
        s = man["shell_width"] or 0.25
        p = man["shell_power"] or 1.0
        KA, KB = man["KA"], man["KB"]
        sigma = 2 * KA / KB
        th = theory(sigma, s, p)
        f = os.path.join(CTL, tag, "summary.csv")
        if not os.path.exists(f):
            continue
        runs = [r for r in csv.DictReader(open(f)) if r["state"] == "done"]
        runs.sort(key=lambda r: (float(r["KI"]), r["init"]))
        print("==== %s  s=%.3f p=%g sigma~=%.4f   continuum S1/E/S2 = %.3f %.3f %.3f   Deserno = %.3f %.3f %.3f" % (
            tag, s, p, sigma, th["S1"], th["E"], th["S2"], *deserno_lines(sigma)))
        print(" init    w~   G_sim  E~_sim | continuum stable states (G, E~)            | dG(near)  dE(near) | conv faces steps")
        flat, full = [], []
        for r in runs:
            KI = float(r["KI"]); w = 4 * KI / KB
            G = -float(r["E_bead"]) / (2 * math.pi * KI)
            E = 2 * float(r["Total_E"]) / (math.pi * KB)
            st = [x for x in states_at(th, w) if x[2]]
            near = min(st, key=lambda x: abs(x[0] - G)) if st else None
            lo = min(st, key=lambda x: x[1]) if st else None
            sts = " ".join("(%.3f,%.3f)" % (x[0], x[1]) for x in st)
            dG = G - near[0] if near else float("nan")
            dE = E - near[1] if near else float("nan")
            print(" %-5s %5.2f %6.3f %8.3f | %-45s | %+7.3f %+8.3f | %s %s %s" % (
                r["init"], w, G, E, sts, dG, dE, r["converged"], r["faces"], r["steps"]))
            (flat if r["init"] == "flat" else full).append((w, G, E))
            rows.append(dict(tag=tag, s=s, p=p, sigma_t=sigma, init=r["init"], KI=KI, w_t=w, G_sim=G, E_sim=E,
                             G_cont=near[0] if near else "", E_cont=near[1] if near else "",
                             E_cont_lowest=lo[1] if lo else "", dG=dG, dE=dE, steps=r["steps"], faces=r["faces"],
                             converged=r["converged"]))
        b1 = bracket([(w, G) for w, G, E in flat], lambda g: g > GJUMP)
        # S2: smallest w at which the full start is still enveloped: look from high w downwards
        fs = sorted((w, G) for w, G, E in full)
        b2 = None
        for (w0, g0), (w1, g1) in zip(fs[:-1], fs[1:]):
            if g0 <= GJUMP < g1:
                b2 = (w0, w1)
        # E: sign change of E_flat - E_full at equal w
        common = sorted(set(w for w, _, _ in flat) & set(w for w, _, _ in full))
        diffs = [(w, [E for ww, _, E in flat if ww == w][0] - [E for ww, _, E in full if ww == w][0]) for w in common]
        bE = None
        for (w0, d0), (w1, d1) in zip(diffs[:-1], diffs[1:]):
            if d0 * d1 < 0:
                t = d0 / (d0 - d1)
                bE = (w0, w1, w0 + t * (w1 - w0))
        print(" -> S1 sim in %s (jump of the flat start above G=%.1f), continuum %.3f" % (b1, GJUMP, th["S1"]))
        print(" -> S2 sim in %s (full start enveloped above), continuum %.3f" % (b2, th["S2"]))
        print(" -> E  sim in %s (flat-full energy crossing, linear), continuum %.3f" % (bE, th["E"]))
        print(" -> E_flat - E_full:", " ".join("w=%.1f:%+.3f" % d for d in diffs))
    if args.out and rows:
        with open(args.out, "w", newline="") as fh:
            wr = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
            wr.writeheader(); wr.writerows(rows)


if __name__ == "__main__":
    main()
