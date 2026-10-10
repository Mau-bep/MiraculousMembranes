#!/usr/bin/env python
"""s -> 0 convergence of the soft-shell lines to Deserno's (p = 1, R = 10): shifts, local exponents d ln(shift)/d ln s and
power-law / (a sqrt(s) + b s) fits.  Writes ../data/softshell_s0_exponents.csv and prints a markdown table.

    micromamba run -n mir_mem python analyze_s0.py
"""
import os
import csv
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")


def fl(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return float("nan")


rows = list(csv.DictReader(open(os.path.join(DATA, "softshell_lines.csv"))))
out = []
print("| sigma~ | line | s | line(soft) | Deserno | shift | local exponent |")
print("|---|---|---|---|---|---|---|")
for sg in (0.13, 0.27, 0.53, 1.0):
    for key, dkey in (("S1", "deserno_S1"), ("E", "deserno_E"), ("S2", "deserno_S2")):
        d = {}
        for r in rows:
            if abs(fl(r["sigma_t"]) - sg) < 1e-6 and fl(r["p"]) == 1 and fl(r["R"]) == 10.0 and r["set"] in ("s0", "main"):
                v = fl(r[key])
                if np.isfinite(v):
                    d[fl(r["s"])] = (v, fl(r[dkey]))
        ss = sorted(d, reverse=True)
        shifts = {s: d[s][0] - d[s][1] for s in ss}
        for i, s in enumerate(ss):
            ex = float("nan")
            if i > 0 and shifts[s] > 0 and shifts[ss[i - 1]] > 0:
                ex = np.log(shifts[ss[i - 1]] / shifts[s]) / np.log(ss[i - 1] / s)
            print("| %.2f | %s | %g | %.4f | %.4f | %+.4f | %s |" % (sg, key, s, d[s][0], d[s][1], shifts[s], "%.2f" % ex if np.isfinite(ex) else ""))
            out.append(dict(sigma_t=sg, line=key, s=s, soft=d[s][0], deserno=d[s][1], shift=shifts[s], local_exponent=ex))
        # global fits on s <= 0.05
        sm = [s for s in ss if s <= 0.05 and shifts[s] > 0]
        if len(sm) >= 3:
            x = np.log(sm)
            y = np.log([shifts[s] for s in sm])
            slope, icpt = np.polyfit(x, y, 1)
            A = np.vstack([np.sqrt(sm), sm]).T
            coef, res, *_ = np.linalg.lstsq(A, [shifts[s] for s in sm], rcond=None)
            print("|  | %s fit s<=0.05 (%d pts) | | power law exponent %.2f, prefactor %.2f | | a sqrt(s) + b s: a=%.2f b=%.2f | |" % (key, len(sm), slope, np.exp(icpt), coef[0], coef[1]))
            out.append(dict(sigma_t=sg, line=key, s=-1, soft=float("nan"), deserno=float("nan"), shift=float("nan"), local_exponent=slope))
with open(os.path.join(DATA, "softshell_s0_exponents.csv"), "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(out[0].keys()))
    w.writeheader()
    w.writerows(out)
