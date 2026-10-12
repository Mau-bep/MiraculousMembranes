#!/usr/bin/env python
"""Continuum soft-shell prediction at zero tension, for the planar potentials scan (Scripts/Mem_planar_potentials_up.sh).

    OMP_NUM_THREADS=2 micromamba run -n mir_mem python softshell_sigma0.py [--sigma 0] [--s 0.25] [--p 1] [--R 10]

Traces the constrained minimum e(G) of softshell_direct.py (bead + cosine shell of half width s and power p + bending + excess
tension sigma~, disk of radius R clamped at the edge) in the coverage G, from the flat plate to G = 1.985, and writes
../data/softshell_sigma0.csv with columns
    G  wt  e  Ebend  Aex          (units a = 1, kappa = KB/2, energies in pi kappa; wt = de/dG = w~ at which the state is stationary)
so that at w~ the state has E~ = e - w~ G and bending Ebend, the simulation's Total_E = (pi KB/2) E~ (KB = 1), w~ = 4 KI.
Only the cosine shell has a continuum model (cosine, s = 0.25 is the 'cos' condition of the scan); Frenkel, LJ and Morse profiles do not.
Scripts/Plots.py: Softshell_prediction() reads the file.
"""
import os
import sys
import csv
import argparse
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import softshell_direct as sd

ap = argparse.ArgumentParser()
ap.add_argument("--sigma", type=float, default=0.0, help="sigma~ = 2 KA a^2/KB (0: the zero tension of the scan)")
ap.add_argument("--s", type=float, default=0.25)
ap.add_argument("--p", type=float, default=1.0)
ap.add_argument("--R", type=float, default=10.0)
ap.add_argument("--grid", default="mid")
ap.add_argument("--out", default=os.path.join(HERE, "..", "data", "softshell_sigma0.csv"))
args = ap.parse_args()

model = sd.Model(args.sigma, s=args.s, p=args.p, R=args.R)
kw = sd.GRIDS[args.grid]
G_end, wt0 = 1.985, 2.0
try:
    st = sd.initial_state(model, wt0, args.grid)
    G0 = sd.energies(model, st, 0.0)[3]
    br = sd.trace_branch(model, st, wt0, G0, G_end, dG0=0.02, dG_max=0.03, grid_kw=kw, verbose=True)
    ok = br.G[-1] > G_end - 0.05
except sd.SolveFail as exc:
    print("direct start failed:", exc)
    ok = False
if not ok:
    # narrow shell: the plate stays stable up to w~ ~ 3-4, start on the bound branch
    st0, wt0 = sd.initial_state_bound(model, grid=args.grid)
    G0 = sd.energies(model, st0, 0.0)[3]
    br = sd.trace_branch(model, st0, wt0, G0, G_end, dG0=0.02, dG_max=0.03, grid_kw=kw, verbose=True)
a = br.arrays()
print("stopped:", getattr(br, "stopped", None), " G range %.3f .. %.3f, w~ range %.3f .. %.3f" % (a["G"].min(), a["G"].max(), a["wt"].min(), a["wt"].max()))
os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
with open(args.out, "w", newline="") as f:
    f.write("# softshell_sigma0.py: sigma~=%g s=%g p=%g R=%g; G wt e Ebend Aex (a=1, kappa=KB/2, energies in pi kappa)\n" % (args.sigma, args.s, args.p, args.R))
    w = csv.writer(f, delimiter=" ")
    for row in zip(a["G"], a["wt"], a["e"], a["Eb"], a["Aex"]):
        w.writerow(["%.10g" % x for x in row])
print("wrote", args.out)
