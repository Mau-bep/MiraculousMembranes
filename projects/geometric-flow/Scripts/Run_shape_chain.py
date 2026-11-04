#!/usr/bin/env python3
"""One BRANCH of vesicle shapes followed in the reduced volume nu at a fixed reduced spontaneous curvature c0: every nu is a relaxation warm-started from the previous one.

usage:  python3 Run_shape_chain.py CHAIN.json          (CHAIN.json is written by Create_subjob_shape_chain.py)

The shapes of Seifert, Berndl and Lipowsky (Phys. Rev. A 44, 1182 (1991)): the bending energy F_b = (kappa/2) int (C1 + C2 - C0)^2 dA of a vesicle of fixed area and volume, in
the reduced units nu = V / (4 pi/3 R0^3), R0 = sqrt(A / 4 pi), and c0 = C0 R0. Their Figs. 11 and 13 are F_b / (8 pi kappa) against nu for c0 = 2.4 and 3: the energy of each BRANCH of
stationary shapes (prolate-1, prolate-2, pear), not only the lowest one. A minimiser follows one branch for as long as it is a local minimum, so the energy of a branch is one
chain, started from a shape of that branch and walked in nu until it stops being a minimum.

Units of the code: the Bending energy is K_B int (H - H0)^2 dA with H = mean_curv / dual_area = (C1 + C2)/2 (a sphere has H = 1/R: checked, E = 4 pi K_B for H0 = 0), so
    kappa = K_B / 2,   C0 = 2 H0,   c0 = 2 H0 R0   ->   H0 = c0 / (2 R0),   F_b / (8 pi kappa) = E_bend / (4 pi K_B).
R0 follows from the target nu of the link (V fixed, A = (36 pi V^2 / nu^2)^(1/3)), so c0 stays the same along the chain and H0 changes from link to link. The constraints are penalties,
so the nu that is reached is a little off the target (0.0005 to 0.004); the energy of a penalty minimum is exactly the constrained one at the nu and c0 that were reached, so the
chain records nu and c0_eff = 2 H0 sqrt(A / 4 pi) of the final mesh and the curves are drawn against that nu.

0. the SEED (CHAIN.json "seed"): a mesh of the branch, made here as a surface of revolution of volume V0 (so every chain has the length scale the remesher settings are for):
     {"capsule": {"nu": 0.55}}                   a cylinder with hemispherical caps of that reduced volume (the long cylinder that shrinks towards its shape: the target area then
                                                 decreases, the stable direction of the flow). Prolate branches.
     {"spheres": {"q": 1.3, "neck": 0.4}}        two spheres of radii ratio q >= 1 joined by a smooth throat (a rolling-ball fillet of radius "fillet", default 0.5) whose smallest
                                                 radius is neck, both in units of the radius of the small sphere: q = 1 a dumbbell, q > 1 a pear (the starting shape of the pear
                                                 branch, which the symmetric cylinders never reach). The seed must be smooth: a crease is refined by the remesher to edges of a
                                                 thousandth of the neck and the relaxation stalls in sliver polishes
     {"mesh": "path.obj"}                        any closed triangle mesh (scaled to the volume V0)
1. link 0: the seed relaxed at nu0 ("nu0"). Every further link starts from the final mesh of the one before.
2. the LEGS ("legs": name, nu_to, step, stop_when): each starts from link 0 and walks nu to nu_to in steps of step (the sign of nu_to - nu0 sets the direction). A link is checked (the run
   finished, finite energies, the volume and nu close to their targets) and compared with the one before: when the SHAPE changed class (the neck ratio, the lobe ratio or the skewness
   jumped) the branch has ended there, so the step is halved and the link repeated from the last good one, down to min_step; the link that still jumps at min_step is the landing of the
   branch's end, saved with kind "jump" (the state it fell to), and the leg stops. A leg also stops after the first link that satisfies "stop_when": {"skew_below": s} (a pear that
   became symmetric: the continuous bifurcation C^pear) or {"skew_above": s} (a symmetric branch that became asymmetric). A relaxation that fails is repeated with half the step.

Everything goes to batch_dir: seed.obj, <tag>.json (config of a link), <tag>.log and <tag>/1/ (the run folder: Input_file.json, Output_data.txt, Convergence_log.txt, Final_state.obj),
chain_summary.csv (one line per relaxation: kind seed / point / jump / rejected / failed, nu, F = F_b / 8 pi kappa, c0_eff, the shape numbers, steps, status converged or capped) and
chain.log. RESUME: the links are cached by name and by the link they started from, so after a time limit or a failure sbatch the same subjob again and the chain replays its decisions
from the finished links (delete chain_summary.csv, or the folder, to start over). Analyze_shape_chains.py turns the folders into the curves with snapshots of the shapes.
Needs only json and numpy, so it runs inside the sbatch job.
"""
import copy
import csv
import hashlib
import json
import math
import os
import re
import shutil
import subprocess
import sys
import time

import numpy as np

FIELDS = ["tag", "from_tag", "leg", "kind", "nu_target", "nu", "F", "E_total", "c0_eff", "H0", "steps", "status", "seconds", "nV", "length", "rmax", "r_big", "r_small", "rneck",
          "neck_ratio", "lobe_ratio", "skew", "nlobes", "mesh", "note"]
NUMERIC = {"nu_target", "nu", "F", "E_total", "c0_eff", "H0", "steps", "seconds", "nV", "length", "rmax", "r_big", "r_small", "rneck", "neck_ratio", "lobe_ratio", "skew", "nlobes"}


# ----------------------------------------------------------------------------------------------------------------------
# meshes
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


def write_obj(path, V, F):
    with open(path, "w") as f:
        for p in V:
            f.write("v %.12g %.12g %.12g\n" % tuple(p))
        for t in F:
            f.write("f %d %d %d\n" % (t[0] + 1, t[1] + 1, t[2] + 1))


def volume_area(V, F):
    a, b, c = V[F[:, 0]], V[F[:, 1]], V[F[:, 2]]
    return np.einsum("ij,ij->i", a, np.cross(b, c)).sum() / 6.0, 0.5 * np.linalg.norm(np.cross(b - a, c - a), axis=1).sum()


def reduced_volume(vol, area):
    return 6.0 * math.sqrt(math.pi) * vol / area ** 1.5


# ----------------------------------------------------------------------------------------------------------------------
# seeds: surfaces of revolution about the y axis (the axis of the cylinders in input/)
def capsule_profile(rho, length):
    """(z, r) from the north pole: a cylinder of radius rho and length `length` with hemispherical caps."""
    th = np.linspace(0.0, 0.5 * math.pi, 200)
    z = np.concatenate([0.5 * length + rho * np.cos(th), -0.5 * length - rho * np.sin(th)])
    r = np.concatenate([rho * np.sin(th), rho * np.cos(th)])
    return z, r


def capsule_for_nu(nu, vol):
    """Radius and cylinder length of the capsule of volume vol and reduced volume nu (bisection on the radius: nu grows with it, 1 at the sphere)."""
    rmax = (3.0 * vol / (4.0 * math.pi)) ** (1.0 / 3.0)

    def nu_of(rho):
        length = (vol - 4.0 * math.pi * rho ** 3 / 3.0) / (math.pi * rho ** 2)
        return reduced_volume(vol, 2.0 * math.pi * rho * length + 4.0 * math.pi * rho ** 2), length

    lo, hi = 1e-3 * rmax, rmax * (1.0 - 1e-9)
    if not 0.0 < nu < 1.0:
        raise SystemExit("capsule: nu must be in (0, 1), not %g" % nu)
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if nu_of(mid)[0] < nu else (lo, mid)
    rho = 0.5 * (lo + hi)
    return rho, nu_of(rho)[1]


def two_sphere_profile(q, neck, fillet=0.5):
    """(z, r) from the north pole of a sphere of radius q (up) and one of radius 1 (down) joined by a rolling-ball fillet: a toroidal throat of meridian radius `fillet` and smallest radius
    `neck` (< 1, in units of the small radius) tangent to both spheres. A plain union has a crease (a circle of enormous bending energy and curvature): the remesher refines it to edges
    of a thousandth of the neck radius and the run stalls there, so the seed must be smooth (the radii of curvature of the mesh are what the remesher resolves)."""
    if not (q >= 1.0 and 0.0 < neck < 1.0 and fillet > 0.0):
        raise SystemExit("spheres: need q >= 1, 0 < neck < 1 and fillet > 0 (lengths in units of the small radius), got q = %g, neck = %g, fillet = %g" % (q, neck, fillet))
    rf = neck + fillet                                          # radius of the centre of the fillet circle (at z = 0)
    z1, z2 = math.sqrt((q + fillet) ** 2 - rf ** 2), -math.sqrt((1.0 + fillet) ** 2 - rf ** 2)    # centres of the big and the small sphere
    th1 = math.atan2(rf, -z1)                                   # the big sphere from its north pole to the tangent point
    a1, a2 = math.atan2(-rf, z1), math.atan2(-rf, z2)           # the fillet circle between the two tangent points (through its lowest point, the neck)
    b0 = math.atan2(rf, -z2)                                    # the small sphere from its tangent point to its south pole
    th, al, be = np.linspace(0.0, th1, 400), np.linspace(a1, a2, 200)[1:], np.linspace(b0, math.pi, 400)[1:]
    z = np.concatenate([z1 + q * np.cos(th), fillet * np.cos(al), z2 + np.cos(be)])
    r = np.concatenate([q * np.sin(th), rf + fillet * np.sin(al), np.sin(be)])
    r[0] = r[-1] = 0.0
    return z, r


def revolve(z, r, spacing):
    """Closed triangle mesh of the surface of revolution of the profile (z, r) (pole to pole) with edges of about `spacing`, axis y, faces outwards."""
    s = np.concatenate([[0.0], np.cumsum(np.hypot(np.diff(z), np.diff(r)))])
    n_around = max(int(round(2.0 * math.pi * r.max() / spacing)), 8)
    m = max(int(round(s[-1] / spacing)), 8)
    t = np.linspace(0.0, s[-1], m + 1)
    z, r = np.interp(t, s, z), np.interp(t, s, r)
    ang = 2.0 * math.pi * np.arange(n_around) / n_around
    V = [[0.0, z[0], 0.0]]
    for j in range(1, m):
        V += [[r[j] * math.cos(a), z[j], r[j] * math.sin(a)] for a in ang]
    V.append([0.0, z[m], 0.0])
    ring = lambda j: 1 + (j - 1) * n_around + np.arange(n_around)
    F = []
    for k in range(n_around):
        k1 = (k + 1) % n_around
        F.append([0, ring(1)[k], ring(1)[k1]])
        for j in range(1, m - 1):
            a, a1, b, b1 = ring(j)[k], ring(j)[k1], ring(j + 1)[k], ring(j + 1)[k1]
            F += [[a1, a, b], [a1, b, b1]]
        F.append([len(V) - 1, ring(m - 1)[k1], ring(m - 1)[k]])
    V, F = np.array(V), np.array(F, dtype=int)
    if volume_area(V, F)[0] < 0.0:
        F = F[:, ::-1]
    return V, F


def make_seed(seed, vol, path, spacing=0.05):
    """Write the seed mesh (volume vol) to path. Returns a description."""
    if "capsule" in seed:
        rho, length = capsule_for_nu(float(seed["capsule"]["nu"]), vol)
        z, r = capsule_profile(rho, length)
        what = "capsule of nu %g (radius %.4f, cylinder length %.4f)" % (seed["capsule"]["nu"], rho, length)
    elif "spheres" in seed:
        fillet = float(seed["spheres"].get("fillet", 0.5))
        z, r = two_sphere_profile(float(seed["spheres"]["q"]), float(seed["spheres"]["neck"]), fillet)
        what = "two spheres, radii ratio %g, neck %g and fillet %g of the small radius" % (seed["spheres"]["q"], seed["spheres"]["neck"], fillet)
    elif "mesh" in seed:
        V, F = read_obj(seed["mesh"])
        v0 = abs(volume_area(V, F)[0])
        write_obj(path, V * (vol / v0) ** (1.0 / 3.0), F)
        return "mesh " + seed["mesh"]
    else:
        raise SystemExit('"seed" must have one of "capsule", "spheres", "mesh"')
    V, F = revolve(z, r, 1.0)                    # unit scale first, to scale to the volume
    scale = (vol / volume_area(V, F)[0]) ** (1.0 / 3.0)
    V, F = revolve(z * scale, r * scale, spacing)
    V *= (vol / volume_area(V, F)[0]) ** (1.0 / 3.0)
    write_obj(path, V, F)
    return what


# ----------------------------------------------------------------------------------------------------------------------
# the shape of a relaxed mesh
def shape_frame(V, F):
    """Area weighted centroid, axis e (largest principal axis, signed so that the long tail of the area distribution points up: a pear stands on its big lobe), a second axis u, and
    the skewness of the area along the axis."""
    a, b, c = V[F[:, 0]], V[F[:, 1]], V[F[:, 2]]
    cen = (a + b + c) / 3.0
    area = 0.5 * np.linalg.norm(np.cross(b - a, c - a), axis=1)
    A = area.sum()
    c0 = (area[:, None] * cen).sum(0) / A
    d = cen - c0
    C = (area[:, None, None] * d[:, :, None] * d[:, None, :]).sum(0) / A
    _, U = np.linalg.eigh(C)
    e, u = U[:, 2], U[:, 1]
    s = d @ e
    m3 = (area * s ** 3).sum() / A
    if m3 < 0.0:
        e, m3 = -e, -m3
    sigma = math.sqrt((area * s ** 2).sum() / A)
    return c0, e, u, m3 / sigma ** 3


def shape_numbers(V, F, bins=60):
    """The numbers that say which branch a mesh is on: the reduced volume, the length and largest radius, the profile r(s) along the axis, its lobes (local maxima that rise at
    least 3 percent above the valley between them) and neck, and |skewness|. neck_ratio is the radius of the neck over the smaller outer lobe (1 without a neck), lobe_ratio the
    smaller of the two largest lobes over the larger (1 for one lobe), skew the skewness of the area along the axis (0 for an up-down symmetric shape)."""
    vol, area = volume_area(V, F)
    c0, e, u, skew = shape_frame(V, F)
    s = (V - c0) @ e
    r = np.linalg.norm((V - c0) - np.outer(s, e), axis=1)
    edges = np.linspace(s.min(), s.max() + 1e-12, bins + 1)
    idx = np.clip(np.digitize(s, edges) - 1, 0, bins - 1)
    prof = np.full(bins, np.nan)
    for j in range(bins):
        sel = r[idx == j]
        if len(sel):
            prof[j] = sel.max()
    good = ~np.isnan(prof)
    prof = np.interp(np.arange(bins), np.arange(bins)[good], prof[good])
    prof = np.convolve(np.pad(prof, 1, mode="edge"), np.ones(3) / 3.0, mode="valid")
    peaks = [i for i in range(1, bins - 1) if prof[i] >= prof[i - 1] and prof[i] > prof[i + 1] and prof[i] > 0.3 * prof.max()]
    changed = True
    while changed and len(peaks) > 1:
        changed = False
        for p, q in zip(peaks[:-1], peaks[1:]):
            if prof[p:q + 1].min() > 0.97 * min(prof[p], prof[q]):
                peaks.remove(p if prof[p] < prof[q] else q)
                changed = True
                break
    if len(peaks) >= 2:
        values = sorted((prof[p] for p in peaks), reverse=True)
        rneck = prof[peaks[0]:peaks[-1] + 1].min()
        neck_ratio = rneck / min(prof[peaks[0]], prof[peaks[-1]])
        r_big, r_small = values[0], values[1]
    else:
        rneck = r_big = r_small = float(prof.max())
        neck_ratio = 1.0
    return {"nu": reduced_volume(vol, area), "length": float(s.max() - s.min()), "rmax": float(r.max()), "r_big": float(r_big), "r_small": float(r_small), "rneck": float(rneck),
            "neck_ratio": float(neck_ratio), "lobe_ratio": float(r_small / r_big), "skew": float(skew), "nlobes": max(len(peaks), 1), "nV": len(V)}


def meridian_outline(V, F):
    """The cut of the mesh by the plane through its axis: the outline of the shape as segments ((w0, s0), (w1, s1)) (w across, s along the axis, the centroid at the origin, the big
    lobe at the bottom). Both sides of the shape are in it."""
    c0, e, u, _ = shape_frame(V, F)
    n = np.cross(e, u)
    d = (V - c0) @ n
    segs = []
    for tri in F:
        dd = d[tri]
        if not (dd.max() > 0.0 > dd.min()):
            continue
        pts = []
        for i, j in ((0, 1), (1, 2), (2, 0)):
            if (dd[i] > 0.0) != (dd[j] > 0.0):
                t = dd[i] / (dd[i] - dd[j])
                p = V[tri[i]] + t * (V[tri[j]] - V[tri[i]]) - c0
                pts.append((p @ u, p @ e))
        if len(pts) == 2:
            segs.append(pts)
    return np.array(segs)


# ----------------------------------------------------------------------------------------------------------------------
# main_cluster
def run_main_cluster(binary, config_path, log_path, threads, nice):
    cmd = [binary, config_path]
    if nice:
        cmd = ["nice", "-n", str(nice)] + cmd
    env = dict(os.environ, OMP_NUM_THREADS=str(threads))
    with open(log_path, "w") as log:
        return subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, env=env).returncode


def find_run_folder(first_dir):
    nums = sorted(int(d) for d in os.listdir(first_dir) if d.isdigit() and os.path.isdir(os.path.join(first_dir, d))) if os.path.isdir(first_dir) else []
    return os.path.join(first_dir, str(nums[-1])) + "/" if nums else None


def read_last_row(run):
    with open(os.path.join(run, "Output_data.txt")) as f:
        header = f.readline().split()
        rows = [l.split() for l in f if l.strip()]
    return dict(zip(header, [float(x) for x in rows[-1]]))


def run_problem(rc, text, run):
    """Why a finished main_cluster run is not usable, or None. main_cluster returns 0 and writes a Final_state.obj also when the run broke."""
    if rc != 0:
        return "main_cluster returned %d" % rc
    if re.search(r"The line search hit a nan", text):
        return "the line search hit a nan"
    m = re.search(r"still stalled at t = (\d+), ending the run", text)
    if m:
        return "BFGS-Normal stalled and was ended at step %s" % m.group(1)
    if "Ending sim due to small TS" in text:
        return "ended due to a small timestep"
    if not re.search(r"Final state saved as step (\d+), total energy (\S+)", text):
        return "no 'Final state saved' line in the log"
    if run is None or not os.path.exists(os.path.join(run, "Final_state.obj")):
        return "no Final_state.obj"
    return None


# ----------------------------------------------------------------------------------------------------------------------
class Chain:
    def __init__(self, spec, spec_path):
        self.spec = spec
        self.batch = os.path.abspath(spec["batch_dir"])
        os.makedirs(self.batch, exist_ok=True)
        self.binary = os.path.abspath(spec.get("binary", "../build/bin/main_cluster"))
        self.template = json.load(open(spec["config"]))
        self.c0 = float(spec["c0"])
        self.KB, self.KV, self.KA = float(spec.get("KB", 1.0)), float(spec.get("KV", 1e4)), float(spec.get("KA", 1e4))
        self.dA = float(spec.get("dA", 4e-4))
        self.V0 = float(spec.get("V0", 0.5035))
        self.switch = spec.get("switch", [1000, 4000])
        self.timesteps = int(spec.get("timesteps", 60000))
        self.threads, self.nice = int(spec.get("threads", 1)), int(spec.get("nice", 0))
        self.nu_tol = float(spec.get("nu_tol", 0.01))
        self.jump = dict({"neck_ratio": 0.15, "lobe_ratio": 0.08, "skew": 0.08}, **spec.get("jump", {}))
        self.log = open(os.path.join(self.batch, "chain.log"), "a")
        self.summary_path = os.path.join(self.batch, "chain_summary.csv")
        self.rows = []
        self.cache = {}
        shutil.copy(spec_path, os.path.join(self.batch, "chain.json"))
        if os.path.exists(self.summary_path):
            with open(self.summary_path) as f:
                for r in csv.DictReader(f):
                    row = {k: (float(v) if (k in NUMERIC and v != "") else v) for k, v in r.items()}
                    self.rows.append(row)
                    self.cache[row["tag"]] = row

    # ---- bookkeeping
    def say(self, msg):
        print(msg, flush=True)
        self.log.write("%s %s\n" % (time.strftime("%Y-%m-%d %H:%M:%S"), msg))
        self.log.flush()

    def write_summary(self):
        tmp = self.summary_path + ".tmp"
        with open(tmp, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(FIELDS)
            for r in self.rows:
                w.writerow([r.get(k, "") for k in FIELDS])
            f.flush()
            os.fsync(f.fileno())
        os.replace(tmp, self.summary_path)

    def record(self, row):
        self.rows = [r for r in self.rows if r["tag"] != row["tag"]] + [row]
        self.cache[row["tag"]] = row
        self.write_summary()

    # ---- one relaxation
    def link_config(self, init_mesh, first_dir, nu):
        a_bar = (36.0 * math.pi * self.V0 ** 2 / nu ** 2) ** (1.0 / 3.0)
        r0 = math.sqrt(a_bar / (4.0 * math.pi))
        h0 = self.c0 / (2.0 * r0)
        cfg = copy.deepcopy(self.template)
        cfg["init_file"] = os.path.abspath(init_mesh)
        cfg["first_dir"] = first_dir
        cfg["timesteps"] = self.timesteps
        cfg["Energies"] = [{"Name": "Bending", "constants": [self.KB, h0]}, {"Name": "Volume_constraint", "constants": [self.KV, self.V0]},
                           {"Name": "Area_constraint", "constants": [self.KA, -1.0, nu, self.dA]}]
        # BFGS must not start while the area target is still moving
        _, area0 = volume_area(*read_obj(init_mesh))
        ramp = int(abs(a_bar - area0) / self.dA) + 300
        sw = [max(self.switch[0], ramp), max(self.switch[1], max(self.switch[0], ramp) + 3000)]
        cfg["Switch_times"] = [sw[0], -2500, -1, sw[1], -1, -1]
        if not self.spec.get("bfgs_normal", True):
            cfg["Switch_times"][3] = -1                          # stay in BFGS (it remeshes), the "stopping" block then ends the run
        cfg["remesher"].update(self.spec.get("remesher", {}))
        cfg["stopping"].update(self.spec.get("stopping", {}))
        cfg.update(self.spec.get("config_overrides", {}))        # any other top level key of the config (remesh_every, adapt_remesh, ...)
        return cfg, h0

    def relax(self, tag, leg, kind, from_row, nu):
        """Relax from the mesh of from_row at the target nu: the row of the result (a cached one when this very relaxation was done before)."""
        old = self.cache.get(tag)
        if old is not None and old.get("from_tag") == from_row["tag"] and abs(old["nu_target"] - nu) < 1e-9 and (old["kind"] == "failed" or os.path.exists(old["mesh"])):
            self.say("%s: cached (nu %.4f, F %s)" % (tag, nu, old["F"] if old["F"] != "" else "-"))
            return dict(old, kind=kind if old["kind"] != "failed" else "failed", leg=leg)
        first_dir = os.path.join(self.batch, tag) + "/"
        if os.path.isdir(first_dir):
            shutil.rmtree(first_dir)
        cfg, h0 = self.link_config(from_row["mesh"], first_dir, nu)
        cfg_path = os.path.join(self.batch, tag + ".json")
        json.dump(cfg, open(cfg_path, "w"), indent=4)
        log_path = os.path.join(self.batch, tag + ".log")
        t0 = time.time()
        rc = run_main_cluster(self.binary, cfg_path, log_path, self.threads, self.nice)
        seconds = time.time() - t0
        text = open(log_path, errors="replace").read()
        run = find_run_folder(first_dir)
        row = {"tag": tag, "from_tag": from_row["tag"], "leg": leg, "kind": kind, "nu_target": nu, "H0": h0, "seconds": round(seconds, 1)}
        problem = run_problem(rc, text, run)
        if problem is None:
            try:
                last = read_last_row(run)
                V, F = read_obj(run + "Final_state.obj")
                vol, area = volume_area(V, F)
                nums = shape_numbers(V, F)
                steps = int(re.search(r"Final state saved as step (\d+)", text).group(1))
                converged = re.search(r"converged at t = \d+", text) is not None
                row.update(nums)
                row.update({"F": last["Bending"] / (4.0 * math.pi * self.KB), "E_total": last["Total_E"], "c0_eff": 2.0 * h0 * math.sqrt(area / (4.0 * math.pi)), "steps": steps,
                            "status": "converged" if converged else ("capped" if steps >= self.timesteps else "ended"), "mesh": run + "Final_state.obj"})
                if not all(math.isfinite(row[k]) for k in ("F", "E_total", "nu")):
                    problem = "non-finite energy or reduced volume"
                elif abs(vol - self.V0) > 5e-3 * self.V0:
                    problem = "the volume %.5g is not the target %.5g" % (vol, self.V0)
                elif abs(nums["nu"] - nu) > self.nu_tol:
                    problem = "the reduced volume %.4f is not the target %.4f" % (nums["nu"], nu)
                elif abs(last["Volume"] - vol) > 1e-4 * self.V0 or abs(last["Area"] - area) > 1e-4 * area:
                    problem = "Final_state.obj is not the state of the last line of Output_data.txt"
            except (OSError, ValueError, KeyError, IndexError, AttributeError, np.linalg.LinAlgError) as exc:
                problem = "unreadable output: %r" % (exc,)
        if problem is not None:
            row.update({"kind": "failed", "note": problem})
            self.say("%s: FAILED at nu %.4f: %s" % (tag, nu, problem))
        else:
            self.say("%s: nu %.4f -> %.4f  F %.5f  neck %.2f lobes %.2f skew %.3f  %d steps %s  %.0f s" % (tag, nu, row["nu"], row["F"], row["neck_ratio"], row["lobe_ratio"], row["skew"],
                                                                                                               row["steps"], row["status"], seconds))
        self.record(row)
        return row

    # ---- the decisions
    def jumped(self, prev, cur):
        return (abs(cur["neck_ratio"] - prev["neck_ratio"]) > self.jump["neck_ratio"] or abs(cur["lobe_ratio"] - prev["lobe_ratio"]) > self.jump["lobe_ratio"]
                or abs(cur["skew"] - prev["skew"]) > self.jump["skew"])

    def run_leg(self, leg, origin):
        name, nu_to = leg["name"], float(leg["nu_to"])
        nominal = abs(float(leg["step"]))
        min_step = float(leg.get("min_step", nominal / 8.0))
        stop = leg.get("stop_when", {})
        cur, cur_nu, step, k = origin, float(origin["nu_target"]), nominal, 0
        sign = 1.0 if nu_to > cur_nu else -1.0
        while (nu_to - cur_nu) * sign > 1e-9:
            s = min(step, (nu_to - cur_nu) * sign)
            target = round(cur_nu + sign * s, 6)
            k += 1
            tag = "%s_%03d_nu%.4f" % (name, k, target)
            row = self.relax(tag, name, "point", cur, target)
            if row["kind"] == "failed":
                if s / 2.0 >= min_step - 1e-12:
                    step = s / 2.0
                    continue
                self.say("leg %s: stops, the relaxation at nu %.4f keeps failing" % (name, target))
                return
            if self.jumped(cur, row):
                if s / 2.0 >= min_step - 1e-12:
                    row["kind"] = "rejected"
                    self.record(row)
                    step = s / 2.0
                    self.say("leg %s: the shape changed class between nu %.4f and %.4f, the step is halved to %g" % (name, cur_nu, target, step))
                    continue
                row["kind"] = "jump"
                self.record(row)
                self.say("leg %s: the branch ends between nu %.4f and %.4f, it falls to F %.5f (neck %.2f, lobes %.2f, skew %.3f)" % (name, cur_nu, target, row["F"], row["neck_ratio"],
                                                                                                                                     row["lobe_ratio"], row["skew"]))
                return
            cur, cur_nu, step = row, target, min(nominal, 2.0 * s)
            if "skew_below" in stop and row["skew"] < float(stop["skew_below"]):
                self.say("leg %s: stops at nu %.4f, the skewness %.3f is below %g (the shape is symmetric)" % (name, target, row["skew"], stop["skew_below"]))
                return
            if "skew_above" in stop and row["skew"] > float(stop["skew_above"]):
                self.say("leg %s: stops at nu %.4f, the skewness %.3f is above %g (the shape lost its symmetry)" % (name, target, row["skew"], stop["skew_above"]))
                return
        self.say("leg %s: reached nu %.4f" % (name, cur_nu))

    def run(self):
        spec = self.spec
        with open(self.binary, "rb") as f:
            self.say("binary %s md5 %s" % (self.binary, hashlib.md5(f.read()).hexdigest()))
        seed_path = os.path.join(self.batch, "seed.obj")
        if not os.path.exists(seed_path):
            self.say("seed: " + make_seed(spec["seed"], self.V0, seed_path))
        nu0 = float(spec["nu0"])
        origin = self.relax("link0", "seed", "seed", {"tag": "seed", "mesh": seed_path}, nu0)
        if origin["kind"] == "failed":
            self.say("chain stops: link 0 at nu %.4f failed" % nu0)
            return False
        for leg in spec["legs"]:
            self.run_leg(leg, origin)
        self.say("chain done")
        return True


def main(path):
    spec = json.load(open(path))
    return 0 if Chain(spec, path).run() else 1


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise SystemExit(__doc__)
    sys.exit(main(sys.argv[1]))
