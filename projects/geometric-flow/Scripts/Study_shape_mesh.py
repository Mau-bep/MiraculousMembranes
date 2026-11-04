#!/usr/bin/env python3
"""How the remesher settings change the energy and the cost of a relaxed vesicle shape (the thin necks of the pear and prolate-2 branches in particular).

usage:  python3 Study_shape_mesh.py run  [--workers 6] [--cases sphere pear ...] [--groups angle size_max size_min aspect] [--cap 100000]
        micromamba run -n mir_mem python3 Study_shape_mesh.py plot --out ../Studies/shape_chains/mesh

Every run is ONE link of Run_shape_chain.py: a seed relaxed at a fixed nu and c0 with other remesher settings (the same seed, the same stopping rules), so the runs of a case differ only
in the mesh. The remesher (deps/geometric-central remeshing.cpp, RemeshMonitor.cpp) sizes the mesh from the curvature: the shape operator S of every face is estimated from the dihedral
angles of its edges, the sizing field is S^2 / refine_angle^2 clamped to [1/size_max^2, 1/size_min^2], and an edge is split when its length times the square root of the sizing exceeds 1,
so the edge length is about refine_angle / (the largest principal curvature) between size_min and size_max: refine_angle is the angle the normal turns across one edge. At a neck of radius
r_n the largest curvature is 1 / r_n, the edges there are refine_angle * r_n long (refine_angle = 0.2: 31 vertices around the neck), and size_min only matters once it is above that.
Cases (c0 = 2.4 unless said): sphere (nu = 1, the exact energy (1 - c0/2)^2 is known, stiff constraints), prolate1 (nu = 0.80, smooth), prolate2 (nu = 0.69, neck 0.09), pear (nu = 0.72,
neck 0.05) and pear3 (c0 = 3, nu = 0.78: a small bud on a big sphere, the case whose default run loops in sliver polishes), the last three the thin necks. Groups: angle (refine_angle 0.5 ... 0.07, the resolution), size_max, size_min (up to the 0.05 of the wrapping templates) and aspect (aspect_min).
Results: ../Results/Shape_mesh_study/<case>/<group>_<value>/ (a chain folder each) and mesh_study.csv; `plot` writes the figures and a table.
"""
import argparse
import csv
import json
import math
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import Run_shape_chain as R  # noqa: E402

ROOT = "../Results/Shape_mesh_study"
TEMPLATE = "../Config_files/Shape_c0_2.4_prolate1_ConfigFile.json"     # any config written by Create_subjob_shape_chain.py: the runner overrides what it needs


def ensure_template():
    """TEMPLATE when it exists, else the config rendered from ../Templates/Shape_chain.txt with its defaults (needs jinja2, like the Create_subjob scripts): the study then runs without a chain."""
    global TEMPLATE
    if os.path.exists(TEMPLATE):
        return
    from jinja2 import Environment, FileSystemLoader
    os.makedirs(ROOT, exist_ok=True)
    TEMPLATE = os.path.join(ROOT, "study_template_ConfigFile.json")
    data = json.loads(Environment(loader=FileSystemLoader("../Templates/")).get_template("Shape_chain.txt").render(first_dir="../Results/Shape_chains/"))
    json.dump(data, open(TEMPLATE, "w"), indent=4)
    print("template config", TEMPLATE)
BASE = {"size_max": 0.15, "size_min": 0.005, "refine_angle": 0.2, "aspect_min": 0.2}

CASES = {
    "sphere": dict(c0=2.4, seed={"capsule": {"nu": 0.9}}, nu=1.0, extra={"KA": 1e6, "KV": 1e6}, exact=(1.0 - 2.4 / 2.0) ** 2),
    "prolate1": dict(c0=2.4, seed={"capsule": {"nu": 0.70}}, nu=0.80, extra={}),
    "prolate2": dict(c0=2.4, seed={"capsule": {"nu": 0.55}}, nu=0.69, extra={}),
    "pear": dict(c0=2.4, seed={"spheres": {"q": 1.3, "neck": 0.2}}, nu=0.72, extra={}),
    "pear3": dict(c0=3.0, seed={"spheres": {"q": 2.0, "neck": 0.3}}, nu=0.78, extra={}),
}
GROUPS = {
    "angle": [dict(refine_angle=a) for a in (0.5, 0.4, 0.3, 0.2, 0.15, 0.1, 0.07)],
    "size_max": [dict(size_max=s) for s in (0.15, 0.07)],
    "size_min": [dict(size_min=s) for s in (0.002, 0.01, 0.05)],
    "aspect": [dict(aspect_min=a) for a in (0.1, 0.4)],
    # what runs after the first BFGS phase: BFGS-Normal (the default: fixed mesh, polished by remesh + gradient descent when slivers appear) or BFGS with the remesher on
    "integrator": [dict(name="bfgs_only", bfgs_normal=False, stopping={"bfgs_switch": False, "bfgs_tol_E": 1e-5}),
                   dict(name="bfgs_every1", bfgs_normal=False, stopping={"bfgs_switch": False, "bfgs_tol_E": 1e-5}, config_overrides={"remesh_every": 1})],
}
GROUP_CASES = {"angle": ["sphere", "prolate1", "prolate2", "pear", "pear3"], "size_max": ["prolate1", "pear3"], "size_min": ["prolate2", "pear", "pear3"], "aspect": ["pear3"],
               "integrator": ["prolate2", "pear", "pear3"]}
# the finest angle is only run for the cheap cases (the neck meshes have 4 times the vertices at half the angle)
EXPENSIVE_MIN_ANGLE = {"prolate2": 0.1, "pear": 0.1, "pear3": 0.1}
SPECIAL = ("name", "bfgs_normal", "stopping", "config_overrides")


def run_name(group, over):
    if "name" in over:
        return "%s_%s" % (group, over["name"])
    k, v = next(iter(over.items()))
    return "%s_%g" % (group, v)


BINARY = "../build/bin/main_cluster"


def make_spec(case, group, over, cap):
    c = CASES[case]
    folder = os.path.abspath(os.path.join(ROOT, case, run_name(group, over)))
    spec = {"config": os.path.abspath(TEMPLATE), "binary": os.path.abspath(BINARY), "batch_dir": folder, "c0": c["c0"], "KB": 1.0, "V0": 0.5035, "dA": 4e-4,
            "seed": c["seed"], "nu0": c["nu"], "legs": [], "switch": [1000, 4000], "timesteps": cap, "threads": 1, "nice": 0,
            "remesher": dict(BASE, **{k: v for k, v in over.items() if k not in SPECIAL}), **c["extra"], **{k: v for k, v in over.items() if k in SPECIAL and k != "name"}}
    return folder, spec


def one_run(args):
    folder, spec = args
    os.makedirs(folder, exist_ok=True)
    path = os.path.join(folder, "study_spec.json")
    json.dump(spec, open(path, "w"), indent=1)
    t0 = time.time()
    ch = R.Chain(spec, path)
    ok = ch.run()
    return folder, ok, time.time() - t0


def run(args):
    jobs = []
    for group in args.groups:
        for case in GROUP_CASES[group]:
            if args.cases and case not in args.cases:
                continue
            for over in GROUPS[group]:
                if over.get("refine_angle", 1.0) < EXPENSIVE_MIN_ANGLE.get(case, 0.0):
                    continue
                folder, spec = make_spec(case, group, over, args.cap)
                done = os.path.exists(os.path.join(folder, "chain_summary.csv"))
                if done and not args.redo:
                    continue
                if group != "angle" and "name" not in over and over == {k: BASE[k] for k in over} and os.path.exists(os.path.join(ROOT, case, "angle_%g" % BASE["refine_angle"], "chain_summary.csv")):
                    continue    # the base setting is the angle_0.2 run
                jobs.append((folder, spec))
    # the long fine-mesh runs first, so the pool does not end on them
    jobs.sort(key=lambda j: (j[1]["remesher"]["refine_angle"], j[1]["remesher"]["size_max"]))
    print("%d runs, %d workers" % (len(jobs), args.workers), flush=True)
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        for folder, ok, sec in pool.map(one_run, jobs):
            print("%-60s %s %.0f s" % (os.path.relpath(folder, ROOT), "ok" if ok else "FAILED", sec), flush=True)
    collect()


# ---- the remesh stress test: meshes on which the remesher of the DEFAULT settings (size_min 0.0005) blows up, and which settings survive them
STRESS_DIR = "../Studies/shape_chains/stress_meshes"
STRESS = {
    # name: c0, nu of the link, seed, steps to run before the mesh is taken (None: the converged final state), what it is
    "pear_c2.4": dict(c0=2.4, nu=0.72, seed={"spheres": {"q": 1.3, "neck": 0.2}}, steps=None, note="pear at c0 = 2.4, v = 0.72 after BFGS-Normal (the end state of link 0 of the pear chain: its neck mesh has corners down to 12 deg)"),
    "pear_c3": dict(c0=3.0, nu=0.78, seed={"spheres": {"q": 2.0, "neck": 0.3}}, steps=8300, note="pear at c0 = 3, v = 0.78 in BFGS-Normal, the last saved mesh before the first sliver guard (step 8305)"),
    "prolate2_c2.4": dict(c0=2.4, nu=0.70, seed={"capsule": {"nu": 0.55}}, steps=None, note="control: prolate-2 at c0 = 2.4, v = 0.70 (neck ratio 0.17), a healthy neck mesh"),
}
STRESS_SETTINGS = ([("angle", dict(refine_angle=a)) for a in (0.4, 0.3, 0.2, 0.15, 0.1)] + [("size_min", dict(size_min=m)) for m in (0.0005, 0.002, 0.005, 0.01, 0.02, 0.05)]
                   + [("aspect", dict(aspect_min=a)) for a in (0.1, 0.2, 0.3, 0.4, 0.5)] + [("size_max", dict(size_max=m)) for m in (0.3, 0.15, 0.1, 0.07)])
OLD_DEFAULT = {"size_max": 0.15, "size_min": 0.0005, "refine_angle": 0.2, "aspect_min": 0.2}


def stress_make(args):
    """Make the stress meshes with the old default remesher (size_min 0.0005), which still reaches them from a clean seed, and save them in STRESS_DIR."""
    os.makedirs(STRESS_DIR, exist_ok=True)
    for name, st in STRESS.items():
        if args.only and name not in args.only:
            continue
        dest = os.path.join(STRESS_DIR, name + ".obj")
        if os.path.exists(dest) and not args.redo:
            continue
        folder = os.path.abspath(os.path.join(ROOT, "stress", name))
        spec = {"config": os.path.abspath(TEMPLATE), "binary": os.path.abspath(BINARY), "batch_dir": folder, "c0": st["c0"], "KB": 1.0, "V0": 0.5035, "dA": 4e-4, "seed": st["seed"],
                "nu0": st["nu"], "legs": [], "switch": [1000, 4000], "timesteps": 60000, "threads": 1, "nice": 0, "remesher": OLD_DEFAULT}
        os.makedirs(folder, exist_ok=True)
        path = os.path.join(folder, "study_spec.json")
        json.dump(spec, open(path, "w"), indent=1)
        ch = R.Chain(spec, path)
        if st["steps"] is None:
            ok = ch.run()
            row = ch.cache.get("link0")
            if not ok or row is None or row["kind"] == "failed":
                print(name, "FAILED")
                continue
            shutil_copy(row["mesh"], dest)
        else:
            # run st["steps"] steps with a mesh every 250 and take the last one before the first sliver polish
            if not os.path.exists(os.path.join(folder, "seed.obj")):
                R.make_seed(st["seed"], 0.5035, os.path.join(folder, "seed.obj"))
            cfg, _ = ch.link_config(os.path.join(folder, "seed.obj"), os.path.join(folder, "run") + "/", st["nu"])
            cfg["timesteps"], cfg["save_interval"] = st["steps"], 250
            cfg_path = os.path.join(folder, "run.json")
            json.dump(cfg, open(cfg_path, "w"), indent=1)
            R.run_main_cluster(ch.binary, cfg_path, os.path.join(folder, "run.log"), 1, 0)
            run = R.find_run_folder(os.path.join(folder, "run"))
            steps = sorted(int(f[9:-4]) for f in os.listdir(run) if f.startswith("membrane_") and f.endswith(".obj"))
            first = None
            if os.path.exists(run + "Polish_log.txt"):
                for line in open(run + "Polish_log.txt"):
                    if "slivers" in line:
                        first = int(line.split()[0])
                        break
            use = max([t for t in steps if first is None or t < first])
            print(name, "mesh of step", use, "(first sliver polish at", first, ")")
            shutil_copy(run + "membrane_%d.obj" % use, dest)
        open(os.path.join(STRESS_DIR, name + ".txt"), "w").write("c0 %g, v %g: %s\n" % (st["c0"], st["nu"], st["note"]))
        print("wrote", dest)


def shutil_copy(src, dst):
    import shutil
    shutil.copy(src, dst)


def stress_one(job):
    name, (group, over), cfg, tag_dir = job
    c = json.loads(json.dumps(cfg))
    d = tag_dir + "/"
    os.makedirs(d, exist_ok=True)
    c["first_dir"] = d
    c["remesher"].update(over)
    json.dump(c, open(d + "cfg.json", "w"))
    import re
    import subprocess
    p = subprocess.run([os.path.abspath(BINARY), d + "cfg.json"], capture_output=True, text=True, env=dict(os.environ, OMP_NUM_THREADS="1"))
    text = p.stdout
    undone = len(re.findall(r"undone", text))
    first = re.search(r"Remesh at t = 0 undone \(Bending: (\S+) -> (\S+)\)", text)
    ok = re.search(r"Final state saved as step (\d+), total energy (\S+)", text)
    v = re.findall(r"vertices (\d+) mean edge length (\S+)", text)
    return {"mesh": name, "group": group, "setting": ",".join("%s=%g" % kv for kv in over.items()), "undone": undone, "first_ratio": (float(first.group(2)) / float(first.group(1))) if first else 1.0,
            "E": float(ok.group(2)) if ok else float("nan"), "nV": int(v[-1][0]) if v else 0, "edge": float(v[-1][1]) if v else float("nan")}


def stress_robust(args):
    """150 steps of BFGS (remeshing every 5 steps) from every stress mesh for every setting of the remesher: how often the remesh is undone (rollback) because it makes the energy explode."""
    from concurrent.futures import ThreadPoolExecutor
    jobs = []
    for name, st in STRESS.items():
        mesh = os.path.abspath(os.path.join(STRESS_DIR, name + ".obj"))
        if not os.path.exists(mesh):
            print("missing", mesh, "(python3 Study_shape_mesh.py stress_make)")
            continue
        folder = os.path.abspath(os.path.join(ROOT, "stress", name))
        os.makedirs(folder, exist_ok=True)
        spec = {"config": os.path.abspath(TEMPLATE), "binary": os.path.abspath(BINARY), "batch_dir": folder, "c0": st["c0"], "KB": 1.0, "V0": 0.5035, "dA": 4e-4, "seed": st["seed"],
                "nu0": st["nu"], "legs": [], "switch": [1, 4000], "timesteps": 150, "threads": 1, "nice": 0, "remesher": OLD_DEFAULT}
        path = os.path.join(folder, "robust_spec.json")
        json.dump(spec, open(path, "w"), indent=1)
        ch = R.Chain(spec, path)
        for group, over in STRESS_SETTINGS:
            tag = "robust_%s_%s" % (group, "_".join("%s%g" % kv for kv in over.items()))
            cfg, _ = ch.link_config(mesh, os.path.join(folder, tag) + "/", st["nu"])
            cfg["timesteps"], cfg["save_interval"], cfg["Switch_times"] = 150, 1000000, [1, -2500, -1, -1, -1, -1]
            cfg["remesher"] = dict(OLD_DEFAULT)
            jobs.append((name, (group, over), cfg, os.path.join(folder, tag)))
    with ThreadPoolExecutor(max_workers=args.workers) as ex:
        res = list(ex.map(stress_one, jobs))
    path = os.path.join(ROOT, "stress_robust.csv")
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        keys = ["mesh", "group", "setting", "undone", "first_ratio", "E", "nV", "edge"]
        w.writerow(keys)
        for r in res:
            w.writerow([("%.6g" % r[k] if isinstance(r[k], float) else r[k]) for k in keys])
    print("wrote", path)
    try:
        stress_plot(args)
    except ImportError:
        print("no matplotlib here: run `Study_shape_mesh.py stress_plot --out %s` where there is one" % args.out)


def stress_plot(args):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    rows = list(csv.DictReader(open(os.path.join(ROOT, "stress_robust.csv"))))
    settings = []
    for r in rows:
        if r["setting"] not in settings:
            settings.append(r["setting"])
    meshes = [m for m in STRESS if any(r["mesh"] == m for r in rows)]
    fig, ax = plt.subplots(figsize=(1.6 * len(meshes) + 3.2, 0.26 * len(settings) + 1.4))
    for j, m in enumerate(meshes):
        for i, st in enumerate(settings):
            rr = [r for r in rows if r["mesh"] == m and r["setting"] == st]
            if not rr:
                continue
            u = int(rr[0]["undone"])
            ax.plot([j], [i], "s", color="#2ca02c" if u == 0 else "#d62728", ms=11)
            ax.text(j + 0.22, i, "%s" % ("ok, %s vertices" % rr[0]["nV"] if u == 0 else "undone %d, x%.0f" % (u, float(rr[0]["first_ratio"]))), va="center", fontsize=6.5)
    ax.set_xticks(range(len(meshes)))
    ax.set_xticklabels(meshes, fontsize=8)
    ax.set_yticks(range(len(settings)))
    ax.set_yticklabels(settings, fontsize=7)
    ax.set_xlim(-0.4, len(meshes) - 0.1 + 0.9)
    ax.invert_yaxis()
    ax.set_title("remesh of a thin-neck mesh: green = it works, red = undone (Bending x factor at t = 0)", fontsize=8)
    fig.tight_layout()
    fig.savefig(args.out + "_stress.png", dpi=170)
    plt.close(fig)
    print("wrote", args.out + "_stress.png")


def mesh_stats(path):
    V, F = R.read_obj(path)
    a, b, c = V[F[:, 0]], V[F[:, 1]], V[F[:, 2]]
    e = np.concatenate([np.linalg.norm(a - b, axis=1), np.linalg.norm(b - c, axis=1), np.linalg.norm(c - a, axis=1)])
    ab, ac, bc = b - a, c - a, c - b
    cosA = (ab * ac).sum(1) / (np.linalg.norm(ab, axis=1) * np.linalg.norm(ac, axis=1))
    cosB = (-ab * bc).sum(1) / (np.linalg.norm(ab, axis=1) * np.linalg.norm(bc, axis=1))
    cosC = (-ac * -bc).sum(1) / (np.linalg.norm(ac, axis=1) * np.linalg.norm(bc, axis=1))
    ang = np.degrees(np.arccos(np.clip(np.concatenate([cosA, cosB, cosC]), -1, 1)))
    return {"edge_mean": float(e.mean()), "edge_min": float(e.min()), "edge_max": float(e.max()), "angle_min": float(ang.min()), "angle_max": float(ang.max())}


def collect():
    out = []
    for case in CASES:
        base = os.path.join(ROOT, case)
        if not os.path.isdir(base):
            continue
        for name in sorted(os.listdir(base)):
            folder = os.path.join(base, name)
            summ = os.path.join(folder, "chain_summary.csv")
            if not os.path.exists(summ):
                continue
            rows = list(csv.DictReader(open(summ)))
            row = rows[0]
            spec = json.load(open(os.path.join(folder, "study_spec.json")))
            rem = dict(BASE, **spec["remesher"])
            rec_extra = {"integrator": "bfgs_normal" if spec.get("bfgs_normal", True) else "bfgs", "remesh_every": spec.get("config_overrides", {}).get("remesh_every", 5)}
            rec = {"case": case, "run": name, "kind": row["kind"], "note": row["note"], "size_max": rem["size_max"], "size_min": rem["size_min"], "refine_angle": rem["refine_angle"],
                   "aspect_min": rem["aspect_min"], **rec_extra}
            log = os.path.join(folder, "link0.log")
            text = open(log, errors="replace").read() if os.path.exists(log) else ""
            rec["rollbacks"] = text.count("undone")
            rec["polishes"] = text.count(": slivers at")
            if row["kind"] != "failed" and row["F"] != "":
                rec.update({k: float(row[k]) for k in ("nu", "F", "c0_eff", "steps", "seconds", "nV", "rneck", "neck_ratio", "lobe_ratio", "skew", "length")})
                rec["status"] = row["status"]
                rec.update(mesh_stats(row["mesh"]))
            out.append(rec)
    keys = ["case", "run", "kind", "status", "note", "integrator", "remesh_every", "refine_angle", "size_max", "size_min", "aspect_min", "nu", "c0_eff", "F", "steps", "seconds", "nV", "edge_mean", "edge_min", "edge_max", "angle_min",
            "angle_max", "rneck", "neck_ratio", "lobe_ratio", "skew", "length", "rollbacks", "polishes"]
    with open(os.path.join(ROOT, "mesh_study.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(keys)
        for r in out:
            w.writerow([("%.6g" % r[k] if isinstance(r.get(k), float) else r.get(k, "")) for k in keys])
    print("wrote", os.path.join(ROOT, "mesh_study.csv"), len(out), "runs")


def load():
    rows = list(csv.DictReader(open(os.path.join(ROOT, "mesh_study.csv"))))
    for r in rows:
        for k, v in r.items():
            if k not in ("case", "run", "kind", "status", "note", "integrator"):
                r[k] = float(v) if v != "" else float("nan")
    return rows


def plot(args):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    rows = load()
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    colors = {"sphere": "k", "prolate1": "#1f77b4", "prolate2": "#2ca02c", "pear": "#d62728"}
    ok = lambda r: r["kind"] != "failed" and not math.isnan(r["F"])

    # 1. the energy against the resolution (angle group)
    fig, axs = plt.subplots(2, 2, figsize=(9, 6.6))
    ref = {}
    for case in CASES:
        pts = sorted((r for r in rows if r["case"] == case and r["run"].startswith("angle_") and ok(r)), key=lambda r: r["refine_angle"])
        if not pts:
            continue
        th = np.array([r["refine_angle"] for r in pts])
        F = np.array([r["F"] for r in pts])
        if "exact" in CASES[case]:
            ref[case] = CASES[case]["exact"]
        else:
            # Richardson: F = F0 + c th^2 through the three finest
            fine = np.argsort(th)[:3]
            c, f0 = np.polyfit(th[fine] ** 2, F[fine], 1)
            ref[case] = f0
        axs[0, 0].plot(th ** 2, (F - ref[case]) / ref[case] * 100, "o-", color=colors[case], label=case + " (F_ref %.4f)" % ref[case])
        axs[0, 1].plot(th, [r["nV"] for r in pts], "o-", color=colors[case], label=case)
        axs[1, 0].plot(th, [r["seconds"] for r in pts], "o-", color=colors[case])
        axs[1, 1].plot(th, [r["rneck"] for r in pts], "o-", color=colors[case]) if case in ("prolate2", "pear") else None
    axs[0, 0].set_xlabel(r"refine_angle$^2$")
    axs[0, 0].set_ylabel("(F - F_ref) / F_ref  [%]")
    axs[0, 0].axhline(0, color="0.7", lw=0.6)
    axs[0, 0].legend(fontsize=7)
    axs[0, 1].set_xlabel("refine_angle")
    axs[0, 1].set_ylabel("vertices")
    axs[0, 1].set_yscale("log")
    axs[1, 0].set_xlabel("refine_angle")
    axs[1, 0].set_ylabel("seconds (one core)")
    axs[1, 0].set_yscale("log")
    axs[1, 1].set_xlabel("refine_angle")
    axs[1, 1].set_ylabel("neck radius of the relaxed shape")
    for a in axs.ravel():
        a.set_xscale("linear")
    fig.tight_layout()
    fig.savefig(args.out + "_angle.png", dpi=170)
    plt.close(fig)

    # 2. the other groups: F and the neck against the setting
    fig, axs = plt.subplots(2, 3, figsize=(11, 6))
    for j, (group, key) in enumerate((("size_max", "size_max"), ("size_min", "size_min"), ("aspect", "aspect_min"))):
        for case in GROUP_CASES[group]:
            pts = sorted((r for r in rows if r["case"] == case and r["run"].startswith(group + "_") and ok(r)), key=lambda r: r[key])
            base = [r for r in rows if r["case"] == case and r["run"] == "angle_0.2" and ok(r)]
            pts = sorted(pts + base, key=lambda r: r[key]) if base else pts
            pts = [r for i, r in enumerate(pts) if i == 0 or r[key] != pts[i - 1][key]]
            if not pts:
                continue
            axs[0, j].plot([r[key] for r in pts], [(r["F"] - ref.get(case, r["F"])) / ref.get(case, r["F"]) * 100 for r in pts], "o-", color=colors[case], label=case)
            axs[1, j].plot([r[key] for r in pts], [r["rneck"] for r in pts], "o-", color=colors[case])
            bad = [r for r in rows if r["case"] == case and r["run"].startswith(group + "_") and not ok(r)]
            for r in bad:
                axs[0, j].plot([r[key]], [0], "x", color=colors[case], ms=9)
        axs[0, j].set_xlabel(key)
        axs[0, j].set_ylabel("(F - F_ref) / F_ref  [%] (x: the run failed)")
        axs[1, j].set_xlabel(key)
        axs[1, j].set_ylabel("neck radius")
        if group == "size_min":
            axs[0, j].set_xscale("log")
            axs[1, j].set_xscale("log")
        axs[0, j].legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(args.out + "_settings.png", dpi=170)
    plt.close(fig)

    # 3. robustness: which settings finish, which are capped (the step cap) or fail, with the remesh rollbacks and sliver polishes of the run
    names = []
    for r in rows:
        if r["run"] not in names:
            names.append(r["run"])
    cases = [c for c in CASES if any(r["case"] == c for r in rows)]
    fig, ax = plt.subplots(figsize=(max(7, 0.55 * len(names)), 0.9 * len(cases) + 1.6))
    for j, name in enumerate(names):
        for i, case in enumerate(cases):
            rr = [r for r in rows if r["case"] == case and r["run"] == name]
            if not rr:
                continue
            r = rr[0]
            color = "#d62728" if r["kind"] == "failed" else ("#ff7f0e" if r.get("status") == "capped" else "#2ca02c")
            ax.plot([j], [i], "s", color=color, ms=13)
            if r["rollbacks"] > 0 or r.get("polishes", 0) > 0:
                ax.text(j, i, "%d\n%d" % (r["rollbacks"], r.get("polishes", 0)), ha="center", va="center", fontsize=5.5, color="white")
    ax.set_xticks(range(len(names)))
    ax.set_xticklabels(names, rotation=70, ha="right", fontsize=7)
    ax.set_yticks(range(len(cases)))
    ax.set_yticklabels(cases)
    ax.set_title("green: converged, orange: step cap, red: failed (numbers: rollbacks / sliver polishes)", fontsize=8)
    ax.margins(0.04)
    fig.tight_layout()
    fig.savefig(args.out + "_robustness.png", dpi=170)
    plt.close(fig)

    # the table
    with open(args.out + "_table.md", "w") as f:
        f.write("| case | run | F | nu | steps | seconds | vertices | mean edge | neck radius | rollbacks | polishes | status |\n|---|---|---|---|---|---|---|---|---|---|---|---|\n")
        for r in rows:
            if ok(r):
                f.write("| %s | %s | %.5f | %.4f | %d | %.0f | %d | %.4f | %.4f | %d | %d | %s |\n" % (r["case"], r["run"], r["F"], r["nu"], r["steps"], r["seconds"], r["nV"], r["edge_mean"],
                                                                                                           r["rneck"], r["rollbacks"], r.get("polishes", 0), r["status"]))
            else:
                f.write("| %s | %s | FAILED | | | | | | | %d | %d | %s |\n" % (r["case"], r["run"], r["rollbacks"], r.get("polishes", 0), r["note"]))
    print("wrote", args.out + "_angle.png,", args.out + "_settings.png,", args.out + "_robustness.png,", args.out + "_table.md")


def main():
    global BINARY
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run")
    r.add_argument("--workers", type=int, default=6)
    r.add_argument("--cases", nargs="*", default=None, choices=list(CASES))
    r.add_argument("--groups", nargs="*", default=list(GROUPS), choices=list(GROUPS))
    r.add_argument("--cap", type=int, default=30000, help="largest number of steps of a run (a run that reaches it is marked capped)")
    r.add_argument("--binary", default=BINARY, help="the main_cluster to run (a copy keeps the study on one binary while the source changes)")
    r.add_argument("--redo", action="store_true", help="repeat the runs that have a chain_summary.csv")
    p = sub.add_parser("plot")
    p.add_argument("--out", required=True)
    c = sub.add_parser("collect")
    m = sub.add_parser("stress_make", help="make the stress meshes of the remesh test (STRESS_DIR)")
    m.add_argument("--only", nargs="*", default=None)
    m.add_argument("--redo", action="store_true")
    m.add_argument("--binary", default=BINARY)
    rb = sub.add_parser("stress_robust", help="150 steps from every stress mesh for every remesher setting: is the remesh undone?")
    rb.add_argument("--workers", type=int, default=4)
    rb.add_argument("--out", required=True)
    rb.add_argument("--binary", default=BINARY)
    sp = sub.add_parser("stress_plot")
    sp.add_argument("--out", required=True)
    args = ap.parse_args()
    if getattr(args, "binary", None):
        BINARY = args.binary
    if args.cmd in ("run", "stress_make", "stress_robust"):
        ensure_template()
    if args.cmd == "stress_make":
        stress_make(args)
    elif args.cmd == "stress_robust":
        stress_robust(args)
    elif args.cmd == "stress_plot":
        stress_plot(args)
    elif args.cmd == "run":
        run(args)
    elif args.cmd == "plot":
        plot(args)
    else:
        collect()


if __name__ == "__main__":
    main()
