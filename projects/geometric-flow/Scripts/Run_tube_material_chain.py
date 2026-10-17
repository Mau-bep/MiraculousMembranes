#!/usr/bin/env python3
"""A bead walked through the big planar membrane along x, one frozen position at a time: every position is a relaxation warm-started from the previous one.

usage:  python3 Run_tube_material_chain.py CHAIN.json          (CHAIN.json is written by Create_subjob_tube_material.py)

The material-property sister of Mem_tube_up.sh (there every final x is a separate job that pulls a bead there, here one chain walks the bead out and back):

0. the START STATE (CHAIN.json "start_config"): the old method of Mem_tube_up.sh, a bead pulled by a spring through the flat disk that freezes by itself at x_start
   (default 0.2); its Final_state.obj is the mesh of the chain. It is checked (the bead froze, the run did not break) and saved in <batch_dir>/start_x<x>/, it is not a data
   point (the spring is still in its energy). "start_run" takes the mesh of a finished run instead; with neither the chain starts from the flat disk (bead at -radius).
1. link 0: the frozen bead (the config of CHAIN.json, no spring) at positions[0] relaxed to convergence (the "stopping" block of the config). It is a data point.
2. the way out: for every further position the previous membrane is the initial mesh, the bead (frozen) is moved to [x, 0, 0] and the membrane relaxes. A move larger than
   "max_step" is split into waypoints (saved, kind "waypoint", left out of the curves). Every result is checked: no nan line search, no stall, the total energy did not
   rise during the relaxation, the bead is not through a face of the mesh. A relaxation that fails is retried through the middle position, then the chain stops.
3. "return": the way back from the last position down to the first one and on to "back_positions"[-1] (0.2 -> 6.0 -> 0.2 -> 0.0), each position relaxed from the state of the way
   back (not the way out): the hysteresis of the tube. The runs of the way back are saved in their own folder, <batch_dir>/way_back/ (leg "out" in chain_summary.csv).

Everything goes to batch_dir: start_x0.2/, link000_x0.2/, link001_x0.3/, ... and way_back/link0NN_x.../ for the way back (each is a normal run folder: Input_file.json, Output_data.txt, Final_state.obj, Bead_0_data.txt,
...), chain_summary.csv (one line per run: kind, leg, x, the folder of the run, the force on the bead, energies, ...) and chain.log. Analyze_tube_material.py <batch_dir> turns them
into the curves. RESUME: the chain continues after the last line of chain_summary.csv when it finds one, so after a time limit or a failure sbatch the same subjob again
(delete chain_summary.csv, or the folder, to start over). Needs only json and numpy, so it runs inside the sbatch job.
"""
import copy
import csv
import json
import math
import os
import re
import subprocess
import sys
import time

import numpy as np

FIELDS = ["link", "kind", "leg", "x", "run", "rundir", "step", "seconds", "nV", "x_tip", "d_bead", "Area", "E_total", "E_bend", "E_tension", "E_bead", "Fx", "Fy", "Fz", "capped"]


# ----------------------------------------------------------------------------------------------------------------------
def read_obj(path):
    V, F = [], []
    with open(path) as f:
        for line in f:
            p = line.split()
            if p and p[0] == "v":
                V.append([float(x) for x in p[1:4]])
            elif p and p[0] == "f":
                F.append([int(x.split("/")[0]) - 1 for x in p[1:4]])
    return np.array(V), np.array(F, dtype=int)


def last_mesh(run):
    final = os.path.join(run, "Final_state.obj")
    if os.path.exists(final):
        return final
    meshes = sorted((f for f in os.listdir(run) if f.startswith("membrane_") and f.endswith(".obj")), key=lambda s: int(s[9:-4]))
    if not meshes:
        raise RuntimeError("no mesh in " + run)
    return os.path.join(run, meshes[-1])


def last_bead_row(run):
    rows = [l.split() for l in open(os.path.join(run, "Bead_0_data.txt")) if l.strip() and not l.startswith("#")]
    return [float(x) for x in rows[-1]]


def closest_points(p, A, B, C):
    """Closest point to p on every triangle (A, B, C): Ericson, Real-Time Collision Detection 5.1.5, vectorised (the regions are applied in reverse order of precedence)."""
    ab, ac, ap = B - A, C - A, p - A
    d1, d2 = (ab * ap).sum(1), (ac * ap).sum(1)
    bp = p - B
    d3, d4 = (ab * bp).sum(1), (ac * bp).sum(1)
    cp = p - C
    d5, d6 = (ab * cp).sum(1), (ac * cp).sum(1)
    vc, vb, va = d1 * d4 - d3 * d2, d5 * d2 - d1 * d6, d3 * d6 - d5 * d4
    tiny = 1e-300
    denom = np.where(np.abs(va + vb + vc) < tiny, tiny, va + vb + vc)
    Q = A + ab * (vb / denom)[:, None] + ac * (vc / denom)[:, None]
    m = (va <= 0) & (d4 - d3 >= 0) & (d5 - d6 >= 0)
    Q = np.where(m[:, None], B + (C - B) * ((d4 - d3) / np.maximum((d4 - d3) + (d5 - d6), tiny))[:, None], Q)
    m = (vb <= 0) & (d2 >= 0) & (d6 <= 0)
    Q = np.where(m[:, None], A + ac * (d2 / np.maximum(d2 - d6, tiny))[:, None], Q)
    m = (d6 >= 0) & (d5 <= d6)
    Q = np.where(m[:, None], C, Q)
    m = (vc <= 0) & (d1 >= 0) & (d3 <= 0)
    Q = np.where(m[:, None], A + ab * (d1 / np.maximum(d1 - d3, tiny))[:, None], Q)
    m = (d3 >= 0) & (d4 <= d3)
    Q = np.where(m[:, None], B, Q)
    m = (d1 <= 0) & (d2 <= 0)
    Q = np.where(m[:, None], A, Q)
    # the regions above assume a proper triangle: a face with two coincident corners (a sliver) would give a point that is too far; the nearest of these and the three edges is exact
    best, best_d = Q, ((Q - p) ** 2).sum(1)
    for P0, P1 in ((A, B), (B, C), (C, A)):
        e = P1 - P0
        ee = (e * e).sum(1)
        t = np.clip(np.where(ee > 0, ((p - P0) * e).sum(1) / np.where(ee > 0, ee, 1.0), 0.0), 0.0, 1.0)
        S = P0 + e * t[:, None]
        sd = ((S - p) ** 2).sum(1)
        better = sd < best_d
        best, best_d = np.where(better[:, None], S, best), np.where(better, sd, best_d)
    return best


def bead_distances(V, F, centre):
    """(distance from the centre of the bead to the nearest point of the membrane, faces included, distance to the nearest vertex). The LJ interaction acts on the VERTICES, so a bead
    can sit inside a face of a coarse mesh without any vertex feeling it: then the first is much smaller than the second. A healthy relaxed state has the vertex distance 0.91 to 1.1 of the
    radius and the face distance 0.94 to 1.08 (0.78 at the lowest, on a coarse flat mesh whose faces are as large as the bead). nan for a non-finite mesh or position."""
    c = np.asarray(centre, float)
    if not (np.isfinite(V).all() and np.isfinite(c).all()) or len(F) == 0:
        return float("nan"), float("nan")
    Q = closest_points(c, V[F[:, 0]], V[F[:, 1]], V[F[:, 2]])
    return float(np.sqrt(((Q - c) ** 2).sum(1).min())), float(np.sqrt(((V - c) ** 2).sum(1).min()))


def bead_distance(V, F, centre):
    """Distance from the centre of the bead to the nearest point of the membrane, faces included (bead_distances)."""
    return bead_distances(V, F, centre)[0]


def bead_problem(V, F, centre, radius):
    """Why the bead is not resting on a healthy membrane, or None. Limits from 78 healthy relaxed states of the old pulling runs (README): the nearest vertex is not in the LJ core
    (at least 0.85 radius), the centre is not much closer to a face than to any vertex (a bead inside a face of the mesh: at least 0.7 of the vertex distance; a flat mesh
    with faces as large as the bead has 0.84 at the lowest), the membrane is within 1.4 radius (the bead has not lost it) and lies in front of the bead (the bead is not on the wrong
    side of a flat membrane: the tip of the membrane is at least 0.5 radius beyond the centre). Every comparison is written so that a nan fails it."""
    d, dv = bead_distances(V, F, centre)
    if not (dv >= 0.85 * radius):
        return "a vertex is inside the LJ core of the bead (%.3f from its centre, its radius is %g)" % (dv, radius)
    if not (d >= 0.7 * dv):
        return "the bead is inside a face of the mesh (%.3f from its centre to the membrane, %.3f to the nearest vertex, its radius is %g)" % (d, dv, radius)
    if not (d <= 1.4 * radius):
        return "the bead has lost the membrane (%.3f from its centre, its radius is %g)" % (d, radius)
    if not (np.asarray(V)[:, 0].max() - centre[0] >= 0.5 * radius):
        return "the membrane is not in front of the bead (its tip is %.3f beyond the centre, the radius is %g)" % (np.asarray(V)[:, 0].max() - centre[0], radius)
    return None


def read_energy_rows(run):
    with open(os.path.join(run, "Output_data.txt")) as f:
        header = f.readline().split()
        rows = [[float(x) for x in l.split()] for l in f if l.strip()]
    col = {}
    for i, name in enumerate(header):
        col.setdefault(name, i)
    return header, col, rows


def summary(run, radius=None):
    """The numbers of one run (energies of the last row of Output_data.txt in the energy unit of the config, the force on the bead from the last row of Bead_0_data.txt: step,
    x, y, z, Fx, Fy, Fz), the position of the tip of the tube (largest x of the mesh) and the distance of the bead to the mesh."""
    header, col, rows = read_energy_rows(run)
    last = rows[-1]
    bead = last_bead_row(run)
    V, F = read_obj(last_mesh(run))
    area = 0.5 * np.linalg.norm(np.cross(V[F[:, 1]] - V[F[:, 0]], V[F[:, 2]] - V[F[:, 0]]), axis=1).sum()
    beads = [i for i, name in enumerate(header) if name == "Bead"]
    get = lambda name: last[col[name]] if name in col else float("nan")
    return {"step": last[col["step"]], "nV": len(V), "x_tip": float(V[:, 0].max()), "d_bead": bead_distance(V, F, bead[1:4]), "Area": area,
            "E_total": get("Total_E"), "E_bend": get("Bending"), "E_tension": get("Surface_tension"), "E_bead": sum(last[i] for i in beads),
            "Fx": bead[4], "Fy": bead[5], "Fz": bead[6]}


# ----------------------------------------------------------------------------------------------------------------------
def run_main_cluster(binary, config_path, log_path, nice):
    cmd = [binary, config_path]
    if nice:
        cmd = ["nice", "-n", str(nice)] + cmd
    with open(log_path, "w") as log:
        proc = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT)
    return proc.returncode


def check_run(log_path, run, radius, check_bead=True, check_energy=True):
    """None when the relaxation is usable, else why not. main_cluster returns 0 and writes a Final_state.obj also when the run broke, so the log, the energies and the
    geometry are read (check_energy=False for the pull of the start state, whose energy rises on purpose): no nan line search, no stalled BFGS-Normal that was ended, the total energy did not RISE during the run (a minimisation that ends above where it
    started, a vertex inside the bead after a remesh, is not a result), the energy and the force are finite and sane, and the bead rests on the membrane (bead_problem).
    'Sim broke or timestep very small' alone is not a failure (a membrane already at its minimum ends gradient descent that way). Unreadable output is a failure too."""
    try:
        text = open(log_path, errors="replace").read()
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
        header, col, rows = read_energy_rows(run)
        e0, e1 = rows[0][col["Total_E"]], rows[-1][col["Total_E"]]
        if not math.isfinite(e1):
            return "total energy %g" % e1
        if check_energy and math.isfinite(e0) and e1 > e0 + max(5.0, 3e-3 * abs(e0)):
            return "the total energy rose from %.6g to %.6g during the run" % (e0, e1)
        bead = last_bead_row(run)
        if not all(math.isfinite(v) for v in bead) or abs(bead[4]) > 1e4:
            return "the force on the bead is %s" % bead[4:7]
        if check_bead:
            V, F = read_obj(last_mesh(run))
            return bead_problem(V, F, bead[1:4], radius)
    except (OSError, ValueError, IndexError, KeyError, RuntimeError) as exc:
        return "unreadable output: %r" % (exc,)
    return None


def find_run_folder(first_dir):
    """The folder main_cluster made inside first_dir (the first free number: 1 in a fresh directory, the highest one after a retry)."""
    if not os.path.isdir(first_dir):
        raise RuntimeError("no folder " + first_dir)
    nums = sorted(int(d) for d in os.listdir(first_dir) if d.isdigit() and os.path.isdir(os.path.join(first_dir, d)))
    if not nums:
        raise RuntimeError("main_cluster made no run folder in " + first_dir)
    return os.path.join(first_dir, str(nums[-1])) + "/"


def numbered(first_dir):
    return {int(d) for d in os.listdir(first_dir) if d.isdigit() and os.path.isdir(os.path.join(first_dir, d))} if os.path.isdir(first_dir) else set()


def new_run_folder(first_dir, before):
    """The folder the last main_cluster call made in first_dir (it takes the first free number, which is not always the highest: a stale folder may sit above a gap). None when it made none."""
    new = numbered(first_dir) - before
    return os.path.join(first_dir, str(max(new))) + "/" if new else None


def try_find_run_folder(first_dir):
    try:
        return find_run_folder(first_dir)
    except RuntimeError:
        return None


# ----------------------------------------------------------------------------------------------------------------------
def main(chain_path):
    spec = json.load(open(chain_path))
    binary = os.path.abspath(spec.get("binary", "../build/bin/main_cluster"))
    batch = os.path.abspath(spec["batch_dir"])
    os.makedirs(batch, exist_ok=True)
    nice = spec.get("nice", 0)
    switch = spec.get("switch", [100, 1500])
    max_step = spec.get("max_step")
    link_stopping, link_timesteps = spec.get("link_stopping"), spec.get("link_timesteps")
    positions = [float(x) for x in spec["positions"]]
    # the way back: the list backwards to positions[0] and on below it to x_back (0.2 -> 6.0 -> 0.2 -> 0.0); its runs are saved in their own folder, back_dir inside batch_dir
    back_positions = [float(x) for x in spec["back_positions"]] if "back_positions" in spec else list(reversed(positions[:-1]))
    back_dir = os.path.join(batch, spec.get("back_dir", "way_back"))
    template = json.load(open(spec["config"]))      # the frozen bead: every run of the chain is built from it (the spring of the start run never gets in)
    radius = float(spec.get("radius", template["Beads"][0]["radius"]))
    log = open(os.path.join(batch, "chain.log"), "a")

    def say(msg):
        print(msg, flush=True)
        log.write("%s %s\n" % (time.strftime("%Y-%m-%d %H:%M:%S"), msg))
        log.flush()

    summary_path = os.path.join(batch, "chain_summary.csv")
    rows = []

    def write_summary():
        tmp = summary_path + ".tmp"
        with open(tmp, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(FIELDS)
            for r in rows:
                w.writerow([r.get(k, "") for k in FIELDS])
            f.flush()
            os.fsync(f.fileno())
        os.replace(tmp, summary_path)

    # the schedule: link 0 at positions[0], then the way out, then the way back (ends at positions[0] again)
    sched = [(0, positions[0], "in")] + [(k, x, "in") for k, x in enumerate(positions[1:], start=1)]
    if spec.get("return", True):
        sched += [(len(positions) + j, x, "out") for j, x in enumerate(back_positions)]
    final_tag = lambda k, x: "link%03d_x%g" % (k, x)

    def relax_config(init_mesh, first_dir, x, cold):
        cfg = copy.deepcopy(template)
        cfg["init_file"] = os.path.abspath(init_mesh)
        cfg["rescale"] = 1.0
        cfg["first_dir"] = first_dir
        if not cold:
            cfg["Switch_times"] = switch
            if link_stopping:
                cfg["stopping"].update(link_stopping)
            if link_timesteps:
                cfg["timesteps"] = int(link_timesteps)
        cfg["Beads"][0]["Pos"] = [float(x), 0.0, 0.0]
        cfg["Beads"][0]["state"] = "froze"
        return cfg

    def run_one(k, kind, tag, x, leg, init_mesh, cold=False):
        """One relaxation of the frozen bead at [x, 0, 0] from init_mesh. Returns the row (with 'mesh' = its Final_state.obj) or None when main_cluster failed or the result is not usable."""
        base = back_dir if leg == "out" else batch                     # the way back has its own folder: runs, configs and logs
        os.makedirs(base, exist_ok=True)
        fdir = os.path.join(base, tag) + "/"
        os.makedirs(fdir, exist_ok=True)
        cfg = relax_config(init_mesh, fdir, x, cold)
        cfg_path = os.path.join(base, tag + ".json")
        json.dump(cfg, open(cfg_path, "w"), indent=4)
        say("link %03d: x = %g (%s%s)" % (k, x, "way back" if leg == "out" else "way out", ", " + kind if kind != "point" else ""))
        t0 = time.time()
        log_path = os.path.join(base, tag + ".log")
        before = numbered(fdir)
        rc = run_main_cluster(binary, cfg_path, log_path, nice)
        run = new_run_folder(fdir, before)
        why = "main_cluster returned %d" % rc if rc != 0 else check_run(log_path, run, radius)
        if why:
            say("link %03d failed at x = %g: %s, see %s.log" % (k, x, why, tag))
            return None
        try:
            r = summary(run)
        except (OSError, ValueError, IndexError, KeyError, RuntimeError) as exc:
            say("link %03d failed at x = %g: unreadable output: %r, see %s.log" % (k, x, exc, tag))
            return None
        r.update({"link": k, "kind": kind, "leg": leg, "x": x, "run": tag, "rundir": os.path.relpath(run, batch), "seconds": round(time.time() - t0),
                  "capped": int(r["step"] >= cfg["timesteps"] - 1), "mesh": os.path.join(run, "Final_state.obj")})
        rows.append(r)
        write_summary()
        say("link %03d done: x = %g, E_total = %.4f, Fx = %.4f, x_tip = %.3f, bead %.3f from the membrane, %d steps, %d s%s" % (
            k, x, r["E_total"], r["Fx"], r["x_tip"], r["d_bead"], r["step"], r["seconds"], "  (hit the cap of %d steps: not converged)" % cfg["timesteps"] if r["capped"] else ""))
        return r

    def relax_with_retry(k, kind, tag, x, leg, prev):
        """run_one, and when it fails once more in two half moves through the middle position (the bead moved too far for the membrane)."""
        r = run_one(k, kind, tag, x, leg, prev["mesh"])
        if r is not None:
            return r
        x_mid = round(0.5 * (prev["x"] + x), 9)
        say("link %03d: retrying in two moves, through x = %g" % (k, x_mid))
        mid = run_one(k, "mid", "link%03d_m_x%g" % (k, x_mid), x_mid, leg, prev["mesh"])
        if mid is None:
            return None
        return run_one(k, kind, tag + "_retry", x, leg, mid["mesh"])

    # --- resume after the last line of chain_summary.csv ---------------------------------------------------------
    # the definition of the chain is stored with the first run: the same batch folder with another list of positions (the generator run again with other options) would resume
    # the old chain with the link numbers of the new one
    chain_def = {"positions": positions, "back_positions": back_positions, "max_step": max_step, "return": bool(spec.get("return", True)), "radius": radius, "KA_KB": [e["constants"][0] for e in template["Energies"]],
                 "start": "run" if spec.get("start_run") else ("pull" if spec.get("start_config") else "flat")}
    def_path = os.path.join(batch, "chain_spec.json")
    prev = None
    first = 0
    if os.path.exists(summary_path):
        if os.path.exists(def_path) and json.load(open(def_path)) != chain_def:
            say("chain_summary.csv in %s belongs to another chain (positions, max_step, return, radius, KA, KB or start differ from chain_spec.json): use another batch folder or delete the old one" % batch)
            return 1
        complete = lambda r: all(r.get(k) not in (None, "") for k in FIELDS)     # a row cut by a kill is incomplete
        old = [r for r in csv.DictReader(open(summary_path)) if complete(r)]
        have_mesh = [r for r in old if os.path.exists(os.path.join(batch, r["rundir"], "Final_state.obj"))]
        if old and not have_mesh:
            say("chain_summary.csv has %d finished runs but none of their Final_state.obj exists: delete the summary to start over" % len(old))
            return 1
        if old:
            rows = [dict(r) for r in old]
            last = have_mesh[-1]
            if last is not old[-1]:
                say("warning: the meshes of the last %d runs are gone, resuming from %s" % (len(old) - 1 - old.index(last), last["run"]))
                rows = [dict(r) for r in old[:old.index(last) + 1]]
            prev = {"x": float(last["x"]), "mesh": os.path.join(batch, last["rundir"], "Final_state.obj")}
            k_done = int(float(last["link"]))
            first = k_done + 1 if last["kind"] == "point" else k_done
            say("resuming: %d runs in chain_summary.csv, the last is %s (link %d, %s), continuing with link %d" % (len(rows), last["run"], k_done, last["kind"], first))
    json.dump(chain_def, open(def_path, "w"))

    # --- the start state and link 0 ---------------------------------------------------------------------------------
    if prev is None:
        x0 = positions[0]
        if spec.get("start_run"):
            src = spec["start_run"]
            try:
                mesh = src if src.endswith(".obj") else last_mesh(src)
            except (RuntimeError, OSError):
                mesh = None
            if mesh is None or not os.path.exists(mesh):
                say("start_run %s has no mesh (a run folder with Final_state.obj or membrane_<step>.obj, or an .obj file)" % src)
                return 1
            say("start state: the mesh %s" % mesh)
        elif spec.get("start_config"):
            start_cfg = json.load(open(spec["start_config"]))
            sdir = start_cfg["first_dir"]
            sdir = sdir if os.path.isabs(sdir) else os.path.normpath(os.path.join(os.getcwd(), sdir))
            os.makedirs(sdir, exist_ok=True)
            slog = os.path.join(batch, "start_x%g.log" % x0)
            srun = None
            for n in sorted(numbered(sdir)):                      # a pull that was validated before (marker file) is reused, nothing else
                cand = os.path.join(sdir, str(n)) + "/"
                if os.path.exists(os.path.join(cand, "start_ok")) and os.path.exists(os.path.join(cand, "Final_state.obj")):
                    srun = cand
            if srun is not None:
                say("start state: reusing the validated pull %s" % srun)
            else:
                say("start state: the bead is pulled by the spring from x = %g to x = %g and frozen there -> %s" % (start_cfg["Beads"][0]["Pos"][0], x0, sdir))
                before = numbered(sdir)
                rc = run_main_cluster(binary, spec["start_config"], slog, nice)
                srun = new_run_folder(sdir, before)
                text = open(slog, errors="replace").read()
                if rc != 0 or srun is None or not os.path.exists(os.path.join(srun, "Final_state.obj")):
                    say("the start run failed (return code %d), see %s" % (rc, slog))
                    return 1
                if "Bead has reached the final position" not in text:
                    try:
                        reached = "the pull stopped at x = %.3f" % last_bead_row(srun)[1]
                    except (OSError, ValueError, IndexError):
                        reached = "no bead position"
                    say("the bead did not reach x = %g in %d steps (%s): raise --start_timesteps or check %s" % (x0, start_cfg["timesteps"], reached, slog))
                    return 1
                why = check_run(slog, srun, radius, check_bead=False, check_energy=False)   # the pull raises the energy of the membrane on purpose
                if why:
                    say("the start run is not usable: %s, see %s" % (why, slog))
                    return 1
            try:
                bx = last_bead_row(srun)[1]
                V, F = read_obj(last_mesh(srun))
            except (OSError, ValueError, IndexError, RuntimeError) as exc:
                say("the start run is unreadable: %r" % (exc,))
                return 1
            why = bead_problem(V, F, [bx, 0.0, 0.0], radius)
            say("start state: bead at x = %.4f (the chain starts at %g), %d vertices, %.3f from its centre to the membrane" % (bx, x0, len(V), bead_distance(V, F, [bx, 0.0, 0.0])))
            if not (abs(bx - x0) <= 0.02) or why:
                say("the start state is not usable (bead at %.4f%s)" % (bx, ": " + why if why else ""))
                return 1
            open(os.path.join(srun, "start_ok"), "w").write("validated\n")
            mesh = last_mesh(srun)
        else:
            mesh = template["init_file"]
            say("start state: the flat membrane %s with the bead at x = %g" % (mesh, x0))
        r = run_one(0, "point", final_tag(0, x0), x0, "in", mesh, cold=True)
        if r is None:
            say("giving up at link 000: the first relaxation failed")
            return 1
        prev = r
        first = 1

    # --- the walk ----------------------------------------------------------------------------------------------------
    for k, x, leg in sched:
        if k < first:
            continue
        cur = prev["x"]
        n_steps = max(1, int(math.ceil(abs(x - cur) / max_step - 1e-9))) if max_step else 1
        for i in range(1, n_steps):
            x_way = round(cur + (x - cur) * i / n_steps, 9)
            r = relax_with_retry(k, "waypoint", "link%03d_w_x%g" % (k, x_way), x_way, leg, prev)
            if r is None:
                say("giving up at a waypoint of link %03d (x = %g): fix the cause or sbatch the job again to retry from the last finished run" % (k, x_way))
                return 1
            prev = r
        r = relax_with_retry(k, "point", final_tag(k, x), x, leg, prev)
        if r is None:
            say("giving up at link %03d (x = %g): fix the cause or sbatch the job again to retry from the last finished run" % (k, x))
            return 1
        prev = r
    say("chain finished, %d runs" % len(rows))
    return 0


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise SystemExit(__doc__)
    sys.exit(main(sys.argv[1]))
