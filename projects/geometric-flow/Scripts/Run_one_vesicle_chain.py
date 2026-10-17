#!/usr/bin/env python3
"""One bead wrapped by a vesicle (or a flat disk), scanned in a control parameter: every value is a relaxation warm-started from the previous one.

usage:  python3 Run_one_vesicle_chain.py CHAIN.json          (CHAIN.json is written by Create_subjob_one_vesicle.py)

The sister of Run_two_vesicle_chain.py: there the two beads move along the vesicle, here the bead stays where it is and one number of the config changes
("scan": {"var", "list", "cold", "max_step"} of CHAIN.json): f (wrapped area fraction of the patch), theta (patch angle), W (adhesion), sigma (tension),
excess (relative excess area of a fixed-area vesicle), volume (fraction of the sphere of radius R_v) or cov (target of the Coverage energy);
wt and st are the adhesion and the tension in the dimensionless w~ = 4 W a^2/KB and sigma~ = 2 sigma a^2/KB of the Deserno mapping.

1. link 0 (cold start): optionally a quick wrap with the stiff Coverage energy (the "prewrap_config" of CHAIN.json), then main_cluster on the config of CHAIN.json
   (the tangent start of Templates/Wrapping_one_vesicle.txt at the cold value of the parameter) to convergence (the "stopping" block of the config). With
   "start_run": "<run folder>" in CHAIN.json a finished run is the link 0 instead (set scan.cold to the value of that run).
2. the links: the previous membrane (Final_state.obj) is the initial mesh, the bead stays frozen where it was (where a free cold start left it), the patch axis
   (if any) points to the centre of the vesicle again (away from it for a bead inside, down on the disk), the volume and area targets are absolute, and the
   parameter takes the next value of the list. A change larger than "max_step" is split into waypoints (the contact line of a patch advances by a * dtheta
   per relaxation: a large step leaves the membrane behind the edge of the patch). The first link repeats the cold value when the list starts with it:
   it is a relaxed point like the others, the cold start is not one.
3. "return": the list is run backwards too (a hysteresis check: the two branches must agree within the noise).

Everything goes to batch_dir: link00_cold/, link01_f0.3/, ... (each is a normal run folder: Input_file.json, Output_data.txt, Final_state.obj, Bead_0_data.txt, ...),
chain_summary.csv (one line per run) and chain.log. Analyze_one_vesicle.py <batch_dir> turns the folders into the curves. Needs only json and numpy, so it runs inside
the sbatch job.
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


def volume_centre(V, F):
    a, b, c = V[F[:, 0]], V[F[:, 1]], V[F[:, 2]]
    trip = (a * np.cross(b, c)).sum(1)
    vol = trip.sum() / 6.0
    return (((a + b + c) / 4.0) * (trip / 6.0)[:, None]).sum(0) / vol, vol


def mesh_area(V, F):
    return 0.5 * np.linalg.norm(np.cross(V[F[:, 1]] - V[F[:, 0]], V[F[:, 2]] - V[F[:, 0]]), axis=1).sum()


def last_mesh(run):
    final = os.path.join(run, "Final_state.obj")
    if os.path.exists(final):
        return final
    meshes = sorted((f for f in os.listdir(run) if f.startswith("membrane_") and f.endswith(".obj")), key=lambda s: int(s[9:-4]))
    if not meshes:
        raise RuntimeError("no mesh in " + run)
    return os.path.join(run, meshes[-1])


def last_bead_positions(run, n):
    out = []
    for i in range(n):
        rows = [l.split() for l in open(os.path.join(run, "Bead_%d_data.txt" % i)) if l.strip() and not l.startswith("#")]
        out.append(np.array([float(x) for x in rows[-1][1:4]]))
    return out


def summary(run, cfg):
    with open(os.path.join(run, "Output_data.txt")) as f:
        header = f.readline().split()
        rows = [[float(x) for x in l.split()] for l in f if l.strip()]
    col = {}
    for i, name in enumerate(header):
        col.setdefault(name, i)
    last = rows[-1]
    bend = [e for e in cfg["Energies"] if e["Name"] in ("Bending", "Bending_tan")][0]
    kappa = bend["constants"][0] / 2.0
    beads = [i for i, name in enumerate(header) if name == "Bead"]
    return {"run": os.path.basename(os.path.dirname(run.rstrip("/"))), "step": last[col["step"]], "E_bend_kappa": last[col[bend["Name"]]] / kappa,
            "E_total_kappa": last[col["Total_E"]] / kappa, "E_adh_kappa": sum(last[i] for i in beads) / kappa, "Volume": last[col["Volume"]], "Area": last[col["Area"]]}


# ----------------------------------------------------------------------------------------------------------------------
def run_main_cluster(binary, config_path, log_path, nice):
    cmd = [binary, config_path]
    if nice:
        cmd = ["nice", "-n", str(nice)] + cmd
    with open(log_path, "w") as log:
        proc = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT)
    return proc.returncode


def run_failed(log_path, run, min_steps=50):
    """The run ended badly: the integrator broke ("Sim broke or timestep very small"), finished after a handful of steps, or the energy blew up.
    main_cluster still returns 0 and writes a Final_state.obj then, so the log is read."""
    text = open(log_path, errors="replace").read()
    if "Sim broke" in text:
        return "the integrator broke (timestep very small)"
    m = re.search(r"Final state saved as step (\d+), total energy (\S+)", text)
    if not m:
        return "no 'Final state saved' line in the log"
    step, energy = int(m.group(1)), float(m.group(2))
    if step < min_steps:
        return "finished after %d steps" % step
    if not (abs(energy) < 1e6):
        return "total energy %g" % energy
    return None


def find_run_folder(first_dir):
    """The folder main_cluster made inside first_dir (the first free number: 1 in a fresh directory)."""
    nums = sorted(int(d) for d in os.listdir(first_dir) if d.isdigit() and os.path.isdir(os.path.join(first_dir, d)))
    if not nums:
        raise RuntimeError("main_cluster made no run folder in " + first_dir)
    return os.path.join(first_dir, str(nums[-1])) + "/"


def apply_scan(cfg, var, value, init_area, excess_fixed):
    """Put the value of the scanned parameter into a config (the same places the generator puts the cold value in). init_area: area of the initial mesh of the run,
    Excess_tension reads its target as a fraction of it."""
    Rv = cfg["Rv"]
    A_sphere = 4.0 * math.pi * Rv ** 2
    V_ref = 4.0 / 3.0 * math.pi * Rv ** 3
    if var in ("theta", "f"):
        theta = value if var == "theta" else math.degrees(math.acos(1.0 - 2.0 * value))
        for b in cfg["Beads"]:
            if "patch_angle" in b:
                b["patch_angle"] = theta
    elif var == "W":
        for b in cfg["Beads"]:
            b["inter_str"] = value
    elif var == "wt":
        # w~ = 4 W a^2/KB of the Deserno mapping
        KB = [e for e in cfg["Energies"] if e["Name"] in ("Bending", "Bending_tan")][0]["constants"][0]
        for b in cfg["Beads"]:
            b["inter_str"] = value * KB / (4.0 * b["radius"] ** 2)
    elif var == "cov":
        for e in cfg["Energies"]:
            if e["Name"] == "Coverage":
                e["constants"] = [e["constants"][0]] + [value] * len(cfg["Beads"])
    elif var in ("sigma", "st"):
        sigma = value
        if var == "st":
            # sigma~ = 2 sigma a^2/KB of the Deserno mapping
            KB = [e for e in cfg["Energies"] if e["Name"] in ("Bending", "Bending_tan")][0]["constants"][0]
            sigma = value * KB / (2.0 * cfg["Beads"][0]["radius"] ** 2)
        for e in cfg["Energies"]:
            if e["Name"] in ("Surface_tension", "Excess_tension"):
                e["constants"][0] = sigma
    elif var == "volume":
        for e in cfg["Energies"]:
            if e["Name"] == "Volume_constraint":
                e["constants"][1] = V_ref * value
    excess = value if var == "excess" else excess_fixed
    for e in cfg["Energies"]:
        if e["Name"] == "Area_constraint" and var == "excess":
            e["constants"][1] = A_sphere * (1.0 + excess)
        elif e["Name"] == "Excess_tension":
            # the loader multiplies the second constant by the area of the initial mesh of the run: keep the target absolute
            e["constants"][1] = A_sphere * (1.0 + excess) / init_area


# ----------------------------------------------------------------------------------------------------------------------
def main(chain_path):
    spec = json.load(open(chain_path))
    binary = os.path.abspath(spec.get("binary", "../build/bin/main_cluster"))
    batch = os.path.abspath(spec["batch_dir"])
    os.makedirs(batch, exist_ok=True)
    nice = spec.get("nice", 0)
    switch = spec.get("switch", [100, 1500])
    log = open(os.path.join(batch, "chain.log"), "a")

    def say(msg):
        stamp = time.strftime("%Y-%m-%d %H:%M:%S")
        print(msg, flush=True)
        log.write("%s %s\n" % (stamp, msg))
        log.flush()

    scan = spec["scan"]
    var, values, cold_value, max_step = scan["var"], [float(v) for v in scan["list"]], float(scan["cold"]), scan.get("max_step")
    link_stopping, link_timesteps = spec.get("link_stopping"), spec.get("link_timesteps")
    excess_fixed = float(spec.get("excess", 0.05))
    # a bead inside the vesicle ("inside": true in the chain file, "outside": -1 in the config): the patch axis points away from the volume centre
    side = -1.0 if spec.get("inside") else 1.0
    # a flat disk with a clamped edge ("planar": true): no volume centre, the patch axis stays pointing down; the bead sinks freely in the cold start
    # (state "default") and is frozen at that height from the first link on
    planar = bool(spec.get("planar"))
    cold = json.load(open(spec["config"]))
    n_beads = len(cold["Beads"])
    rows = []
    summary_path = os.path.join(batch, "chain_summary.csv")

    def write_summary():
        with open(summary_path, "w", newline="") as f:
            w = csv.writer(f)
            fields = ["link", "leg", "scan_var", "scan_value", "run", "step", "E_bend_kappa", "E_total_kappa", "E_adh_kappa", "Volume", "Area", "seconds"]
            w.writerow(fields)
            for r in rows:
                w.writerow([r.get(k, "") for k in fields])

    # --- link 0: the cold start (or an existing converged run, "start_run") ----------------------------------------
    t0 = time.time()
    if spec.get("start_run"):
        run = os.path.abspath(spec["start_run"]) + "/"
        cfg = json.load(open(os.path.join(run, "Input_file.json")))
        say("link 00: starting from the finished run %s" % run)
    else:
        cold_dir = os.path.join(batch, "link00_cold") + "/"
        cfg = copy.deepcopy(cold)
        cfg["first_dir"] = cold_dir
        if spec.get("prewrap_config"):
            # a quick wrap with the stiff Coverage energy first; the patch run then starts from that mesh, the bead where it was
            pre_dir = os.path.join(batch, "link00a_prewrap") + "/"
            pre = json.load(open(spec["prewrap_config"]))
            pre["first_dir"] = pre_dir
            os.makedirs(pre_dir, exist_ok=True)
            pre_path = os.path.join(batch, "link00a_prewrap.json")
            json.dump(pre, open(pre_path, "w"), indent=4)
            say("link 00a (pre-wrap with Coverage, %d steps) -> %s" % (pre["timesteps"], pre_dir))
            rc = run_main_cluster(binary, pre_path, os.path.join(batch, "link00a_prewrap.log"), nice)
            pre_run = find_run_folder(pre_dir)
            why = run_failed(os.path.join(batch, "link00a_prewrap.log"), pre_run, min_steps=100)
            if rc != 0 or why or not os.path.exists(os.path.join(pre_run, "Final_state.obj")):
                say("link 00a failed (return code %d, %s), see link00a_prewrap.log" % (rc, why))
                return 1
            cfg["init_file"] = os.path.abspath(last_mesh(pre_run))
            cfg["rescale"] = 1.0
            cfg.pop("Displacement", None)
            cfg["Switch_times"] = switch
            # the vesicle has moved up to the bead: the patch axis points to its centre now
            Vp, Fp = read_obj(last_mesh(pre_run))
            vcp, _ = volume_centre(Vp, Fp)
            for b in cfg["Beads"]:
                if "patch_axis" in b:
                    ax = side * (vcp - np.array(b["Pos"]))
                    b["patch_axis"] = [float(x) for x in ax / np.linalg.norm(ax)]
        os.makedirs(cold_dir, exist_ok=True)
        cfg_path = os.path.join(batch, "link00_cold.json")
        json.dump(cfg, open(cfg_path, "w"), indent=4)
        say("link 00 (cold start at %s = %g) -> %s" % (var, cold_value, cold_dir))
        rc = run_main_cluster(binary, cfg_path, os.path.join(batch, "link00_cold.log"), nice)
        run = find_run_folder(cold_dir)
        why = run_failed(os.path.join(batch, "link00_cold.log"), run, min_steps=100)
        if rc != 0 or why or not os.path.exists(os.path.join(run, "Final_state.obj")):
            say("link 00 failed (return code %d, %s), see link00_cold.log" % (rc, why))
            return 1
    r = summary(run, cfg)
    r.update({"link": 0, "leg": "cold", "scan_var": var, "scan_value": cold_value, "seconds": round(time.time() - t0)})
    rows.append(r)
    write_summary()
    say("link 00 done: E_total = %.4f kappa, E_bend = %.4f kappa, E_adhesion = %.4f kappa after %d steps, %d s" % (r["E_total_kappa"], r["E_bend_kappa"], r["E_adh_kappa"], r["step"], r["seconds"]))
    prev_run, prev_cfg = run, cfg

    # --- the schedule: the list, and back along it (ends at the cold value) -------------------------------------------
    sched = [(v, "in") for v in values]
    if spec.get("return", False):
        back = [(v, "out") for v in reversed(values[:-1])]
        if abs(cold_value - values[0]) > 1e-12:
            back.append((cold_value, "out"))  # the cold start is not on the list: the way back ends there
        sched += back

    def attempt(k, tag, value, leg, prev_run, prev_cfg):
        """One link: put the value into the config, start from the previous membrane, relax. Returns (run, cfg) or None if main_cluster failed."""
        p = last_bead_positions(prev_run, n_beads)
        V, F = read_obj(last_mesh(prev_run))
        vc = None if planar else volume_centre(V, F)[0]
        fdir = os.path.join(batch, tag) + "/"
        os.makedirs(fdir, exist_ok=True)
        cfg = copy.deepcopy(prev_cfg)
        cfg["init_file"] = os.path.abspath(last_mesh(prev_run))
        cfg["rescale"] = 1.0
        cfg.pop("Displacement", None)
        cfg["first_dir"] = fdir
        cfg["Switch_times"] = switch
        if link_stopping:
            cfg["stopping"].update(link_stopping)
        if link_timesteps:
            cfg["timesteps"] = int(link_timesteps)
        apply_scan(cfg, var, value, mesh_area(V, F), excess_fixed)
        for i, b in enumerate(cfg["Beads"]):
            b["Pos"] = [float(x) for x in p[i]]
            b["state"] = "froze"  # the bead is frozen from the first link on (a free cold start, --free_cold and --planar, lets it settle first)
            if planar:
                pass  # the patch axis stays pointing down
            elif "patch_axis" in b:
                ax = side * (vc - p[i])
                b["patch_axis"] = [float(x) for x in ax / np.linalg.norm(ax)]
        cfg_path = os.path.join(batch, tag + ".json")
        json.dump(cfg, open(cfg_path, "w"), indent=4)
        say("link %02d: %s = %g (%s leg)" % (k, var, value, leg))
        t0 = time.time()
        rc = run_main_cluster(binary, cfg_path, os.path.join(batch, tag + ".log"), nice)
        run = find_run_folder(fdir)
        why = run_failed(os.path.join(batch, tag + ".log"), run)
        if rc != 0 or why or not os.path.exists(os.path.join(run, "Final_state.obj")):
            say("link %02d failed (return code %d, %s), see %s.log" % (k, rc, why, tag))
            return None
        r = summary(run, cfg)
        r.update({"link": k, "leg": leg, "scan_var": var, "scan_value": value, "seconds": round(time.time() - t0)})
        rows.append(r)
        write_summary()
        say("link %02d done: %s = %g, E_total = %.4f kappa, E_bend = %.4f kappa, E_adhesion = %.4f kappa, %d steps, %d s" % (
            k, var, value, r["E_total_kappa"], r["E_bend_kappa"], r["E_adh_kappa"], r["step"], r["seconds"]))
        return run, cfg

    cur = cold_value
    for k, (value, leg) in enumerate(sched, start=1):
        tag = "link%02d_%s%g" % (k, var, value)
        # a patch edge that advances by more than about a tenth of a bead radius per relaxation leaves the membrane behind: go in steps of at most max_step
        n_steps = max(1, int(math.ceil(abs(value - cur) / max_step - 1e-9))) if max_step else 1
        v_from = cur
        for i in range(1, n_steps):
            v_way = v_from + (value - v_from) * i / n_steps
            say("link %02d: waypoint %d of %d at %s = %g" % (k, i, n_steps, var, v_way))
            res = attempt(k, "%s_w%d" % (tag, i), v_way, leg, prev_run, prev_cfg)
            if res is None:
                say("giving up at the waypoint of link %02d" % k)
                return 1
            prev_run, prev_cfg = res
            cur = v_way
        res = attempt(k, tag, value, leg, prev_run, prev_cfg)
        if res is None:
            # one retry in two half steps (the parameter moved too far for the membrane): through the middle value first
            v_mid = 0.5 * (cur + value)
            say("link %02d: retrying in two steps, through %s = %g" % (k, var, v_mid))
            res = attempt(k, tag + "_mid", v_mid, leg, prev_run, prev_cfg)
            if res is not None:
                res = attempt(k, tag + "_retry", value, leg, res[0], res[1])
        if res is None:
            say("giving up at link %02d" % k)
            return 1
        prev_run, prev_cfg = res
        cur = value
    say("chain finished, %d runs" % len(rows))
    return 0


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise SystemExit(__doc__)
    sys.exit(main(sys.argv[1]))
