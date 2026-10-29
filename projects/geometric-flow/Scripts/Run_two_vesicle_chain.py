#!/usr/bin/env python3
"""Two wrapped beads on a vesicle, scanned in the separation by moving the beads along the vesicle in small steps.

usage:  python3 Run_two_vesicle_chain.py CHAIN.json          (CHAIN.json is written by Create_subjob_two_vesicle.py)

1. link 0 (cold start): main_cluster on the config of CHAIN.json (the tangent start of Templates/Wrapping_two_vesicle.txt at the first
   separation); the beads wrap, the run goes on to convergence (the "stopping" block of the config). With "start_run": "<run folder>"
   in CHAIN.json a finished run is the link 0 instead (to reuse a cold start, to extend a chain, or to try another list of separations).
2. the next links: the two beads are rotated about the volume centre of the vesicle by +-(phi_new - phi_old)/2, so that the distance s between
   them takes the next value of "s_list" (in D_p). "distance" in CHAIN.json says which distance: "arc" (the old default) = R_v phi, the arc of the
   sphere; "chord" = the straight distance between the bead centres d/D_p, d^2 = R_0^2 + R_1^2 - 2 R_0 R_1 cos(phi) with R_i the distance of bead i
   from the volume centre. With "chord" the beads reach contact at s = 1 and every value of "s_list" below 1 is refused (the beads never overlap).
   The membrane of the previous link
   is the initial mesh (Final_state.obj), the beads stay frozen, the patch axes (if any) point to the vesicle centre again, the
   absolute volume (and area) targets are the same. A link starts close to its minimum and takes minutes, not the hours of a cold
   start, and the meshes of neighbouring separations are related, which removes most of the noise a fresh mesh brings into the
   energy (re-meshing and re-relaxing a converged state changes E_bend by 0.01-0.04 kappa).
   A step larger than "max_step" (default 0.5 D_p) is split into waypoints (the cold start of a floppy vesicle ends 5 D_p apart where its beads sit deep in
   their pockets, and a jump to 3 D_p unwrapped both beads).
3. "return": the list is run backwards too (a hysteresis check: the two branches must agree within the noise).

Everything goes to batch_dir: link00_cold/, link01_s2.50/, ... (each is a normal run folder: Input_file.json, Output_data.txt, Final_state.obj,
Bead_*_data.txt, ...), chain_summary.csv (one line per link) and chain.log. Analyze_two_vesicle.py <batch_dir> turns the folders into u(s).
Needs only json and numpy, so it runs inside the sbatch job.
"""
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


def rotate(v, axis, angle):
    """Rodrigues rotation of v about the unit vector axis."""
    c, s = math.cos(angle), math.sin(angle)
    return v * c + np.cross(axis, v) * s + axis * np.dot(axis, v) * (1.0 - c)


def arc_s(run, Rv, a, n=2, planar=False):
    """s = R_v phi / D_p of the last state of a run: phi is the angle between the bead centres seen from the volume centre.
    Planar disk (no volume centre): s is the distance between the bead centres in D_p."""
    p = last_bead_positions(run, n)
    if planar:
        return float(np.linalg.norm(p[0] - p[1])) / (2.0 * a)
    V, F = read_obj(last_mesh(run))
    vc, _ = volume_centre(V, F)
    u = [(q - vc) / np.linalg.norm(q - vc) for q in p]
    return Rv * math.acos(max(-1.0, min(1.0, float(u[0] @ u[1])))) / (2.0 * a)


def summary(run, cfg, s_target, chord=False):
    with open(os.path.join(run, "Output_data.txt")) as f:
        header = f.readline().split()
        rows = [[float(x) for x in l.split()] for l in f if l.strip()]
    col = {}
    for i, name in enumerate(header):
        col.setdefault(name, i)
    last = rows[-1]
    bend = [e for e in cfg["Energies"] if e["Name"] in ("Bending", "Bending_tan")][0]  # Bending since 2026-10-12, Bending_tan before
    kappa = bend["constants"][0] / 2.0
    out = {"s_target": s_target, "s_actual": arc_s(run, cfg["Rv"], cfg["Beads"][0]["radius"], planar=bool(cfg.get("boundary")) or chord), "run": os.path.basename(os.path.dirname(run.rstrip("/"))), "step": last[col["step"]],
           "E_bend_kappa": last[col[bend["Name"]]] / kappa, "E_total_kappa": last[col["Total_E"]] / kappa,
           "Volume": last[col["Volume"]], "Area": last[col["Area"]]}
    return out


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
    main_cluster still returns 0 and writes a Final_state.obj then, so the log is read. A "Sim broke" that follows a vanishing gradient is a converged run."""
    text = open(log_path, errors="replace").read()
    # A line search that stops on a vanishing gradient ("The energy diff is quite small and so is the gradient", then the projection) returns -1 like a collapsed one and
    # main_cluster prints "Sim broke" for both: the run is converged, not broken, and short because it started close to its minimum
    converged = re.search(r"quite small and so is the gradient\s+The energy diff is\S+\s+The projection is(\S+)\s+Sim broke", text)
    if converged and float(converged.group(1)) < 1e-4:
        min_steps = 0
    elif "Sim broke" in text:
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


def main(chain_path):
    spec = json.load(open(chain_path))
    here = os.getcwd()  # the Scripts folder: the paths of the configs are relative to it, like the init_file
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

    cold = json.load(open(spec["config"]))
    # beads inside the vesicle ("inside": true in the chain file, "outside": -1 in the config): the patch axis points away from the volume centre
    side = -1.0 if spec.get("inside") else 1.0
    # a flat disk with a clamped edge ("planar": true): no volume centre, the beads are translated along the line joining them; they sink freely in
    # the cold start (state "default") and are frozen at that height from the first link on
    planar = bool(spec.get("planar"))
    dist_mode = "chord" if (planar or spec.get("distance", "arc") == "chord") else "arc"  # the disk has no volume centre: its s is always the straight distance
    chord = dist_mode == "chord"
    a = cold["Beads"][0]["radius"]
    Rv = cold["Rv"]
    n_beads = len(cold["Beads"])
    if n_beads != 2:
        raise SystemExit("two beads expected")
    rows = []
    summary_path = os.path.join(batch, "chain_summary.csv")

    def write_summary():
        with open(summary_path, "w", newline="") as f:
            w = csv.writer(f)
            fields = ["link", "s_target", "s_actual", "run", "step", "E_bend_kappa", "E_total_kappa", "Volume", "Area", "seconds", "distance"]
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
        cfg = json.loads(json.dumps(cold))
        cfg["first_dir"] = cold_dir
        if spec.get("prewrap_config"):
            # a quick wrap with the stiff Coverage energy first; the patch run then starts from that mesh, beads where they were
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
            # the vesicle has moved up to the beads: the patch axes point to its centre now
            Vp, Fp = read_obj(last_mesh(pre_run))
            vcp, _ = volume_centre(Vp, Fp)
            for b in cfg["Beads"]:
                if "patch_axis" in b:
                    ax = side * (vcp - np.array(b["Pos"]))
                    b["patch_axis"] = [float(x) for x in ax / np.linalg.norm(ax)]
        os.makedirs(cold_dir, exist_ok=True)
        cfg_path = os.path.join(batch, "link00_cold.json")
        json.dump(cfg, open(cfg_path, "w"), indent=4)
        say("link 00 (cold start) -> %s" % cold_dir)
        rc = run_main_cluster(binary, cfg_path, os.path.join(batch, "link00_cold.log"), nice)
        run = find_run_folder(cold_dir)
        why = run_failed(os.path.join(batch, "link00_cold.log"), run, min_steps=100)
        if rc != 0 or why or not os.path.exists(os.path.join(run, "Final_state.obj")):
            say("link 00 failed (return code %d, %s), see link00_cold.log" % (rc, why))
            return 1
    r = summary(run, cfg, None, chord)
    r.update({"link": 0, "seconds": round(time.time() - t0), "distance": dist_mode})
    rows.append(r)
    write_summary()
    say("link 00 done: s = %.3f D_p, E_bend = %.4f kappa after %d steps, %d s" % (r["s_actual"], r["E_bend_kappa"], r["step"], r["seconds"]))
    s_cold = r["s_actual"]
    prev_run, prev_cfg = run, cfg

    # --- the other links ------------------------------------------------------------------------------------------
    s_list = list(spec["s_list"])
    if chord and min(s_list) < 1.0 - 1e-9:
        raise SystemExit("distance chord: the beads would overlap (s = %g D_p < 1), every s of s_list must be >= 1" % min(s_list))
    if spec.get("return", False):
        # back out along the same values (ends at the separation of the cold start)
        s_list = s_list + list(reversed(s_list[:-1])) + [s_cold]
    def attempt(k, tag, s, prev_run, prev_cfg):
        """One link: move the beads to the arc s, relax. Returns (run, cfg) or None if main_cluster failed."""
        p = last_bead_positions(prev_run, n_beads)
        if planar:
            # the beads move along the line joining them, in the plane of the disk (x is the normal), about their midpoint
            n01 = p[1] - p[0]
            n01[0] = 0.0
            n01 /= np.linalg.norm(n01)
            mid = 0.5 * (p[0] + p[1])
            d_new = s * 2.0 * a
            phi_old = phi_new = float("nan")
            pos = [mid - 0.5 * d_new * n01, mid + 0.5 * d_new * n01]
            pos[0][0] = pos[1][0] = mid[0]
        else:
            V, F = read_obj(last_mesh(prev_run))
            vc, vol = volume_centre(V, F)
            u = [(q - vc) / np.linalg.norm(q - vc) for q in p]
            phi_old = math.acos(max(-1.0, min(1.0, float(u[0] @ u[1]))))
            if chord:
                # the straight distance d = 2 a s between the centres: d^2 = R0^2 + R1^2 - 2 R0 R1 cos(phi), R_i the distance of bead i from the volume centre
                R0, R1 = (float(np.linalg.norm(q - vc)) for q in p)
                cos_phi = (R0 * R0 + R1 * R1 - (s * 2.0 * a) ** 2) / (2.0 * R0 * R1)
                if abs(cos_phi) > 1.0:
                    raise RuntimeError("no angle puts the beads %.3f apart (R_0 = %.3f, R_1 = %.3f)" % (s * 2.0 * a, R0, R1))
                phi_new = math.acos(cos_phi)
            else:
                phi_new = s * 2.0 * a / Rv
            axis = np.cross(u[0], u[1])
            axis /= np.linalg.norm(axis)
            dphi = 0.5 * (phi_new - phi_old)
            # bead 0 turns away from (towards) bead 1: the rotation axis u0 x u1 turns u0 towards u1, so rotate u0 by -dphi and u1 by +dphi
            pos = [vc + rotate(p[0] - vc, axis, -dphi), vc + rotate(p[1] - vc, axis, +dphi)]
        fdir = os.path.join(batch, tag) + "/"
        os.makedirs(fdir, exist_ok=True)
        cfg = json.loads(json.dumps(prev_cfg))
        cfg["init_file"] = os.path.abspath(last_mesh(prev_run))
        cfg["rescale"] = 1.0
        cfg.pop("Displacement", None)
        cfg["first_dir"] = fdir
        cfg["theta"] = float(np.linalg.norm(pos[0] - pos[1]))
        if cfg["theta"] < 2.0 * a * (1.0 - 1e-9):
            raise RuntimeError("the beads would overlap: distance %.4f < 2 a = %.4f" % (cfg["theta"], 2.0 * a))
        cfg["Switch_times"] = switch
        for i, b in enumerate(cfg["Beads"]):
            b["Pos"] = [float(x) for x in pos[i]]
            b["state"] = "froze"  # the beads are frozen from the first link on (a free cold start, "--free_cold" and --planar, lets them settle first)
            if planar:
                pass  # the patch axes stay pointing down
            elif "patch_axis" in b:
                ax = side * (vc - pos[i])
                b["patch_axis"] = [float(x) for x in ax / np.linalg.norm(ax)]
        cfg_path = os.path.join(batch, tag + ".json")
        json.dump(cfg, open(cfg_path, "w"), indent=4)
        say("link %02d: s = %.3f D_p (%s; phi %.4f -> %.4f rad, bead distance %.3f)" % (k, s, dist_mode, phi_old, phi_new, cfg["theta"]) if not planar else
            "link %02d: s = %.3f D_p (bead distance %.3f, height %.3f)" % (k, s, cfg["theta"], pos[0][0]))
        t0 = time.time()
        rc = run_main_cluster(binary, cfg_path, os.path.join(batch, tag + ".log"), nice)
        run = find_run_folder(fdir)
        why = run_failed(os.path.join(batch, tag + ".log"), run)
        if rc != 0 or why or not os.path.exists(os.path.join(run, "Final_state.obj")):
            say("link %02d failed (return code %d, %s), see %s.log" % (k, rc, why, tag))
            return None
        r = summary(run, cfg, s, chord)
        r.update({"link": k, "seconds": round(time.time() - t0), "distance": dist_mode})
        rows.append(r)
        write_summary()
        say("link %02d done: s = %.3f D_p, E_bend = %.4f kappa, %d steps, %d s" % (k, r["s_actual"], r["E_bend_kappa"], r["step"], r["seconds"]))
        return run, cfg

    max_step = float(spec.get("max_step", 0.5))
    for k, s in enumerate(s_list, start=1):
        tag = "link%02d_s%.2f" % (k, s)
        # a bead moved by more than about its own radius leaves its pocket (a step of 2.2 D_p unwrapped both beads in a test): go in steps of at most max_step
        s_from = rows[-1]["s_actual"]
        n_steps = max(1, int(math.ceil(abs(s - s_from) / max_step - 1e-9)))
        for i in range(1, n_steps):
            s_way = s_from + (s - s_from) * i / n_steps
            say("link %02d: waypoint %d of %d at s = %.3f D_p" % (k, i, n_steps, s_way))
            res = attempt(k, "%s_w%d" % (tag, i), s_way, prev_run, prev_cfg)
            if res is None:
                say("giving up at the waypoint of link %02d" % k)
                return 1
            prev_run, prev_cfg = res
        res = attempt(k, tag, s, prev_run, prev_cfg)
        if res is None:
            # one retry in two half steps (the beads moved too far for the mesh): through the middle value first
            s_mid = 0.5 * (rows[-1]["s_actual"] + s)
            say("link %02d: retrying in two steps, through s = %.3f D_p" % (k, s_mid))
            res = attempt(k, tag + "_mid", s_mid, prev_run, prev_cfg)
            if res is not None:
                res = attempt(k, tag + "_retry", s, res[0], res[1])
        if res is None:
            say("giving up at link %02d" % k)
            return 1
        prev_run, prev_cfg = res
    say("chain finished, %d links" % len(rows))
    return 0


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise SystemExit(__doc__)
    sys.exit(main(sys.argv[1]))
