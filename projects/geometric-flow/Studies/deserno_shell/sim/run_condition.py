#!/usr/bin/env python
"""Run one 'condition' = a batch of planar-bead simulations on this machine, and summarise it.

    micromamba run -n mir_mem python run_condition.py launch  --tag T --KA 0.2667 --KIs 1.4,1.6 --inits flat,full \
          --shell_width 0.1 --size_min 0.02 --size_max 0.3 --jobs 3 [--extra "--shell_power 2"] [--set remesher.refine_angle=0.3]
    micromamba run -n mir_mem python run_condition.py status  --tag T
    micromamba run -n mir_mem python run_condition.py collect --tag T [--out data/T.csv]
    micromamba run -n mir_mem python run_condition.py plot    --tag T[,T2] --out ../figures/conv_T.png

What it does (see protocol.md): every run gets its own config (Scripts/Create_subjob_planar.py, always --rescale 1.0 and --tension
excess, then the json is post-edited: Switch_times[1], timesteps, any --set key=value), its own batch folder
Results/<tag>_<init>_KI<KI>/ (so the numbered-folder creation never races), and is started through slot_run.sh (global CPU
semaphore) by a detached worker that keeps at most --jobs of them waiting/running. `launch` returns at once.
Control files (run.sh, stdout.log, start/end stamps, manifest.json, summary.csv) live in Results/_sim_ctl/<tag>/ (gitignored).
`collect` deletes Coverage_data.txt of each batch folder first (PostProcessing APPENDS), runs the binary set's PostProcessing,
and writes summary.csv with one row per run.
"""
import argparse, csv, json, math, os, shlex, subprocess, sys, time

HERE = os.path.dirname(os.path.abspath(__file__))
STUDY = os.path.dirname(HERE)
P = os.path.dirname(os.path.dirname(STUDY))          # projects/geometric-flow
SCRIPTS = os.path.join(P, "Scripts")
RESULTS = os.path.join(P, "Results")
CTL_ROOT = os.path.join(RESULTS, "_sim_ctl")
SLOT_RUN = os.path.join(STUDY, "slot_run.sh")
DEFAULT_BINARY = os.path.join(P, "build", "shellpower", "main_cluster")


# ----------------------------------------------------------------------------- helpers
def ctl_dir(tag):
    return os.path.join(CTL_ROOT, tag)


def fnum(x):
    return ("%g" % float(x))


def run_id(init, KI):
    return "%s_KI%s" % (init, fnum(KI))


def parse_KIs(s):
    """'1.4,1.6' or 'a:b:step' (inclusive range)."""
    if ":" in s:
        a, b, st = (float(x) for x in s.split(":"))
        n = int(round((b - a) / st))
        return [round(a + i * st, 10) for i in range(n + 1)]
    return [float(x) for x in s.split(",") if x]


def set_path(d, path, value):
    """d['a'][1]['b'] = value for path 'a.1.b' (integer keys index lists)."""
    keys = path.split(".")
    cur = d
    for k in keys[:-1]:
        cur = cur[int(k)] if isinstance(cur, list) else cur[k]
    k = keys[-1]
    if isinstance(cur, list):
        cur[int(k)] = value
    else:
        cur[k] = value


def read_json(p):
    with open(p) as f:
        return json.load(f)


def write_json(p, d):
    with open(p, "w") as f:
        json.dump(d, f, indent=4)


# ----------------------------------------------------------------------------- launch
def make_config(args, KI, init):
    rid = run_id(init, KI)
    uniq = "%s_%s" % (args.tag, rid)
    batch = uniq
    cmd = [sys.executable, "Create_subjob_planar.py", "--inter_str", fnum(KI), "--radius", fnum(args.radius),
           "--KA", fnum(args.KA), "--KB", fnum(args.KB), "--init", init, "--tension", "excess", "--rescale", "1.0",
           "--batch_tag", batch, "--unique_tag", uniq]
    for name in ("shell_width", "shell_power", "size_min", "size_max"):
        v = getattr(args, name)
        if v is not None:
            cmd += ["--" + name, fnum(v)]
    if args.extra:
        cmd += shlex.split(args.extra)
    subprocess.run(cmd, cwd=SCRIPTS, check=True, stdout=subprocess.DEVNULL)
    cfg_path = os.path.join(P, "Config_files", uniq + "_ConfigFile.json")
    cfg = read_json(cfg_path)
    if args.switch_normal is not None:
        cfg["Switch_times"][1] = int(args.switch_normal)
    if args.timesteps is not None:
        cfg["timesteps"] = int(args.timesteps)
    for kv in args.set or []:
        k, v = kv.split("=", 1)
        try:
            v = json.loads(v)
        except ValueError:
            pass
        set_path(cfg, k, v)
    write_json(cfg_path, cfg)
    return dict(id=rid, init=init, KI=float(KI), batch=batch, config=cfg_path, unique=uniq)


def cmd_launch(args):
    d = ctl_dir(args.tag)
    mp = os.path.join(d, "manifest.json")
    if os.path.exists(mp) and not args.dry_run and not read_json(mp).get("dry_run"):
        sys.exit("tag %s already exists (%s): use a new tag or delete the folder" % (args.tag, d))
    binary = os.path.abspath(args.binary)
    if not os.access(binary, os.X_OK):
        sys.exit("binary not executable: %s" % binary)
    if args.snapshot_binary:   # copy so that a rebuild cannot change a running campaign
        os.makedirs(d, exist_ok=True)
        snap = os.path.join(d, "bin")
        os.makedirs(snap, exist_ok=True)
        for b in ("main_cluster", "PostProcessing"):
            src = os.path.join(os.path.dirname(binary), b)
            if os.path.exists(src):
                subprocess.run(["cp", "-p", src, os.path.join(snap, b)], check=True)
        binary = os.path.join(snap, "main_cluster")
    KIs = parse_KIs(args.KIs)
    inits = args.inits.split(",")
    runs = [make_config(args, KI, init) for KI in KIs for init in inits]
    os.makedirs(d, exist_ok=True)
    for r in runs:
        rd = os.path.join(d, r["id"])
        os.makedirs(rd, exist_ok=True)
        r["ctl"] = rd
        cfg_rel = os.path.relpath(r["config"], SCRIPTS)
        tmo = ("timeout -s TERM %d " % args.timeout) if args.timeout else ""
        with open(os.path.join(rd, "run.sh"), "w") as f:
            f.write("#!/usr/bin/env bash\n# runs inside a slot_run.sh slot: stamps the real start/end (slot waiting excluded)\n"
                    "cd %s\ndate +%%s.%%N > %s/start\n%s%s %s > %s/stdout.log 2>&1\necho $? > %s/exit_code\ndate +%%s.%%N > %s/end\n"
                    % (shlex.quote(SCRIPTS), shlex.quote(rd), tmo, shlex.quote(binary), shlex.quote(cfg_rel),
                       shlex.quote(rd), shlex.quote(rd), shlex.quote(rd)))
        os.chmod(os.path.join(rd, "run.sh"), 0o755)
    manifest = dict(tag=args.tag, binary=binary, KA=args.KA, KB=args.KB, radius=args.radius,
                    shell_width=args.shell_width, shell_power=args.shell_power, size_min=args.size_min,
                    size_max=args.size_max, extra=args.extra, set=args.set, switch_normal=args.switch_normal,
                    timesteps=args.timesteps, jobs=args.jobs, dry_run=args.dry_run, launched=time.time(), runs=runs)
    write_json(os.path.join(d, "manifest.json"), manifest)
    print("%d runs, configs in %s" % (len(runs), os.path.join(P, "Config_files")))
    for r in runs:
        print("  %-22s %s" % (r["id"], os.path.relpath(r["config"], P)))
    if args.dry_run:
        print("dry run: nothing started (ctl folder %s written)" % d)
        return
    log = open(os.path.join(d, "worker.log"), "w")
    subprocess.Popen([sys.executable, os.path.abspath(__file__), "_worker", "--tag", args.tag, "--jobs", str(args.jobs)],
                     stdout=log, stderr=log, stdin=subprocess.DEVNULL, start_new_session=True, cwd=SCRIPTS)
    print("worker started (detached, <= %d concurrent submissions through slot_run.sh); poll with: status --tag %s"
          % (args.jobs, args.tag))


def cmd_worker(args):
    d = ctl_dir(args.tag)
    m = read_json(os.path.join(d, "manifest.json"))
    pending = list(m["runs"])
    active = []
    while pending or active:
        active = [(r, p) for r, p in active if p.poll() is None]
        while pending and len(active) < args.jobs:
            r = pending.pop(0)
            p = subprocess.Popen([SLOT_RUN, os.path.join(r["ctl"], "run.sh")], stdin=subprocess.DEVNULL,
                                 stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, cwd=SCRIPTS)
            active.append((r, p))
            print(time.strftime("%H:%M:%S"), "submitted", r["id"], flush=True)
        time.sleep(5)
    print(time.strftime("%H:%M:%S"), "all done", flush=True)


# ----------------------------------------------------------------------------- inspection of one run
def result_dir(r):
    base = os.path.join(RESULTS, r["batch"])
    if not os.path.isdir(base):
        return None
    nums = [int(x) for x in os.listdir(base) if x.isdigit() and os.path.isdir(os.path.join(base, x))]
    return os.path.join(base, str(max(nums))) + "/" if nums else None


def stamp(rd, name):
    p = os.path.join(rd, name)
    if os.path.exists(p):
        try:
            return float(open(p).read().split()[0])
        except (ValueError, IndexError):
            return None
    return None


def conv_rows(res):
    """Convergence_log.txt -> list of (step, integ, E, dE_rel, passed)."""
    out = []
    p = os.path.join(res, "Convergence_log.txt")
    if not os.path.exists(p):
        return out
    for line in open(p):
        t = line.split()
        if line.startswith("#") or len(t) < 7:
            continue
        try:
            out.append((int(t[0]), t[1], float(t[2]), float(t[3]), int(t[6])))
        except ValueError:
            pass
    return out


def output_rows(res):
    p = os.path.join(res, "Output_data.txt")
    if not os.path.exists(p):
        return None, []
    lines = open(p).read().splitlines()
    hdr = lines[0].split()
    rows = []
    for l in lines[1:]:
        t = l.split()
        try:
            rows.append([float(x) for x in t])
        except ValueError:
            pass
    return hdr, rows


def state_of(rd, res):
    if os.path.exists(os.path.join(rd, "end")):
        return "done" if open(os.path.join(rd, "exit_code")).read().strip() == "0" else "failed"
    if os.path.exists(os.path.join(rd, "start")):
        return "running"
    return "queued"


def cmd_status(args):
    m = read_json(os.path.join(ctl_dir(args.tag), "manifest.json"))
    print("%-18s %-8s %8s %-11s %10s %4s %8s" % ("run", "state", "step", "integrator", "E", "pass", "elapsed_s"))
    for r in m["runs"]:
        res = result_dir(r)
        st = state_of(r["ctl"], res)
        step, integ, E, ps = "-", "-", "-", "-"
        if res:
            c = conv_rows(res)
            if c:
                step, integ, E, ps = c[-1][0], c[-1][1], "%.6f" % c[-1][2], c[-1][4]
        t0, t1 = stamp(r["ctl"], "start"), stamp(r["ctl"], "end")
        el = "-" if t0 is None else "%.0f" % ((t1 or time.time()) - t0)
        print("%-18s %-8s %8s %-11s %10s %4s %8s" % (r["id"], st, step, integ, E, ps, el))


def cmd_wait(args):
    """Block until every run of the tag has finished (or --max_s seconds passed)."""
    m = read_json(os.path.join(ctl_dir(args.tag), "manifest.json"))
    t0 = time.time()
    while time.time() - t0 < args.max_s:
        if all(os.path.exists(os.path.join(r["ctl"], "end")) for r in m["runs"]):
            print("all runs of %s finished" % args.tag)
            return
        time.sleep(5)
    print("timeout, not all finished")


# ----------------------------------------------------------------------------- collect
def read_obj(path):
    V, F = [], []
    for line in open(path):
        if line.startswith("v "):
            V.append([float(x) for x in line.split()[1:4]])
        elif line.startswith("f "):
            F.append([int(x.split("/")[0]) - 1 for x in line.split()[1:4]])
    return V, F


def mesh_stats(res, step, bead, a, band=0.3):
    import numpy as np
    p = os.path.join(res, "membrane_%d.obj" % step)
    if not os.path.exists(p):
        p = os.path.join(res, "Final_state.obj")
    V, F = read_obj(p)
    V, F = np.array(V), np.array(F)
    cen = V[F].mean(axis=1)
    d = np.linalg.norm(cen - np.array(bead), axis=1)
    near = np.abs(d - a) < band
    out = dict(faces=len(F), faces_near=int(near.sum()))
    if near.any():
        Fn = F[near]
        e = np.concatenate([np.linalg.norm(V[Fn[:, i]] - V[Fn[:, (i + 1) % 3]], axis=1) for i in range(3)])
        out.update(h_min=e.min(), h_mean=e.mean(), h_max=e.max())
    else:
        out.update(h_min=float("nan"), h_mean=float("nan"), h_max=float("nan"))
    return out


COLS = ["run", "KA", "KI", "w_t", "sigma_t", "init", "z", "CoverageUnion", "Total_E", "E_bead", "E_bend", "E_tension",
        "steps", "wall_s", "state", "stop", "converged", "dE_last_window", "dE_last_1000", "integrator", "passed",
        "faces", "faces_near", "h_min", "h_mean", "h_max", "BeadX", "result_dir"]


def cmd_collect(args):
    d = ctl_dir(args.tag)
    m = read_json(os.path.join(d, "manifest.json"))
    pp = os.path.join(os.path.dirname(m["binary"]), "PostProcessing")
    if not os.path.exists(pp):
        pp = os.path.join(os.path.dirname(DEFAULT_BINARY), "PostProcessing")
    rows = []
    for r in m["runs"]:
        res = result_dir(r)
        row = dict.fromkeys(COLS, "")
        row.update(run=r["id"], KA=m["KA"], KI=r["KI"], init=r["init"], w_t=4 * r["KI"] * m["radius"] ** 2 / m["KB"],
                   sigma_t=2 * m["KA"] * m["radius"] ** 2 / m["KB"])
        st = state_of(r["ctl"], res)
        row["state"] = st
        rows.append(row)
        if not res:
            continue
        row["result_dir"] = os.path.relpath(res, P)
        t0, t1 = stamp(r["ctl"], "start"), stamp(r["ctl"], "end")
        if t0 and t1:
            row["wall_s"] = round(t1 - t0, 1)
        hdr, orows = output_rows(res)
        if orows:
            last = orows[-1]
            row["steps"] = int(last[1])
            for name, col in (("Total_E", "Total_E"), ("E_bead", "Bead"), ("E_bend", "Bending_tan"),
                              ("E_tension", "Excess_tension")):
                if col in hdr:
                    row[name] = last[hdr.index(col)]
        c = conv_rows(res)
        patience = read_json(os.path.join(res, "Input_file.json")).get("stopping", {}).get("patience", 4)
        if c:
            row["integrator"], row["passed"] = c[-1][1], c[-1][4]
            row["dE_last_window"] = c[-1][3] * max(abs(c[-1][2]), 1.0)
            ref = [x for x in c if x[0] <= c[-1][0] - 1000]
            if ref:
                row["dE_last_1000"] = abs(c[-1][2] - ref[-1][2])
            row["converged"] = int(c[-1][1] == "BFGS-Normal" and c[-1][4] >= patience)
        log = os.path.join(r["ctl"], "stdout.log")
        txt = open(log, errors="replace").read() if os.path.exists(log) else ""
        row["stop"] = ("converged" if "BFGS-Normal converged" in txt else "stalled" if "ending the run" in txt
                       else "maxsteps" if st == "done" else st)
        if st != "done":
            continue
        # PostProcessing appends to Coverage_data.txt: delete first
        cov = os.path.join(RESULTS, r["batch"], "Coverage_data.txt")
        if os.path.exists(cov):
            os.remove(cov)
        ppcfg = os.path.join(r["ctl"], "pp.json")
        write_json(ppcfg, {"first_dir": "../Results/%s/" % r["batch"]})
        subprocess.run([pp, ppcfg], cwd=SCRIPTS, stdout=open(os.path.join(r["ctl"], "pp.log"), "w"),
                       stderr=subprocess.STDOUT)
        if os.path.exists(cov):
            line = [l for l in open(cov) if not l.startswith("#")]
            if line:
                t = line[-1].split()
                # DIR KA KB KI BeadRadius rc CoveredArea CoverageUnion MultilayerFrac BeadX Area Excess_tension Bending_tan Bead Total_E
                row["CoverageUnion"] = float(t[7])
                row["z"] = 2 * float(t[7])
                row["BeadX"] = float(t[9])
                row["Total_E"] = float(t[14])
        bp = os.path.join(res, "Bead_0_data.txt")
        brows = [l.split() for l in open(bp) if not l.startswith("#") and l.strip()]
        bead = [float(x) for x in brows[-1][1:4]]
        try:
            row.update(mesh_stats(res, row["steps"], bead, m["radius"]))
        except Exception as e:   # keep going, the row is still useful
            print("mesh stats failed for", r["id"], e, file=sys.stderr)
    out = args.out or os.path.join(d, "summary.csv")
    with open(out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=COLS)
        w.writeheader()
        for row in rows:
            w.writerow(row)
    if out != os.path.join(d, "summary.csv"):
        import shutil
        shutil.copy(out, os.path.join(d, "summary.csv"))
    def f(x, fmt):
        return fmt % x if isinstance(x, (int, float)) else str(x)
    print("%-16s %5s %5s %5s %8s %8s %6s %6s %6s %3s %9s %9s %6s %6s %6s %6s" % (
        "run", "KI", "w~", "z", "Total_E", "steps", "wall", "state", "stop", "cv", "dE_win", "dE_1000", "faces", "h_min", "h_mean", "h_max"))
    for x in rows:
        print("%-16s %5s %5s %5s %8s %8s %6s %6s %6s %3s %9s %9s %6s %6s %6s %6s" % (
            x["run"], f(x["KI"], "%.2f"), f(x["w_t"], "%.2f"), f(x["z"], "%.3f"), f(x["Total_E"], "%.4f"),
            x["steps"], f(x["wall_s"], "%.0f"), x["state"], x["stop"], x["converged"],
            f(x["dE_last_window"], "%.1e"), f(x["dE_last_1000"], "%.1e"), x["faces"],
            f(x["h_min"], "%.3f"), f(x["h_mean"], "%.3f"), f(x["h_max"], "%.3f")))
    print("summary:", out)


# ----------------------------------------------------------------------------- plot
def cmd_plot(args):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axs = plt.subplots(1, 2, figsize=(11, 4.2))
    for tag in args.tag.split(","):
        m = read_json(os.path.join(ctl_dir(tag), "manifest.json"))
        for r in m["runs"]:
            if args.runs and not any(s in r["id"] for s in args.runs.split(",")):
                continue
            res = result_dir(r)
            if not res:
                continue
            c = conv_rows(res)
            hdr, orows = output_rows(res)
            lab = "%s %s" % (tag, r["id"])
            if orows and hdr:
                i = hdr.index("Total_E")
                axs[0].plot([x[1] for x in orows], [x[i] for x in orows], lw=0.8, label=lab)
            if c:
                steps = [x[0] for x in c]
                Ef = c[-1][2]
                axs[1].semilogy(steps, [max(abs(x[2] - Ef), 1e-8) for x in c], lw=0.8, label=lab)
    axs[0].set_xlabel("step"); axs[0].set_ylabel("Total_E (every 500 steps)")
    axs[1].set_xlabel("step"); axs[1].set_ylabel("|E(step) - E(final)|  (every 50 steps)")
    axs[1].axhline(1e-3, color="k", ls=":", lw=0.8)
    axs[0].legend(fontsize=6); axs[0].grid(alpha=0.3); axs[1].grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(args.out, dpi=130)
    print("wrote", args.out)


# ----------------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    l = sub.add_parser("launch")
    l.add_argument("--tag", required=True)
    l.add_argument("--KA", type=float, required=True)
    l.add_argument("--KB", type=float, default=1.0)
    l.add_argument("--radius", type=float, default=1.0)
    l.add_argument("--KIs", required=True, help="1.4,1.6 or start:stop:step")
    l.add_argument("--inits", default="flat,full")
    l.add_argument("--shell_width", type=float, default=None)
    l.add_argument("--shell_power", type=float, default=None)
    l.add_argument("--size_min", type=float, default=0.05)
    l.add_argument("--size_max", type=float, default=0.5)
    l.add_argument("--switch_normal", type=int, default=30000, help="Switch_times[1]: step of the BFGS -> BFGS-Normal switch (old runs: 30000)")
    l.add_argument("--timesteps", type=int, default=None, help="hard cap of the number of steps (template: 200000)")
    l.add_argument("--set", action="append", help="json override 'a.b.2=value' (value parsed as json), repeatable, e.g. stopping.window=100 remesher.refine_angle=0.3")
    l.add_argument("--extra", default="", help="extra arguments for Create_subjob_planar.py, e.g. '--shell_power 2'")
    l.add_argument("--binary", default=DEFAULT_BINARY)
    l.add_argument("--snapshot_binary", action="store_true", help="copy the binaries into the control folder first")
    l.add_argument("--jobs", type=int, default=3, help="max concurrent submissions (<= 3 for this study)")
    l.add_argument("--timeout", type=int, default=0, help="kill a run after this many wall seconds (no final state then)")
    l.add_argument("--dry_run", action="store_true")
    w = sub.add_parser("_worker"); w.add_argument("--tag", required=True); w.add_argument("--jobs", type=int, default=3)
    s = sub.add_parser("status"); s.add_argument("--tag", required=True)
    wt = sub.add_parser("wait"); wt.add_argument("--tag", required=True); wt.add_argument("--max_s", type=float, default=540)
    c = sub.add_parser("collect"); c.add_argument("--tag", required=True); c.add_argument("--out", default=None)
    p = sub.add_parser("plot"); p.add_argument("--tag", required=True, help="comma separated tags")
    p.add_argument("--out", required=True); p.add_argument("--runs", default="", help="comma separated substrings of run ids")
    args = ap.parse_args()
    {"launch": cmd_launch, "_worker": cmd_worker, "status": cmd_status, "wait": cmd_wait, "collect": cmd_collect, "plot": cmd_plot}[args.cmd](args)


if __name__ == "__main__":
    main()
