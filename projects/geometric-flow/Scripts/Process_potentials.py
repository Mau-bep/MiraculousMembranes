import os
import re
import sys
import json
import glob
import argparse
import subprocess
import tempfile

# Runs PostProcessing on every batch folder of the potentials scan (Mem_planar_potentials_up.sh, ../Results/Potentials_<cond>_<init>/)
# and gathers all the Coverage_data.txt in one file, like CoverageShellSum/CoverageSum.txt, with the columns
#   Cond Init Profile shell_width shell_alpha shell_power  +  the columns of Coverage_data.txt (DIR KA KB KI BeadRadius rc CoveredArea ...
#   ... BeadX Excess_tension Bending_tan Bead Total_E)
# so the flat and full starts of every condition can be compared (the equilibrium state is the one with the lower Total_E).
#
# PostProcessing needs a config only for first_dir, the bead and the energy constants: it is built from the Input_file.json that
# every run saves in its own folder (first_dir replaced by the batch folder), so the Config_files folder is not needed.
# PostProcessing APPENDS to Coverage_data.txt: the old file of the batch is removed first, otherwise every rerun would double the rows
# (--keep_old skips the PostProcessing step and only gathers what is already there).
# Runs without Input_file.json (unfinished) are skipped by PostProcessing itself.
#
# Usage (from Scripts/):  python3 Process_potentials.py
#         python3 Process_potentials.py --prefix Potentials --results ../Results --out ../Results/Potentials_sum/CoverageSum.txt

parser = argparse.ArgumentParser(description="Runs PostProcessing on the batches of the potentials scan and gathers them in one file.")
parser.add_argument("--results", type=str, default="../Results", help="folder with the batch folders (default ../Results)")
parser.add_argument("--prefix", type=str, default="Potentials", help="batches are <prefix>_<cond>_<flat|full> (default Potentials)")
parser.add_argument("--out", type=str, default=None, help="output file (default <results>/<prefix>_sum/CoverageSum.txt)")
parser.add_argument("--post", type=str, default="../build/bin/PostProcessing", help="PostProcessing executable")
parser.add_argument("--keep_old", action="store_true", help="do not run PostProcessing, gather the existing Coverage_data.txt files")
parser.add_argument("--append", action="store_true", help="append to the output file instead of rewriting it (rows can then be duplicated)")
args = parser.parse_args()

results = os.path.abspath(args.results)
out = os.path.abspath(args.out if args.out else os.path.join(results, args.prefix + "_sum", "CoverageSum.txt"))
post = os.path.abspath(args.post)
tag = re.compile(r"^{}_(?P<cond>.+)_(?P<init>flat|full)$".format(re.escape(args.prefix)))


def numeric_subdirs(batch):
    return sorted((d for d in os.listdir(batch) if d[:1].isdigit() and os.path.isdir(os.path.join(batch, d))), key=lambda d: (len(d), d))


def run_post(batch):
    """PostProcessing on one batch folder. Returns True if it ran and wrote Coverage_data.txt."""
    subdirs = numeric_subdirs(batch)
    configs = [os.path.join(batch, d, "Input_file.json") for d in subdirs if os.path.isfile(os.path.join(batch, d, "Input_file.json"))]
    if not configs:
        print("  no finished runs (no Input_file.json), skipped")
        return False
    with open(configs[0]) as f:
        cfg = json.load(f)
    cfg["first_dir"] = batch + "/"
    old = os.path.join(batch, "Coverage_data.txt")
    if os.path.exists(old):
        os.remove(old)
    with tempfile.NamedTemporaryFile("w", suffix=".json", delete=False) as tmp:
        json.dump(cfg, tmp)
        tmp_path = tmp.name
    try:
        proc = subprocess.run([post, tmp_path], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    finally:
        os.remove(tmp_path)
    # PostProcessing's main returns 1 on success
    if not os.path.exists(old):
        print("  PostProcessing failed (exit {}), last lines:\n    {}".format(proc.returncode, "\n    ".join(proc.stdout.strip().splitlines()[-4:])))
        return False
    return True


def profile_of(batch):
    """Profile, shell_width, shell_alpha and shell_power of the bead of the first finished run of the batch."""
    for d in numeric_subdirs(batch):
        p = os.path.join(batch, d, "Input_file.json")
        if os.path.isfile(p):
            with open(p) as f:
                b = json.load(f)["Beads"][0]
            return (b.get("shell_profile", "cosine"), b.get("shell_width", 0.25), b.get("shell_alpha", "-"), b.get("shell_power", 1))
    return ("?", "?", "?", "?")


batches = []
for name in sorted(os.listdir(results)):
    m = tag.match(name)
    if m and os.path.isdir(os.path.join(results, name)):
        batches.append((name, m.group("cond"), m.group("init")))
if not batches:
    sys.exit("No folders {}_<cond>_<flat|full> in {}".format(args.prefix, results))
if not args.keep_old and not os.path.isfile(post):
    sys.exit("PostProcessing not found at {} (build it, or use --keep_old)".format(post))

os.makedirs(os.path.dirname(out), exist_ok=True)
header = None
rows = 0
mode = "a" if args.append and os.path.exists(out) else "w"
if mode == "a":
    with open(out) as f:
        header = f.readline().rstrip("\n") or None
with open(out, mode) as w:
    for name, cond, init in batches:
        batch = os.path.join(results, name)
        print(name)
        if not args.keep_old and not run_post(batch):
            continue
        path = os.path.join(batch, "Coverage_data.txt")
        if not os.path.isfile(path):
            print("  no Coverage_data.txt, skipped")
            continue
        prefix = "{} {} {} {} {} {}".format(cond, init, *profile_of(batch))
        n = 0
        with open(path) as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                if line.startswith("#"):
                    # the header, repeated in the file when PostProcessing was run more than once
                    cols = "# Cond Init Profile shell_width shell_alpha shell_power " + line.lstrip("#").strip()
                    if header is None:
                        header = cols
                        w.write(header + "\n")
                    elif cols != header:
                        print("  warning: the columns of {} differ from the first batch".format(path))
                    continue
                if header is None:
                    sys.exit("{}: data before the header".format(path))
                w.write(prefix + " " + line + "\n")
                n += 1
        print("  {} rows ({} {} {}, shell_width {})".format(n, cond, init, prefix.split()[2], prefix.split()[3]))
        rows += n
print("{} rows from {} batches in {}".format(rows, len(batches), out))
