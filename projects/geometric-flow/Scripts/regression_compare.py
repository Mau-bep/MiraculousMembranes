#!/usr/bin/env python3
"""Run the regression configs with a reference and a new main_cluster binary and
compare their outputs numerically.

Usage:
    python3 regression_compare.py --ref ../build/bin/main_cluster_ref \
        --new ../build/bin/main_cluster --out /tmp/regression [--rtol 1e-10]

Each config in ../regression/configs/ is run once per binary with
first_dir overridden to <out>/<ref|new>/<config name>/. Every text output
(Output_data.txt, Bead_*_data.txt, *.obj, ...) found in the run folder is then
compared token by token: numeric tokens with a relative tolerance, the rest exactly.
"""
import argparse
import json
import math
import os
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
CONFIG_DIR = os.path.join(HERE, '..', 'regression', 'configs')
# Files whose content depends on wall-clock time
SKIP = ('timing', 'Timing')


def run(binary, config_path, out_dir, log_path):
    with open(config_path) as f:
        data = json.load(f)
    if os.path.exists(out_dir):
        shutil.rmtree(out_dir)
    os.makedirs(out_dir)
    data['first_dir'] = out_dir + '/'
    tmp_config = os.path.join(out_dir, 'config.json')
    with open(tmp_config, 'w') as f:
        json.dump(data, f, indent=4)
    with open(log_path, 'w') as log:
        proc = subprocess.run([os.path.abspath(binary), tmp_config],
                              stdout=log, stderr=subprocess.STDOUT, cwd=out_dir)
    # main_cluster writes into first_dir/1/
    return proc.returncode, os.path.join(out_dir, '1')


def to_float(tok):
    try:
        return float(tok)
    except ValueError:
        return None


def compare_file(a, b, rtol, atol):
    with open(a) as fa, open(b) as fb:
        ta = fa.read().split()
        tb = fb.read().split()
    if len(ta) != len(tb):
        return 'token count %d vs %d' % (len(ta), len(tb)), float('inf')
    worst = 0.0
    for i, (x, y) in enumerate(zip(ta, tb)):
        if x == y:
            continue
        fx, fy = to_float(x), to_float(y)
        if fx is None or fy is None:
            return 'token %d: %r vs %r' % (i, x, y), float('inf')
        if math.isnan(fx) and math.isnan(fy):
            continue
        err = abs(fx - fy) / max(abs(fx), abs(fy), atol)
        worst = max(worst, err)
    if worst > rtol:
        return 'max rel err %.3e' % worst, worst
    return None, worst


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--ref', required=True)
    parser.add_argument('--new', required=True)
    parser.add_argument('--out', required=True)
    parser.add_argument('--rtol', type=float, default=1e-10)
    parser.add_argument('--atol', type=float, default=1e-12)
    parser.add_argument('--only', nargs='*', help='subset of config names')
    parser.add_argument('--skip-ref', action='store_true',
                        help='reuse existing reference outputs')
    args = parser.parse_args()

    configs = sorted(c for c in os.listdir(CONFIG_DIR) if c.endswith('.json'))
    if args.only:
        configs = [c for c in configs if any(o in c for o in args.only)]

    failed = False
    for cfg in configs:
        name = cfg[:-5]
        path = os.path.join(CONFIG_DIR, cfg)
        ref_root = os.path.abspath(os.path.join(args.out, 'ref', name))
        new_root = os.path.abspath(os.path.join(args.out, 'new', name))
        if args.skip_ref and os.path.isdir(os.path.join(ref_root, '1')):
            rc_ref, ref_dir = 0, os.path.join(ref_root, '1')
        else:
            rc_ref, ref_dir = run(args.ref, path, ref_root, ref_root + '.log')
        rc_new, new_dir = run(args.new, path, new_root, new_root + '.log')
        print('== %s (exit ref=%d new=%d)' % (name, rc_ref, rc_new))
        if rc_ref != rc_new:
            failed = True
        ref_files = sorted(f for f in os.listdir(ref_dir) if not any(s in f for s in SKIP))
        new_files = set(os.listdir(new_dir)) if os.path.isdir(new_dir) else set()
        worst_all = 0.0
        for f in ref_files:
            if f == 'Input_file.json':
                continue
            if f not in new_files:
                print('   MISSING  %s' % f)
                failed = True
                continue
            msg, worst = compare_file(os.path.join(ref_dir, f), os.path.join(new_dir, f),
                                      args.rtol, args.atol)
            worst_all = max(worst_all, worst)
            if msg:
                print('   DIFF     %s: %s' % (f, msg))
                failed = True
        print('   %d files compared, max rel err %.3e' % (len(ref_files), worst_all))
    sys.exit(1 if failed else 0)


if __name__ == '__main__':
    main()
