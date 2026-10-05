#!/usr/bin/env python3
"""Run the benchmark configs and show where the time goes.

Usage:
    python3 run_bench.py --bin ../build/bin/main_cluster --out /tmp/bench [--only wrapping] [-j 2]

Each config in ../regression/bench/ runs with first_dir overridden to
<out>/<config name>/ (the runs are independent, -j runs them in parallel).
The table comes from the Timing.txt of each run: wall time, remeshing,
the energy and gradient evaluations (count and time) and the most
expensive energy terms. For exactness use regression_compare.py with
--config-dir ../regression/bench.
"""
import argparse
import json
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
BENCH_DIR = os.path.join(HERE, '..', 'regression', 'bench')


def run(binary, config_path, out_dir):
    with open(config_path) as f:
        data = json.load(f)
    os.makedirs(out_dir, exist_ok=True)
    data['first_dir'] = out_dir + '/'
    tmp_config = os.path.join(out_dir, 'config.json')
    with open(tmp_config, 'w') as f:
        json.dump(data, f, indent=4)
    with open(os.path.join(out_dir, 'log.txt'), 'w') as log:
        proc = subprocess.run([os.path.abspath(binary), tmp_config],
                              stdout=log, stderr=subprocess.STDOUT, cwd=out_dir)
    return proc.returncode


def newest_run(out_dir):
    runs = [d for d in os.listdir(out_dir) if d.isdigit()]
    return os.path.join(out_dir, max(runs, key=int)) if runs else None


def read_timing(run_dir):
    res = {}
    with open(os.path.join(run_dir, 'Timing.txt')) as f:
        for line in f:
            k, v = line.split()
            res[k] = float(v)
    return res


def report(name, t):
    print('== %s: %d vertices at the end, %d remeshes' % (name, t['vertices'], t['remeshes']))
    print('   wall %.1f s = integrate %.1f + remesh %.1f + save %.1f' % (
        t['wall_ms'] / 1e3, t['integrate_ms'] / 1e3, t['remesh_ms'] / 1e3, t['save_ms'] / 1e3))
    print('   energies:  %7d calls %8.1f s (%.3f ms each)' % (
        t['energy_calls'], t['energy_ms'] / 1e3, t['energy_ms'] / max(t['energy_calls'], 1)))
    print('   gradients: %7d calls %8.1f s (%.3f ms each)' % (
        t['gradient_calls'], t['gradient_ms'] / 1e3, t['gradient_ms'] / max(t['gradient_calls'], 1)))
    terms = sorted(((v, k) for k, v in t.items() if k.startswith(('energy_ms_', 'gradient_ms_'))), reverse=True)
    for v, k in terms[:6]:
        print('     %-40s %8.1f s' % (k, v / 1e3))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--bin', required=True)
    parser.add_argument('--out', required=True)
    parser.add_argument('--only', nargs='*', help='subset of config names')
    parser.add_argument('-j', type=int, default=1, help='runs in parallel')
    args = parser.parse_args()

    configs = sorted(c for c in os.listdir(BENCH_DIR) if c.endswith('.json'))
    if args.only:
        configs = [c for c in configs if any(o in c for o in args.only)]
    jobs = [(c[:-5], os.path.join(BENCH_DIR, c), os.path.abspath(os.path.join(args.out, c[:-5]))) for c in configs]
    with ThreadPoolExecutor(max_workers=args.j) as pool:
        codes = list(pool.map(lambda j: run(args.bin, j[1], j[2]), jobs))

    failed = False
    for (name, _, out_dir), code in zip(jobs, codes):
        run_dir = newest_run(out_dir)
        if code != 0 or run_dir is None or not os.path.exists(os.path.join(run_dir, 'Timing.txt')):
            print('== %s: failed (exit %d), see %s/log.txt' % (name, code, out_dir))
            failed = True
            continue
        report(name, read_timing(run_dir))
    sys.exit(1 if failed else 0)


if __name__ == '__main__':
    main()
