#!/usr/bin/env python3
"""Calibrate f_tol of the "quality" adaptive remeshing ("adapt_remesh": "quality").

The idea: take one representative input file and run it
  * ref        remeshing as often as the loop allows (fixed, remesh_every = 1,
               which remeshes every second step): the "always fresh mesh" answer
  * ref_rep<k> the same with rescale changed by k * 1e-6, to measure how much
               the observables move from a perturbation that should not matter
               (chaotic runs drift apart anyway; this is the noise floor)
  * current    the input file's own remeshing settings, for comparison
  * ftol_<x>   the quality mode with f_tol = x
then pick the largest f_tol whose observables stay as close to ref as the
noise floor (or --tol, whichever is larger), and report what it costs.

Steps:
    # 1. write the configs (paths in the input file are resolved against --cwd,
    #    by default projects/geometric-flow/build, where the input files expect to run)
    python3 calibrate_remesh_tol.py prepare ../Config_files/My_run.json --out ../Calibration/My_run \
        --ftol 0.005 0.01 0.02 0.05 0.1 --timesteps 20000

    # 2a. run them here (each job gets an equal share of the cores, so the timings compare)
    python3 calibrate_remesh_tol.py run ../Calibration/My_run --exe ../build/bin/main_cluster -j 4
    # 2b. or on the cluster: every line of <out>/run_all.sh is one job
    #     (MAIN_CLUSTER=/path/to/main_cluster), run them with the same resources

    # 3. compare
    python3 calibrate_remesh_tol.py analyze ../Calibration/My_run [--tol 0.01]

Every run writes Remesh_log.txt (mesh quality and energy before/after each
remesh) and Timing.txt (wall time split into remeshing / integrating / saving).
The logs also give the defect growth rate per step, from which analyze
says which f_tol gives which remesh interval.

What came out of testing it on regression config d (coverage wrapping):
remeshing raises the energy every time (smoothing, splits and collapses move
the discrete energy up), so the frequently remeshed ref relaxes more slowly and
at a fixed step it is simply less relaxed. Compare converged end states: analyze
warns when a run has not converged, and only recommends a value when ref has.
"""
import argparse
import concurrent.futures
import json
import math
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_CWD = os.path.normpath(os.path.join(HERE, '..', 'build'))
REMESH_SWITCHES = ('Remesh_always', 'Restore_remeshing', 'Adapt_remesh', 'No_remesh')


# ----------------------------------------------------------------------------- prepare

def absolute(path, cwd):
    return path if os.path.isabs(path) else os.path.normpath(os.path.join(cwd, path))


def tag(x):
    return ('%g' % x).replace('-', 'm')


def prepare(args):
    with open(args.base) as f:
        base = json.load(f)
    out = os.path.abspath(args.out)
    for sub in ('configs', 'runs', 'logs'):
        os.makedirs(os.path.join(out, sub), exist_ok=True)

    base['init_file'] = absolute(base['init_file'], args.cwd)
    if not os.path.exists(base['init_file']):
        sys.exit('init_file %s does not exist (paths are resolved against --cwd %s)' % (base['init_file'], args.cwd))
    if args.timesteps:
        base['timesteps'] = args.timesteps
    if not base.get('remeshing', True):
        sys.exit('The input file has "remeshing": false, there is nothing to calibrate')
    base['Remesh_log'] = True
    base.pop('Subfolder', None)
    clash = [s for s in base.get('Switches', []) if s in REMESH_SWITCHES]
    if clash:
        print('Warning: the switches %s change the remeshing during the run; the comparison is only '
              'clean before the first of them' % clash)

    runs = {}

    def add(name, data, kind, ftol=None):
        data = json.loads(json.dumps(data))
        data['first_dir'] = os.path.join(out, 'runs', name) + '/'
        path = os.path.join(out, 'configs', name + '.json')
        with open(path, 'w') as f:
            json.dump(data, f, indent=4)
        runs[name] = {'config': path, 'first_dir': data['first_dir'], 'kind': kind, 'ftol': ftol}

    ref = dict(base, adapt_remesh=False, remesh_every=1)
    ref.pop('remesh_quality', None)
    add('ref', ref, 'ref')
    for k in range(1, args.replicas + 1):
        rep = dict(ref, rescale=base.get('rescale', 1.0) * (1.0 + k * args.perturbation))
        add('ref_rep%d' % k, rep, 'replica')
    add('current', base, 'current')
    for x in args.ftol:
        q = dict(base, adapt_remesh='quality',
                 remesh_quality={'f_tol': x, 'min_every': args.min_every, 'max_every': args.max_every})
        q['remesh_every'] = 1
        add('ftol_' + tag(x), q, 'quality', x)

    manifest = {'base': os.path.abspath(args.base), 'cwd': args.cwd, 'timesteps': base['timesteps'],
                'tol': args.tol, 'runs': runs}
    with open(os.path.join(out, 'manifest.json'), 'w') as f:
        json.dump(manifest, f, indent=4)
    with open(os.path.join(out, 'run_all.sh'), 'w') as f:
        f.write('# One job per line, all paths absolute. Set MAIN_CLUSTER to the binary.\n')
        for name, r in runs.items():
            f.write('"${MAIN_CLUSTER:-main_cluster}" %s > %s 2>&1\n'
                    % (r['config'], os.path.join(out, 'logs', name + '.log')))
    print('Wrote %d configs to %s' % (len(runs), os.path.join(out, 'configs')))
    print('Run them with:  %s run %s --exe <main_cluster> -j <N>' % (sys.argv[0], out))
    print('or submit the lines of %s' % os.path.join(out, 'run_all.sh'))


# ----------------------------------------------------------------------------- run

def run_one(exe, name, r, log_dir, threads):
    env = dict(os.environ, OMP_NUM_THREADS=str(threads))
    os.makedirs(r['first_dir'], exist_ok=True)
    with open(os.path.join(log_dir, name + '.log'), 'w') as log:
        proc = subprocess.run([exe, r['config']], stdout=log, stderr=subprocess.STDOUT, cwd=r['first_dir'], env=env)
    return name, proc.returncode


def run(args):
    manifest = load_manifest(args.dir)
    exe = os.path.abspath(args.exe)
    todo = {n: r for n, r in manifest['runs'].items() if args.force or run_dir(r) is None}
    if not todo:
        print('All runs are finished (use --force to run them again)')
        return
    threads = max(1, (os.cpu_count() or 1) // args.jobs)
    print('Running %d jobs, %d at a time, %d threads each' % (len(todo), args.jobs, threads))
    log_dir = os.path.join(args.dir, 'logs')
    with concurrent.futures.ThreadPoolExecutor(args.jobs) as pool:
        futures = [pool.submit(run_one, exe, n, r, log_dir, threads) for n, r in todo.items()]
        for fut in concurrent.futures.as_completed(futures):
            name, rc = fut.result()
            print('  %-14s exit %d' % (name, rc))


# ----------------------------------------------------------------------------- analyze

def load_manifest(d):
    with open(os.path.join(d, 'manifest.json')) as f:
        return json.load(f)


def run_dir(r):
    """Newest numbered folder of a run that finished (has Timing.txt)."""
    root = r['first_dir']
    if not os.path.isdir(root):
        return None
    nums = sorted((int(n) for n in os.listdir(root) if n.isdigit()), reverse=True)
    for n in nums:
        d = os.path.join(root, str(n))
        if os.path.exists(os.path.join(d, 'Timing.txt')):
            return d
    return None


def numeric_rows(path):
    rows = []
    with open(path) as f:
        for line in f:
            if not line.strip() or line.startswith('#'):
                continue
            try:
                rows.append([float(t) for t in line.split()])
            except ValueError:
                continue  # header line
    return rows


def read_run(d):
    res = {}
    with open(os.path.join(d, 'Timing.txt')) as f:
        for line in f:
            k, v = line.split()
            res[k] = float(v)

    with open(os.path.join(d, 'Output_data.txt')) as f:
        header = f.readline().split()
    n_energies = header.index('Total_E') - 4
    out = numeric_rows(os.path.join(d, 'Output_data.txt'))
    # columns: time step Volume Area <energies> Total_E ...
    res['traj'] = {int(r[1]): r[4 + n_energies] for r in out}
    last = out[-1]
    res.update(step=int(last[1]), V=last[2], A=last[3], E=last[4 + n_energies], E0=out[0][4 + n_energies])

    beads = []
    i = 0
    while os.path.exists(os.path.join(d, 'Bead_%d_data.txt' % i)):
        b = numeric_rows(os.path.join(d, 'Bead_%d_data.txt' % i))
        if b:
            beads.append(b[-1][1:4])
        i += 1
    res['beads'] = beads

    log = numeric_rows(os.path.join(d, 'Remesh_log.txt'))
    # t since ops nV bad_b long_b short_b flip_b bad_a long_a short_a flip_a E_b E_a
    res['log'] = log
    # Energy put in by the remeshes (the first one, on the input mesh, is left out)
    res['kick'] = sum(r[13] - r[12] for r in log if r[0] > 0)
    removed = res['E0'] - res['E'] + res['kick']  # what the minimizer took out
    res['undone'] = res['kick'] / removed if removed > 0 else float('nan')
    return res


def convergence(run, conv):
    """(relative energy change over the last 10% of the steps, first saved step
    after which E stays within conv of its final value), both relative to the
    total energy change of the run."""
    steps = sorted(run['traj'])
    E = [run['traj'][t] for t in steps]
    span = abs(E[0] - E[-1]) or 1e-300
    t90 = min(steps, key=lambda t: abs(t - 0.9 * steps[-1]))
    tail = abs(run['traj'][t90] - E[-1]) / span
    t_conv = steps[-1]
    for t, e in zip(reversed(steps), reversed(E)):
        if abs(e - E[-1]) > conv * span:
            break
        t_conv = t
    return tail, t_conv


def quality_triggers(log, ftol):
    """Remeshes where the defect growth since the previous one exceeded f_tol;
    the others were asked for by the integrator or hit max_every."""
    return sum(1 for prev, cur in zip(log, log[1:]) if cur[4] - prev[8] > ftol)


def drift_rate(logs, min_interval=10):
    """Median growth of the bad edge fraction per step, from the remesh
    intervals of at least min_interval steps (shorter ones are dominated by
    what the previous remesh left behind)."""
    rates = []
    for log in logs:
        for prev, cur in zip(log, log[1:]):
            if cur[1] >= min_interval:
                rates.append((cur[4] - prev[8]) / cur[1])
    rates.sort()
    return rates[len(rates) // 2] if rates else float('nan')


def deviations(run, ref):
    """How far a run ends up from the reference, each as a dimensionless number."""
    common = sorted(set(run['traj']) & set(ref['traj']))
    e_ref = [ref['traj'][t] for t in common]
    scale = (max(e_ref) - min(e_ref)) if len(e_ref) > 1 else 0.0
    scale = scale if scale > 0 else max(abs(ref['E']), 1e-300)
    traj = max((abs(run['traj'][t] - ref['traj'][t]) for t in common), default=float('nan')) / scale
    bead = 0.0
    for p, q in zip(run['beads'], ref['beads']):
        bead = max(bead, math.sqrt(sum((a - b) ** 2 for a, b in zip(p, q))))
    size = math.sqrt(ref['A'] / (4 * math.pi))  # radius of the sphere with the same area
    return {
        'E_end': abs(run['E'] - ref['E']) / scale,
        'E_traj': traj,
        'V_end': abs(run['V'] - ref['V']) / abs(ref['V']),
        'A_end': abs(run['A'] - ref['A']) / abs(ref['A']),
        'bead': bead / size,
    }


METRICS = ('E_end', 'V_end', 'A_end', 'bead', 'E_traj')
METRIC_HELP = {
    'E_end': '|E - E_ref| at the last step, over the range of E_ref(t)',
    'V_end': 'relative volume difference at the end',
    'A_end': 'relative area difference at the end',
    'bead': 'largest bead position difference at the end, over the membrane size sqrt(A/4pi)',
    'E_traj': 'largest |E(t) - E_ref(t)| over the saved steps, same scale as E_end',
}
# For a relaxation only the end state matters; add E_traj (--accept) when the
# path itself is the result (pulling, prescribed bead motion)
DEFAULT_ACCEPT = ('E_end', 'V_end', 'A_end', 'bead')


def analyze(args):
    manifest = load_manifest(args.dir)
    tol = args.tol if args.tol is not None else manifest.get('tol', 0.01)
    accept = args.accept or list(DEFAULT_ACCEPT)
    data, missing = {}, []
    for name, r in manifest['runs'].items():
        d = run_dir(r)
        if d is None:
            missing.append(name)
            continue
        data[name] = read_run(d)
        data[name].update(kind=r['kind'], ftol=r['ftol'])
    if missing:
        print('Not finished (skipped): %s' % ', '.join(missing))
    if args.ref not in data:
        sys.exit('The reference run "%s" has not finished' % args.ref)
    ref = data[args.ref]

    for name, r in data.items():
        if r['step'] < ref['step']:
            print('Warning: %s stopped at step %d, the reference at %d' % (name, r['step'], ref['step']))
        r['dev'] = deviations(r, ref) if name != args.ref else {m: 0.0 for m in METRICS}
        r['tail'], r['t_conv'] = convergence(r, args.conv)
        r['by_tol'] = quality_triggers(r['log'], r['ftol']) if r['kind'] == 'quality' else None

    # Spread between ref and its perturbed copies (only meaningful against ref)
    reps = [r for r in data.values() if r['kind'] == 'replica'] if args.ref == 'ref' else []
    noise = {m: max([r['dev'][m] for r in reps], default=0.0) for m in METRICS}
    limit = {m: max(tol, args.noise_factor * noise[m]) for m in METRICS}

    print('\nDeviation from %s (dimensionless):' % args.ref)
    for m in METRICS:
        print('  %-7s %s' % (m, METRIC_HELP[m]))
    print('Accepted when %s <= max(tol = %g, %g x replica spread)\n' % (', '.join(accept), tol, args.noise_factor))

    order = ['ref'] + sorted((n for n in data if data[n]['kind'] == 'replica')) + \
        (['current'] if 'current' in data else []) + \
        sorted((n for n in data if data[n]['kind'] == 'quality'), key=lambda n: data[n]['ftol'])
    order = [n for n in order if n in data]
    head = '%-12s %8s %6s %8s %7s %7s %7s %6s %7s' % ('run', 'remeshes', 'by_tol', 'interval', 'wall_s', 'speedup',
                                                    'undone', 'tail', 't_conv')
    head += ''.join(' %7s' % m for m in METRICS) + '  ok'
    print(head)
    print('-' * len(head))
    rows_csv = []
    best = None
    for name in order:
        r = data[name]
        interval = (r['step'] + 1) / max(1.0, r['remeshes'])
        speed = ref['wall_ms'] / r['wall_ms'] if r['wall_ms'] > 0 else float('nan')
        ok = all(r['dev'][m] <= limit[m] for m in accept)
        flag = '' if name == args.ref else ('yes' if ok else 'no')
        line = '%-12s %8d %6s %8.1f %7.1f %7.2f %6.0f%% %6.3f %7d' % (
            name, r['remeshes'], '-' if r['by_tol'] is None else r['by_tol'], interval, r['wall_ms'] / 1e3, speed,
            100 * r['undone'], r['tail'], r['t_conv'])
        line += ''.join(' %7.2g' % r['dev'][m] for m in METRICS) + '  ' + flag
        print(line)
        rows_csv.append([name, '' if r['ftol'] is None else r['ftol'], r['remeshes'],
                         '' if r['by_tol'] is None else r['by_tol'], interval, r['wall_ms'], r['remesh_ms'], speed,
                         r['kick'], r['undone'], r['tail'], r['t_conv']] + [r['dev'][m] for m in METRICS] + [flag])
        if r['kind'] == 'quality' and ok and (best is None or r['ftol'] > data[best]['ftol']):
            best = name

    print('by_tol  remeshes triggered by f_tol; the rest were asked for by the integrator (BFGS remesh_flag,')
    print('        Newton fallback) or hit max_every')
    print('undone  energy the remeshes added, as a share of what the minimizer removed')
    print('tail    energy change over the last 10%% of the steps, over the total change (converged: < %g)' % args.conv)
    print('t_conv  first saved step after which E stays within %g of its final value (same scale)' % args.conv)
    print('wall_s is only comparable between runs with the same load')

    if reps:
        print('\nReplica spread (noise floor): ' + ', '.join('%s %.2g' % (m, noise[m]) for m in METRICS))
    rate = drift_rate([r['log'] for r in data.values()])
    if rate == rate and rate > 0:
        print('Defect growth: %.3g of the edges per step, so f_tol = %.3g remeshes about every 10 steps, '
              '%.3g about every 50 (max_every caps it)' % (rate, 10 * rate, 50 * rate))

    unconverged = [n for n in order if data[n]['tail'] > args.conv]
    if unconverged:
        print('\nNot converged: %s. Their end states are still moving, so a deviation there measures how' %
              ', '.join(unconverged))
        print('fast each run relaxes as much as where it ends up. Run longer (prepare --timesteps) for a clean answer.')
    if ref['undone'] == ref['undone'] and ref['undone'] > 0.1:
        print('\nRemeshing in %s undoes %.0f%% of the relaxation, so it relaxes more slowly than runs that remesh'
              % (args.ref, 100 * ref['undone']))
        print('less; compare converged end states rather than states at the same step.')

    quality = sorted((n for n in data if data[n]['kind'] == 'quality'), key=lambda n: data[n]['ftol'])
    if ref['tail'] > args.conv:
        print('\nNo recommendation: %s has not converged, so there is no end state to compare with.' % args.ref)
    elif quality and all(data[n]['by_tol'] == 0 for n in quality):
        print('\nNo f_tol in the scan ever triggered a remesh: the mesh stayed within all of them, so this run cannot')
        print('tell them apart. Use a longer run, one where the membrane deforms more, or smaller f_tol values.')
    elif best:
        b = data[best]
        print('\nRecommended: f_tol = %g (largest accepted): %d remeshes instead of %d, %.2fx the speed of %s'
              % (b['ftol'], b['remeshes'], ref['remeshes'], ref['wall_ms'] / b['wall_ms'], args.ref))
        if best == quality[-1]:
            print('  It is the largest value scanned; try larger ones too.')
        if b['by_tol'] < b['remeshes'] / 2:
            print('  Only %d of its remeshes came from f_tol, the rest from max_every or the integrator: at this'
                  % b['by_tol'])
            print('  tolerance the run is limited by max_every (%s), which is worth scanning too.'
                  % 'prepare --max-every')
        bad = [n for n in quality if data[n]['ftol'] < b['ftol'] and not
               all(data[n]['dev'][m] <= limit[m] for m in accept)]
        if bad:
            print('  Smaller values that failed: %s; the acceptance is not monotonic, so treat this with care.'
                  % ', '.join(bad))
    elif quality:
        print('\nNo f_tol passed; scan smaller values (or loosen --tol).')

    if args.csv:
        with open(args.csv, 'w') as f:
            f.write(','.join(['run', 'ftol', 'remeshes', 'by_tol', 'interval', 'wall_ms', 'remesh_ms', 'speedup',
                              'kick', 'undone', 'tail', 't_conv'] + list(METRICS) + ['ok']) + '\n')
            for row in rows_csv:
                f.write(','.join(str(x) for x in row) + '\n')
        print('Wrote %s' % args.csv)


# ----------------------------------------------------------------------------- main

def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest='cmd', required=True)

    p = sub.add_parser('prepare', help='write the calibration configs')
    p.add_argument('base', help='representative input file')
    p.add_argument('--out', required=True)
    p.add_argument('--ftol', type=float, nargs='+', default=[0.005, 0.01, 0.02, 0.05, 0.1])
    p.add_argument('--replicas', type=int, default=2, help='perturbed copies of ref (noise floor)')
    # Output_data.txt has 6 significant digits, so smaller changes are invisible
    p.add_argument('--perturbation', type=float, default=1e-6, help='relative change of rescale per replica')
    p.add_argument('--min-every', type=int, default=1)
    p.add_argument('--max-every', type=int, default=100)
    p.add_argument('--timesteps', type=int, help='override the number of steps')
    p.add_argument('--tol', type=float, default=0.01, help='default acceptance tolerance stored for analyze')
    p.add_argument('--cwd', default=DEFAULT_CWD, help='folder relative paths in the input file refer to')

    p = sub.add_parser('run', help='run the configs locally')
    p.add_argument('dir')
    p.add_argument('--exe', required=True)
    p.add_argument('-j', '--jobs', type=int, default=1)
    p.add_argument('--force', action='store_true', help='also rerun finished runs')

    p = sub.add_parser('analyze', help='compare the runs to the reference')
    p.add_argument('dir')
    p.add_argument('--tol', type=float, help='acceptance tolerance (default: the one given to prepare)')
    p.add_argument('--noise-factor', type=float, default=2.0)
    p.add_argument('--ref', default='ref', help='run to compare against (e.g. the smallest f_tol)')
    p.add_argument('--accept', nargs='+', choices=METRICS, help='metrics that decide (default: the end state)')
    p.add_argument('--conv', type=float, default=0.01, help='convergence tolerance, relative to the energy change')
    p.add_argument('--csv', help='also write the table to this file')

    args = parser.parse_args()
    {'prepare': prepare, 'run': run, 'analyze': analyze}[args.cmd](args)


if __name__ == '__main__':
    main()
