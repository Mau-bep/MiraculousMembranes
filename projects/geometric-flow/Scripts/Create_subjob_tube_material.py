import os
import json
import math
import argparse
from jinja2 import Environment, FileSystemLoader
from cli_args import number_str
from subjob import write_subjob

# A bead walked through the big planar membrane along x, one frozen position at a time, every position relaxed from the previous one (the material-property
# sister of Mem_tube_up.sh: there each position is a separate job that pulls the bead there, here one job walks the bead out and back). The START STATE of the chain is made
# the old way: the bead is pulled by a spring through the flat membrane up to x_start (0.2 by default), where it freezes (../Templates/Tube_start.txt); the chain starts from that
# mesh (a frozen bead that jumps into the flat sheet makes the vertices slide aside and pokes through a face: Studies/tube_material/README.md). Writes the config of the frozen
# bead (../Templates/Tube_material.txt), the config of the pull, the chain file (the positions Run_tube_material_chain.py walks through) and the sbatch subjob.

parser = argparse.ArgumentParser(description="Writes the configs, the chain file and the sbatch subjob of one walk of a bead through the planar membrane (one KA, KB).")
parser.add_argument("--KA", type=number_str, default="16.0", help="surface tension constant KA (the Surface_tension energy is KA * area), default of Mem_tube_up.sh: 32.0 and 16.0")
parser.add_argument("--KB", type=number_str, default="4.0", help="bending constant KB (the Bending energy is KB H^2, kappa = KB/2)")
parser.add_argument("--radius", type=number_str, default="0.3", help="bead radius")
parser.add_argument("--x_start", type=number_str, default="0.2",
                    help="first x of the chain (default 0.2: the start state is the bead pulled by the spring to there and frozen, the force is linear in x below that and needs no chain). "
                         "-radius (or less) means no pull: the chain starts from the flat membrane with the bead below it, which only works with --max_step 0.005 or less and a fine mesh (README)")
parser.add_argument("--x_end", type=number_str, default="6.0", help="last x of the bead on the way out")
parser.add_argument("--step", type=number_str, default="0.1", help="distance between two positions of the list")
parser.add_argument("--max_step", type=number_str, default="0.025",
                    help="largest move of the bead between two relaxations: a bigger move of the list goes through waypoints (default 0.025: the positions of the list are 0.1 apart, each is "
                         "reached in four relaxations, saved as runs too and left out of the curves; the tube snaps from the dimple near x = 1.8-2.6, where small moves matter most); 0 walks the list as given")
parser.add_argument("--no_return", action="store_true", help="do not walk back from x_end (the way back is saved as runs too, in its own folder, it shows the hysteresis)")
parser.add_argument("--x_back", type=number_str, default=None,
                    help="last x of the way back: it goes x_end -> x_start and on below it to here (default 0.0: 0.2 -> 6.0 -> 0.2 -> 0.0; -radius is the flat state and never needs the chain; the start "
                         "state of the way out is the pull, the positions below x_start are reached by the chain itself)")
parser.add_argument("--back_dir", type=str, default="way_back", help="folder inside the results folder for the runs of the way back (their configs and logs too)")
parser.add_argument("--inter_str", type=number_str, default="10", help="strength of the bead-membrane LJ interaction (Tube_eq.txt: 10)")
parser.add_argument("--shift", type=number_str, default=None, help="shift of the LJ interaction (default: equal to --inter_str, which makes the energy continuous at the cutoff; Tube_eq.txt: 10 and 10)")
parser.add_argument("--init_file", type=str, default=None, help="starting membrane (default ../../../input/Big_planar_mem.obj, the disk of radius 10 in the plane x = 0)")
parser.add_argument("--start_run", type=str, default=None,
                    help="skip the pull: a finished run (folder with Final_state.obj, or a .obj file) whose mesh has the bead at x_start, e.g. a frame of your own pulling runs; the chain starts from that mesh")
# --- the pull that makes the start state
parser.add_argument("--pull_force", type=number_str, default="200.0", help="force of the spring that pulls the bead in the start run (the Lineal bond of Tube_eq.txt: 200)")
parser.add_argument("--anchor_x", type=number_str, default="6.0", help="x of the frozen anchor of the spring in the start run (must be well away from x_start)")
parser.add_argument("--start_timesteps", type=number_str, default="6000",
                    help="steps of the start run (a cap: with the default mesh the bead reaches 0.2 after 250-600 steps and freezes, the rest relaxes the frozen state; much longer coarsens the dimple; "
                         "a finer mesh (smaller --refine_angle or --size_max) moves the bead more slowly, 1500-4000 steps). The chain stops with a message if the bead did not freeze")
# --- mesh (the remesher of Tube_eq.txt)
parser.add_argument("--size_min", type=number_str, default=None, help="smallest edge length of the remesher (template 0.001)")
parser.add_argument("--size_max", type=number_str, default=None, help="largest edge length of the remesher (template 0.3)")
parser.add_argument("--refine_angle", type=number_str, default="0.25",
                    help="remesher refine_angle of the pull and of the chain (default 0.25; Tube_eq.txt has 0.5). The vertex-only LJ needs more than the 4-8 vertices that 0.5 leaves under a bead of radius 0.3: "
                         "with 0.5 the first relaxation of a frozen bead stalls, with 0.25 a chain ran smoothly (force monotonic in x, bead 1.0 radius from the membrane), 0.15 moves the pulled bead too slowly")
parser.add_argument("--aspect_min", type=number_str, default=None, help="remesher aspect_min (template 0.3)")
parser.add_argument("--remesh_max_rise", type=number_str, default="1e9",
                    help="remesh_rollback max_rise of the frozen-bead runs: a remesh that raises a term of the energy by more than max_rise times its size is undone (template 10; here 1e9: never). "
                         "A remesh that puts a vertex inside the bead raises the bead energy and must be allowed there, or the run polishes for ever; the runner checks every result (README)")
# --- how long each run is (the numbers of the stopping block of the config; the links copy the cold start unless a --link_ value is given)
parser.add_argument("--switch", type=str, default="100,20000", help="Switch_times of the links (BFGS, BFGS-Normal): two numbers")
parser.add_argument("--timesteps", type=number_str, default=None, help="largest number of steps of the relaxation of the first position (template 200000, the stopping block ends a run long before)")
parser.add_argument("--window", type=number_str, default=None, help="stopping: steps of one window (template 50)")
parser.add_argument("--patience", type=number_str, default=None, help="stopping: windows in a row that must pass the test (template 4)")
parser.add_argument("--bfgs_tol", type=number_str, default=None, help="stopping: relative net energy change of a window below which BFGS hands over to BFGS-Normal (template 5e-05)")
parser.add_argument("--normal_tol", type=number_str, default=None, help="stopping: relative net energy change of a window below which BFGS-Normal ends the run (template 1e-08)")
parser.add_argument("--stall_steps", type=number_str, default=None, help="stopping: line searches of tiny steps in a row that reset or end a run (template 10)")
parser.add_argument("--link_timesteps", type=number_str, default="60000",
                    help="largest number of steps of a link; a link that reaches it is marked capped in chain_summary.csv (not converged) and the chain goes on")
parser.add_argument("--link_normal_tol", type=number_str, default=None, help="normal_tol of the links only")
# --- the job
parser.add_argument("--batch_tag", type=str, required=True, help="tag of the batch (also the results folder name, ../Results/<batch_tag>/)")
parser.add_argument("--unique_tag", type=str, required=True, help="unique tag of this walk (prefix of the config, subjob and output files)")
parser.add_argument("--time", type=str, default=None, help="time limit of the sbatch job (default from the number of relaxations: 3 minutes each, at least 24 h, at most 168 h); if it runs out, "
                                                           "sbatch the same subjob again: the chain resumes after the last finished run")
parser.add_argument("--mem", type=str, default="16G", help="memory of the sbatch job")
args = parser.parse_args()

def finite(name, text, low=None, strictly=False):
    """The number of an option, finite (number_str lets nan, inf, .5 and 5. through, which are not all valid JSON) and above a limit."""
    v = float(text)
    if not math.isfinite(v):
        raise SystemExit("%s must be a finite number, not %s" % (name, text))
    if low is not None and (v <= low if strictly else v < low):
        raise SystemExit("%s must be %s %g, not %s" % (name, "above" if strictly else "at least", low, text))
    return v


radius = finite("--radius", args.radius, 0.0, True)
KA = finite("--KA", args.KA, 0.0)
KB = finite("--KB", args.KB, 0.0, True)
x_start = finite("--x_start", args.x_start)
x_end = finite("--x_end", args.x_end)
step = finite("--step", args.step, 0.0, True)
max_step = finite("--max_step", args.max_step, 0.0)
for name in ("inter_str", "pull_force", "anchor_x", "size_min", "size_max", "refine_angle", "aspect_min", "remesh_max_rise", "bfgs_tol", "normal_tol", "link_normal_tol"):
    if getattr(args, name) is not None:
        finite("--" + name, getattr(args, name), 0.0 if name != "anchor_x" else None)
if args.shift is not None:
    finite("--shift", args.shift)
if (x_end - x_start) / step > 5000:
    raise SystemExit("%d positions: --step %s is too small for the range" % ((x_end - x_start) / step, args.step))
if not step > 0.0:
    raise SystemExit("--step must be positive")
if not x_end > x_start:
    raise SystemExit("--x_end (%g) must be larger than --x_start (%g)" % (x_end, x_start))
if max_step < 0.0:
    raise SystemExit("--max_step must not be negative (0 walks the list as given)")
try:
    switch_times = [int(float(x)) for x in args.switch.split(",")]
except ValueError:
    switch_times = []
if len(switch_times) != 2 or not 0 < switch_times[0] < switch_times[1]:
    raise SystemExit("--switch needs two increasing positive step numbers (BFGS, BFGS-Normal), not %r" % args.switch)
if args.start_run is not None and not os.path.exists(args.start_run):
    raise SystemExit("--start_run %s does not exist" % args.start_run)
x_back = finite("--x_back", args.x_back) if args.x_back is not None else min(0.0, x_start)
if x_back > x_start + 1e-9:
    raise SystemExit("--x_back (%g) must not be above --x_start (%g): the way back ends at or below the start of the way out" % (x_back, x_start))
if x_back < -radius - 1e-9:
    raise SystemExit("--x_back %g is below -radius" % x_back)
if "/" in args.back_dir or args.back_dir in ("", ".", ".."):
    raise SystemExit("--back_dir is a folder name, not a path: %r" % args.back_dir)
if x_start < -radius - 1e-9:
    raise SystemExit("--x_start %g is below -radius: the bead would start inside the flat membrane's other side; use -radius or more" % x_start)
flat_start = x_start <= -radius + 1e-9
if flat_start and args.start_run is not None:
    raise SystemExit("--start_run is for a bead above -radius (the pulled state), not for the flat start")
if not flat_start and args.start_run is None and float(args.anchor_x) < x_start + 1.0:
    raise SystemExit("--anchor_x %s must be at least 1 above x_start %g: the spring pulls the bead towards the anchor, the bead freezes at x_start on the way" % (args.anchor_x, x_start))

# the positions: x_start + i step up to x_end (computed from the index, no accumulated rounding), x_end itself is always the last one
n = int((x_end - x_start) / step + 1e-9)
positions = [round(x_start + i * step, 9) for i in range(n + 1)]
if positions[-1] < x_end - 1e-9:
    positions.append(round(x_end, 9))
# the way back: the list backwards down to positions[0], then x_start - i step down to x_back (x_back itself is always the last one)
back_positions = list(reversed(positions[:-1]))
n_down = int((x_start - x_back) / step + 1e-9)
back_positions += [round(x_start - i * step, 9) for i in range(1, n_down + 1)]
if back_positions and back_positions[-1] > x_back + 1e-9:
    back_positions.append(round(x_back, 9))


def relaxations(path):
    return sum(max(1, int(math.ceil(abs(b - a) / max_step - 1e-9))) if max_step > 0 else 1 for a, b in zip(path[:-1], path[1:]))


n_relax = 1 + relaxations(positions) + (0 if args.no_return else relaxations([positions[-1]] + back_positions))
estimate_h = int(math.ceil(n_relax * 3.0 / 60.0))
time_limit = args.time if args.time is not None else "%d:00:00" % min(168, max(24, estimate_h))
if args.time is None and estimate_h > 168:
    print("warning: about %d relaxations are about %d hours at 3 minutes each, the time limit is clamped at 168 h: sbatch the subjob again when it ends, the chain resumes" % (n_relax, estimate_h))
shift = args.shift if args.shift is not None else args.inter_str

Batch_dir = '../Results/{}/'.format(args.batch_tag)


def create_json():
    os.makedirs("../Config_files/", exist_ok=True)
    env = Environment(loader=FileSystemLoader('../Templates/'))
    num = lambda text: repr(float(text))          # the number as JSON text (what was typed may be .5 or 5. or 1_0)
    mesh = {}
    for key in ("size_min", "size_max", "refine_angle", "aspect_min"):
        if getattr(args, key) is not None:
            mesh[key] = num(getattr(args, key))
    # the frozen bead of the chain (the first position, the start mesh replaces the init file when the chain runs)
    params = dict(first_dir=Batch_dir, KA=num(args.KA), KB=num(args.KB), radius=num(args.radius), x=positions[0], inter_str=num(args.inter_str), shift=num(shift),
                  remesh_max_rise=num(args.remesh_max_rise), **mesh)
    if args.init_file is not None:
        params["init_file"] = args.init_file
    for key, option in (("timesteps", args.timesteps), ("stop_window", args.window), ("stop_patience", args.patience), ("stall_steps", args.stall_steps)):
        if option is not None:
            params[key] = int(float(option))
    for key, option in (("bfgs_tol", args.bfgs_tol), ("normal_tol", args.normal_tol)):
        if option is not None:
            params[key] = num(option)
    data = json.loads(env.get_template('Tube_material.txt').render(**params))
    Config_path = '../Config_files/{0}_ConfigFile.json'.format(args.unique_tag)
    with open(Config_path, 'w') as file:
        json.dump(data, file, indent=4)

    chain = {"config": Config_path, "binary": '../build/bin/main_cluster', "batch_dir": Batch_dir, "positions": positions, "back_positions": back_positions, "back_dir": args.back_dir,
             "max_step": max_step if max_step > 0.0 else None, "return": not args.no_return,
             "switch": switch_times, "nice": 0, "radius": radius, "link_timesteps": int(float(args.link_timesteps))}
    if args.link_normal_tol is not None:
        chain["link_stopping"] = {"normal_tol_E": float(args.link_normal_tol)}
    if args.start_run is not None:
        chain["start_run"] = os.path.abspath(args.start_run)
    elif not flat_start:
        # one step of the manual bead is capped at line_search_bead_step radii: 0.015 in length, below the 0.02 wide ball around FinalPos in which the bead freezes
        start = dict(first_dir=Batch_dir + "start_x%g/" % x_start, KA=num(args.KA), KB=num(args.KB), radius=num(args.radius), x_ini=num(-radius), x_final=num(x_start),
                     inter_str=num(args.inter_str), shift=num(shift), pull_force=num(args.pull_force), anchor_x=num(args.anchor_x), bead_step=num(0.015 / radius),
                     timesteps=int(float(args.start_timesteps)), **mesh)
        if args.init_file is not None:
            start["init_file"] = args.init_file
        sdata = json.loads(env.get_template('Tube_start.txt').render(**start))
        Start_path = '../Config_files/{0}_start.json'.format(args.unique_tag)
        with open(Start_path, 'w') as file:
            json.dump(sdata, file, indent=4)
        chain["start_config"] = Start_path
    Chain_path = '../Config_files/{0}_chain.json'.format(args.unique_tag)
    with open(Chain_path, 'w') as file:
        json.dump(chain, file, indent=4)
    return Config_path, Chain_path


os.makedirs('../Subjobs/', exist_ok=True)
os.makedirs('../Outputs/', exist_ok=True)
if os.path.exists(os.path.join(Batch_dir, "chain_summary.csv")):
    print("note: %schain_summary.csv exists: the job will RESUME that chain after its last finished run (delete the folder to start over; the runner refuses a chain with other positions, "
          "KA, KB, radius or start state)" % Batch_dir)

Config_path, Chain_path = create_json()

Output_name = '{0}_output.output'.format(args.unique_tag)

write_subjob('{0}_subjob'.format(args.unique_tag),
             Output_name,
             'python3 Run_tube_material_chain.py {}'.format(Chain_path),
             job_name='Tube_material', mem=args.mem, copy_output_to=Batch_dir, time=time_limit)

print("chain {}: KA {}, KB {}, bead radius {}, x from {:g} to {:g} in {} positions (step {}), start {}, way back {}, about {} relaxations (waypoints included), time limit {}, results in {}".format(
    args.unique_tag, args.KA, args.KB, args.radius, positions[0], positions[-1], len(positions), args.step,
    "from the run " + args.start_run if args.start_run else ("the flat membrane" if flat_start else "the spring pull to x = %g" % x_start),
    "off" if args.no_return else "to x = %g in %s/ (%d positions)" % (back_positions[-1] if back_positions else positions[0], args.back_dir, len(back_positions)), n_relax, time_limit, Batch_dir))
