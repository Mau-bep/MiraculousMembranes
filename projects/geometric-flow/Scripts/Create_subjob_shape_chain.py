import os
import json
import math
import argparse
from jinja2 import Environment, FileSystemLoader
from cli_args import number_str
from subjob import write_subjob

# One BRANCH of vesicle shapes (prolate-1, prolate-2 or pear) followed in the reduced volume nu at a fixed reduced spontaneous curvature c0, the energy-against-volume curves of
# Seifert, Berndl and Lipowsky (Phys. Rev. A 44, 1182 (1991)), Fig. 11 (c0 = 2.4) and Fig. 13 (c0 = 3). Writes the config every link is built from (../Templates/Shape_chain.txt),
# the chain file (the seed, the first nu and the legs that Run_shape_chain.py walks) and the sbatch subjob. The presets below are the seeds and ranges of the three branches at
# c0 = 2.4 and 3; every number can be overridden, and another c0 needs --seed_capsule or --seed_spheres, --nu0 and --nu_up / --nu_down. See Run_shape_chain.py for what a chain does.

# seed (a capsule of reduced volume nu, or two spheres of radii ratio q and neck radius `neck` times the small radius), nu0 (the first link), the legs (end, step, stop_when)
PRESETS = {
    (2.4, "prolate1"): dict(seed={"capsule": {"nu": 0.70}}, nu0=0.74, up=(0.99, 0.01), down=(0.60, 0.01), up_stop={}, down_stop={}),
    (2.4, "prolate2"): dict(seed={"capsule": {"nu": 0.55}}, nu0=0.60, up=(0.72, 0.01), down=None, up_stop={"skew_above": 0.03}, down_stop={}),
    (2.4, "pear"): dict(seed={"spheres": {"q": 1.3, "neck": 0.2}}, nu0=0.72, up=(0.85, 0.01), down=(0.70, 0.005), up_stop={}, down_stop={"skew_below": 0.02}),
    (3.0, "prolate1"): dict(seed={"capsule": {"nu": 0.70}}, nu0=0.85, up=(0.99, 0.01), down=(0.60, 0.01), up_stop={}, down_stop={}),
    (3.0, "prolate2"): dict(seed={"capsule": {"nu": 0.40}}, nu0=0.50, up=(0.72, 0.005), down=None, up_stop={"skew_above": 0.03}, down_stop={}),
    (3.0, "pear"): dict(seed={"spheres": {"q": 2.0, "neck": 0.3}}, nu0=0.78, up=(0.92, 0.01), down=(0.70, 0.005), up_stop={}, down_stop={"skew_below": 0.02}),
}

parser = argparse.ArgumentParser(description="Writes the config, the chain file and the sbatch subjob of one branch of vesicle shapes followed in the reduced volume at a fixed c0.")
parser.add_argument("--c0", type=number_str, required=True, help="reduced spontaneous curvature c0 = C0 R0 of Seifert et al. (H0 of the config is c0 / (2 R0) with R0 = sqrt(A / 4 pi), it changes from link to link)")
parser.add_argument("--branch", type=str, required=True, choices=["prolate1", "prolate2", "pear"], help="the branch: the preset of the seed and the ranges for c0 = 2.4 and 3 (Figs. 11 and 13)")
parser.add_argument("--nu0", type=number_str, default=None, help="reduced volume of the first link (the seed is relaxed there); preset")
parser.add_argument("--nu_up", type=number_str, default=None, help="end of the leg that goes up in nu from nu0 (the area shrinks: the stable direction of the flow); preset. 'none' drops the leg")
parser.add_argument("--nu_down", type=number_str, default=None, help="end of the leg that goes down in nu from nu0 (the area grows); preset. 'none' drops the leg")
parser.add_argument("--step", type=number_str, default=None, help="nu step of both legs (preset: 0.01, 0.005 for the pear's way down); a step is halved (down to a eighth) at a jump or a failure")
parser.add_argument("--seed_capsule", type=number_str, default=None, help="seed: a capsule (cylinder with hemispherical caps) of this reduced volume")
parser.add_argument("--seed_spheres", type=str, default=None, help="seed: two overlapping spheres, 'q,neck': radii ratio q >= 1 and the radius of the circle of intersection in units of the small radius (q = 1 a dumbbell, q > 1 a pear)")
parser.add_argument("--seed_mesh", type=str, default=None, help="seed: a closed triangle mesh (scaled to the volume --V0)")
parser.add_argument("--binary", type=str, default="../build/bin/main_cluster", help="the main_cluster that runs the links (a copy of a build keeps a chain on one binary while the source changes)")
parser.add_argument("--V0", type=number_str, default="0.5035", help="volume of every link (0.5035 is the volume of input/Asym_cil_regular.obj; the remesher lengths below are for this size, nu and c0 do not depend on it)")
parser.add_argument("--KB", type=number_str, default="1.0", help="bending constant (the energies of the curves are F_b / (8 pi kappa) = E_bend / (4 pi KB), independent of it)")
parser.add_argument("--dA", type=number_str, default="0.0004", help="the target area moves this much per step towards the area of the target nu (a large jump of the target makes the flow unstable)")
# --- mesh (the remesher) and how long a link is
parser.add_argument("--size_min", type=number_str, default=None, help="smallest edge length of the remesher (template 0.001)")
parser.add_argument("--size_max", type=number_str, default=None, help="largest edge length of the remesher (template 0.15)")
parser.add_argument("--refine_angle", type=number_str, default=None, help="remesher refine_angle: an edge turns the surface normal by about this angle, so the edge length is refine_angle / (the largest principal curvature) between size_min and size_max (template 0.2)")
parser.add_argument("--aspect_min", type=number_str, default=None, help="remesher aspect_min (template 0.2)")
parser.add_argument("--timesteps", type=number_str, default="60000", help="largest number of steps of a link (the stopping block ends a link long before on a smooth branch; a link that reaches the cap is marked capped in chain_summary.csv)")
parser.add_argument("--switch", type=str, default="1000,4000", help="Switch_times of a link (BFGS, BFGS-Normal): two numbers; BFGS is delayed while the area target is still moving")
parser.add_argument("--normal_tol", type=number_str, default=None, help="stopping: relative net energy change of a window below which BFGS-Normal ends a link (template 1e-8)")
parser.add_argument("--nu_tol", type=number_str, default="0.01", help="a link whose reduced volume is further than this from its target is a failure")
parser.add_argument("--jump_neck", type=number_str, default="0.15", help="a change of the neck ratio between two links larger than this is a jump of the branch")
parser.add_argument("--jump_lobes", type=number_str, default="0.08", help="a change of the lobe ratio larger than this is a jump")
parser.add_argument("--jump_skew", type=number_str, default="0.08", help="a change of the skewness larger than this is a jump")
# --- the job
parser.add_argument("--batch_tag", type=str, required=True, help="tag of the batch (also the results folder name, ../Results/<batch_tag>/)")
parser.add_argument("--unique_tag", type=str, required=True, help="unique tag of this chain (prefix of the config, subjob and output files)")
parser.add_argument("--time", type=str, default=None, help="time limit of the sbatch job (default from the number of links: 3 minutes each, at least 10 h, at most 168 h); if it runs out, sbatch the same subjob again: the chain resumes")
parser.add_argument("--mem", type=str, default="8G", help="memory of the sbatch job")
parser.add_argument("--cpus", type=int, default=1, help="cpus of the job (a link is a mesh of a few thousand vertices: one is enough)")
args = parser.parse_args()


def number(name, text, low=None, high=None):
    v = float(text)
    if not math.isfinite(v) or (low is not None and v < low) or (high is not None and v > high):
        raise SystemExit("%s = %s is out of range" % (name, text))
    return v


c0 = number("--c0", args.c0, 0.0)
preset = next((p for (c, b), p in PRESETS.items() if b == args.branch and abs(c - c0) < 1e-9), None)
custom = args.seed_capsule is not None or args.seed_spheres is not None or args.seed_mesh is not None
if preset is None and not (custom and args.nu0 is not None):
    raise SystemExit("no preset for c0 = %g and the branch %s (presets: %s): give a seed (--seed_capsule, --seed_spheres or --seed_mesh) and --nu0 with --nu_up and/or --nu_down"
                     % (c0, args.branch, ", ".join("%g/%s" % k for k in PRESETS)))
preset = preset or dict(seed=None, nu0=None, up=None, down=None, up_stop={}, down_stop={})

if args.seed_capsule is not None:
    seed = {"capsule": {"nu": number("--seed_capsule", args.seed_capsule, 0.05, 0.999)}}
elif args.seed_spheres is not None:
    q, neck = (float(x) for x in args.seed_spheres.split(","))
    seed = {"spheres": {"q": q, "neck": neck}}
elif args.seed_mesh is not None:
    if not os.path.exists(args.seed_mesh):
        raise SystemExit("--seed_mesh %s does not exist" % args.seed_mesh)
    seed = {"mesh": os.path.abspath(args.seed_mesh)}
else:
    seed = preset["seed"]
nu0 = number("--nu0", args.nu0, 0.05, 1.0) if args.nu0 is not None else preset["nu0"]


def leg_end(option, preset_leg):
    if option is not None and option.lower() == "none":
        return None
    if option is not None:
        return number("nu end", option, 0.05, 1.0)
    return preset_leg[0] if preset_leg else None


up_end, down_end = leg_end(args.nu_up, preset["up"]), leg_end(args.nu_down, preset["down"])
step_up = number("--step", args.step, 1e-4, 0.2) if args.step is not None else (preset["up"][1] if preset["up"] else 0.01)
step_down = number("--step", args.step, 1e-4, 0.2) if args.step is not None else (preset["down"][1] if preset["down"] else 0.01)
legs = []
if up_end is not None:
    if not up_end > nu0:
        raise SystemExit("--nu_up (%g) must be above --nu0 (%g)" % (up_end, nu0))
    legs.append({"name": "up", "nu_to": up_end, "step": step_up, **({"stop_when": preset["up_stop"]} if preset["up_stop"] else {})})
if down_end is not None:
    if not down_end < nu0:
        raise SystemExit("--nu_down (%g) must be below --nu0 (%g)" % (down_end, nu0))
    legs.append({"name": "down", "nu_to": down_end, "step": step_down, **({"stop_when": preset["down_stop"]} if preset["down_stop"] else {})})
if not legs:
    raise SystemExit("no leg: give --nu_up and/or --nu_down")
try:
    switch_times = [int(float(x)) for x in args.switch.split(",")]
except ValueError:
    switch_times = []
if len(switch_times) != 2 or not 0 < switch_times[0] < switch_times[1]:
    raise SystemExit("--switch needs two increasing positive step numbers (BFGS, BFGS-Normal), not %r" % args.switch)

Batch_dir = '../Results/{}/'.format(args.batch_tag)
n_links = 1 + sum(int(math.ceil(abs(l["nu_to"] - nu0) / l["step"])) for l in legs)
estimate_h = int(math.ceil(n_links * 1.5 * 3.0 / 60.0))
time_limit = args.time if args.time is not None else "%d:00:00" % min(168, max(10, estimate_h))


def create_json():
    os.makedirs("../Config_files/", exist_ok=True)
    env = Environment(loader=FileSystemLoader('../Templates/'))
    num = lambda text: repr(float(text))
    params = dict(first_dir=Batch_dir, KB=num(args.KB), dA=num(args.dA), timesteps=int(float(args.timesteps)), switch_bfgs=switch_times[0], switch_normal=switch_times[1])
    for key in ("size_min", "size_max", "refine_angle", "aspect_min"):
        if getattr(args, key) is not None:
            params[key] = num(getattr(args, key))
    if args.normal_tol is not None:
        params["normal_tol"] = num(args.normal_tol)
    data = json.loads(env.get_template('Shape_chain.txt').render(**params))
    Config_path = '../Config_files/{0}_ConfigFile.json'.format(args.unique_tag)
    with open(Config_path, 'w') as file:
        json.dump(data, file, indent=4)
    chain = {"config": Config_path, "binary": args.binary, "batch_dir": Batch_dir, "c0": c0, "KB": float(args.KB), "V0": float(args.V0), "dA": float(args.dA),
             "seed": seed, "nu0": nu0, "legs": legs, "switch": switch_times, "timesteps": int(float(args.timesteps)), "threads": args.cpus, "nice": 0,
             "nu_tol": float(args.nu_tol), "jump": {"neck_ratio": float(args.jump_neck), "lobe_ratio": float(args.jump_lobes), "skew": float(args.jump_skew)}}
    Chain_path = '../Config_files/{0}_chain.json'.format(args.unique_tag)
    with open(Chain_path, 'w') as file:
        json.dump(chain, file, indent=4)
    return Config_path, Chain_path


os.makedirs('../Subjobs/', exist_ok=True)
os.makedirs('../Outputs/', exist_ok=True)
if os.path.exists(os.path.join(Batch_dir, "chain_summary.csv")):
    print("note: %schain_summary.csv exists: the job will RESUME that chain from its finished links (delete the folder to start over)" % Batch_dir)

Config_path, Chain_path = create_json()

Output_name = '{0}_output.output'.format(args.unique_tag)

write_subjob('{0}_subjob'.format(args.unique_tag),
             Output_name,
             'python3 Run_shape_chain.py {}'.format(Chain_path),
             job_name='Shape_chain', mem=args.mem, cpus=args.cpus, copy_output_to=Batch_dir, time=time_limit)

print("chain {}: c0 {}, branch {}, seed {}, nu0 {:g}, legs {}, about {} links, time limit {}, results in {} (locally: python3 Run_shape_chain.py {})".format(
    args.unique_tag, args.c0, args.branch, json.dumps(seed), nu0, ", ".join("%s to %g (step %g)" % (l["name"], l["nu_to"], l["step"]) for l in legs), n_links, time_limit, Batch_dir, Chain_path))
