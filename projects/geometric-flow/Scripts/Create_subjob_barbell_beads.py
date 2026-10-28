import os
import json
import math
import argparse
from jinja2 import Environment, FileSystemLoader
from cli_args import number_str
from subjob import write_subjob
from Create_barbell_obj import barbell_mesh, check, write_obj

# Barbell with beads of radius 1: two frozen beads inside a membrane that has bending and surface tension but no volume constraint, so the membrane is
# pressed on the beads and a tube is left between them (the old pipeline, Create_subjob_barbell.py, had beads of radius 0.2 and a start mesh that was already a thin tube).
# Writes the config (from ../Templates/Barbell_beads.txt), the start mesh if needed, and the sbatch subjob of ONE run (KA, KB). Loop over KA and KB in Mem_barbell_beads_up.sh.
# The tube of the relaxed barbell has radius R_t = sqrt(KB / (4 KA)) (KA = the Surface_tension constant, kappa = KB/2); measured 0.349 / 0.224 for KA 4 / 10 at KB 2.
#
# --init mesh   (default) start from a barbell mesh of revolution (Create_barbell_obj.py): caps around the beads, a tube of the expected radius between them. Relaxes in minutes.
# --init shrink start from a sphere of radius R0 around both beads, placed at their final positions: the surface tension shrinks the membrane onto the beads and pulls the
#               tube out; slower (tens of thousands of steps) and the tube is only as good as the run gets, but it is the "no volume, let it shrink" protocol.

parser = argparse.ArgumentParser(description="Writes the config, the start mesh and the sbatch subjob of one barbell run.")
parser.add_argument("--KA", type=number_str, required=True, help="surface tension constant (the Surface_tension energy), the tube radius is sqrt(KB/(4 KA))")
parser.add_argument("--KB", type=number_str, required=True, help="bending constant KB = 2 kappa")
parser.add_argument("--Nsim", type=number_str, default="1", help="simulation number (part of the file names, e.g. a replica)")
parser.add_argument("--d", type=number_str, default="10.0", help="distance between the bead centres (the beads sit at x = +-d/2), default 10 = beads at -5 and 5")
parser.add_argument("--radius", type=number_str, default="1.0", help="bead radius (fixed at 1 in the study)")
parser.add_argument("--inter_str", type=number_str, default="10.0", help="strength of the repulsive core (Shifted-LJ) between beads and membrane")
parser.add_argument("--init", choices=["mesh", "shrink"], default="mesh", help="start mesh: mesh = barbell of revolution (default), shrink = sphere around both beads that shrinks")
parser.add_argument("--R0", type=number_str, default=None, help="--init shrink: radius of the start sphere (default d/2 + 1.5 radii)")
parser.add_argument("--tube", type=number_str, default=None, help="--init mesh: radius of the tube of the start mesh (default sqrt(KB/(4 KA)), at most one bead radius)")
parser.add_argument("--timesteps", type=int, default=None, help="step limit (default 100000 for mesh, 200000 for shrink; the run stops earlier when it converges)")
parser.add_argument("--size_min", type=number_str, default=None, help="smallest edge length of the remesher (default min(0.02, R_t/10))")
parser.add_argument("--size_max", type=number_str, default=None, help="largest edge length of the remesher (template 0.5)")
parser.add_argument("--refine_angle", type=number_str, default=None, help="remesher refine_angle (template 0.3)")
parser.add_argument("--batch_tag", type=str, required=True, help="tag of the batch (also the results folder name)")
parser.add_argument("--unique_tag", type=str, default=None, help="prefix of the config, subjob and output files (default Barbell_beads_KA_<KA>_KB_<KB>_Nsim_<Nsim>)")
parser.add_argument("--time", type=str, default=None, help="time limit of the sbatch job (default 12:00:00 for mesh, 48:00:00 for shrink)")
args = parser.parse_args()

a = float(args.radius)
d = float(args.d)
KA = float(args.KA)
KB = float(args.KB)
if d < 2.0 * a + 0.5:
    raise SystemExit("--d %g: the beads (radius %g) are too close, d must be at least %g" % (d, a, 2.0 * a + 0.5))
if args.unique_tag is None:
    args.unique_tag = "Barbell_beads_KA_{}_KB_{}_Nsim_{}".format(args.KA, args.KB, args.Nsim)
Rt = math.sqrt(KB / (4.0 * KA))   # radius of the free tube
Batch_dir = '../Results/{}/'.format(args.batch_tag)
timesteps = args.timesteps if args.timesteps is not None else (100000 if args.init == "mesh" else 200000)
time_limit = args.time if args.time is not None else ("12:00:00" if args.init == "mesh" else "48:00:00")

os.makedirs("../Config_files/", exist_ok=True)
os.makedirs('../Subjobs/', exist_ok=True)
os.makedirs('../Outputs/', exist_ok=True)

params = dict(first_dir=Batch_dir, dist=d, radius=a, KB=args.KB, sigma=args.KA, inter_str=args.inter_str, timesteps=timesteps,
              beads=[{"pos": [d / 2.0, 0.0, 0.0]}, {"pos": [-d / 2.0, 0.0, 0.0]}])
params["size_min"] = args.size_min if args.size_min is not None else round(min(0.02, Rt / 10.0), 4)
for key in ("size_max", "refine_angle"):
    if getattr(args, key) is not None:
        params[key] = getattr(args, key)

if args.init == "mesh":
    tube = float(args.tube) if args.tube is not None else min(Rt, a)
    tube = max(tube, 0.1)
    Mesh_path = '../Config_files/{}_init.obj'.format(args.unique_tag)
    V, F = barbell_mesh(d, a=a, gap=0.15 * a, rt=tube, edge=min(0.15, max(0.6 * tube, 0.05)))
    bad, chi, vol, area = check(V, F)
    if bad or chi != 2 or vol <= 0:
        raise SystemExit("the start mesh is not a closed outward-oriented sphere")
    write_obj(Mesh_path, V, F)
    params.update(init_file=Mesh_path, R0=1.0)
else:
    R0 = float(args.R0) if args.R0 is not None else d / 2.0 + 1.5 * a
    if R0 < d / 2.0 + 1.2 * a:
        raise SystemExit("--R0 %g leaves no room for the beads: it must be larger than d/2 + 1.2 a = %g" % (R0, d / 2.0 + 1.2 * a))
    params.update(R0=R0)   # Simplesphere.obj (radius 1) scaled by R0, the default init_file of the template

env = Environment(loader=FileSystemLoader('../Templates/'))
template = env.get_template('Barbell_beads.txt')
data = json.loads(template.render(**params))
Config_path = '../Config_files/{}_ConfigFile.json'.format(args.unique_tag)
with open(Config_path, 'w') as file:
    json.dump(data, file, indent=4)

print("{}: init {}, beads at +-{:g}, expected tube radius {:.3f}".format(args.unique_tag, args.init, d / 2.0, Rt))

write_subjob('{}_subjob'.format(args.unique_tag),
             '{}_output.output'.format(args.unique_tag),
             '../build/bin/main_cluster {}'.format(Config_path),
             job_name=args.batch_tag, mem='4G', copy_output_to=Batch_dir, time=time_limit)
