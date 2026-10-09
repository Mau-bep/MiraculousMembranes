import os
import json
import argparse
from jinja2 import Environment, FileSystemLoader
from cli_args import number_str
from subjob import write_subjob

# One bead wrapping on a planar membrane (Big_planar_mem.obj, a disk in the yz plane, or Planar_full.obj).
# Writes the config from ../Templates/Wrapping_planar.txt and the sbatch subjob for one run.

parser = argparse.ArgumentParser(description="Writes the config and the sbatch subjob for one run.")
parser.add_argument("--inter_str", type=number_str, required=True, help="bead-membrane adhesion strength")
parser.add_argument("--radius", type=number_str, required=True, help="bead radius")
parser.add_argument("--KA", type=number_str, required=True, help="surface tension constant KA")
parser.add_argument("--KB", type=number_str, required=True, help="bending constant KB")
parser.add_argument("--init", choices=["flat", "full"], default="flat",
                    help="flat: Big_planar_mem.obj with the bead above the plane, "
                         "full: Planar_full.obj, the membrane already wrapped around the bead")
parser.add_argument("--shell_width", type=number_str, default=None,
                    help="half width of the adhesion shell as a fraction of the bead radius (default of the code: 0.25); "
                         "a narrower shell needs a finer mesh")
parser.add_argument("--shell_power", type=number_str, default=None,
                    help="exponent p >= 1 of the adhesion shell, weight = [1/2 (1 + cos(pi x))]^p (default of the code: 1; 2, 4 = cos^4, cos^8 of pi x/2)")
parser.add_argument("--size_min", type=number_str, default=None, help="smallest edge length of the remesher (template default)")
parser.add_argument("--size_max", type=number_str, default=None, help="largest edge length of the remesher (template default)")
parser.add_argument("--rescale", type=number_str, default=None,
                    help="rescale factor of the mesh about its centre of mass (template default 3.0; the old runs used 1.0: "
                         "Big_planar_mem.obj is already a disk of radius 10 and Planar_full.obj is wrapped around a bead of radius 1, "
                         "the bead itself is not rescaled)")
parser.add_argument("--tension", choices=["surface", "excess"], default="surface",
                    help="surface: KA * area, excess: KA * (area - initial area) only while the area is above the initial one")
parser.add_argument("--batch_tag", type=str, required=True, help="tag of the batch (also the results folder name)")
parser.add_argument("--unique_tag", type=str, required=True, help="unique tag of this run (prefix of the config, subjob and output files)")
args = parser.parse_args()

Strength = args.inter_str
radius = args.radius
KA = args.KA
KB = args.KB
Excess = args.tension == "excess"
Init = args.init
# Optional template overrides, left out of the render (so the template defaults apply) when not given
Overrides = {k: v for k, v in (("shell_width", args.shell_width), ("shell_power", args.shell_power), ("rescale", args.rescale), ("size_min", args.size_min), ("size_max", args.size_max)) if v is not None}
Batch_tag = args.batch_tag
Unique_tag = args.unique_tag

# Planar_full.obj is a wrapped state: the bead (radius 1) sits at this position in it
FULL_BEAD_POS = (-1.82316, 0.0744237, 0.00342679)

FLAT_MESH = '../../../input/Big_planar_mem.obj'
FULL_MESH = '../../../input/Planar_full.obj'


def obj_area(path):
    """Total area of the triangles of an obj file, path relative to the Scripts folder like the init_file."""
    here = os.path.dirname(os.path.abspath(__file__))
    vertices, area = [], 0.0
    with open(os.path.join(here, path)) as f:
        for line in f:
            t = line.split()
            if not t:
                continue
            if t[0] == 'v':
                vertices.append([float(x) for x in t[1:4]])
            elif t[0] == 'f':
                a, b, c = (vertices[int(x.split('/')[0]) - 1] for x in t[1:4])
                u = [b[i] - a[i] for i in range(3)]
                v = [c[i] - a[i] for i in range(3)]
                cross = [u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]]
                area += 0.5 * sum(x * x for x in cross) ** 0.5
    return area


# Each batch gets its own results folder
Batch_dir = '../Results/{}/'.format(Batch_tag)


def Create_json_wrapping_planar(ka, kb, r, inter_str):
    os.makedirs("../Config_files/", exist_ok=True)
    env = Environment(loader=FileSystemLoader('../Templates/'))

    template = env.get_template('Wrapping_planar.txt')
    if Init == "full":
        if float(r) != 1.0:
            print("Warning: Planar_full.obj was made with a bead of radius 1, not {}".format(r))
        init_file = FULL_MESH
        xpos, ypos, zpos = FULL_BEAD_POS
        # Excess_tension counts the area above the flat disk's, not above the wrapped mesh's:
        # the target area is this fraction of the area of the initial mesh
        area_factor = obj_area(FLAT_MESH) / obj_area(FULL_MESH)
    else:
        # The adhesion is soft, so the bead starts with its centre one radius above the plane (x = 0)
        init_file = FLAT_MESH
        xpos, ypos, zpos = float(r), 0.0, 0.0
        area_factor = 1.0
    output_from_parsed_template = template.render(KA=ka, KB=kb, radius=r, init_file=init_file, xpos=xpos, ypos=ypos, zpos=zpos,
                                                  interaction=inter_str, excess=Excess, area_factor=area_factor, **Overrides)

    data = json.loads(output_from_parsed_template)

    data['first_dir'] = Batch_dir

    Config_path = '../Config_files/{0}_ConfigFile.json'.format(Unique_tag)

    sim_path = data['first_dir']

    with open(Config_path, 'w') as file:
        json.dump(data, file, indent=4)

    return Config_path, sim_path


os.makedirs('../Subjobs/', exist_ok=True)
os.makedirs('../Outputs/', exist_ok=True)

Config_path, sim_path = Create_json_wrapping_planar(KA, KB, radius, Strength)


Output_name = '{0}_output.output'.format(Unique_tag)

write_subjob('{0}_subjob'.format(Unique_tag),
             Output_name,
             '../build/bin/main_cluster {}'.format(Config_path),
             job_name=Batch_tag, mem='16G', copy_output_to=sim_path)
