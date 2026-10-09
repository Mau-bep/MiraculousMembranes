import os
import json
import argparse
from jinja2 import Environment, FileSystemLoader
from cli_args import number_str
from subjob import write_subjob

# One bead wrapping on a planar membrane (Big_planar_mem.obj, a disk in the yz plane).
# Writes the config from ../Templates/Wrapping_planar.txt and the sbatch subjob for one run.

parser = argparse.ArgumentParser(description="Writes the config and the sbatch subjob for one run.")
parser.add_argument("--inter_str", type=number_str, required=True, help="bead-membrane adhesion strength")
parser.add_argument("--radius", type=number_str, required=True, help="bead radius")
parser.add_argument("--KA", type=number_str, required=True, help="surface tension constant KA")
parser.add_argument("--KB", type=number_str, required=True, help="bending constant KB")
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
Batch_tag = args.batch_tag
Unique_tag = args.unique_tag

# Each batch gets its own results folder
Batch_dir = '../Results/{}/'.format(Batch_tag)


def Create_json_wrapping_planar(ka, kb, r, inter_str):
    os.makedirs("../Config_files/", exist_ok=True)
    env = Environment(loader=FileSystemLoader('../Templates/'))

    template = env.get_template('Wrapping_planar.txt')
    # The adhesion is soft, so the bead starts with its centre one radius above the plane (x = 0)
    output_from_parsed_template = template.render(KA=ka, KB=kb, radius=r, xpos=float(r), interaction=inter_str, excess=Excess)

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
