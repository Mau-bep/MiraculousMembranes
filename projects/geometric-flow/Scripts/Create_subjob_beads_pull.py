import sys 
import os
import argparse
from cli_args import number_str
from subjob import write_subjob
#I want to create a script that ask for the step size and creates a subjob that uses it

#Lets assume we are on the directory where i can store all the data 







v=1.0
parser = argparse.ArgumentParser(description="Writes the config (if any) and the sbatch subjob for one run.")
parser.add_argument("--radius", type=float, required=True, help="bead radius")
parser.add_argument("--inter_str", type=number_str, required=True, help="bead-membrane interaction strength")
parser.add_argument("--init_cond", type=number_str, required=True, help="initial condition")
parser.add_argument("--Nsim", type=number_str, required=True, help="simulation number (part of the file names)")
parser.add_argument("--KB", type=number_str, required=True, help="bending constant KB")
parser.add_argument("--KA", type=number_str, required=True, help="surface tension constant KA")
args = parser.parse_args()

radius = args.radius
Strength = args.inter_str
Init_cond = args.init_cond
Curv = 0.0
Nsim = args.Nsim
KB = args.KB
KA = args.KA
os.makedirs('../Subjobs/',exist_ok=True)
os.makedirs('../Outputs/',exist_ok=True)


write_subjob('subjob_serial_bead_pull_radius_{}_KA_{}_Strg_{}_init_cond_{}_Nsim_{}_KB_{}'.format(radius,KA,Strength,Init_cond,Nsim,KB),
             'output_serial_bead_pull_radius_{}_KA_{}_Strg_{}_Init_cond_{}_Nsim_{}_KB_{}'.format(radius,KA,Strength,Init_cond,Nsim,KB),
             '../build/bin/main_cluster_pulling_beads {} {} {} {} {} {} {}'.format(Curv,Strength,Init_cond,Nsim,KA,radius,KB),
             mem='2G', cpus=2)
