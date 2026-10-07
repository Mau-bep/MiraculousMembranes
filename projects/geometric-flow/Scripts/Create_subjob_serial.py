import sys 
import os
import argparse
from cli_args import number_str
from subjob import write_subjob
#I want to create a script that ask for the step size and creates a subjob that uses it

#Lets assume we are on the directory where i can store all the data 








parser = argparse.ArgumentParser(description="Writes the config (if any) and the sbatch subjob for one run.")
parser.add_argument("--v", type=float, required=True, help="v")
parser.add_argument("--KB", type=number_str, required=True, help="bending constant KB")
parser.add_argument("--init_cond", type=int, required=True, help="initial condition")
parser.add_argument("--Nsim", type=int, required=True, help="simulation number (part of the file names)")
args = parser.parse_args()

v = args.v
KB = args.KB
Init_cond = args.init_cond

Nsim = args.Nsim



os.makedirs('../Subjobs/',exist_ok=True)
os.makedirs('../Outputs/',exist_ok=True)


write_subjob('subjob_serial_correct_v_{}_KB_{}_init_cond_{}_Nsim_{}'.format(v,KB,Init_cond,Nsim),
             'output_Mem3DG_v_{}_KB_{}_evol_Init_cond_{}_Nsim_{}'.format(v,KB,Init_cond,Nsim),
             '../build/bin/main_cluster {} {} {} {}'.format(v,Init_cond,Nsim,KB),
             time='30:00:00', mem='16G')
