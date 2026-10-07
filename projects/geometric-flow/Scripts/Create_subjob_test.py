import sys 
import os
import argparse
from subjob import write_subjob
#I want to create a script that ask for the step size and creates a subjob that uses it

#Lets assume we are on the directory where i can store all the data 








parser = argparse.ArgumentParser(description="Writes the config (if any) and the sbatch subjob for one run.")
parser.add_argument("--v", type=float, required=True, help="v")
parser.add_argument("--c0", type=float, required=True, help="c0")
parser.add_argument("--KA", type=float, required=True, help="surface tension constant KA")
parser.add_argument("--KB", type=float, required=True, help="bending constant KB")
args = parser.parse_args()

v = args.v
c0 = args.c0
KA = args.KA
KB = args.KB



os.makedirs('../Subjobs/',exist_ok=True)
os.makedirs('../Outputs/',exist_ok=True)


write_subjob('subjob_test_th_memshape_v_{}_c0_{}_KA_{}_KB_{}'.format(v,c0,KA,KB),
             'output_test_th_Mem3DG_v_{}_c0_{}_KA_{}_KB_{}'.format(v,c0,KA,KB),
             '../build/bin/Grad_tests_th {} {} {} {}'.format(v,c0,KA,KB),
             job_name='TestMem', time='40:00:00', mem='16G', mail_user='mrojasve@ist.ac.at')
