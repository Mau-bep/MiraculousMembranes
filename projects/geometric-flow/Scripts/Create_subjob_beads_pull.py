import sys 
import os
import argparse
from cli_args import number_str
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


f=open('../Subjobs/subjob_serial_bead_pull_radius_{}_KA_{}_Strg_{}_init_cond_{}_Nsim_{}_KB_{}'.format(radius,KA,Strength,Init_cond,Nsim,KB),'w')

f.write('#!/bin/bash \n')
f.write('# \n')

f.write('#SBATCH --job-name=Mem3DGpa\n')
f.write('#SBATCH --output=../Outputs/output_serial_bead_pull_radius_{}_KA_{}_Strg_{}_Init_cond_{}_Nsim_{}_KB_{}\n'.format(radius,KA,Strength,Init_cond,Nsim,KB))
f.write('#\n')
f.write('#number of CPUs to be used\n')
f.write('#SBATCH --ntasks=1\n')
f.write('#Define the number of hours the job should run. \n')
f.write('#Maximum runtime is limited to 10 days, ie. 240 hours\n')
f.write('#SBATCH --time=10:00:00\n')

f.write('#\n')
f.write('#Define the amount of system RAM used by your job in GigaBytes\n')
f.write('#SBATCH --mem=2G\n')
f.write('#\n')

#f.write('#Send emails when a job starts, it is finished or it exits\n')
#f.write('#SBATCH --mail-user=mrojasve@ist.ac.at\n')
#f.write('#SBATCH --mail-type=ALL\n')
#f.write('#\n')

f.write('#Pick whether you prefer requeue or not. If you use the --requeue\n')
f.write('#option, the requeued job script will start from the beginning, \n')
f.write('#potentially overwriting your previous progress, so be careful.\n')
f.write('#For some people the --requeue option might be desired if their\n')

f.write('#application will continue from the last state.\n')
f.write('#Do not requeue the job in the case it fails.\n')
f.write('#SBATCH --no-requeue\n')
f.write('#\n')

f.write('#Define the "gpu" partition for GPU-accelerated jobs\n')
f.write('#####SBATCH --partition=gpu\n')
f.write('#Define the number of GPUs used by your job\n')
f.write('#######SBATCH --gres=gpu:1\n')
f.write('#Define the GPU architecture (GTX980 in the example, other options are GTX1080Ti, K40)\n')
f.write('########SBATCH --constraint=GTX980\n')

f.write('\n')
f.write('#Do not export the local environment to the compute nodes\n')
f.write('#SBATCH --export=NONE\n')
f.write('\n')

f.write('unset SLURM_EXPORT_ENV\n')
f.write('#for single-CPU jobs make sure that they use a single thread\n')
f.write('export OMP_NUM_THREADS=2\n')
f.write('#SBATCH --nodes=1\n')
f.write('#SBATCH --cpus-per-task=2\n')

f.write('\n')

f.write('#load an CUDA software module\n')
f.write('#module load cuda/11.1.0\n')
f.write('#export XLA_FLAGS=--xla_gpu_cuda_data_dir=/usr/lib/cuda\n')
f.write('#print out the list of GPUs before the job is started\n')
f.write('#srun /usr/bin/nvidia-smi\n')
f.write("#run your CUDA binary through SLURM's srun\n")
f.write("#scontrol -o show nodes | awk '{ print $1, $6, $4, $15, $16}' | sort -n | grep gpu62'\n")
# f.write(" #printf ' \n' \n")

# f.write('source /nfs/scistore16/wojtgrp/mrojasve/.bashrc\n')
f.write('export PATH="/nfs/scistore16/wojtgrp/mrojasve/.local/bin:$PATH"\n')
f.write('echo $PATH\n')

# f.write('source ~/anaconda3/etc/profile.d/conda.sh\n')
# f.write('conda activate: jax_cpu\n')


# f.write("#printf ' \n' \n")
# f.write("#printf '==========================================================================\n'\n")
# f.write("#printf ' \n'\n")



f.write('pwd\n')

# f.write('srun time -v ../build/bin/main_cluster_pulling {} {} {} {} {} {}\n'.format(rc,Strength,Init_cond,Nsim,KB,KA))
f.write('srun time -v ../build/bin/main_cluster_pulling_beads {} {} {} {} {} {} {}\n'.format(Curv,Strength,Init_cond,Nsim,KA,radius,KB))



f.write('\n')
f.write('#sacct --format="JobID, State, AllocGRES, AllocNodes, CPUTime, ReqMem, MaxRSS, AveRSS, Elapsed" --units=G | head -n 1\n')
f.write('#sacct --format="JobID, State, AllocGRES, AllocNodes, CPUTime, ReqMem, MaxRSS, AveRSS, Elapsed" --units=G | tail -n 1\n')
f.write('\n')

f.close()
