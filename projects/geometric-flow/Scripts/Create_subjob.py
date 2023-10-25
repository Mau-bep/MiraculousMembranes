import sys 
import os
#I want to create a script that ask for the step size and creates a subjob that uses it

#Lets assume we are on the directory where i can store all the data 







# Nsim=int(sys.argv[1])
# ini_config=int(sys.argv[2])
# fin_config=int(sys.argv[3])
# Target_val=float(sys.argv[4])

v=float(sys.argv[1])
c0=float(sys.argv[2])
KA=float(sys.argv[3])
KB=float(sys.argv[4])



os.makedirs('../Subjobs/',exist_ok=True)
os.makedirs('../Outputs/',exist_ok=True)


f=open('../Subjobs/subjob_parallel_memshape_v_{}_c0_{}_KA_{}_KB_{}'.format(v,c0,KA,KB),'w')

f.write('#!/bin/bash \n')
f.write('# \n')

f.write('#SBATCH --job-name=Mem3DGpa\n')
f.write('#SBATCH --output=../Outputs/output_parallel_Mem3DG_v_{}_c0_{}_KA_{}_KB_{}\n'.format(v,c0,KA,KB))
f.write('#\n')
f.write('#number of CPUs to be used\n')
f.write('#SBATCH --ntasks=1\n')
f.write('#Define the number of hours the job should run. \n')
f.write('#Maximum runtime is limited to 10 days, ie. 240 hours\n')
f.write('#SBATCH --time=60:00:00\n')

f.write('#\n')
f.write('#Define the amount of system RAM used by your job in GigaBytes\n')
f.write('#SBATCH --mem=16G\n')
f.write('#\n')

f.write('#Send emails when a job starts, it is finished or it exits\n')
f.write('#SBATCH --mail-user=mrojasve@ist.ac.at\n')
f.write('#SBATCH --mail-type=ALL\n')
f.write('#\n')

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


f.write('srun time -v ../build/bin/main_cluster_parallel {} {} {} {}\n'.format(v,c0,KA,KB))



f.write('\n')
f.write('#sacct --format="JobID, State, AllocGRES, AllocNodes, CPUTime, ReqMem, MaxRSS, AveRSS, Elapsed" --units=G | head -n 1\n')
f.write('#sacct --format="JobID, State, AllocGRES, AllocNodes, CPUTime, ReqMem, MaxRSS, AveRSS, Elapsed" --units=G | tail -n 1\n')
f.write('\n')

f.close()
