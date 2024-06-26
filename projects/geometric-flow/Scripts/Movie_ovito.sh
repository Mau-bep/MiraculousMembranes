#!/bin/bash
#$ -cwd
#$ -j y

# source /mnt/nfs/clustersw/Debian/bullseye/anaconda3/2023.04/activate_anaconda3_2023.04.txt
#module load anaconda3/2023.04
#module load ffmpeg/6.0
#conda init 


#module load python/3.9.7
#python3 -m venv MauEnv 
#source MauEnv/bin/activate
#pip install ffmpeg-python
#pip install ffmpeg
#pip install opencv-python --upgrade --force-reinstall --with-ffmepg
#pip install opencv-python-headless --upgrade --force-reinstall --with-ffmpeg

# module load python/3.9.7
#source /mnt/nfs/clustersw/Debian/bullseye/anaconda3/2023.04/activate_anaconda3_2023.04.txt
#EXPORT PATH=/nfs/scistore16/wojtgrp/mrojasve/.local/bin:/nfs/scistore16/wojtgrp/mrojasve/anaconda3/bin/:$PATH

#source activate /nfs/scistore16/wojtgrp/mrojasve/Cluster_Folders/projects/geometric-flow/build/Mauenv
#conda install ffmpeg
#foreach v ( 0.1 0.3 0.58 0.6 0.65 0.75 0.8 1.0 )
# foreach v ( 0.2 0.3 0.4 0.59 0.65 0.8 1.0  )
# foreach KB ( 0.005 )

set Nsim=1

# conda activate ovito

#!/bin/tcsh
#$ -cwd
#$ -j y


for nu in 0.5 0.55 0.6 0.65 0.7 0.75 0.85 0.9 0.95 1.0
do
for KB in 0.1
do
for Init_cond in 1
do 
    python Movie_ovito.py ${nu} ${Init_cond} 1 ${KB}
done 
done
done




for nu in 0.45 0.5 0.55 0.6 0.65 0.7 0.75
do
for KB in 0.1
do
for Init_cond in 2
do 
    python Movie_ovito.py ${nu} ${Init_cond} 1 ${KB}
done 
done
done



for nu in 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.55 0.6
do
for KB in 0.1
do
for Init_cond in 3
do 
    python Movie_ovito.py ${nu} ${Init_cond} 1 ${KB}
done 
done
done




#deactivate

# conda deactivate ovito
