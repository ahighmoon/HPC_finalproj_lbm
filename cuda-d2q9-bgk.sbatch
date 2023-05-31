#!/bin/bash

#SBATCH --job-name cuda-d2q9-bgk
#SBATCH --time 00:10:00
#SBATCH --exclusive 
#SBATCH --gres=gpu:1
echo Running on host `hostname`
echo Slurm job ID is $SLURM_JOB_ID
#echo This job runs on the following machines:
#echo `echo $SLURM_JOB_NODELIST | uniq`

./cuda-d2q9-bgk input_256x128.params obstacles_256x128.dat
./cuda-d2q9-bgk input_512x256.params obstacles_512x256.dat
./cuda-d2q9-bgk input_1024x512.params obstacles_1024x512.dat