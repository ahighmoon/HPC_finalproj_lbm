#!/bin/bash

#SBATCH --job-name d2q9-bgk
#SBATCH --nodes 1
#SBATCH --time 00:10:00
#SBATCH --exclusive 

echo Running on host `hostname`
echo Directory is `pwd`
echo Slurm job ID is $SLURM_JOB_ID

./d2q9-bgk input_256x128.params obstacles_256x128.dat
./d2q9-bgk input_512x256.params obstacles_512x256.dat
./d2q9-bgk input_1024x512.params obstacles_1024x512.dat
