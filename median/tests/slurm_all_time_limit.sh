#!/bin/bash
#SBATCH --job-name=time_limit
#SBATCH --output=../bin/time_limit.txt
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=6
#SBATCH --mail-type=ALL
#SBATCH --mem-per-cpu=20000
#SBATCH --mail-user=david.blumenthal@wzw.tum.de

export OMP_NUM_THREADS=6

srun --ntasks 1 --cpus-per-task $OMP_NUM_THREADS ../bin/time_limit_tests Letter
srun --ntasks 1 --cpus-per-task $OMP_NUM_THREADS ../bin/time_limit_tests AIDS
srun --ntasks 1 --cpus-per-task $OMP_NUM_THREADS ../bin/time_limit_tests Mutagenicity