#!/bin/bash
#SBATCH --job-name=median_aids_edit
#SBATCH --output=../bin/median_aids_edit.txt
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=6
#SBATCH --mail-type=ALL
#SBATCH --mem-per-cpu=50000
#SBATCH --mail-user=david.blumenthal@wzw.tum.de

export OMP_NUM_THREADS=6

srun --ntasks 1 --cpus-per-task $OMP_NUM_THREADS ../bin/median_tests AIDS-EDIT
