#!/bin/bash
#SBATCH --job-name=median_baseline
#SBATCH --output=../bin/median_baseline.txt
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=6
#SBATCH --mail-type=ALL
#SBATCH --mem-per-cpu=20000
#SBATCH --mail-user=david.blumenthal@wzw.tum.de

export OMP_NUM_THREADS=6

srun --ntasks 1 --cpus-per-task $OMP_NUM_THREADS ../bin/median_tests AIDS -b
srun --ntasks 1 --cpus-per-task $OMP_NUM_THREADS ../bin/median_tests Mutagenicity -b
srun --ntasks 1 --cpus-per-task $OMP_NUM_THREADS ../bin/median_tests Letter-1 -b
srun --ntasks 1 --cpus-per-task $OMP_NUM_THREADS ../bin/median_tests Letter-2 -b
srun --ntasks 1 --cpus-per-task $OMP_NUM_THREADS ../bin/median_tests AIDS-EDIT -b
srun --ntasks 1 --cpus-per-task $OMP_NUM_THREADS ../bin/median_tests S-MOL-5 -b
