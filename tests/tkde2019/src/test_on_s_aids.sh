#!/bin/bash
#SBATCH --job-name=S-AIDS
#SBATCH --output=../bin/S-AIDS.txt
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mail-type=ALL
#SBATCH --mem-per-cpu=20000
#SBATCH --mail-user=david.blumenthal@wzw.tum.de

for NL in 01 04 07 10 13 16 19
do
	srun --ntasks 1 --cpus-per-task 8 ../bin/train_walks S-AIDS_NL${NL}
	srun --ntasks 1 --cpus-per-task 8 ../bin/train_subgraph S-AIDS_NL${NL}
	srun --ntasks 1 --cpus-per-task 8 ../bin/train_ring S-AIDS_NL${NL}
	srun --ntasks 1 --cpus-per-task 8 ../bin/train_ml --no-svm S-AIDS_NL${NL}
	srun --ntasks 1 --cpus-per-task 8 ../bin/run_tests --quick S-AIDS_NL${NL}
done