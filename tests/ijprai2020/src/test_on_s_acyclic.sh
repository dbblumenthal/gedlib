#!/bin/bash
#SBATCH --job-name=S-acyclic
#SBATCH --output=../bin/S-acyclic.txt
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mail-type=ALL
#SBATCH --mem-per-cpu=20000
#SBATCH --mail-user=david.blumenthal@wzw.tum.de

for NL in 03 05 07 09
do
	srun --ntasks 1 --cpus-per-task 8 ../bin/train_walks S-acyclic_NL${NL}
	srun --ntasks 1 --cpus-per-task 8 ../bin/train_subgraph S-acyclic_NL${NL}
	srun --ntasks 1 --cpus-per-task 8 ../bin/train_ring S-acyclic_NL${NL}
	srun --ntasks 1 --cpus-per-task 8 ../bin/train_ml --no-svm S-acyclic_NL${NL}
	srun --ntasks 1 --cpus-per-task 8 ../bin/run_tests --quick S-acyclic_NL${NL}
done