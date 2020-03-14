#!/bin/bash
#SBATCH --job-name=soa_ibd
#SBATCH --output=soa_ibd.txt
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=30000
#SBATCH --mail-type=ALL
#SBATCH --mail-user=david.blumenthal@wzw.tum.de

export OMP_NUM_THREADS=6

for CLASS in C NO YES
do
	srun --ntasks 1 --cpus-per-task $OMP_NUM_THREADS ./test_VLDBJ_final/chemical-edit-distances ../../../data/datasets/IBD/data/IBD_${CLASS}.ds -q ../../../src/edit_costs/otu_distances.csv -m lsape_multi_bunke -i lsape_multi_bunke -p 1 -l 60000 > results/res_IBD_${CLASS}.csv
done