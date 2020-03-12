#!/bin/bash
#SBATCH --job-name=soa_aids
#SBATCH --output=soa_aids.txt
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=20000
#SBATCH --mail-type=ALL
#SBATCH --mail-user=david.blumenthal@wzw.tum.de

export OMP_NUM_THREADS=6

for ID in {0..4}
do
	for PERCENT in {1..9}0
	do
		srun --ntasks 1 --cpus-per-task $OMP_NUM_THREADS ./test_VLDBJ_final/chemical-edit-distances ../../../data/datasets/AIDS/data/AIDS -I ${ID} -P ${PERCENT} -m lsape_multi_bunke -i lsape_multi_bunke -p 1  > results/res_AIDS-${PERCENT}-${ID}.csv
	done
done
