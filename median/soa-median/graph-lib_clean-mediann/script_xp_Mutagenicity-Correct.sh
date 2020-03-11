#!/bin/bash
#SBATCH --job-name=soa_muta
#SBATCH --output=soa_muta.txt
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=20000
#SBATCH --mail-type=ALL
#SBATCH --mail-user=david.blumenthal@wzw.tum.de

export OMP_NUM_THREADS=6

srun --ntasks 1 --cpus-per-task $OMP_NUM_THREADS ./test_VLDBJ_final/chemical-edit-distances ../datasets/MUTA/Mutagenicity-Correct -I 0 -P 10 -m lsape_multi_bunke -i lsape_multi_bunke -p 1  > results/res_Mutagenicity-Correct-10-0.csv 
srun --ntasks 1 --cpus-per-task $OMP_NUM_THREADS ./test_VLDBJ_final/chemical-edit-distances ../datasets/MUTA/Mutagenicity-Correct -I 0 -P 20 -m lsape_multi_bunke -i lsape_multi_bunke -p 1  > results/res_Mutagenicity-Correct-20-0.csv 
srun --ntasks 1 --cpus-per-task $OMP_NUM_THREADS ./test_VLDBJ_final/chemical-edit-distances ../datasets/MUTA/Mutagenicity-Correct -I 0 -P 30 -m lsape_multi_bunke -i lsape_multi_bunke -p 1  > results/res_Mutagenicity-Correct-30-0.csv 
srun --ntasks 1 --cpus-per-task $OMP_NUM_THREADS ./test_VLDBJ_final/chemical-edit-distances ../datasets/MUTA/Mutagenicity-Correct -I 0 -P 40 -m lsape_multi_bunke -i lsape_multi_bunke -p 1  > results/res_Mutagenicity-Correct-40-0.csv 
srun --ntasks 1 --cpus-per-task $OMP_NUM_THREADS ./test_VLDBJ_final/chemical-edit-distances ../datasets/MUTA/Mutagenicity-Correct -I 0 -P 50 -m lsape_multi_bunke -i lsape_multi_bunke -p 1  > results/res_Mutagenicity-Correct-50-0.csv 
srun --ntasks 1 --cpus-per-task $OMP_NUM_THREADS ./test_VLDBJ_final/chemical-edit-distances ../datasets/MUTA/Mutagenicity-Correct -I 0 -P 60 -m lsape_multi_bunke -i lsape_multi_bunke -p 1  > results/res_Mutagenicity-Correct-60-0.csv 
srun --ntasks 1 --cpus-per-task $OMP_NUM_THREADS ./test_VLDBJ_final/chemical-edit-distances ../datasets/MUTA/Mutagenicity-Correct -I 0 -P 70 -m lsape_multi_bunke -i lsape_multi_bunke -p 1  > results/res_Mutagenicity-Correct-70-0.csv 
srun --ntasks 1 --cpus-per-task $OMP_NUM_THREADS ./test_VLDBJ_final/chemical-edit-distances ../datasets/MUTA/Mutagenicity-Correct -I 0 -P 80 -m lsape_multi_bunke -i lsape_multi_bunke -p 1  > results/res_Mutagenicity-Correct-80-0.csv 
srun --ntasks 1 --cpus-per-task $OMP_NUM_THREADS ./test_VLDBJ_final/chemical-edit-distances ../datasets/MUTA/Mutagenicity-Correct -I 0 -P 90 -m lsape_multi_bunke -i lsape_multi_bunke -p 1  > results/res_Mutagenicity-Correct-90-0.csv 
srun --ntasks 1 --cpus-per-task $OMP_NUM_THREADS ./test_VLDBJ_final/chemical-edit-distances ../datasets/MUTA/Mutagenicity-Correct -I 1 -P 10 -m lsape_multi_bunke -i lsape_multi_bunke -p 1  > results/res_Mutagenicity-Correct-10-1.csv 
srun --ntasks 1 --cpus-per-task $OMP_NUM_THREADS ./test_VLDBJ_final/chemical-edit-distances ../datasets/MUTA/Mutagenicity-Correct -I 1 -P 20 -m lsape_multi_bunke -i lsape_multi_bunke -p 1  > results/res_Mutagenicity-Correct-20-1.csv 
srun --ntasks 1 --cpus-per-task $OMP_NUM_THREADS ./test_VLDBJ_final/chemical-edit-distances ../datasets/MUTA/Mutagenicity-Correct -I 1 -P 30 -m lsape_multi_bunke -i lsape_multi_bunke -p 1  > results/res_Mutagenicity-Correct-30-1.csv 
srun --ntasks 1 --cpus-per-task $OMP_NUM_THREADS ./test_VLDBJ_final/chemical-edit-distances ../datasets/MUTA/Mutagenicity-Correct -I 1 -P 40 -m lsape_multi_bunke -i lsape_multi_bunke -p 1  > results/res_Mutagenicity-Correct-40-1.csv 
srun --ntasks 1 --cpus-per-task $OMP_NUM_THREADS ./test_VLDBJ_final/chemical-edit-distances ../datasets/MUTA/Mutagenicity-Correct -I 1 -P 50 -m lsape_multi_bunke -i lsape_multi_bunke -p 1  > results/res_Mutagenicity-Correct-50-1.csv 
srun --ntasks 1 --cpus-per-task $OMP_NUM_THREADS ./test_VLDBJ_final/chemical-edit-distances ../datasets/MUTA/Mutagenicity-Correct -I 1 -P 60 -m lsape_multi_bunke -i lsape_multi_bunke -p 1  > results/res_Mutagenicity-Correct-60-1.csv 
srun --ntasks 1 --cpus-per-task $OMP_NUM_THREADS ./test_VLDBJ_final/chemical-edit-distances ../datasets/MUTA/Mutagenicity-Correct -I 1 -P 70 -m lsape_multi_bunke -i lsape_multi_bunke -p 1  > results/res_Mutagenicity-Correct-70-1.csv 
srun --ntasks 1 --cpus-per-task $OMP_NUM_THREADS ./test_VLDBJ_final/chemical-edit-distances ../datasets/MUTA/Mutagenicity-Correct -I 1 -P 80 -m lsape_multi_bunke -i lsape_multi_bunke -p 1  > results/res_Mutagenicity-Correct-80-1.csv 
srun --ntasks 1 --cpus-per-task $OMP_NUM_THREADS ./test_VLDBJ_final/chemical-edit-distances ../datasets/MUTA/Mutagenicity-Correct -I 1 -P 90 -m lsape_multi_bunke -i lsape_multi_bunke -p 1  > results/res_Mutagenicity-Correct-90-1.csv 
srun --ntasks 1 --cpus-per-task $OMP_NUM_THREADS ./test_VLDBJ_final/chemical-edit-distances ../datasets/MUTA/Mutagenicity-Correct -I 2 -P 10 -m lsape_multi_bunke -i lsape_multi_bunke -p 1  > results/res_Mutagenicity-Correct-10-2.csv 
srun --ntasks 1 --cpus-per-task $OMP_NUM_THREADS ./test_VLDBJ_final/chemical-edit-distances ../datasets/MUTA/Mutagenicity-Correct -I 2 -P 20 -m lsape_multi_bunke -i lsape_multi_bunke -p 1  > results/res_Mutagenicity-Correct-20-2.csv 
srun --ntasks 1 --cpus-per-task $OMP_NUM_THREADS ./test_VLDBJ_final/chemical-edit-distances ../datasets/MUTA/Mutagenicity-Correct -I 2 -P 30 -m lsape_multi_bunke -i lsape_multi_bunke -p 1  > results/res_Mutagenicity-Correct-30-2.csv 
srun --ntasks 1 --cpus-per-task $OMP_NUM_THREADS ./test_VLDBJ_final/chemical-edit-distances ../datasets/MUTA/Mutagenicity-Correct -I 2 -P 40 -m lsape_multi_bunke -i lsape_multi_bunke -p 1  > results/res_Mutagenicity-Correct-40-2.csv 
srun --ntasks 1 --cpus-per-task $OMP_NUM_THREADS ./test_VLDBJ_final/chemical-edit-distances ../datasets/MUTA/Mutagenicity-Correct -I 2 -P 50 -m lsape_multi_bunke -i lsape_multi_bunke -p 1  > results/res_Mutagenicity-Correct-50-2.csv 
srun --ntasks 1 --cpus-per-task $OMP_NUM_THREADS ./test_VLDBJ_final/chemical-edit-distances ../datasets/MUTA/Mutagenicity-Correct -I 2 -P 60 -m lsape_multi_bunke -i lsape_multi_bunke -p 1  > results/res_Mutagenicity-Correct-60-2.csv 
srun --ntasks 1 --cpus-per-task $OMP_NUM_THREADS ./test_VLDBJ_final/chemical-edit-distances ../datasets/MUTA/Mutagenicity-Correct -I 2 -P 70 -m lsape_multi_bunke -i lsape_multi_bunke -p 1  > results/res_Mutagenicity-Correct-70-2.csv 
srun --ntasks 1 --cpus-per-task $OMP_NUM_THREADS ./test_VLDBJ_final/chemical-edit-distances ../datasets/MUTA/Mutagenicity-Correct -I 2 -P 80 -m lsape_multi_bunke -i lsape_multi_bunke -p 1  > results/res_Mutagenicity-Correct-80-2.csv 
srun --ntasks 1 --cpus-per-task $OMP_NUM_THREADS ./test_VLDBJ_final/chemical-edit-distances ../datasets/MUTA/Mutagenicity-Correct -I 2 -P 90 -m lsape_multi_bunke -i lsape_multi_bunke -p 1  > results/res_Mutagenicity-Correct-90-2.csv 
srun --ntasks 1 --cpus-per-task $OMP_NUM_THREADS ./test_VLDBJ_final/chemical-edit-distances ../datasets/MUTA/Mutagenicity-Correct -I 3 -P 10 -m lsape_multi_bunke -i lsape_multi_bunke -p 1  > results/res_Mutagenicity-Correct-10-3.csv 
srun --ntasks 1 --cpus-per-task $OMP_NUM_THREADS ./test_VLDBJ_final/chemical-edit-distances ../datasets/MUTA/Mutagenicity-Correct -I 3 -P 20 -m lsape_multi_bunke -i lsape_multi_bunke -p 1  > results/res_Mutagenicity-Correct-20-3.csv 
srun --ntasks 1 --cpus-per-task $OMP_NUM_THREADS ./test_VLDBJ_final/chemical-edit-distances ../datasets/MUTA/Mutagenicity-Correct -I 3 -P 30 -m lsape_multi_bunke -i lsape_multi_bunke -p 1  > results/res_Mutagenicity-Correct-30-3.csv 
srun --ntasks 1 --cpus-per-task $OMP_NUM_THREADS ./test_VLDBJ_final/chemical-edit-distances ../datasets/MUTA/Mutagenicity-Correct -I 3 -P 40 -m lsape_multi_bunke -i lsape_multi_bunke -p 1  > results/res_Mutagenicity-Correct-40-3.csv 
srun --ntasks 1 --cpus-per-task $OMP_NUM_THREADS ./test_VLDBJ_final/chemical-edit-distances ../datasets/MUTA/Mutagenicity-Correct -I 3 -P 50 -m lsape_multi_bunke -i lsape_multi_bunke -p 1  > results/res_Mutagenicity-Correct-50-3.csv 
srun --ntasks 1 --cpus-per-task $OMP_NUM_THREADS ./test_VLDBJ_final/chemical-edit-distances ../datasets/MUTA/Mutagenicity-Correct -I 3 -P 60 -m lsape_multi_bunke -i lsape_multi_bunke -p 1  > results/res_Mutagenicity-Correct-60-3.csv 
srun --ntasks 1 --cpus-per-task $OMP_NUM_THREADS ./test_VLDBJ_final/chemical-edit-distances ../datasets/MUTA/Mutagenicity-Correct -I 3 -P 70 -m lsape_multi_bunke -i lsape_multi_bunke -p 1  > results/res_Mutagenicity-Correct-70-3.csv 
srun --ntasks 1 --cpus-per-task $OMP_NUM_THREADS ./test_VLDBJ_final/chemical-edit-distances ../datasets/MUTA/Mutagenicity-Correct -I 3 -P 80 -m lsape_multi_bunke -i lsape_multi_bunke -p 1  > results/res_Mutagenicity-Correct-80-3.csv 
srun --ntasks 1 --cpus-per-task $OMP_NUM_THREADS ./test_VLDBJ_final/chemical-edit-distances ../datasets/MUTA/Mutagenicity-Correct -I 3 -P 90 -m lsape_multi_bunke -i lsape_multi_bunke -p 1  > results/res_Mutagenicity-Correct-90-3.csv 
srun --ntasks 1 --cpus-per-task $OMP_NUM_THREADS ./test_VLDBJ_final/chemical-edit-distances ../datasets/MUTA/Mutagenicity-Correct -I 4 -P 10 -m lsape_multi_bunke -i lsape_multi_bunke -p 1  > results/res_Mutagenicity-Correct-10-4.csv 
srun --ntasks 1 --cpus-per-task $OMP_NUM_THREADS ./test_VLDBJ_final/chemical-edit-distances ../datasets/MUTA/Mutagenicity-Correct -I 4 -P 20 -m lsape_multi_bunke -i lsape_multi_bunke -p 1  > results/res_Mutagenicity-Correct-20-4.csv 
srun --ntasks 1 --cpus-per-task $OMP_NUM_THREADS ./test_VLDBJ_final/chemical-edit-distances ../datasets/MUTA/Mutagenicity-Correct -I 4 -P 30 -m lsape_multi_bunke -i lsape_multi_bunke -p 1  > results/res_Mutagenicity-Correct-30-4.csv 
srun --ntasks 1 --cpus-per-task $OMP_NUM_THREADS ./test_VLDBJ_final/chemical-edit-distances ../datasets/MUTA/Mutagenicity-Correct -I 4 -P 40 -m lsape_multi_bunke -i lsape_multi_bunke -p 1  > results/res_Mutagenicity-Correct-40-4.csv 
srun --ntasks 1 --cpus-per-task $OMP_NUM_THREADS ./test_VLDBJ_final/chemical-edit-distances ../datasets/MUTA/Mutagenicity-Correct -I 4 -P 50 -m lsape_multi_bunke -i lsape_multi_bunke -p 1  > results/res_Mutagenicity-Correct-50-4.csv 
srun --ntasks 1 --cpus-per-task $OMP_NUM_THREADS ./test_VLDBJ_final/chemical-edit-distances ../datasets/MUTA/Mutagenicity-Correct -I 4 -P 60 -m lsape_multi_bunke -i lsape_multi_bunke -p 1  > results/res_Mutagenicity-Correct-60-4.csv 
srun --ntasks 1 --cpus-per-task $OMP_NUM_THREADS ./test_VLDBJ_final/chemical-edit-distances ../datasets/MUTA/Mutagenicity-Correct -I 4 -P 70 -m lsape_multi_bunke -i lsape_multi_bunke -p 1  > results/res_Mutagenicity-Correct-70-4.csv 
srun --ntasks 1 --cpus-per-task $OMP_NUM_THREADS ./test_VLDBJ_final/chemical-edit-distances ../datasets/MUTA/Mutagenicity-Correct -I 4 -P 80 -m lsape_multi_bunke -i lsape_multi_bunke -p 1  > results/res_Mutagenicity-Correct-80-4.csv 
srun --ntasks 1 --cpus-per-task $OMP_NUM_THREADS ./test_VLDBJ_final/chemical-edit-distances ../datasets/MUTA/Mutagenicity-Correct -I 4 -P 90 -m lsape_multi_bunke -i lsape_multi_bunke -p 1  > results/res_Mutagenicity-Correct-90-4.csv 
