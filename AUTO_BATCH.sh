#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem-per-cpu=2G
#SBATCH --time=0-5:00:00
#SBATCH --output=Mixed_test_second_gen_6.stdout
#SBATCH --job-name="Mixed_test_second_gen_6"
#SBATCH -p batch 
export OMP_NUM_THREADS 12
mkdir Animate_Mixed_test_second_gen_6
./program Animate_Mixed_test_second_gen_6 -ADH 0
