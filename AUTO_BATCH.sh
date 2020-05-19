#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem-per-cpu=2G
#SBATCH --time=0-12:00:00
#SBATCH --output=first_nutrient_tests_long_8.stdout
#SBATCH --job-name="first_nutrient_tests_long_8"
#SBATCH -p long.q 
export OMP_NUM_THREADS 12
mkdir Animate_first_nutrient_tests_long_8
./program Animate_first_nutrient_tests_long_8 -ADH 0 -HAPLOID 1
