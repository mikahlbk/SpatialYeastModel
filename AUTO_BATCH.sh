#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem-per-cpu=2G
#SBATCH --time=0-4:00:00
#SBATCH --output=final_nutrient_sims_longer_50.stdout
#SBATCH --job-name="final_nutrient_sims_longer_50"
#SBATCH -p fast.q 
export OMP_NUM_THREADS 12
mkdir Animate_final_nutrient_sims_longer_50
./program Animate_final_nutrient_sims_longer_50 -Budding 0 -nutrient_depletion 1 -start_from_four 0 -division 2
