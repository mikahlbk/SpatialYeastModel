#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem-per-cpu=2G
#SBATCH --time=0-4:00:00
#SBATCH --output=hypothesis_sims_3.stdout
#SBATCH --job-name="hypothesis_sims_3"
#SBATCH -p fast.q 
export OMP_NUM_THREADS 12
mkdir Animate_hypothesis_sims_3
./program Animate_hypothesis_sims_3 -Budding 1 -division 2 -nutrient_decay .00027
