#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem-per-cpu=2G
#SBATCH --time=0-12:00:00
#SBATCH --output=Protein_Sims_2_40.stdout
#SBATCH --job-name="Protein_Sims_2_40"
#SBATCH -p batch 
export OMP_NUM_THREADS 12
mkdir Animate_Protein_Sims_2_40
./program Animate_Protein_Sims_2_40 -ADH 0 -Initial_Protein 40 -growth_rate .06 -competition_term 35
