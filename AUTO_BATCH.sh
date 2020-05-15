#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem-per-cpu=2G
#SBATCH --time=0-1:00:00
#SBATCH --output=Haploid_sims_20.stdout
#SBATCH --job-name="Haploid_sims_20"
#SBATCH -p short 
export OMP_NUM_THREADS 12
mkdir Animate_Haploid_sims_20
./program Animate_Haploid_sims_20 -ADH 0 -division 0
