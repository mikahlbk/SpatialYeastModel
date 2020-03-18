#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem-per-cpu=2G
#SBATCH --time=01:15:00
#SBATCH --output=ADH_Final_608.stdout
#SBATCH --job-name="ADH_Final_608"
#SBATCH -p short

set OMP_NUM_THREADS=12
mkdir ADH_Final_608
./program ADH_Final_608
