#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem-per-cpu=2G
#SBATCH --time=01:15:00
#SBATCH --output=Mixed_init_5.stdout
#SBATCH --job-name="Mixed_init_5"
#SBATCH -p short

set OMP_NUM_THREADS=12
mkdir Mixed_init_5
./program Mixed_init_5
