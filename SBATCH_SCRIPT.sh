#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem-per-cpu=2G
#SBATCH --time=5:00:00
#SBATCH --output=NO_ADH_Final_5.stdout
#SBATCH --job-name="NO_ADH_Final_5"
#SBATCH -p batch 

set OMP_NUM_THREADS=12
mkdir NO_ADH_Final_5
./program NO_ADH_Final_5
              
