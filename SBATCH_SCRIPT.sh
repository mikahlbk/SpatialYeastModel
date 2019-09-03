#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G
#SBATCH --time=24:15:00
#SBATCH --output=Animate_Jonsson_bucket_test_day_2.stdout
#SBATCH --job-name="Animate_Jonsson_bucket_test_day_2"
#SBATCH -p batch	

set OMP_NUM_THREADS=12
mkdir Animate_Jonsson_bucket_test_day_2
./program Animate_Jonsson_bucket_test_day_2
                   
