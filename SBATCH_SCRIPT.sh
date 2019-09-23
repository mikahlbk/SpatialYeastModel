#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G
#SBATCH --time=00:15:00
#SBATCH --output=Animate_Hertz_test_6.stdout
#SBATCH --job-name="Animate_Hertz_test_6"
#SBATCH -p short	

set OMP_NUM_THREADS=12
mkdir Animate_Hertz_test_6
./program Animate_Hertz_test_6
                   
