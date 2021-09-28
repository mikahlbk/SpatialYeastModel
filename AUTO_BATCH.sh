#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=2G
#SBATCH --time=0-12:00:00
#SBATCH --output=first_official_runs_4.stdout
#SBATCH --job-name="first_official_runs_4"
#SBATCH -p nodes 
export OMP_NUM_THREADS=16
mkdir Animate_first_official_runs_4
./program Animate_first_official_runs_4 -Budding 1 -Nutrient_Condition 1 -Substrate_Consumption 0 -Bistable_protein_dyn .005 -Initial_Protein 100 -Replication_rate .018 -Division_bias .36
