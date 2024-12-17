#!/bin/bash -e
#SBATCH --account massey03345
#SBATCH -J roary-merged
#SBATCH --time 20:00:00
#SBATCH --mem 6GB
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 24
#SBATCH -e roary-merged_%j.err
#SBATCH -o roary-merged_%j.log


module purge
module load Roary/3.13.0-gimkl-2020a

# Get the current date in YYYY-MM-DD format
current_date=$(date +%F)

# Define the output directory
output_dir="./Paraburks_all_roary_output_${current_date}"

# Run Roary with the combined options
roary -r -s -p 24 -i 70 -f "$output_dir" -e -n -v ./Paraburks/*.gff
roary -p 24 -f ./roary_output --mafft -i 95 -cd 99 --group_limit 200000 ./Paraburks/*.gff
