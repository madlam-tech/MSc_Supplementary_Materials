#!/bin/bash -e
#SBATCH --account massey03345
#SBATCH -J roary-merged
#SBATCH --time 20:00:00
#SBATCH --mem 4GB
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 24
#SBATCH -e roary-merged_%j.err
#SBATCH -o roary-merged_%j.log

module purge
module load Roary/3.13.0-gimkl-2020a

# Get the current date in YYYY-MM-DD format
current_date=$(date +%F)

# Define the output directory
output_dir="./sml_Pb_ags_roary_output_${current_date}"

# Run Roary with the combined options and increased group limit
roary -r -s -p24 -I 70 -f "$output_dir" -e -n -v --group_limit 200000 ./sml_set_Pb_ags/*.gff
