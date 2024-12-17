#!/bin/bash -e
#SBATCH --account=massey03345
#SBATCH --job-name=roary-merged
#SBATCH --time=1:00:00
#SBATCH --mem=6GB
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --error=roary-merged_%j.err
#SBATCH --output=roary-merged_%j.log

module purge
module load Roary/3.13.0-gimkl-2020a

# Get the current date in YYYY-MM-DD format
current_date=$(date +%F)

# Define the output directory
output_dir="./Acinetobacter_roary_output_${current_date}"

# Run Roary with the combined options
roary -r -s -p 24 -i 95 -cd 99 -f "$output_dir" -e -n -v --mafft --group_limit 300000 ./gffs_4_roary/*.gff

echo "Roary analysis completed."
