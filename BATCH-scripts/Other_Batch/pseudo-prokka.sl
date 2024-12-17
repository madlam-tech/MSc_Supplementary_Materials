#!/bin/bash -e
#SBATCH --account massey03345
#SBATCH -J pseudo-prokka
#SBATCH --time 04:00:00
#SBATCH --mem 4GB
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 12
#SBATCH -e pseudo-prokka%j.log

module purge
module load prokka/1.14.5-GCC-9.2.0
module load barrnap/0.9-GCC-9.2.0

# Define the directory containing your .ffn files
input_dir="/nesi/nobackup/massey03345/fnas_4_annotation/pseudofinder/pseudogene_output/"

# Define the directory where you want Prokka to save the annotation results
output_dir="/nesi/nobackup/massey03345/fnas_4_annotation/pseudofinder/pseudogene_output/prokka_reannotate"

# Create the output directory if it doesn't exist
mkdir -p "$output_dir"

# Loop through each .ffn file in the input directory
for ffn_file in "$input_dir"/*.ffn; do
    # Check if the file exists and is readable
    if [ -r "$ffn_file" ]; then
        # Extract the filename without the path and extension
        filename=$(basename "$ffn_file" .ffn)
        
        # Run Prokka to annotate the .ffn file
        prokka --outdir "$output_dir/$filename" --prefix "$filename" "$ffn_file"
    else
        echo "Error: Unable to read file: $ffn_file"
    fi
done

