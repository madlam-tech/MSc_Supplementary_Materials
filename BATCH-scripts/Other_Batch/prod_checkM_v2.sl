#!/bin/bash -e
#SBATCH --account       massey03345
#SBATCH --job-name      prodigal_checkm
#SBATCH --time          10:00:00
#SBATCH --mem           10GB
#SBATCH --cpus-per-task 8
#SBATCH --error         prod%j.err
#SBATCH --output        prod%j.out

# Load modules
module purge
module load GCCcore/11.3.0
module load prodigal/2.6.3-GCC-11.3.0
module load CheckM/1.2.1-gimkl-2022a-Python-3.10.5

# Output directory
output_dir="prod_checkM"
mkdir -p "$output_dir"

# Function to process directories recursively
process_directory() {
    local dir="$1"

    # Loop through .fna files in the specified directory
    for filename in $(find "$dir" -name "*.fna"); do
        # Check if the file exists
        if [ -f "$filename" ]; then
            # Extract basename without extension
            basename=$(basename "$filename" .fna)

            # Create base directory for this file's output
            base_dir="$output_dir/$basename"
            mkdir -p "$base_dir/predictions"
            mkdir -p "$base_dir/checkM"

            # Run Prodigal on the .fna file
            prodigal -i "$filename" \
                     -o "$base_dir/predictions/${basename}_genes.gbk" \
                     -a "$base_dir/predictions/${basename}_proteins.faa" \
                     -d "$base_dir/predictions/${basename}_nucl.fnn" \
                     -p meta

            # Run CheckM analysis using your usual CLI command
            checkm lineage_wf --genes -t 8 --pplacer_threads 8 -x faa --tab_table -f "$base_dir/checkM/checkm.txt" "$base_dir/predictions/" "$base_dir/checkM/"

            # Capture the output of the CheckM quality assessment
            qa_output=$(cat "$base_dir/checkM/checkm.txt")
            echo -e "QA Output for $basename:\n$qa_output\n" >> "$base_dir/checkM/summary_qa.txt"
        fi
    done

    # Recursively process subdirectories
    for subdir in "$dir"/*/; do
        if [ -d "$subdir" ]; then
            process_directory "$subdir"
        fi
    done
}

# Start processing from the fnas directory
process_directory "$(pwd)/fnas"
