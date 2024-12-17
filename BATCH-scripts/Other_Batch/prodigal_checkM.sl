#!/bin/bash -e
#SBATCH --account       massey03345
#SBATCH --job-name      prodigal_checkm
#SBATCH --time          10:00:00
#SBATCH --mem           10GB
#SBATCH --cpus-per-task 4
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

    # Loop through .fna files in the current directory
    for filename in "$dir"/*.fna; do
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

            # Run CheckM analysis on the .fna file
            checkm analyze -x fna --threads 4 -o 2 "$filename" "$base_dir/checkM"
            checkm tree -x fna --threads 4 "$filename" "$base_dir/checkM"

            # CheckM quality assessment
            checkm qa "$base_dir/checkM/lineage.ms" "$base_dir/checkM" --threads 4 -o 2 > "$base_dir/checkM/${basename}_qa.txt"
            
            # Capture the output of the QA assessment
            qa_output=$(cat "$base_dir/checkM/${basename}_qa.txt")
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

# Start processing from the current directory
process_directory "$(pwd)"
