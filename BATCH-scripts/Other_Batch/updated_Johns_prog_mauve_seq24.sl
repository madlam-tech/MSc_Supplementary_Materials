#!/bin/bash -e

#SBATCH --account       massey03345
#SBATCH --job-name      progressiveMauve_24
#SBATCH --time          04:00:00        # Increase time as needed
#SBATCH --mem           20G             # Increase memory
#SBATCH --cpus-per-task 8               # Increase CPUs
#SBATCH --output        ./job_%j/slurm_%j.out
#SBATCH --error         ./job_%j/slurm_%j.err

module purge 2>/dev/null
module load Java

# Define the main sequence file
seq_file="sequence_24.fna"

# Check if the main sequence file exists
if [[ ! -f $seq_file ]]; then
    echo "Sequence file $seq_file does not exist."
    exit 1
fi

# Loop through all sequence files and compare to the main sequence file
for x in sequence_*.fna; do
    [[ $x == $seq_file ]] && continue
    basename="${x%%.*}"
    index2=${basename##*_}
    output_file="alignment_24_${index2}.xmfa"
    output_guide_file="guide_tree_24_${index2}.txt"
    
    echo "Running progressiveMauve for ${seq_file} and ${x}"
    
    # Run progressiveMauve with additional logging
    set -x  # Enable debugging
    ./progressiveMauve --output=${output_file} --output-guide-tree=${output_guide_file} --seed-family ${seq_file} ${x} > "${output_file}.log" 2>&1
    set +x  # Disable debugging
    
    if [ $? -ne 0 ]; then
        echo "Error in running progressiveMauve for ${seq_file} and ${x}. See ${output_file}.log for details."
        exit 1
    fi
done

echo -e "\n-------------finished all runs-------------"

