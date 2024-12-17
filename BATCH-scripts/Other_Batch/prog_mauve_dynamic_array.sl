#!/bin/bash -e

#SBATCH --account       massey03345
#SBATCH --job-name      progressiveMauve
#SBATCH --time          01:00:00
#SBATCH --mem           8G
#SBATCH --cpus-per-task 2
#SBATCH --output        ./job_%j/slurm_%j.out
#SBATCH --error         ./job_%j/slurm_%j.err

module purge 2>/dev/null
module load Java

# Get the number of sequence files
num_seq_files=$(find . -maxdepth 1 -name "sequence_*.fna" | wc -l)

# Dynamically set the array size based on the number of sequence files
if (( ${num_seq_files} < 2 )); then
    echo "ERROR: Not enough sequence files. At least 2 are required."
    exit 1
fi

# Export the dynamic array size for SLURM
export SLURM_ARRAY_TASK_COUNT=${num_seq_files}

echo "Comparing '${num_seq_files}' sequence files in the current directory."

index1=${SLURM_ARRAY_TASK_ID}
seq_file=$(find . -maxdepth 1 -name "sequence_*.fasta" | sed -n "${index1}p")

for x in ./sequence_*.fasta; do
    [[ $x == $seq_file ]] && continue
    basename="${x%%.*}"
    index2=${basename##*_}
    output_file="alignment_${index1}_${index2}.xmfa"
    output_guide_file="guide_tree_${index1}_${index2}.txt"
    ./progressiveMauve --output=${output_file} --output-guide-tree=${output_guide_file} --seed-family ${seq_file} ${x}
done

echo "Finished all runs."
