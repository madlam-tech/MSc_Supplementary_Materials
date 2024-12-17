#!/bin/bash -e

#SBATCH --account       massey03345
#SBATCH --job-name      picard
#SBATCH --time          00:20:00
#SBATCH --mem           12GB
#SBATCH --cpus-per-task 12
#SBATCH --error         %x_%j.err
#SBATCH --output        %x_%j.out


module purge
module load picard/2.26.10-Java-11.0.4

for filename in *.sam

do


# Replace 'your_input.sam' with the actual name of your input SAM file
input_sam=${filename}

# Picard tools JAR file path. Adjust the path if needed.
#picard_jar="/path/to/picard.jar"

# Output file names
sorted_bam=${filename}.bam
dedup_bam=${filename}_deduplicated_output.bam
index_metrics= ${filename}_index_metrics.txt

# Step 1: Sort SAM file
SortSam \
    INPUT="$input_sam" \
    OUTPUT="$sorted_bam" \
    SORT_ORDER=coordinate

# Step 2: Mark duplicates
MarkDuplicates \
    INPUT="$sorted_bam" \
    OUTPUT="$dedup_bam" \
    METRICS_FILE="$index_metrics"

# Step 3: Build BAM index
BuildBamIndex \
    INPUT="$dedup_bam"

done
