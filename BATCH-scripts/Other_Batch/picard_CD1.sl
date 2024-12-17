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


SartSam INPUT=CD1_alignments.sam OUTPUT=CD1_alignments.bam SORT_ORDER=coordinate \
MarkDuplicates INPUT=CD1_alignments.bam OUTPUT=CD1_alignments_dedup.bam METRICS_FILE=CD1_alignments_index_metrics.txt
BuildBamIndex INPUT=CD1_alignments_dedup.bam≈
