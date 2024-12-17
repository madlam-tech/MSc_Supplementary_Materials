#!/bin/bash -e
#SBATCH --account massey03345
#SBATCH -J samtools
#SBATCH --time 00:05:00
#SBATCH --mem 1GB
#SBATCH --ntasks 2
#SBATCH --cpus-per-task 2
#SBATCH -e st.err
#SBATCH -o st.out


module purge
module load SAMtools/1.9-GCC-7.4.0


for filename in *.bam

do

samtools collate -o ${filename}col ${filename}

# Add ms and MC tags for markdup to use later
samtools fixmate -m ${filename}col ${filename}fm
# Markdup needs position order
samtools sort -n -o ${filename}ps ${filename}fm
# Finally mark duplicates. Use -r flag to remove duplicates, and -s to print stats.
samtools markdup -l 150 -r -s -T -S-f ${filename}stat -c -m s ${filename}ps ${filename}MD


done
