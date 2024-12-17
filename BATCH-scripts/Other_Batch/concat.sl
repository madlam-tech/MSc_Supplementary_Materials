#!/bin/bash -e
#SBATCH --account        massey03345
#SBATCH --job-name       cat07_220922
#SBATCH --mem            1G
#SBATCH --cpus-per-task  2
#SBATCH --time           00:10:00
#SBATCH --output         slurmout.%j.out
#SBATCH -e cat07.err
#SBATCH -o cat07.out



for name in *.fastq
do
    if [[ "$name" =~ ([0-9-]+)_.*(..)\.fastq ]]; then
        outfile="${BASH_REMATCH[1]}_${BASH_REMATCH[2]}.fastq"

        cat "$name" >>"$outfile"
    fi
done

