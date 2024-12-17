#!/bin/bash -e
#SBATCH --account       massey03345
#SBATCH --job-name      metaxa2
#SBATCH --time          00:15:00
#SBATCH --mem           1GB
#SBATCH --cpus-per-task 4
#SBATCH --error         10-metaxa%j.err

module purge
module load Metaxa2/2.2.3-gimkl-2022a
# Output directory
mkdir -p ribosomes/

# Run Metaxa2
for ribosome_type in ssu lsu; do
  metaxa2  -g ${ribosome_type} --mode genome \
          -i 8_medaka.fasta -o ribosomes/assembly.fa.${ribosome_type}
done
