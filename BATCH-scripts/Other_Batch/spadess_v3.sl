#!/bin/bash -e

#SBATCH --account       massey03345
#SBATCH --job-name      spades_hybrid_assembly
#SBATCH --time          02:00:00
#SBATCH --mem           12GB
#SBATCH --cpus-per-task 12
#SBATCH --output=spades_hybrid_assembly_%j.out
#SBATCH --error=spades_hybrid_assembly_%j.err

module purge
module load SPAdes/3.15.4-gimkl-2022a-Python-3.10.5

mkdir -p ./CD_1_spades_output_deduped

# Run SPAdes from a BAM file
spades.py --only-assembler --careful --cov-cutoff auto --phred-offset 33 \
  	--pe1-1 CD1_alignments_deduplicated_output.bam --pe1-2 CD1_alignments_deduplicated_output.bam \
	--trusted-contigs AX4.fasta \
  	-o CD_1_spades_output_deduped \
 	-t 8

echo "SPAdes assembly completed successfully."


