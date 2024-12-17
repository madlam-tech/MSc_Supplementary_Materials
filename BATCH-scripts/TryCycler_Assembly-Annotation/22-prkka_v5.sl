#!/bin/bash -e
#SBATCH --account massey03345
#SBATCH -J prokka_v5
#SBATCH --time 02:30:00
#SBATCH --mem 5GB
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 12
#SBATCH -e 22-prokka_multi_v5%j.log
mkdir -p prokka_v5

module purge
module load prokka/1.14.5-GCC-9.2.0
module load RNAmmer/1.2-GCC-9.2.0-Perl-5.30.1


for filename in *.fasta

do


BASENAME=$(basename ${filename} .fasta)



prokka ${filename} \
	--force \
	--outdir prokka-2023-05-05_v5/prokka.${BASENAME} \
	--prefix ${BASENAME} --addgenes --locustag PAGRI \
	--genus Paraburkholderia --species agricolaris --strain BaQS159 \
	--kingdom Bacteria. --gcode 11 --usegenus
	--proteins GCF_000756125.1_ASM75612v1_genomic.gbff 
	--gram neg \
	--addgenes \
	--rnammer
 
done

