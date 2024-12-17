#!/bin/bash -e
#SBATCH --account massey03345
#SBATCH -J prokka_v12
#SBATCH --time 06:30:00
#SBATCH --mem 4GB
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 12
#SBATCH -e 22-prokka_multi_v12%j.log
mkdir -p prokka_v12

module purge
module load prokka/1.14.5-GCC-9.2.0
module load barrnap/0.9-GCC-9.2.0


for filename in *.fna

do


BASENAME=$(basename ${filename} _1.fna)


prokka ${filename} \
	--force \
	--outdir prokka-2023-05-12_v13/prokka.${BASENAME} \
	--prefix ${BASENAME} --addgenes --locustag ${BASENAME} \
	--proteins GCF_000756125.1_ASM75612v1_genomic.gbff \
	--genus Paraburkholderia \
	--kingdom Bacteria --gcode 11 --usegenus \
	--addgenes \
	--compliant --rfam 	
done
