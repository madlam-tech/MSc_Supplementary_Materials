#!/bin/bash -e
#SBATCH --account massey03345
#SBATCH -J prokka
#SBATCH --time 02:30:00
#SBATCH --mem 1GB
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 12
#SBATCH -e 21-prokka_multi_%j.log
#SBATCH -o 21-prokka_multi%j.log

mkdir -p prokka_v2

module purge
module load prokka/1.14.5-GCC-9.2.0
module load RNAmmer/1.2-GCC-9.2.0-Perl-5.30.1

for filename in *.fasta

do


BASENAME=$(basename ${filename} .fasta)

prokka --compliant --rfam  ${filename} \
	--dbdir /scale_wlg_persistent/filesets/opt_nesi/CS400_centos7_bdw/prokka/1.14.5-GCC-9.2.0/db \
	--outdir ./prokka_v2/prokka.${BASENAME}
	--prefix ${BASENAME} \
	--addgenes \
	--addmrna \
	--coverage 70 

done
