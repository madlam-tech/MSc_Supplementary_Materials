#!/bin/bash -e
#SBATCH --account massey03345
#SBATCH -J 2-mm-racon.sl
#SBATCH --time 8:00:00
#SBATCH --mem 4GB
#SBATCH --ntasks 4
#SBATCH --cpus-per-task 8
#SBATCH -e 2-mm-racon_%j.err
#SBATCH -o 2-mm-racon_%j.out
#SBATCH --export NONE

module purge
module load Porechop/0.2.4-gimkl-2020a-Python-3.8.2
module load nanofilt/2.6.0-gimkl-2020a-Python-3.8.2


#clean up reads
porechop -i ./*fq -o ./porechopped.fastq

NanoFilt –l 500 --headcrop 10 < ./porechop/porechopped.fastq \
 > ./nanofilt_trimmed.fastq


#assemble, convert gfa -> fasta, polish with racon

module purge
module load minimap2/2.24-GCC-9.2.0
module load miniasm/0.3-20191007-GCC-11.3.0
module load Racon/1.5.0-GCC-11.3.0

minimap2 –x ava-ont \
 ./nanofilt_trimmed.fastq \
| gzip -1 > ./minimap.paf.gz

miniasm -f \
./nanofilt_trimmed.fastq \
./minimap.paf.gz > ./miniasm.gfa

#convert gfa to fasta

awk ’/^S/{print “>”$2”\n”$3}’ ./miniasm.gfa > ./miniasm.fasta
assembly-stats ./miniasm.fasta > ./miniasm.stats

#racon polish

racon ./nanofilt_trimmed.fastq \
./minimap.paf.gz \
 ./miniasm.fasta \ 
> ./miniasm.racon.consensus.fasta





