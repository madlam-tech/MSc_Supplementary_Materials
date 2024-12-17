#!/bin/sh -e
#SBATCH --account massey03345
#SBATCH -J minipolish
#SBATCH --time 00:15:00
#SBATCH --mem 5GB
#SBATCH --ntasks 2
#SBATCH --cpus-per-task 2
#SBATCH -e mp.err
#SBATCH -o mp.out


module purge 
module load Python/3.10.5-gimkl-2022a
module load miniasm/0.3-20191007-GCC-11.3.0
module load minimap2/2.24-GCC-9.2.0


minimap2 -t 8 -x ava-ont sample_01.fastq.gz  > overlaps.paf
miniasm -f sample_01.fastq.gz  overlaps.paf > assembly.gfa
minipolish.py -t 8 sample_01.fastq.gz assembly.gfa > polished.gfa
