#!/bin/bash -e
#SBATCH --account massey03345
#SBATCH -J tryf_v2  
#SBATCH --time 00:15:00
#SBATCH --mem 2GB
#SBATCH --ntasks 4
#SBATCH --cpus-per-task 8
#SBATCH -e tryFRR.err
#SBATCH -o tryFRR.out
#SBATCH --mail-user=madlam@me.com
#SBATCH --mail-type=ALL
#SBATCH --export NONE

threads=32
mkdir -p assemblies_v3



module purge
module load Python/3.10.5-gimkl-2022a
module load any2fasta/0.4.2-GCC-9.2.0
module load minimap2/2.24-GCC-9.2.0
module load Raven/1.5.0-GCC-9.2.0
module load Racon/1.5.0-GCC-11.3.0
module load miniasm/0.3-20191007-GCC-11.3.0

./miniasm_and_minipolish.sh ./sample_02.fastq "$threads" > assemblies/assembly_02.gfa && any2fasta assemblies/assembly_02.gfa > assemblies/assembly_02.fasta


