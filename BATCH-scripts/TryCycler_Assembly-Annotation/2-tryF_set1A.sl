#!/bin/bash -e
#SBATCH --account massey03345
#SBATCH -J 2-tryF01-04.sl
#SBATCH --time 3:00:00
#SBATCH --mem 8GB
#SBATCH --ntasks 4
#SBATCH --cpus-per-task 8
#SBATCH -e 2-tryF_01-04%j.err
#SBATCH --export NONE

threads=32

module purge
module load Flye/2.9.1-gimkl-2022a-Python-3.10.5

flye --nano-hq ./sample_01.fastq --threads "$threads" --out-dir assembly_01 && cp assembly_01/assembly.fasta assemblies/assembly_01.fasta && cp assembly_01/assembly_graph.gfa assemblies/assembly_01.gfa && rm -r assembly_01

flye --nano-hq ./sample_02.fastq --threads "$threads" --out-dir assembly_02 && cp assembly_02/assembly.fasta assemblies/assembly_02.fasta && cp assembly_02/assembly_graph.gfa assemblies/assembly_02.gfa && rm -r assembly_02

flye --nano-hq ./sample_03.fastq --threads "$threads" --out-dir assembly_03 && cp assembly_03/assembly.fasta assemblies/assembly_03.fasta && cp assembly_03/assembly_graph.gfa assemblies/assembly_03.gfa && rm -r assembly_03

flye --nano-hq ./sample_04.fastq --threads "$threads" --out-dir assembly_04 && cp assembly_04/assembly.fasta assemblies/assembly_04.fasta && cp assembly_04/assembly_graph.gfa assemblies/assembly_04.gfa && rm -r assembly_04
