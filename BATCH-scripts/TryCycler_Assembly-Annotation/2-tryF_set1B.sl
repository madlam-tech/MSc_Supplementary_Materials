#!/bin/bash -e
#SBATCH --account massey03345
#SBATCH -J 2-tryF05-08.sl
#SBATCH --time 3:00:00
#SBATCH --mem 8GB
#SBATCH --ntasks 4
#SBATCH --cpus-per-task 8
#SBATCH -e 2-tryF_05-08%j.err
#SBATCH --export NONE

threads=32

module purge
module load Flye/2.9.1-gimkl-2022a-Python-3.10.5

flye --nano-hq ./sample_05.fastq --threads "$threads" --out-dir assembly_05 && cp assembly_05/assembly.fasta assemblies/assembly_05.fasta && cp assembly_05/assembly_graph.gfa assemblies/assembly_05.gfa && rm -r assembly_05

flye --nano-hq ./sample_06.fastq --threads "$threads" --out-dir assembly_06 && cp assembly_06/assembly.fasta assemblies/assembly_06.fasta && cp assembly_06/assembly_graph.gfa assemblies/assembly_06.gfa && rm -r assembly_06

flye --nano-hq ./sample_07.fastq --threads "$threads" --out-dir assembly_07 && cp assembly_07/assembly.fasta assemblies/assembly_07.fasta && cp assembly_07/assembly_graph.gfa assemblies/assembly_07.gfa && rm -r assembly_07

flye --nano-hq ./sample_08.fastq --threads "$threads" --out-dir assembly_08 && cp assembly_08/assembly.fasta assemblies/assembly_08.fasta && cp assembly_08/assembly_graph.gfa assemblies/assembly_08.gfa && rm -r assembly_08
