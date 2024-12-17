#!/bin/bash -e
#SBATCH --account massey03345
#SBATCH -J tryfRR_v4  
#SBATCH --time 4:00:00
#SBATCH --mem 4GB
#SBATCH --ntasks 4
#SBATCH --cpus-per-task 8
#SBATCH -e tryFRR.err
#SBATCH -o tryFRR.out
#SBATCH --mail-user=madlam@me.com
#SBATCH --mail-type=FAIL
#SBATCH --export NONE

threads=32
rm -r assemblies
mkdir assemblies


module purge
module load Flye/2.9.1-gimkl-2022a-Python-3.10.5


flye --nano-hq ./sample_01.fastq --threads "$threads" --out-dir assembly_01 -i 5 && cp assembly_01/assembly.fasta assemblies/assembly_01.fasta && cp assembly_01/assembly_graph.gfa assemblies/assembly_01.gfa && rm -r assembly_01

flye --nano-hq ./sample_04.fastq --threads "$threads" --out-dir assembly_04 -i 5 && cp assembly_04/assembly.fasta assemblies/assembly_04.fasta && cp assembly_04/assembly_graph.gfa assemblies/assembly_04.gfa && rm -r assembly_04

flye --nano-hq ./sample_07.fastq --threads "$threads" --out-dir assembly_07 -i 5 && cp assembly_07/assembly.fasta assemblies/assembly_07.fasta && cp assembly_07/assembly_graph.gfa assemblies/assembly_07.gfa && rm -r assembly_07

flye --nano-hq ./sample_10.fastq --threads "$threads" --out-dir assembly_10 -i 5 && cp assembly_10/assembly.fasta assemblies/assembly_10.fasta && cp assembly_10/assembly_graph.gfa assemblies/assembly_10.gfa && rm -r assembly_10

module purge
module load Python/3.10.5-gimkl-2022a
module load any2fasta/0.4.2-GCC-9.2.0
module load minimap2/2.24-GCC-9.2.0
module load Raven/1.5.0-GCC-9.2.0
module load Racon/1.5.0-GCC-11.3.0
module load miniasm/0.3-20191007-GCC-11.3.0

./miniasm_and_minipolish.sh ./sample_02.fastq "$threads" > assemblies/assembly_02.gfa && any2fasta assemblies/assembly_02.gfa > assemblies/assembly_02.fasta

./miniasm_and_minipolish.sh ./sample_05.fastq "$threads" > assemblies/assembly_05.gfa && any2fasta assemblies/assembly_05.gfa > assemblies/assembly_05.fasta

./miniasm_and_minipolish.sh ./sample_08.fastq "$threads" > assemblies/assembly_08.gfa && any2fasta assemblies/assembly_08.gfa > assemblies/assembly_08.fasta

./miniasm_and_minipolish.sh ./sample_11.fastq "$threads" > assemblies/assembly_11.gfa && any2fasta assemblies/assembly_11.gfa > assemblies/assembly_11.fasta

raven --threads "$threads" --disable-checkpoints --graphical-fragment-assembly assemblies/assembly_03.gfa ./sample_03.fastq > assemblies/assembly_03.fasta

raven --threads "$threads" --disable-checkpoints --graphical-fragment-assembly assemblies/assembly_06.gfa ./sample_06.fastq > assemblies/assembly_06.fasta

raven --threads "$threads" --disable-checkpoints --graphical-fragment-assembly assemblies/assembly_09.gfa ./sample_09.fastq > assemblies/assembly_09.fasta

raven --threads "$threads" --disable-checkpoints --graphical-fragment-assembly assemblies/assembly_12.gfa ./sample_12.fastq > assemblies/assembly_12.fasta



