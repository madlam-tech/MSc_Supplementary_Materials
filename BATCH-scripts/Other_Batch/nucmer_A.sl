#!/bin/bash
#SBATCH --account massey03345
#SBATCH -J nucmer
#SBATCH --time 4:00:00
#SBATCH --mem 4GB
#SBATCH --ntasks 4
#SBATCH --cpus-per-task 8
#SBATCH -e Nucmer_%j.log

module purge
module load MUMmer/4.0.0rc1-GCCcore-9.2.0


for filename in *fasta
#mkdir -p ${filename}_out
do 

# run nucmer
nucmer -p ${filename}_out/${filename}nucmer_output  CP102429.1.fna ${filename}

# run show-coords to generate a coordinate file
show-coords -rcl ${filename}_out/nucmer_output.delta > ${filename}_out/nucmer_coords.txt

done
