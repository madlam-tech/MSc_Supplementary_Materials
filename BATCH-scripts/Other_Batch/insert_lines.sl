#!/bin/bash -e
#SBATCH –account massey03345
#SBATCH -J insert_top_line
#SBATCH –time 00:05:00
#SBATCH –mem 512mB
#SBATCH –ntasks 1
#SBATCH –cpus-per-task 8
#SBATCH -e top_line.err
#SBATCH -o top_line.out

for filename in *.breport
do

    sed ‘1 i\100\t\0\t\0\t\U\t\0\t\unclassified’ -I ${filename}

done

