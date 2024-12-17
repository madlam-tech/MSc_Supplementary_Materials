#!/bin/bash -e

#SBATCH --account	massey03345
#SBATCH --job-name      spades_BC01
#SBATCH --time          02:00:00
#SBATCH --mem           20GB
#SBATCH --cpus-per-task 12
#SBATCH --error         spades_assembly.err
#SBATCH --output        spades_assembly.out

module purge
module load SPAdes/3.15.4-gimkl-2022a-Python-3.10.5

spades.py -k 33,55,77,99,121 -t 12 *fastq_runid  -o spades_assembly/

