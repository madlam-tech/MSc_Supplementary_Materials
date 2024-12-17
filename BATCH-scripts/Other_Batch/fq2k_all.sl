#!/bin/bash -e
#SBATCH --account massey03345
#SBATCH -J kraken2-2022-10-11
#SBATCH --time 00:10:00
#SBATCH --mem 65GB
#SBATCH --ntasks 2
#SBATCH --cpus-per-task 2
#SBATCH -e k2EO1219_all.err
#SBATCH -o k2EO1219_all.out
#SBATCH --export NONE

module purge
module load Kraken2/2.1.1-GCC-9.2.0

	kraken2 --db $KRAKEN2_DEFAULT_DB \
        --report EO1219_all_k2.tax \
        --report-minimizer-data \
        --minimum-hit-groups 3  \
        --paired EO1219unmapped.bamR1.fq EO1219unmapped.bamR1.fq\
        --output $EO1219_all_kraken2.txt
