#!/bin/bash -e
#SBATCH --account massey03345
#SBATCH -J bracken-combine
#SBATCH --time 0:30:00
#SBATCH --mem 1GB
#SBATCH --ntasks 2
#SBATCH --cpus-per-task 2
#SBATCH -e bc.err
#SBATCH -o bc.out

module purge
module load Python/3.7.3-gimkl-2018b


python kraken-combine.py --display-headers -k *kr -o kraken-combine_2022-11-09.tsv

