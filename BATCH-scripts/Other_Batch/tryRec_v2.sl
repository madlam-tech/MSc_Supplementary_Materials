#!/bin/bash -e
#SBATCH --account massey03345
#SBATCH -J tryRec_2
#SBATCH --time 4:00:00
#SBATCH --mem 4GB
#SBATCH --ntasks 4
#SBATCH --cpus-per-task 8
#SBATCH -e tryRec_v2.err
#SBATCH -o tryRec_v2.out
#SBATCH --export NONE

threads=32


module purge
module load Trycycler/0.5.3-gimkl-2022a-Python-3.10.5


trycycler reconcile --reads ./ma1265_4.fq  --cluster_dir ./cluster_001/1_contigs/*fasta


