#!/bin/bash -e
#SBATCH --account massey03345
#SBATCH -J tryRec_v3
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

for dir in cluster_*
do

trycycler reconcile --reads ../*.fq  --cluster_dir ./${dir}/1_contigs/*fasta

done



