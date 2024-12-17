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


module purge

module load Trycycler/0.5.3-gimkl-2022a-Python-3.10.5



trycycler reconcile --reads ./ma1265_3_BC06_07_09/ma1265_3.fq  --cluster_dir ./ma1265_3_BC06_07_09/tryclust/cluster_001
trycycler reconcile --reads ./ma1265_3_BC06_07_09/ma1265_3.fq  --cluster_dir ./ma1265_3_BC06_07_09/tryclust/cluster_002
trycycler reconcile --reads ./ma1265_3_BC06_07_09/ma1265_3.fq  --cluster_dir ./ma1265_3_BC06_07_09/tryclust/cluster_003


