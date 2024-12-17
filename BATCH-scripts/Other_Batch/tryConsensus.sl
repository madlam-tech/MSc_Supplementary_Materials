#!/bin/bash -e
#SBATCH --account massey03345
#SBATCH -J tryCons_v1
#SBATCH --time 4:00:00
#SBATCH --mem 6GB
#SBATCH --ntasks 4
#SBATCH --cpus-per-task 2
#SBATCH -e tryCons_v1.err
#SBATCH -o tryCons_v1.out
#SBATCH --export NONE

threads=32

module purge
module load Trycycler/0.5.3-gimkl-2022a-Python-3.10.5

trycycler consensus  --cluster_dir ./cluster_001
trycycler consensus  --cluster_dir ./cluster_002
trycycler consensus  --cluster_dir ./cluster_003
trycycler consensus  --cluster_dir ./cluster_004
trycycler consensus  --cluster_dir ./cluster_005
