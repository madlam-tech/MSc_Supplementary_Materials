#!/bin/bash -e
#SBATCH --account massey03345
#SBATCH -J tryMSA_v1
#SBATCH --time 03:00:00
#SBATCH --mem 12GB
#SBATCH --cpus-per-task 32
#SBATCH -e 5-tryMSA_%j.err

threads=32


module purge
module load Trycycler/0.5.3-gimkl-2022a-Python-3.10.5

trycycler msa --cluster_dir ./tryCluster/cluster_001
#trycycler msa --cluster_dir ./tryCluster/cluster_002
#trycycler msa --cluster_dir ./tryCluster/cluster_003
#trycycler msa --cluster_dir ./tryCluster/cluster_004
#trycycler msa --cluster_dir ./tryCluster/cluster_005
#trycycler msa --cluster_dir ./tryCluster/cluster_006
#trycycler msa --cluster_dir ./tryCluster/cluster_007
#trycycler msa --cluster_dir ./tryCluster/cluster_008
#trycycler msa --cluster_dir ./tryCluster/cluster_009
#trycycler msa --cluster_dir ./tryCluster/cluster_010
#trycycler msa --cluster_dir ./tryCluster/cluster_011
#trycycler msa --cluster_dir ./tryCluster/cluster_012