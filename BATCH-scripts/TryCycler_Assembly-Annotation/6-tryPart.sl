#!/bin/bash -e
#SBATCH --account massey03345
#SBATCH -J tryPart06
#SBATCH --time 4:00:00
#SBATCH --mem 2GB
#SBATCH --ntasks 4
#SBATCH --cpus-per-task 2
#SBATCH -e 6-tryPart%j.log
#SBATCH --export NONE

threads=32

module purge
module load Trycycler/0.5.3-gimkl-2022a-Python-3.10.5

trycycler partition --reads ../filtered/BC06.fq --cluster_dirs ./clusters/cluster_001
#trycycler partition --reads ../filtered/BC06.fq --cluster_dirs ./clusters/cluster_002
#trycycler partition --reads ../filtered/BC06.fq --cluster_dirs ./clusters/cluster_003
#trycycler partition --reads ../filtered/BC06.fq --cluster_dirs ./clusters/cluster_004
#trycycler partition --reads ../filtered/BC06.fq --cluster_dirs ./clusters/cluster_005
#trycycler partition --reads ../filtered/BC06.fq --cluster_dirs ./clusters/cluster_006
#trycycler partition --reads ../filtered/BC06.fq --cluster_dirs ./clusters/cluster_007
#trycycler partition --reads ../filtered/BC06.fq --cluster_dirs ./clusters/cluster_008
#trycycler partition --reads ../filtered/BC06.fq --cluster_dirs ./clusters/cluster_009
#trycycler partition --reads ../filtered/BC06.fq --cluster_dirs ./clusters/cluster_010
#trycycler partition --reads ../filtered/BC06.fq --cluster_dirs ./clusters/cluster_011
#trycycler partition --reads ../filtered/BC06.fq --cluster_dirs ./clusters/cluster_012

