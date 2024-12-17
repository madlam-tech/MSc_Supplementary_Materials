#!/bin/bash -e
#SBATCH --account massey03345
#SBATCH -J tryPart
#SBATCH --time 4:00:00
#SBATCH --mem 12GB
#SBATCH --cpus-per-task 32
#SBATCH -e 6-tryPart%j.err
#SBATCH -o 6-tryPart%j.out
#SBATCH --export NONE

threads=32

module purge
module load Trycycler/0.5.3-gimkl-2022a-Python-3.10.5


trycycler partition --reads ../../filtered/*fq --cluster_dir ./cluster_001
#trycycler partition --reads ../../filtered/*fq --cluster_dir ./cluster_002
#trycycler partition --reads ../../filtered/*fq --cluster_dir ./cluster_003
#trycycler partition --reads ../../filtered/*fq --cluster_dir ./cluster_004
#trycycler partition --reads ../../filtered/*fq --cluster_dir ./cluster_005
#trycycler partition --reads ../../filtered/*fq --cluster_dir ./cluster_006
#trycycler partition --reads ../../filtered/*fq --cluster_dir ./cluster_007
#trycycler partition --reads ../../filtered/*fq --cluster_dir ./cluster_008
#trycycler partition --reads ../../filtered/*fq --cluster_dir ./cluster_009
#trycycler partition --reads ../../filtered/*fq --cluster_dir ./cluster_010
#trycycler partition --reads ../../filtered/*fq --cluster_dir ./cluster_011
#trycycler partition --reads ../../filtered/*fq --cluster_dir ./cluster_012
