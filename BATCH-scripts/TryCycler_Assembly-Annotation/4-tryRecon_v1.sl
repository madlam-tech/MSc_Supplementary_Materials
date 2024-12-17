#!/bin/bash -e
#SBATCH --account massey03345
#SBATCH -J tryRec_v3
#SBATCH --time 4:00:00
#SBATCH --mem 4GB
#SBATCH --ntasks 4
#SBATCH --cpus-per-task 2
#SBATCH -e 4-tryRec_%j.log
#SBATCH --export NONE

threads=32

module purge
module load Trycycler/0.5.3-gimkl-2022a-Python-3.10.5


trycycler reconcile --reads ./../filtered/*fq --cluster_dir ./clusters/cluster_001
trycycler reconcile --reads ./../filtered/*fq --cluster_dir ./clusters/cluster_001B
trycycler reconcile --reads ./../filtered/*fq --cluster_dir ./clusters/cluster_002
trycycler reconcile --reads ./../filtered/*fq --cluster_dir ./clusters/cluster_003
trycycler reconcile --reads ./../filtered/*fq --cluster_dir ./clusters/cluster_003B
trycycler reconcile --reads ./../filtered/*fq --cluster_dir ./clusters/cluster_004
trycycler reconcile --reads ./../filtered/*fq --cluster_dir ./clusters/cluster_005
trycycler reconcile --reads ./../filtered/*fq --cluster_dir ./clusters/cluster_006
trycycler reconcile --reads ./../filtered/*fq --cluster_dir ./clusters/cluster_008
trycycler reconcile --reads ./../filtered/*fq --cluster_dir ./clusters/cluster_010
trycycler reconcile --reads ./../filtered/*fq --cluster_dir ./clusters/cluster_011
trycycler reconcile --reads ./../filtered/*fq --cluster_dir ./clusters/cluster_017
