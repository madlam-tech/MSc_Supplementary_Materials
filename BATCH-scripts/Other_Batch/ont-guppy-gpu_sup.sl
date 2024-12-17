#!/bin/sh
#SBATCH --account        massey03345
#SBATCH --job-name       ont-guppy-gpu_job
#SBATCH --gpus-per-node  P100:1
#SBATCH --mem            6G
#SBATCH --cpus-per-task  4
#SBATCH --time           26:00:00
#SBATCH --output         slurmout.%j.out
#SBATCH -e guppy%j.err

module purge
module load ont-guppy-gpu/6.2.1


guppy_basecaller -i  /nesi/nobackup/massey03345/nanopore/2023-03-10/2023-03-10/20230310_2205_MN24479_FAQ90885_71171b83/2023-03-10/fast5 \
-s /nesi/nobackup/massey03345/nanopore/2023-03-10 \
    --config /opt/nesi/CS400_centos7_bdw/ont-guppy-gpu/6.2.1/data/dna_r9.4.1_450bps_sup.cfg \
    --detect_barcodes \
    --barcode_kits SQK-RBK004 \
    --disable_qscore_filtering \
    --device auto  \
    --recursive \
    --records_per_fastq 4000
