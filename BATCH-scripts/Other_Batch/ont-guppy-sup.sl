#!/bin/sh
#SBATCH --account        massey03345
#SBATCH --job-name       ont-guppy-gpu_job
#SBATCH --gpus-per-node  P100:1
#SBATCH --mem            6G
#SBATCH --cpus-per-task  4
#SBATCH --time           24:00:00
#SBATCH --output         slurmout.%j.out
#SBATCH -e guppy.err
#SBATCH -o guppy.out


module purge
module load ont-guppy-gpu/6.2.1


guppy_basecaller -i  /nesi/nobackup/massey03345/nanopore/2022-12-01/no_sample/20221201_1833_MN24479_FAU42483_a7b602bc/fast5 \
	-s  /nesi/nobackup/massey03345/nanopore/2022-12-01/no_sample/20221201_1833_MN24479_FAU42483_a7b602bc/fastq/files \
	--config /opt/nesi/CS400_centos7_bdw/ont-guppy-gpu/6.2.1/data/dna_r9.4.1_450bps_sup.cfg \
	--detect_barcodes \
	--barcode_kits SQK-RBK004 \
	--disable_qscore_filtering \
	--device auto  \
	--recursive \
	--records_per_fastq 4000
