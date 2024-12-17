#!/bin/sh
#!/bin/bash -e
#SBATCH --account massey03345
#SBATCH -J guppy-basecaller
#SBATCH --time 04:00:00
#SBATCH --mem 1GB
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH -e guppy.err
#SBATCH -o guppy.out
#SBATCH --export NONE

module purge
module load ont-guppy-gpu/6.2.1


for filename in *fast5
do 
    

guppy_basecaller -i fast5 -s fastq --disable_qscore_filtering --barcode_kits SQK-RBK004 --detect_barcodes --verbose_logs --compress_fastq --fast5_out -c dna_r9.4.1_450bps_sup.cfg

done
