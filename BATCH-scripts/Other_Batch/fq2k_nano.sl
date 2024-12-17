#!/bin/sh
#!/bin/bash -e
#SBATCH --account massey03345
#SBATCH -J kraken2_09
#SBATCH --time 01:00:00
#SBATCH --mem 65GB
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH -e k2_nano_BC09.err
#SBATCH -o k2_nano_BC09.out
#SBATCH --export NONE

module purge
module load Kraken2/2.1.1-GCC-9.2.0

for filename in *.fastq

do 
    kraken2 --db $KRAKEN2_DEFAULT_DB --report ${filename}_kraken2.tax -i ${filename}.fastq --output ${filename}_kraken2.txt
    
done 


