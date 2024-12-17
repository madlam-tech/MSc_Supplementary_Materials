#!/bin/sh
#!/bin/bash -e
#SBATCH --account massey03345
#SBATCH -J kraken2
#SBATCH --time 01:00:00
#SBATCH --mem 100GB
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH -e k11.err
#SBATCH -o k11.out
#SBATCH --export NONE

module purge
module load Kraken2/2.1.1-GCC-9.2.0

for filename in $(ls *.bamR1.fq | sed 's/.bamR1.fq//')

do 
    kraken2 --db $KRAKEN2_DEFAULT_DB --threads 16 --fastq --report ${filename}k2.report --paired ${filename}.bamR1.fq ${filename}.bamR2.fq --output ${filename}k22.txt
    
done 
