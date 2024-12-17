#!/bin/sh
#!/bin/bash -e
#SBATCH --account massey03345
#SBATCH -J kraken2
#SBATCH --time 00:20:00
#SBATCH --mem 50GB
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH -e k300821.err
#SBATCH -o k300821.out
#SBATCH --export NONE

module purge
module load Kraken2/2.1.1-GCC-9.2.0

for filename in $(ls *.R1.fq | sed 's/.R1.fq//')

do 
    kraken2 --db $KRAKEN2_DEFAULT_DB --report ${filename}_kraken2.report --paired ${filename}.R1.fq ${filename}.R2.fq --output ${filename}_kraken2.txt
    
done 
