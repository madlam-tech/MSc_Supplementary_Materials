#!/bin/sh
#!/bin/bash -e
#SBATCH --account massey03345
#SBATCH -J kraken2-2022-08-01
#SBATCH --time 00:10:00
#SBATCH --mem 65GB
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH -e k2.err
#SBATCH -o k12.out
#SBATCH --export NONE

module purge
module load Kraken2/2.1.1-GCC-9.2.0
for filename in $(ls *.fq | sed 's/.fq//')

do 
    kraken2 --db $KRAKEN2_DEFAULT_DB --report ${filename}_kraken2.tax  ${filename}.fq  --output ${filename}_kraken2.txt
    
done 

