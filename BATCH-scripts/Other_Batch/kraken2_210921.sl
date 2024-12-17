#!/bin/sh
#!/bin/bash -e
#SBATCH --account massey03345
#SBATCH -J kraken2
#SBATCH --time 01:30:00
#SBATCH --mem 40GB
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH -e k210921.err
#SBATCH -o k210921.out
#SBATCH --export NONE

module purge
module load Kraken2/2.1.1-GCC-9.2.0

for filename in $(ls *.bamR1.fq | sed 's/.bamR1.fq//')

do 
    kraken2 --use-mpa-style --db $KRAKEN2_DEFAULT_DB --use-names --report ${filename}_k2.report --paired ${filename}.bamR1.fq ${filename}.bamR2.fq --output ${filename}k2.tax
    
done 
