#SBATCH --account massey03345
#SBATCH -J kraken2
#SBATCH --time 00:05:00
#SBATCH --mem 65GB
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH -e k2.err
#SBATCH -o k12.out
#SBATCH --export NONE

module purge
module load Kraken2/2.1.1-GCC-9.2.0
for filename in $(ls *.bamR1.fq | sed 's/.bamR1.fq//')

do 
    kraken2 --db $KRAKEN2_DEFAULT_DB --report ${filename}_kraken2.tax --paired ${filename}.bamR1.fq ${filename}.bamR2.fq --output ${filename}_kraken2.txt
    
done 

