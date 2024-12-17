#!/bin/sh
#!/bin/bash -e
#SBATCH --account massey03345
#SBATCH -J fq2fa
#SBATCH --time 01:00:00
#SBATCH --mem 10GB
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH -e fq2fa.err
#SBATCH -o fq2fa.out
#SBATCH --export NONE

module purge
module load SAMtools/1.9-GCC-7.4.0
module load seqtk/1.3-gimkl-2018b

for filename in *.fq 
do 
    
    fq2fa --merge -filter -1 ${filename}R1.fq -2 ${filename}R2.fq ${filename}_merged.fna
      
done


