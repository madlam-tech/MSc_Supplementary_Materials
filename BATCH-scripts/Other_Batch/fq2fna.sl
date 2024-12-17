#!/bin/sh
#!/bin/bash -e
#SBATCH --account massey03345
#SBATCH -J fq2fa
#SBATCH --time 01:00:00
#SBATCH --mem 10GB
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH -e fq2fna.err
#SBATCH -o fq2fna.out
#SBATCH --export NONE
module purge
module load IDBA-UD/1.1.3-gimkl-2018b

for filename in *.fq

do

  fq2fa --merge ${filename}bamR1.fq ${filename}bamR2.fq ${filename}.fna
  
done
