#!/bin/sh
#!/bin/bash -e
#SBATCH --account massey03345
#SBATCH -J fastqc
#SBATCH --time 01:00:00
#SBATCH --mem 20GB
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8


module purge
module load FastQC/0.11.7

for filename in *.fq

do


	fastqc ${filename}


done
