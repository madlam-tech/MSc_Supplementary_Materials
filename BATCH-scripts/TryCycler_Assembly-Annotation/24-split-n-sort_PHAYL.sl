#!/bin/bash -e
#SBATCH --account massey03345
#SBATCH -J split-n-sort
#SBATCH --time 00:45:00
#SBATCH --mem 1GB
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH -e 25-split-n-sort_%j.log

for filename in *.gfff; 
do BASENAME=$(basename ${filename} _pseudos.gfff); grep PHAYL ${filename} | sort -k1,1 -k4,4n | awk -F "\t" '{split($9,a,";"); split(a[3],b,"="); split(a[1],c,":"); print b[2]}' OFS="\t" > ${BASENAME}.psuedogenes.txt;  sed 's/,/\n/g' ${BASENAME}.psuedogenes.txt | grep PHAYL | sort > temp; mv temp ${BASENAME}.psuedogenes.txt;  done
