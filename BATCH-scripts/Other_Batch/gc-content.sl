#!/bin/bash -e
#SBATCH --account massey03345
#SBATCH -J GC-content
#SBATCH --time 00:15:00
#SBATCH --mem 1GB
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH -e GC_%j.log



for FILE in *.fna; do
        echo "$FILE"
        var1=$(grep -v ">" `echo "$FILE"` | tr -d -c GCgc | wc -c)
        var2=$(grep -v ">" `echo "$FILE"` | tr -d -c ATGCatgc | wc -c)
        echo "$var1, $var2"
        awk -v v1=$var1 -v v2=$var2 'BEGIN {print v1/v2}'
done
