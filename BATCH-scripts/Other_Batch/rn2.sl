#!/bin/bash -e
#SBATCH --account massey03345
#SBATCH -J rename
#SBATCH --time 00:02:00
#SBATCH --mem 512mB
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 2
#SBATCH -e r.err
#SBATCH -o r.out
#SBATCH --export NONE

for filename in *_clas.tax”  


 do
   mv -- "$filename" "${filename/unmapped.bamR1_kneaddata__k2_m3_clas.tax”}/k2.tax" ./rename
 done
