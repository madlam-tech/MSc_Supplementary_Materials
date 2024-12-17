#!/bin/bash -e
#SBATCH --account massey03345
#SBATCH -J rename
#SBATCH --time 00:02:00
#SBATCH --mem 512mB
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH -e r.err
#SBATCH -o r.out
#SBATCH --export NONE




for file_name in *R2fq
do 
  new_file_name=$(sed 's/.R2fq/_R2.fq/g' <<< "$file_name");
  mv "$file_name" "$new_file_name";
done
