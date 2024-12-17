#!/bin/sh
#!/bin/bash -e
#SBATCH --account massey03345
#SBATCH -J bracken-species
#SBATCH --time 00:10:00
#SBATCH --mem 5GB
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH -e bracken_s.err
#SBATCH -o bracken_s.out
#SBATCH --export NONE

module purge
module load  Bracken/2.7-GCC-11.3.0

for filename in *.tax

do 

    bracken -d  /opt/nesi/db/Kraken2/standard-2022-07/ -i ${filename} -l S -o ${filename}.bracken_species
    
done 