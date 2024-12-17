#!/bin/sh
#!/bin/bash -e
#SBATCH --account massey03345
#SBATCH -J bracken
#SBATCH --time 00:04:00
#SBATCH --mem 10GB
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH -e bracken.err
#SBATCH -o bracken.out
#SBATCH --export NONE

module purge
module load Bracken/2.6.2-GCCcore-9.2.0

for filename in *.report

do 

    bracken -d /opt/nesi/db/Kraken2/standard-2018-09 -i ${filename} -l S -o ${filename}.bracken
    
done 

