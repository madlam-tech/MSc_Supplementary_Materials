#!/bin/sh
#!/bin/bash -e
#SBATCH --account massey03345
#SBATCH -J basename
#SBATCH --time 00:05:00
#SBATCH --mem 2GB
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8

for filename in *.fq
do
    name=$(basename ${filename} .fq)
    echo ${name}
done
