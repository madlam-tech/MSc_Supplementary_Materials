#!/bin/sh
#!/bin/bash -e
#SBATCH --account massey03345
#SBATCH -J fqc2mqc
#SBATCH --time 01:00:00
#SBATCH --mem 2gB
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH -e fqc2mqc.err
#SBATCH -o fqc2mqc.out

module purge
module load MultiQC/1.9-gimkl-2020a-Python-3.8.2

multiqc .