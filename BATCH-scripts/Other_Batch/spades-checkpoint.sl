#!/bin/sh
#!/bin/bash -e
#SBATCH --account massey03345
#SBATCH -J spades_assembly
#SBATCH --time 00:50:00
#SBATCH --mem 10GB
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 12
#SBATCH -e spades_assembly.err
#SBATCH -o spades_assembly.out
#SBATCH --export NONE

export SLURM_EXPORT_ENV=ALL

module purge
module load SPAdes/3.13.1-gimkl-2018b

spades.py --meta -k 33,55,77,99,121 -t 12 -1 V14.unmapped.R1.fastq -2 V14.unmapped.R2.fastq -o EO1834ass/

