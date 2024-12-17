#!/bin/sh
#!/bin/bash -e
#SBATCH --account massey03345
#SBATCH -J spades1
#SBATCH --time 2:00:00
#SBATCH --mem 10GB
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 12
#SBATCH -e spades1.err
#SBATCH -o spades1.out
#SBATCH --export NONE

export SLURM_EXPORT_ENV=ALL

module purge
module load SPAdes/3.13.1-gimkl-2018b

for filename in $(ls *.bamR1.fq | sed 's/.bamR1.fq//')

do  

    spades.py --meta -k 33,55,77,99,121 -t 12 -1 ${filename}.bamR1.fq -2 ${filename}.bamR2.fq -o ${filename}_spades_assembly/

done 

