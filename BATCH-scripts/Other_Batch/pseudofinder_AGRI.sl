#!/bin/bash -e
#SBATCH --account massey03345
#SBATCH -J diamondDB
#SBATCH --time 02:00:00
#SBATCH --mem 10GB
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH -e 23-pseudogenes.log

module load DIAMOND/2.1.6-GCC-11.3.0
module load BLASTDB/2024-04
module load BLAST/2.13.0-GCC-11.3.0
module load Python/3.10.5-gimkl-2022a

conda activate pseudofinder

for filename in *gbk

do

basename=$(basename ${filename} .gbk)
mkdir ${basename}

pseudofinder.py annotate --diamond --skip_makedb -g ${filename} -db  uniprot_sprot.dmnd -op ${basename}

./pseudofinder/pseudofinder.py annotate --diamond --skip_makedb -g ${filename} -db  /nesi/nobackup/massey03345/uniprot/reference.dmnd -op ${basename}

done
