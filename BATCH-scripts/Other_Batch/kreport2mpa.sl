#!/bin/bash -e
#SBATCH --account massey03345
#SBATCH -J llu_kreport2krone.py
#SBATCH --time 00:15:00
#SBATCH --mem 2GB
#SBATCH --ntasks 2
#SBATCH --cpus-per-task 2
#SBATCH -e k2mpa_tax.err
#SBATCH -o k2mpa_tax.out


module purge
module load Python/3.10.5-gimkl-2022a

for filename in *txt
do

python  kreport2mpa.py -r  ${filename} -o ${filename}.mpa

done
