#!/bin/sh
#!/bin/bash -e
#SBATCH --account massey03345
#SBATCH -J multiqc
#SBATCH --time 00:10:00
#SBATCH --mem 1GB
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8


module purge
module load MultiQC/1.9-gimkl-2020a-Python-3.8.2

for filename in *.kraken

do
	multiqc ${filename}
done


