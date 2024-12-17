#!/bin/bash -e
#SBATCH --account       massey03345
#SBATCH --job-name      CheckM
#SBATCH --time          01:00:00
#SBATCH --mem           40GB
#SBATCH --cpus-per-task 10
#SBATCH --error         11-checkM%x_%j.log
#SBATCH --output        11-checkM%x_%j.out

# Load modules
module purge
module load CheckM/1.2.1-gimkl-2022a-Python-3.10.5


mkdir -p ./checkM


# Run CheckM
checkm lineage_wf --genes -t $SLURM_CPUS_PER_TASK --pplacer_threads $SLURM_CPUS_PER_TASK \
                  -x faa --tab_table -f ./checkM/checkm.txt \
                  ./predictions/ ./checkM/


done
