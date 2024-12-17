#!/bin/bash -e
#SBATCH --account massey03345
#SBATCH -J rename
#SBATCH --time 00:02:00
#SBATCH --mem 512mB
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH -e r.err
#SBATCH -o r.out
#SBATCH --export NONE

#!/bin/bash
for f in *.br; do
    mv -- "$f" "${f%.br}".k2
done
