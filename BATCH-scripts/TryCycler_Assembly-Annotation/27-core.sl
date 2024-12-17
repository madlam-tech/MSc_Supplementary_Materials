#!/bin/bash -e
#SBATCH --account massey03345
#SBATCH -J 27-core
#SBATCH --time 00:10:00
#SBATCH --mem 1GB
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 4
#SBATCH -e 27-core%j.log

for filename in *.fai

do

BASENAME=$(basename ${filename} .ffn.fai) 

awk '{print $1}' ${filename} > ${BASENAME}_genes.txt

awk 'NR==FNR {a[$1]; next} !($1 in a) {print}' ${BASENAME}.psuedogenes.txt ${BASENAME}_genes.txt > ${BASENAME}_intact.txt

grep -f ${BASENAME}_intact.txt ${BASENAME}.gff > ${BASENAME}.intact.gff

# add headers; ##FASTA; and fastas, e.g.
grep '##' {BASENAME}.gff | sed '$d' > head

echo '##FASTA' >> {BASENAME}.gff

cat head ${BASENAME}.intact.gff {BASENAME}.fsa > temp
mv temp ${BASENAME}.intact.gff

done

