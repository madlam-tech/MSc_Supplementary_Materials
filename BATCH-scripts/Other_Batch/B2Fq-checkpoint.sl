#!/bin/sh
#!/bin/bash -e
#SBATCH --account massey03345
#SBATCH -J bam-2-fastq
#SBATCH --time 01:00:00
#SBATCH --mem 20GB
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH -e b2s2f.err
#SBATCH -o b2s2f.out
#SBATCH --export NONE

module purge
module load SAMtools/1.9-GCC-7.4.0
module load seqtk/1.3-gimkl-2018b

for filename in *.bam 
do 
    
    samtools bam2fq -1 ${filename}R1.fq -2 ${filename}R2.fq -0 /dev/null -s /dev/null -n -F 0x900 ${filename}
      
done



