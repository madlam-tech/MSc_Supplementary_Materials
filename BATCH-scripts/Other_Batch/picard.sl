#!/bin/bash -e
#SBATCH --account massey03345
#SBATCH -J picard_MDedup
#SBATCH --time 00:02:00
#SBATCH --mem 4gB
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 2
#SBATCH -e pMD_%j.err
module purge
module load picard/2.26.10-Java-11.0.4

for filename in *.bam

do

export _JAVA_OPTIONS=-Djava.io.tmpdir=${TMPDIR}
java -Xmx3g -Djava.io.tmpdir=${TMPDIR} -jar ./picard.jar MarkDuplicates I=${filename} O=${filename}pMD M=metricsMD_jvm.tx

done
