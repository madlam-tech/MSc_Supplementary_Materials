#!/bin/sh
#!/bin/bash -e
#SBATCH --account massey03345
#SBATCH -J Braken2
#SBATCH --time 01:00:00
#SBATCH --mem 10GB
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH -e MB.err
#SBATCH -o MB.out
#SBATCH --export NONE

module purge
module load Bracken/2.6.2-GCCcore-9.2.0

metabat2 -t 2 -m 1500 \
         -i spades_assembly/spades_assembly.m1000.fna \
         -a metabat.txt \
         -o metabat/metabat

module purge
module load MetaBAT/2.13-GCC-7.4.0

cd /nesi/nobackup/nesi02659/MGSS_U/<YOUR FOLDER>/5.binning/

# Manual specification of files
jgi_summarize_bam_contig_depths --outputDepth metabat.txt sample1.bam sample2.bam sample3.bam sample4.bam

# Wildcard
jgi_summarize_bam_contig_depths --outputDepth metabat.txt sample*.bam

MetaBAT: Metagenome Binning based on Abundance and Tetranucleotide frequency (version 2:GIT-NOTFOUND; 2019-08-28T11:37:17)
by Don Kang (ddkang@lbl.gov), Feng Li, Jeff Froula, Rob Egan, and Zhong Wang (zhongwang@lbl.gov) 

Allowed options:
  -h [ --help ]                     produce help message
  -i [ --inFile ] arg               Contigs in (gzipped) fasta file format [Mandatory]
  -o [ --outFile ] arg              Base file name and path for each bin. The default output is fasta format.
                                    Use -l option to output only contig names [Mandatory].
  -a [ --abdFile ] arg              A file having mean and variance of base coverage depth (tab delimited; 
                                    the first column should be contig names, and the first row will be 
                                    considered as the header and be skipped) [Optional].
  -m [ --minContig ] arg (=2500)    Minimum size of a contig for binning (should be >=1500).
  --maxP arg (=95)                  Percentage of 'good' contigs considered for binning decided by connection
                                    among contigs. The greater, the more sensitive.
  --minS arg (=60)                  Minimum score of a edge for binning (should be between 1 and 99). The 
                                    greater, the more specific.
  --maxEdges arg (=200)             Maximum number of edges per node. The greater, the more sensitive.
  --pTNF arg (=0)                   TNF probability cutoff for building TNF graph. Use it to skip the 
                                    preparation step. (0: auto).
  --noAdd                           Turning off additional binning for lost or small contigs.
  --cvExt                           When a coverage file without variance (from third party tools) is used 
                                    instead of abdFile from jgi_summarize_bam_contig_depths.
  -x [ --minCV ] arg (=1)           Minimum mean coverage of a contig in each library for binning.
  --minCVSum arg (=1)               Minimum total effective mean coverage of a contig (sum of depth over 
                                    minCV) for binning.
  -s [ --minClsSize ] arg (=200000) Minimum size of a bin as the output.
  -t [ --numThreads ] arg (=0)      Number of threads to use (0: use all cores).
  -l [ --onlyLabel ]                Output only sequence labels as a list in a column without sequences.
  --saveCls                         Save cluster memberships as a matrix format
  --unbinned                        Generate [outFile].unbinned.fa file for unbinned contigs
  --noBinOut                        No bin output. Usually combined with --saveCls to check only contig 
                                    memberships
  --seed arg (=0)                   For exact reproducibility. (0: use random seed)
  -d [ --debug ]                    Debug output
  -v [ --verbose ]                  Verbose output

