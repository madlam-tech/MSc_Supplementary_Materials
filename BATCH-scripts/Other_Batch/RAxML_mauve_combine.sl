#!/bin/bash -e

#SBATCH --account       massey03345
#SBATCH --job-name      RaxML_combine_mauve
#SBATCH --time          1:00:00
#SBATCH --mem           2GB
#SBATCH --cpus-per-task 4
#SBATCH --output 	RaxML_combine_mauve_%j.out
#SBATCH --error		RaxML_combine_mauve_%j.err

module purge
module load RAxML-NG/1.1.0-gimkl-2022a

#!/bin/bash

# Step 1: Combine XMFA files into a single alignment file
cat alignment_1_2.xmfa alignment_1_3.xmfa alignment_1_4.xmfa > combined_alignment.xmfa

# Step 2: (Optional) Clean up the alignment
# Example: Remove poorly aligned regions using Gblocks
# Replace the command with any alignment cleanup tool you prefer
# gblocks combined_alignment.xmfa -t=d -b5=h

# Step 3: Infer a phylogenetic tree from the alignment
# Example: Use RAxML to infer a maximum likelihood tree
raxmlHPC -T 4 -f a -m GTRGAMMA -p 12345 -x 12345 -# 100 -s combined_alignment.xmfa -n combined_tree


