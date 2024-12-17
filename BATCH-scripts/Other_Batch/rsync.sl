#!/bin/bash -e
#SBATCH --account massey03345
#SBATCH -J rsync
#SBATCH --time 6:00:00
#SBATCH --mem 512m
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 2
#SBATCH -e rsync%j
#SBATCH --export NONE


#!/bin/bash

# set variables for local and remote directories
local_dir="/Volumes/Believer/Nanopore_working"
remote_dir=" mahuika:/nesi/nobackup/massey03345/nanopore"

# use rsync to copy files from local to remote
rsync -avz --progress $local_dir $remote_dir
