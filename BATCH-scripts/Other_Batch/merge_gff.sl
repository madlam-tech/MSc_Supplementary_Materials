#!/bin/bash -e
#SBATCH --account massey03345
#SBATCH --job-name merge_annotations
#SBATCH --time 1:00:00
#SBATCH --mem 8GB
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 4
#SBATCH -e merge_annotations_%j.err
#SBATCH -o merge_annotations_%j.log

# Load required modules
module load R/4.0.0

# Create an R script to perform the merge
cat << 'EOF' > merge_annotations.R
# Load required packages
library(dplyr)
library(readr)
library(stringr)

# Define file paths for Prokka, Pseudofinder, and Barrnap rRNA data
prokka_dir <- "/nesi/nobackup/massey03345/fnas_complete_set/prokka/"
pseudofinder_dir <- "/nesi/nobackup/massey03345/fnas_complete_set/pseudogene_output_2024-06-13/"
barrnap_dir <- "/nesi/nobackup/massey03345/fnas_complete_set/barrnap_output_20240605/"

# Function to read Prokka data
read_prokka <- function(file_path) {
  tryCatch({
    read_tsv(file_path, show_col_types = FALSE) %>%
      mutate(Accession = str_extract(basename(file_path), "^[^_]+"))
  }, error = function(e) {
    cat("Error reading Prokka file:", file_path, "\n", e, "\n")
    return(NULL)
  })
}

# Function to read Pseudofinder data
read_pseudofinder <- function(file_path) {
  tryCatch({
    read_tsv(file_path, show_col_types = FALSE) %>%
      mutate(Accession = str_extract(basename(file_path), "^[^_]+"))
  }, error = function(e) {
    cat("Error reading Pseudofinder file:", file_path, "\n", e, "\n")
    return(NULL)
  })
}

# Function to read Barrnap rRNA data
read_barrnap <- function(file_path) {
  tryCatch({
    read_tsv(file_path, show_col_types = FALSE) %>%
      mutate(Accession = str_extract(basename(file_path), "^[^_]+"))
  }, error = function(e) {
    cat("Error reading Barrnap file:", file_path, "\n", e, "\n")
    return(NULL)
  })
}

# List Prokka, Pseudofinder, and Barrnap files
prokka_files <- list.files(prokka_dir, pattern = "\\.tsv$", full.names = TRUE)
pseudofinder_files <- list.files(pseudofinder_dir, pattern = "\\.tsv$", full.names = TRUE)
barrnap_files <- list.files(barrnap_dir, pattern = "\\.tsv$", full.names = TRUE)

# Read all Prokka, Pseudofinder, and Barrnap data
prokka_data <- bind_rows(lapply(prokka_files, read_prokka))
pseudofinder_data <- bind_rows(lapply(pseudofinder_files, read_pseudofinder))
barrnap_data <- bind_rows(lapply(barrnap_files, read_barrnap))

# Merge Prokka, Pseudofinder, and Barrnap data
merged_data <- prokka_data %>%
  full_join(pseudofinder_data, by = "Accession") %>%
  full_join(barrnap_data, by = "Accession")

# View merged data
print(merged_data)

# Save the merged data to a file
write_tsv(merged_data, "/nesi/nobackup/massey03345/fnas_complete_set/merged_prokka_pseudofinder_barrnap_data.tsv")

# Example: How to perform further analysis on the merged data
# Calculate the number of pseudogenes for each genome
pseudogene_counts <- merged_data %>%
  group_by(Accession) %>%
  summarise(pseudogene_count = n())

print(pseudogene_counts)
EOF

# Run the R script
Rscript merge_annotations.R
