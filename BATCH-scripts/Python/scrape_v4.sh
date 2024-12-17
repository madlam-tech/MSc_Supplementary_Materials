#!/bin/bash

# Directory containing log .txt files
directory="./pseudogene_output/"

# Output table file
output_table="pseudofinder_statistics_summary.tsv"
echo -e "File\tInitial ORFs\tInitial Pseudogenes\tNumber of Contigs\tInitial ORFs Joined\tPseudogenes Total\tPseudogenes Too Short\tPseudogenes Too Long\tPseudogenes Fragmented\tPseudogenes No Predicted ORF\tPseudogenes High dN/dS\tPseudogenes Frameshift\tPseudogenes Missing Start Codon\tPseudogenes Missing Stop Codon\tPseudogenes Internal Stop Codon\tPseudogenes Multiple Issues\tIntact Genes" > "$output_table"

# Loop through all log txt files in the specified directory
for log_file in "${directory}"*_log.txt; do
    # Skip if no files are found
    [[ -e "$log_file" ]] || continue

    # Extract the base name of the log file for display in the table
    base_name=$(basename "$log_file" _log.txt)

    # Use awk to extract lines between '####### Statistics #######' and the next '#######'
    stats_block=$(awk '/####### Statistics #######/,/####### [^S]/ {if (!/####### Statistics #######/ && !/####### [^S]/) print}' "$log_file")

    # Debug: Output the stats block to verify correct extraction
    echo "Debug - Stats block for $base_name:"
    echo "$stats_block"

    # Extract each needed value
    initial_orfs=$(echo "$stats_block" | awk '/Initial ORFs:/ {print $3}')
    initial_pseudogenes=$(echo "$stats_block" | awk '/Initial pseudogenes:/ {print $3}')
    num_contigs=$(echo "$stats_block" | awk '/Number of contigs:/ {print $4}')
    initial_orfs_joined=$(echo "$stats_block" | awk '/Inital ORFs joined:/ {print $4}')
    pseudogenes_total=$(echo "$stats_block" | awk '/Pseudogenes \(total\):/ {print $3}')
    pseudogenes_too_short=$(echo "$stats_block" | awk '/Pseudogenes \(too short\):/ {print $4}')
    pseudogenes_too_long=$(echo "$stats_block" | awk '/Pseudogenes \(too long\):/ {print $4}')
    pseudogenes_fragmented=$(echo "$stats_block" | awk '/Pseudogenes \(fragmented\):/ {print $3}')
    pseudogenes_no_predicted_orf=$(echo "$stats_block" | awk '/Pseudogenes \(no predicted ORF\):/ {print $6}')
    pseudogenes_high_dn_ds=$(echo "$stats_block" | awk '/Pseudogenes \(high dN\/dS\):/ {print $4}')
    pseudogenes_frameshift=$(echo "$stats_block" | awk '/Pseudogenes \(frameshift\):/ {print $3}')
    pseudogenes_missing_start_codon=$(echo "$stats_block" | awk '/Pseudogenes \(missing start codon\):/ {print $5}')
    pseudogenes_missing_stop_codon=$(echo "$stats_block" | awk '/Pseudogenes \(missing stop codon\):/ {print $5}')
    pseudogenes_internal_stop_codon=$(echo "$stats_block" | awk '/Pseudogenes \(internal stop codon\):/ {print $5}')
    pseudogenes_multiple_issues=$(echo "$stats_block" | awk '/Pseudogenes \(multiple issues\):/ {print $5}')
    intact_genes=$(echo "$stats_block" | awk '/Intact genes:/ {print $3}')

    # Write the extracted data to the output table
    echo -e "$base_name\t$initial_orfs\t$initial_pseudogenes\t$num_contigs\t$initial_orfs_joined\t$pseudogenes_total\t$pseudogenes_too_short\t$pseudogenes_too_long\t$pseudogenes_fragmented\t$pseudogenes_no_predicted_orf\t$pseudogenes_high_dn_ds\t$pseudogenes_frameshift\t$pseudogenes_missing_start_codon\t$pseudogenes_missing_stop_codon\t$pseudogenes_internal_stop_codon\t$pseudogenes_multiple_issues\t$intact_genes" >> "$output_table"
done

# Feedback
echo "Summary table created at ${output_table}"
