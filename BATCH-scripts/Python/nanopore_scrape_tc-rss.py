import os
import re
import csv

# Define the directory to search
directory = './'

# Define the pattern to search for files
file_pattern = re.compile(r'^1-trysrs.*')

# Define the pattern to scrape the required lines
input_reads_pattern = re.compile(r'^Input reads: (.+)$')
reads_pattern = re.compile(r'^\s+(\d+) reads \([\d,]+ bp\)$')
n50_pattern = re.compile(r'^\s+N50 = ([\d,]+) bp$')
total_read_depth_pattern = re.compile(r'^Total read depth: ([\d.]+)x$')
mean_read_length_pattern = re.compile(r'^Mean read length: ([\d,]+) bp$')

# Prepare to collect data
data = []

# Function to scrape required information from file
def scrape_file(file_path):
    input_reads = ""
    reads = ""
    n50 = ""
    total_read_depth = ""
    mean_read_length = ""
    
    with open(file_path, 'r') as file:
        lines = file.readlines()
        for i, line in enumerate(lines):
            if input_reads_pattern.match(line):
                input_reads = input_reads_pattern.match(line).group(1)
            elif reads_pattern.match(line):
                reads = reads_pattern.match(line).group(1)
            elif n50_pattern.match(line):
                n50 = n50_pattern.match(line).group(1).replace(",", "")
            elif total_read_depth_pattern.match(line):
                total_read_depth = total_read_depth_pattern.match(line).group(1)
            elif mean_read_length_pattern.match(line):
                mean_read_length = mean_read_length_pattern.match(line).group(1).replace(",", "")
    
    data.append([file_path, input_reads, reads, n50, total_read_depth, mean_read_length])

# Walk through directory to find matching files
for root, dirs, files in os.walk(directory):
    for file in files:
        if file_pattern.match(file):
            file_path = os.path.join(root, file)
            scrape_file(file_path)

# Save the collected data to a TSV file
output_file = './nanopore/scraped_tc-rss_data.tsv'
with open(output_file, 'w', newline='') as tsvfile:
    tsvwriter = csv.writer(tsvfile, delimiter='\t')
    tsvwriter.writerow(['File Path', 'Input Reads', 'Reads', 'N50', 'Total Read Depth', 'Mean Read Length'])
    tsvwriter.writerows(data)

print(f"Scraped data saved to {output_file}")
