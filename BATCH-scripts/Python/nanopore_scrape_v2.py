import os
import re
import csv

# Define the directory to search
directory = '/Volumes/Believer'

# Define the patterns to search for files
file_patterns = [re.compile(r'^1-trysrs.*', re.IGNORECASE), re.compile(r'^tc-rss.*', re.IGNORECASE)]

# Define the pattern to scrape the required lines
input_reads_pattern = re.compile(r'^Input reads: (.+)$')
reads_pattern = re.compile(r'^\s*(\d+,\d+) reads \([\d,]+ bp\)$')
n50_pattern = re.compile(r'^\s*N50 = ([\d,]+) bp$')
total_read_depth_pattern = re.compile(r'^Total read depth: ([\d.]+)x$')
mean_read_length_pattern = re.compile(r'^Mean read length: ([\d,]+) bp$')
target_pattern = re.compile(r'^  target: ([\d,]+ bp)$')
keeping_pattern = re.compile(r'^  keeping ([\d,]+ bp)$')

# Prepare to collect data
data = []

# Function to scrape required information from file
def scrape_file(file_path):
    target = ""
    keeping = ""
    input_reads = ""
    reads = ""
    n50 = ""
    total_read_depth = ""
    mean_read_length = ""
    
    with open(file_path, 'r') as file:
        lines = file.readlines()
        for line in lines:
            if target_pattern.match(line) and not target:
                target = target_pattern.match(line).group(1)
            elif keeping_pattern.match(line) and not keeping:
                keeping = keeping_pattern.match(line).group(1)
            elif input_reads_pattern.match(line) and not input_reads:
                input_reads = input_reads_pattern.match(line).group(1)
            elif reads_pattern.match(line) and not reads:
                reads = reads_pattern.match(line).group(1)
            elif n50_pattern.match(line) and not n50:
                n50 = n50_pattern.match(line).group(1)
            elif total_read_depth_pattern.match(line) and not total_read_depth:
                total_read_depth = total_read_depth_pattern.match(line).group(1)
            elif mean_read_length_pattern.match(line) and not mean_read_length:
                mean_read_length = mean_read_length_pattern.match(line).group(1)
    
    data.append([file_path, target, keeping, input_reads, reads, n50, total_read_depth, mean_read_length])

# Walk through directory to find matching files
for root, dirs, files in os.walk(directory):
    for file in files:
        for pattern in file_patterns:
            if pattern.match(file):
                file_path = os.path.join(root, file)
                scrape_file(file_path)

# Save the collected data to a CSV file
output_file = '/Volumes/Believer/scraped_tc-rss_data.csv'
with open(output_file, 'w', newline='') as csvfile:
    csvwriter = csv.writer(csvfile)
    csvwriter.writerow(['File Path', 'Target', 'Keeping', 'Input Reads', 'Reads', 'N50', 'Total Read Depth', 'Mean Read Length'])
    csvwriter.writerows(data)

print(f"Scraped data saved to {output_file}")
