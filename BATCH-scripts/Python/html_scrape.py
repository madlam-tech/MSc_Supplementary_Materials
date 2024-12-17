import os
import pandas as pd
from bs4 import BeautifulSoup
import matplotlib.pyplot as plt

# Function to extract tables from an HTML file and save them as CSV
def extract_tables(html_file, output_dir):
    # Read the HTML file
    with open(html_file, 'r') as file:
        soup = BeautifulSoup(file, 'html.parser')
    
    # Find all tables
    tables = soup.find_all('table')
    for i, table in enumerate(tables):
        # Extract table headers
        headers = [header.text.strip() for header in table.find_all('th')]
        rows = table.find_all('tr')
        
        # Extract table rows
        data = []
        for row in rows:
            columns = row.find_all('td')
            if columns:
                data.append([column.text.strip() for column in columns])
        
        # Create a DataFrame
        df = pd.DataFrame(data)
        
        # Add headers if they exist and the number of headers matches the number of columns
        if headers and len(headers) == df.shape[1]:
            df.columns = headers
        
        # Save DataFrame to CSV
        output_csv = os.path.join(output_dir, f"{os.path.basename(html_file).split('.')[0]}_table_{i+1}.csv")
        df.to_csv(output_csv, index=False)
        print(f"Table {i+1} saved to {output_csv}")

        # Check if the DataFrame contains numeric data
        numeric_df = df.apply(pd.to_numeric, errors='coerce')
        if numeric_df.dropna().shape[1] > 0:
            # Generate and save a plot for the table
            numeric_df.plot(kind='bar', figsize=(10, 6))
            plt.title(f"Table {i+1} from {os.path.basename(html_file)}")
            plt.xlabel('Index' if headers == [] else headers[0])
            plt.ylabel('Values')
            plt.xticks(rotation=45)
            plt.tight_layout()
            output_plot = os.path.join(output_dir, f"{os.path.basename(html_file).split('.')[0]}_table_{i+1}.png")
            plt.savefig(output_plot)
            plt.close()
            print(f"Plot for Table {i+1} saved to {output_plot}")
        else:
            print(f"Table {i+1} does not contain numeric data and will not be plotted.")

# List of HTML files to process
html_files = [
    'ma1263_1.html',
    'ma1265.html',
    'ma1263_3.html',
    'ma102.html'
]

# Define the output directory
output_dir = '/Applications/software/html_scrape'

# Ensure the output directory exists
os.makedirs(output_dir, exist_ok=True)

# Process each HTML file
for html_file in html_files:
    extract_tables(html_file, output_dir)
