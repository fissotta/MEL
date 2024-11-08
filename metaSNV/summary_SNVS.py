# generated by ChatGPT 4.0, prompted by Francisco Issotta
# usage example: Run the script from any directory, and it will output a consolidated TSV file in the calling directory.

import os
import pandas as pd
import glob
import re

# Define the pattern to find your files
folder_path = "*mSNV/filtered/LB/*.freq"

# Collect all relevant files
files = glob.glob(folder_path, recursive=True)

# Initialize a dictionary to store results and a set for unique headers
results = {}
all_headers = set()

if not files:
    print("No files found. Please check the folder path and file pattern.")
else:
    for file_path in files:
        print(f"Processing file: {file_path}")  # Debug print
        # Get the identifier for the file from its path
        file_name = os.path.basename(file_path).replace(".freq", "")
        
        with open(file_path) as f:
            # Read headers (second line)
            headers_line = f.readline().strip()
            headers = headers_line.split("\t")
            headers = [header for header in headers if header.endswith(".bam") or header.endswith(".sorted.bam")]
            
            # Extract unique site IDs (everything before `__`) for consolidated headers
            site_ids = [re.match(r"^(.*?)__", header).group(1) for header in headers if re.match(r"^(.*?)__", header)]
            all_headers.update(site_ids)
            
            # Load data and skip header rows
            data = pd.read_csv(f, sep="\t", skiprows=1, header=None)
            print(f"Data shape: {data.shape}")  # Debug print
            
            if data.empty:
                print(f"No data found in {file_path}. Skipping file.")
                continue
            
            # Calculate sum of values greater than 0.0 for each column
            sums = [data[col].apply(lambda x: x if x > 0 else 0).sum() for col in range(2, data.shape[1])]
            
            # Store sums in a dictionary with site IDs as keys
            results[file_name] = dict(zip(site_ids, sums))

# Create a sorted list of all unique site IDs
sorted_headers = sorted(all_headers)

# Write results to TSV format
output_file_path = os.path.join(os.getcwd(), "consolidated_summary_output.tsv")
with open(output_file_path, "w") as out_file:
    # Write header row
    out_file.write("Reference\t" + "\t".join(sorted_headers) + "\n")
    
    # Write each file's results, filling missing headers with 0
    for file_name, file_data in results.items():
        row = [file_name]
        for header in sorted_headers:
            row.append(str(file_data.get(header, 0)))  # Use 0 if header is missing in file_data
        out_file.write("\t".join(row) + "\n")

print(f"Consolidated TSV output saved to {output_file_path}")