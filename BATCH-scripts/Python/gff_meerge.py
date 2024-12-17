import os
import gfftools.io as gio
from collections import defaultdict

def merge_gff_files_by_basename(input_dir, output_dir):
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)

    # Dictionary to hold lists of files with the same basename
    basename_dict = defaultdict(list)

    # Populate the dictionary
    for filename in os.listdir(input_dir):
        if filename.endswith('.gff'):
            basename = '_'.join(filename.split('_')[:-1])
            basename_dict[basename].append(os.path.join(input_dir, filename))

    # Merge files with the same basename
    for basename, file_list in basename_dict.items():
        if len(file_list) > 1:
            merged = gio.read(file_list[0])
            for gff_file in file_list[1:]:
                merged = merged.merge(gio.read(gff_file))
            output_file = os.path.join(output_dir, f"{basename}_merged.gff")
            merged.write(output_file)
            print(f"Merged {len(file_list)} files into {output_file}")
        else:
            print(f"No files to merge for {basename}")

# Define input and output directories
input_dir = '/nesi/nobackup/massey03345/fnas_complete_set/gffs_2_merge'
output_dir = '/nesi/nobackup/massey03345/fnas_complete_set/gffs_merged'

# Call the function to merge files
merge_gff_files_by_basename(input_dir, output_dir)
