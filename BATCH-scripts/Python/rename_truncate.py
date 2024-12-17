import os
import sys

def rename_files(directory):
    for filename in os.listdir(directory):
        if filename.endswith(".gff"):
            parts = filename.split('_')
            if len(parts) > 2:
                new_name = parts[0] + "_" + parts[1] + ".gff"
                os.rename(os.path.join(directory, filename), os.path.join(directory, new_name))
                print(f"Renamed: {filename} -> {new_name}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python rename_files.py <directory>")
    else:
        directory = sys.argv[1]
        if os.path.isdir(directory):
            rename_files(directory)
        else:
            print(f"Error: {directory} is not a valid directory.")
