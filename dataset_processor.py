import os
import glob
import re

def process_datasets(directory):
    # Find all years
    files = glob.glob(os.path.join(directory, "*_branch_*.txt"))
    years = set()
    for f in files:
        basename = os.path.basename(f)
        match = re.match(r'^(\d{4})_branch_', basename)
        if match:
            years.add(match.group(1))
            
    for year in sorted(list(years)):
        print(f"Processing year {year}...")
        year_files = glob.glob(os.path.join(directory, f"{year}_branch_*.txt"))
        
        # Sort by branch number
        def get_branch_num(filename):
            match = re.search(r'_branch_(\d+)_', filename)
            return int(match.group(1)) if match else 999
            
        year_files.sort(key=get_branch_num)
        
        merged_file = os.path.join(directory, f"{year}_merged.txt")
        with open(merged_file, 'w') as outfile:
            for f in year_files:
                with open(f, 'r') as infile:
                    outfile.write(infile.read())
                    # Ensure newline at end of file to prevent merging issues
                    outfile.write("\n")
        print(f"  Merged {len(year_files)} files into {merged_file}")

if __name__ == "__main__":
    process_datasets("dataset/bhutan_evolution")
