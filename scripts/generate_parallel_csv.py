import pandas as pd
import os
import glob
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--dataset_dir', required=True, help="Original workload dir (e.g. dataset/ca_workload)")
parser.add_argument('--batch_size', type=int, default=1000)
args = parser.parse_args()

in_commits_dir = os.path.join(args.dataset_dir, "01_commits")
out_commits_dir = os.path.join(args.dataset_dir, f"03_commits_par_{args.batch_size}")
os.makedirs(out_commits_dir, exist_ok=True)

files = sorted(glob.glob(f"{in_commits_dir}/commits_*.csv"))

global_cs_map = {}
current_batch_id = None
current_batch_count = 0

for file in files:
    print(f"Processing {file}...")
    df = pd.read_csv(file, dtype={'changeset': str})
    
    unique_css = df['changeset'].unique()
    
    for cs in unique_css:
        if current_batch_count == 0:
            current_batch_id = cs
            
        global_cs_map[cs] = current_batch_id
        current_batch_count += 1
        
        if current_batch_count == args.batch_size:
            current_batch_count = 0
            
    df['changeset'] = df['changeset'].map(global_cs_map)
    
    out_file = os.path.join(out_commits_dir, os.path.basename(file))
    df.to_csv(out_file, index=False)

print(f"Parallel commits generated in {out_commits_dir} with batch size {args.batch_size}.")
