import os
import glob
import csv
from collections import defaultdict

import argparse

parser = argparse.ArgumentParser(description='Prepare OSM Workload.')
parser.add_argument('--indir', required=True, help='Directory containing yearly CSVs from extraction')
parser.add_argument('--outdir', required=True, help='Directory to store the final workload')
parser.add_argument('--base_year', type=int, default=2017, help='Year to construct base snapshot')
args = parser.parse_args()

INPUT_DIR = args.indir
OUTPUT_BASE = args.outdir
BASE_YEAR = args.base_year
os.makedirs(f"{OUTPUT_BASE}/00_build", exist_ok=True)
os.makedirs(f"{OUTPUT_BASE}/01_commits", exist_ok=True)
os.makedirs(f"{OUTPUT_BASE}/02_meta", exist_ok=True)

print("Loading data...")
all_rows = []
for file in glob.glob(f"{INPUT_DIR}/*.csv"):
    with open(file, 'r', encoding='utf-8') as f:
        reader = csv.reader(f)
        next(reader) 
        for r in reader:
            if len(r) >= 8:
                all_rows.append((r[0], r[1], r[2]=='True', r[3], r[4], r[6], r[7]))

print("Sorting globally by timestamp...")
all_rows.sort(key=lambda x: x[4])

print("Computing changeset bounds...")
changeset_meta = {}
for r in all_rows:
    cs = r[3]
    ts = r[4]
    if cs not in changeset_meta:
        changeset_meta[cs] = [ts, ts]
    else:
        if ts < changeset_meta[cs][0]: changeset_meta[cs][0] = ts
        if ts > changeset_meta[cs][1]: changeset_meta[cs][1] = ts

print("Running state machine...")
active_nodes = {}
yearly_stats = defaultdict(lambda: {'I': 0, 'U': 0, 'D': 0, 'Total_Valid': 0})
yearly_commits = defaultdict(list)
build_year = BASE_YEAR

for i, r in enumerate(all_rows):
    node_id, version, visible, cs, ts, lon, lat = r
    year = int(ts[:4])

    op = None
    if node_id not in active_nodes:
        if visible:
            op = 'I'
            # Filter negative coordinates right at the source!
            if float(lon) >= 0 and float(lat) >= 0:
                active_nodes[node_id] = (version, lon, lat, ts)
            else:
                op = None # Discard entirely
    else:
        if visible:
            op = 'U'
            if float(lon) >= 0 and float(lat) >= 0:
                active_nodes[node_id] = (version, lon, lat, ts)
            else:
                op = 'D'
                del active_nodes[node_id]
        else:
            op = 'D'
            del active_nodes[node_id]

    if op:
        yearly_stats[year][op] += 1
        if year > build_year:
            yearly_commits[year].append((
                cs,
                changeset_meta[cs][0],
                changeset_meta[cs][1],
                op,
                node_id,
                version,
                lon if op != 'D' else '',
                lat if op != 'D' else '',
                ts
            ))

    is_last_of_year = (i == len(all_rows)-1) or (int(all_rows[i+1][4][:4]) != year)
    if is_last_of_year:
        yearly_stats[year]['Total_Valid'] = len(active_nodes)
        if year == build_year:
            with open(f"{OUTPUT_BASE}/00_build/base_snapshot_{build_year}.csv", 'w', newline='') as f:
                writer = csv.writer(f)
                writer.writerow(['node_id', 'version', 'lon', 'lat', 'timestamp'])
                for nid, val in active_nodes.items():
                    writer.writerow([nid, val[0], val[1], val[2], val[3]])

print("Writing commits...")
# Remove old commit files to prevent leftover 2007 commits
os.system(f"rm -f {OUTPUT_BASE}/01_commits/commits_*.csv")
for year, commits in yearly_commits.items():
    commits.sort(key=lambda x: (x[1], int(x[0]), x[8]))
    with open(f"{OUTPUT_BASE}/01_commits/commits_{year}.csv", 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['changeset', 'cs_start', 'cs_end', 'op_type', 'node_id', 'version', 'lon', 'lat', 'timestamp'])
        writer.writerows(commits)

print("Data rebuilt successfully with 2007 base!")
