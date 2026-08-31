import os
import glob
import pandas as pd
from collections import defaultdict
import argparse
import time

parser = argparse.ArgumentParser(description='Prepare OSM Workload (Optimized for Large Memory).')
parser.add_argument('--indir', required=True, help='Directory containing yearly CSVs from extraction')
parser.add_argument('--outdir', required=True, help='Directory to store the final workload')
parser.add_argument('--base_year', type=int, default=2017, help='Year to construct base snapshot')
args = parser.parse_args()

INPUT_DIR = args.indir
OUTPUT_BASE = args.outdir
BASE_YEAR = args.base_year
os.makedirs(f"{OUTPUT_BASE}/00_build", exist_ok=True)
os.makedirs(f"{OUTPUT_BASE}/01_commits", exist_ok=True)

print("1. Loading CSVs in parallel via Pandas...")
start_t = time.time()
files = glob.glob(f"{INPUT_DIR}/*.csv")
dfs = []
for file in files:
    print(f"   Reading {file}...")
    df = pd.read_csv(file, dtype={
        'node_id': 'uint64',
        'version': 'uint32',
        'visible': 'bool',
        'changeset': 'uint64',
        'lon': 'float32',
        'lat': 'float32',
        'timestamp': 'string'
    })
    dfs.append(df)

df = pd.concat(dfs, ignore_index=True)
del dfs # Free memory
print(f"   Loaded {len(df)} rows in {time.time() - start_t:.1f}s.")

print("2. Sorting globally by timestamp...")
start_t = time.time()
df.sort_values('timestamp', inplace=True)
print(f"   Sorted in {time.time() - start_t:.1f}s.")

print("3. Computing changeset bounds...")
start_t = time.time()
cs_bounds = df.groupby('changeset')['timestamp'].agg(['min', 'max']).to_dict(orient='index')
print(f"   Computed bounds for {len(cs_bounds)} changesets in {time.time() - start_t:.1f}s.")

print("4. Running high-speed state machine...")
start_t = time.time()

active_nodes = {}
yearly_commits = defaultdict(list)

# Pre-extract numpy arrays for blazing fast iteration
nodes = df['node_id'].values
versions = df['version'].values
visibles = df['visible'].values
css = df['changeset'].values
lons = df['lon'].values
lats = df['lat'].values
timestamps = df['timestamp'].values

# To determine year boundaries fast
years = df['timestamp'].str.slice(0, 4).astype('uint16').values

n_rows = len(df)
for i in range(n_rows):
    if i % 10000000 == 0 and i > 0:
        print(f"   Processed {i} rows...")
        
    nid = nodes[i]
    ver = versions[i]
    vis = visibles[i]
    cs = css[i]
    ts = timestamps[i]
    lon = lons[i]
    lat = lats[i]
    year = years[i]

    op = None
    if nid not in active_nodes:
        if vis:
            if lon >= 0 and lat >= 0:
                op = 'I'
                active_nodes[nid] = (ver, lon, lat, ts)
    else:
        if vis:
            if lon >= 0 and lat >= 0:
                op = 'U'
                active_nodes[nid] = (ver, lon, lat, ts)
            else:
                op = 'D'
                del active_nodes[nid]
        else:
            op = 'D'
            del active_nodes[nid]

    if op and year > BASE_YEAR:
        bounds = cs_bounds[cs]
        yearly_commits[year].append((
            cs, bounds['min'], bounds['max'], op, nid, ver,
            lon if op != 'D' else '',
            lat if op != 'D' else '',
            ts
        ))

    # Check year boundary to save base snapshot
    is_last_of_year = (i == n_rows - 1) or (years[i+1] != year)
    if is_last_of_year and year == BASE_YEAR:
        print(f"   Writing base snapshot for year {BASE_YEAR}...")
        snap_df = pd.DataFrame.from_dict(active_nodes, orient='index', columns=list(['version', 'lon', 'lat', 'timestamp']))
        snap_df.index.name = 'node_id'
        snap_df.to_csv(f"{OUTPUT_BASE}/00_build/base_snapshot_{BASE_YEAR}.csv")

print(f"   State machine completed in {time.time() - start_t:.1f}s.")

print("5. Writing commits...")
start_t = time.time()
os.system(f"rm -f {OUTPUT_BASE}/01_commits/commits_*.csv")
for year, commits in yearly_commits.items():
    print(f"   Writing year {year} ({len(commits)} ops)...")
    cdf = pd.DataFrame(commits, columns=list(['changeset', 'cs_start', 'cs_end', 'op_type', 'node_id', 'version', 'lon', 'lat', 'timestamp']))
    cdf.sort_values(['cs_start', 'changeset', 'timestamp'], inplace=True)
    cdf.to_csv(f"{OUTPUT_BASE}/01_commits/commits_{year}.csv", index=False)

print(f"   All commits written in {time.time() - start_t:.1f}s.")
print("Data rebuilt successfully (Optimized for High Memory Server)!")
