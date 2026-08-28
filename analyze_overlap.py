import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.dates import DateFormatter
import os

sns.set_theme(style="whitegrid")

# 1. Load Data
data_dir = "/data/bhuan102/SILVA-dataset/bhutan_workload/01_commits/"
dfs = []
for year in range(2018, 2027):
    f = f"{data_dir}/commits_{year}.csv"
    if os.path.exists(f):
        df = pd.read_csv(f, usecols=['changeset', 'cs_start', 'cs_end'])
        dfs.append(df)

if not dfs:
    print("No data found.")
    exit(1)

full_df = pd.concat(dfs, ignore_index=True)
cs_df = full_df.drop_duplicates(subset=['changeset']).copy()

# 2. Parse Datetimes
cs_df['cs_start'] = pd.to_datetime(cs_df['cs_start'], utc=True)
cs_df['cs_end'] = pd.to_datetime(cs_df['cs_end'], utc=True)
cs_df['duration_sec'] = (cs_df['cs_end'] - cs_df['cs_start']).dt.total_seconds()

# 3. Determine "Master Branch" Base
# The "master timeline" is defined by the commit end time (cs_end).
cs_df = cs_df.sort_values('cs_end').reset_index(drop=True)
ref_df = cs_df[['cs_end', 'changeset']].rename(columns={'changeset': 'base_commit_id', 'cs_end': 'ref_time'})

cs_df = cs_df.sort_values('cs_start')
merged = pd.merge_asof(cs_df, ref_df, left_on='cs_start', right_on='ref_time', direction='backward')
merged['base_commit_id'] = merged['base_commit_id'].fillna('Base_2017_Snapshot')

final_df = merged.sort_values('cs_end').reset_index(drop=True)

# 4. Detect Overlaps (Sweeping Line)
overlaps_dict = {row['changeset']: [] for _, row in final_df.iterrows()}
events = []
for _, row in final_df.iterrows():
    events.append((row['cs_start'], 'start', row['changeset']))
    events.append((row['cs_end'], 'end', row['changeset']))

# Tie-breaker: process start before end to catch 0-second overlaps
events.sort(key=lambda x: (x[0], 0 if x[1] == 'start' else 1))

active = set()
concurrency_timeline = []

for t, e_type, cid in events:
    if e_type == 'start':
        for active_id in active:
            overlaps_dict[cid].append(active_id)
            overlaps_dict[active_id].append(cid)
        active.add(cid)
    else:
        if cid in active:
            active.remove(cid)
    concurrency_timeline.append({'time': t, 'active_count': len(active)})

final_df['overlapping_count'] = final_df['changeset'].map(lambda x: len(overlaps_dict[x]))

def truncate_list(lst):
    lst = list(set(lst))
    if not lst: return "None"
    if len(lst) > 3: return f"{lst[:3]} + {len(lst)-3} more"
    return str(lst)

final_df['overlapping_ids'] = final_df['changeset'].map(lambda x: truncate_list(overlaps_dict[x]))

# Output CSV
os.makedirs('verification_results', exist_ok=True)
out_csv = 'verification_results/changeset_overlap_analysis.csv'
final_df[['changeset', 'cs_start', 'cs_end', 'duration_sec', 'base_commit_id', 'overlapping_count', 'overlapping_ids']].to_csv(out_csv, index=False)
print(f"Analysis written to {out_csv}")

# Plotting
out_fig_dir = 'results_real_world/bhutan/step_1000/figures'
os.makedirs(out_fig_dir, exist_ok=True)

# Plot 1: Concurrency Timeline
print("Plotting Concurrency Timeline...")
timeline_df = pd.DataFrame(concurrency_timeline)
plt.figure(figsize=(14, 6))
plt.plot(timeline_df['time'], timeline_df['active_count'], linewidth=1.5, color='purple')
plt.title('OSM Bhutan: Number of Active Concurrent Changesets Over Time')
plt.xlabel('Date')
plt.ylabel('Concurrent Transactions')
plt.tight_layout()
plt.savefig(f'{out_fig_dir}/Concurrency_Timeline.png', dpi=300)
plt.close()

# Plot 2: Gantt Chart for peak window
print("Plotting Gantt Chart...")
max_active_row = timeline_df.loc[timeline_df['active_count'].idxmax()]
peak_time = max_active_row['time']

# Filter window around peak
window_start = peak_time - pd.Timedelta(hours=4)
window_end = peak_time + pd.Timedelta(hours=4)
zoom_df = final_df[(final_df['cs_end'] >= window_start) & (final_df['cs_start'] <= window_end)].copy()
zoom_df = zoom_df.sort_values('cs_start')

if len(zoom_df) > 0:
    plt.figure(figsize=(16, max(6, len(zoom_df)*0.4)))
    for i, (_, row) in enumerate(zoom_df.iterrows()):
        start = row['cs_start']
        end = row['cs_end']
        plt.hlines(y=i, xmin=start, xmax=end, linewidth=4, color='teal')
        plt.scatter([start, end], [i, i], color='black', s=20, zorder=5)
        # Annotate
        plt.text(end, i, f" ID:{row['changeset']} (Base: {row['base_commit_id']})", fontsize=10, va='center')

    plt.title(f'Gantt Chart of Intersecting Changesets (Peak Concurrency near {peak_time.date()})')
    plt.xlabel('Timeline (UTC)')
    plt.ylabel('Changesets (Ordered by Start Time)')
    plt.yticks([])
    plt.grid(axis='x', linestyle='--', alpha=0.6)
    plt.tight_layout()
    plt.savefig(f'{out_fig_dir}/Gantt_Peak_Overlap.png', dpi=300)
    plt.close()

print(f"\n--- OVERLAP STATISTICS ---")
print(f"Total Unique Changesets: {len(final_df)}")
print(f"Max Concurrent Overlaps at a single moment: {timeline_df['active_count'].max()}")
print(f"Average Duration of a Changeset: {final_df['duration_sec'].mean():.2f} seconds")
print(f"Max Duration of a Changeset: {final_df['duration_sec'].max():.2f} seconds")
has_overlap = final_df[final_df['overlapping_count'] > 0]
print(f"Changesets experiencing >= 1 overlap: {len(has_overlap)} ({(len(has_overlap) / len(final_df)) * 100:.2f}%)")
