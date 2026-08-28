import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import glob
import os

# Configuration
sns.set_theme(style="whitegrid", palette="tab10")
plt.rcParams.update({'font.size': 12})

dataset = "bhutan"
step_size = 1000
base_dir = f"results_real_world/{dataset}/step_{step_size}"
out_dir = os.path.join(base_dir, "figures")
os.makedirs(out_dir, exist_ok=True)

algos = ["MVZD", "CPAMBB", "Rlog_1yr"]

print("Loading Commit Logs...")
commit_dfs = []
for algo in algos:
    f = os.path.join(base_dir, algo, f"CommitLog_{algo}.txt")
    if os.path.exists(f):
        df = pd.read_csv(f, sep="|", skipinitialspace=True)
        df.columns = [c.strip() for c in df.columns]
        df['Algo'] = algo
        df['Commit_Seq'] = np.arange(len(df))
        commit_dfs.append(df)

if commit_dfs:
    commit_df = pd.concat(commit_dfs, ignore_index=True)
else:
    print("No commit logs found!")
    exit(1)

print("Loading Checkpoint Logs...")
query_dfs = []
for algo in algos:
    files = glob.glob(os.path.join(base_dir, algo, "Checkpoint_*.txt"))
    for f in files:
        df = pd.read_csv(f, sep="|", skipinitialspace=True)
        df.columns = [c.strip() for c in df.columns]
        df['Algo'] = algo
        query_dfs.append(df)

if query_dfs:
    query_df = pd.concat(query_dfs, ignore_index=True)
    # Remove Commit_Avg if it accidentally exists (it shouldn't in the new version)
    query_df = query_df[query_df['QType'] != 'Commit_Avg']
else:
    print("No query logs found!")
    exit(1)

# ==========================================
# Plot 1: Commit Time (Spikes & Time Series)
# ==========================================
print("Plotting Commit Times...")
plt.figure(figsize=(14, 6))
# To make it readable, we can plot every 10th commit or use alpha, and log scale
for algo in algos:
    sub = commit_df[commit_df['Algo'] == algo]
    plt.plot(sub['Commit_Seq'], sub['Time_ms'], alpha=0.6, label=algo, linewidth=1)

plt.yscale('log')
plt.ylabel('Commit Time (ms) - Log Scale')
plt.xlabel('Commit Sequence (Timeline)')
plt.title('Commit Latency Evolution (Highlighting Compaction Spikes in RlogTree)')
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(out_dir, "Commit_Latency_Timeline.png"), dpi=300)
plt.close()

# ==========================================
# Plot 2: Commit Time Distribution (Boxplot)
# ==========================================
plt.figure(figsize=(10, 6))
sns.boxplot(data=commit_df, x='Algo', y='Time_ms', showfliers=True)
plt.yscale('log')
plt.ylabel('Commit Time (ms) - Log Scale')
plt.xlabel('Algorithm')
plt.title('Distribution of Commit Latencies (Median vs Tail)')
plt.tight_layout()
plt.savefig(os.path.join(out_dir, "Commit_Latency_Boxplot.png"), dpi=300)
plt.close()

# ==========================================
# Plot 3: Memory Evolution (Total Memory)
# ==========================================
print("Plotting Memory Evolution...")
plt.figure(figsize=(14, 6))
for algo in algos:
    sub = commit_df[commit_df['Algo'] == algo]
    plt.plot(sub['Commit_Seq'], sub['Total_Index_MB'], label=algo, linewidth=2)

plt.ylabel('Total Memory (MB)')
plt.xlabel('Commit Sequence')
plt.title('Total Memory Footprint Over Time (Historical Version Retention)')
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(out_dir, "Memory_Evolution.png"), dpi=300)
plt.close()

# ==========================================
print("Plotting Overhead Ratio...")
plt.figure(figsize=(14, 6))
for algo in algos:
    sub = commit_df[(commit_df['Algo'] == algo) & (commit_df['Delta_Data_MB'] > 0)].copy()
    if not sub.empty:
        sub['Ratio'] = sub['Delta_Index_MB'] / sub['Delta_Data_MB']
        # smoothing to make it readable
        sub['Smooth_Ratio'] = sub['Ratio'].rolling(window=100).mean()
        plt.plot(sub['Commit_Seq'], sub['Smooth_Ratio'], label=algo, linewidth=2)

plt.ylabel('Space Overhead Ratio (Index_Delta / Data_Delta)')
plt.xlabel('Commit Sequence')
plt.title('Memory Overhead of Multi-versioning (100-Commit Moving Average)')
plt.yscale('log')
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(out_dir, "Memory_Overhead_Ratio.png"), dpi=300)
plt.close()

# ==========================================
# ==========================================
plt.close()

# ==========================================
# Plot 5: Query Latency by Type
# ==========================================
print("Plotting Query Latencies...")
# Group by Algorithm and QType, calculate mean and 95th percentile
query_grouped = query_df.groupby(['Algo', 'QType'])['Time_ms'].agg(['mean', lambda x: np.percentile(x, 95)]).reset_index()
query_grouped.rename(columns={'<lambda_0>': 'p95'}, inplace=True)

plt.figure(figsize=(14, 7))
sns.barplot(data=query_df, x='QType', y='Time_ms', hue='Algo', errorbar=('pi', 95), capsize=0.05)
plt.yscale('log')
plt.ylabel('Query Time (ms) - Log Scale')
plt.xlabel('Query Type')
plt.title('Query Latency Comparison (Mean with 95% Confidence Interval)')
plt.tight_layout()
plt.savefig(os.path.join(out_dir, "Query_Latency_Barplot.png"), dpi=300)
plt.close()

# ==========================================
# Plot 6: Query Nodes Touched
# ==========================================
print("Plotting Nodes Touched...")
plt.figure(figsize=(14, 7))
sns.barplot(data=query_df, x='QType', y='Nodes', hue='Algo', errorbar=None)
plt.yscale('log')
plt.ylabel('Nodes Touched / Evaluated (Log Scale)')
plt.xlabel('Query Type')
plt.title('Query Pruning Efficiency (Internal Traversal + Evaluation Overhead)')
plt.tight_layout()
plt.savefig(os.path.join(out_dir, "Query_Nodes_Touched.png"), dpi=300)
plt.close()

print(f"All plots successfully generated in: {out_dir}")
