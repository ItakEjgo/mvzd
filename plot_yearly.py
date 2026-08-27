import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import argparse

def clean_df(df):
    df.columns = [c.strip() for c in df.columns]
    df['Algo'] = df['Algo'].str.strip()
    return df

def run(year, res_dir):
    out_dir = "plots_and_tables"
    os.makedirs(out_dir, exist_ok=True)
    
    fwd_file = f"{res_dir}/{year}_Forward_Batch.txt"
    bwd_file = f"{res_dir}/{year}_Backward_Batch.txt"
    mid_file = f"{res_dir}/{year}_MidQuery.txt"
    end_file = f"{res_dir}/{year}_EndQuery.txt"
    
    has_fwd = os.path.exists(fwd_file)
    has_bwd = os.path.exists(bwd_file)
    has_mid = os.path.exists(mid_file)
    has_end = os.path.exists(end_file)
    
    # 1. Memory Evolution
    if has_fwd and has_bwd:
        df_fwd = clean_df(pd.read_csv(fwd_file, sep="|", skipinitialspace=True))
        df_bwd = clean_df(pd.read_csv(bwd_file, sep="|", skipinitialspace=True))
        
        plt.figure(figsize=(10, 5))
        algos = df_fwd['Algo'].unique()
        for algo in algos:
            fwd_mem = df_fwd[df_fwd['Algo'] == algo]['Mem_MB'].values
            bwd_mem = df_bwd[df_bwd['Algo'] == algo]['Mem_MB'].values
            mem_series = np.concatenate([fwd_mem, bwd_mem])
            
            linestyle = '--' if 'Rlog' in algo else '-'
            linewidth = 1.5 if 'Rlog' in algo else 2.5
            plt.plot(mem_series, label=algo, linestyle=linestyle, linewidth=linewidth)

        mid_point = len(df_fwd[df_fwd['Algo'] == algos[0]])
        plt.axvline(x=mid_point - 1, color='red', linestyle=':', label='Start Deletion')
        
        plt.title(f"{year} Memory Evolution (Forward Inserts -> Backward Deletes)")
        plt.xlabel("Batch Sequence")
        plt.ylabel("Logical Memory Size (MB)")
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig(f"{out_dir}/{year}_Fig1_Memory_Evolution.png", dpi=300)
        plt.close()

    # 2. Update Breakdown
    if has_fwd:
        df = clean_df(pd.read_csv(fwd_file, sep="|", skipinitialspace=True))
        grouped = df.groupby('Algo')[['Fork_ms', 'Commit_ms', 'Merge_ms']].mean()
        
        fig, ax = plt.subplots(figsize=(10, 6))
        algos = grouped.index
        forks = grouped['Fork_ms'].values
        commits = grouped['Commit_ms'].values
        merges = grouped['Merge_ms'].values
        
        ax.bar(algos, forks, label='Fork', color='#4c72b0')
        ax.bar(algos, commits, bottom=forks, label='Insert/Commit', color='#dd8452')
        ax.bar(algos, merges, bottom=forks+commits, label='Merge', color='#55a868')
        
        plt.title(f"{year} Average Latency Breakdown per Branch (Forward Phase)")
        plt.ylabel("Time (ms)")
        plt.xticks(rotation=45)
        plt.legend()
        plt.grid(axis='y', alpha=0.3)
        plt.tight_layout()
        plt.savefig(f"{out_dir}/{year}_Fig2_Update_Breakdown.png", dpi=300)
        plt.close()

    # 3. Query Latency
    if has_mid:
        df = clean_df(pd.read_csv(mid_file, sep="|", skipinitialspace=True))
        df['QType'] = df['QType'].str.strip()
        grouped = df.groupby(['Algo', 'QType'])['Time_ms'].mean().unstack()
        
        fig, ax = plt.subplots(figsize=(10, 6))
        grouped.plot(kind='bar', ax=ax, width=0.8, colormap='viridis')
        ax.set_yscale('log')
        plt.title(f"{year} Average Query Latency (Log Scale) after Insertions")
        plt.ylabel("Average Time (ms) - Log Scale")
        plt.xticks(rotation=45)
        plt.legend(title="Query Type")
        plt.grid(axis='y', alpha=0.3)
        plt.tight_layout()
        plt.savefig(f"{out_dir}/{year}_Fig3_Query_Latency.png", dpi=300)
        plt.close()

    # 4. Summary Table
    dfs = []
    if has_fwd:
        df_fwd = clean_df(pd.read_csv(fwd_file, sep="|", skipinitialspace=True))
        df_fwd['Inc_Ins_ms'] = df_fwd['Fork_ms'] + df_fwd['Commit_ms'] + df_fwd['Merge_ms']
        dfs.append(df_fwd.groupby('Algo')['Inc_Ins_ms'].mean())
    if has_mid:
        df_mid = clean_df(pd.read_csv(mid_file, sep="|", skipinitialspace=True))
        df_mid['QType'] = df_mid['QType'].str.strip()
        mid_q = df_mid.groupby(['Algo', 'QType'])['Time_ms'].mean().unstack().add_prefix('MidQ_')
        dfs.append(mid_q)
    if has_bwd:
        df_bwd = clean_df(pd.read_csv(bwd_file, sep="|", skipinitialspace=True))
        df_bwd['Inc_Del_ms'] = df_bwd['Fork_ms'] + df_bwd['Commit_ms'] + df_bwd['Merge_ms']
        dfs.append(df_bwd.groupby('Algo')['Inc_Del_ms'].mean())
    if has_end:
        df_end = clean_df(pd.read_csv(end_file, sep="|", skipinitialspace=True))
        df_end['QType'] = df_end['QType'].str.strip()
        end_q = df_end.groupby(['Algo', 'QType'])['Time_ms'].mean().unstack().add_prefix('EndQ_')
        dfs.append(end_q)
        
    if dfs:
        summary = pd.concat(dfs, axis=1).round(3)
        summary.index.name = "Method"
        
        # Format as Markdown
        md_path = f"{out_dir}/{year}_Summary_Table.md"
        with open(md_path, "w") as f:
            f.write(f"# Year {year} Performance Summary\n\n")
            f.write(summary.to_csv(sep="|"))
            
        csv_path = f"{out_dir}/{year}_Summary_Table.csv"
        summary.to_csv(csv_path)
        print(f"Generated successfully in {out_dir}/!")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--year", type=int, required=True)
    parser.add_argument("--res_dir", type=str, default="verification_results")
    args = parser.parse_args()
    run(args.year, args.res_dir)
