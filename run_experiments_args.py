import glob
import re
import os
import shutil
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import logging
import time
import subprocess
from matplotlib.ticker import MaxNLocator

# Removed global PARLAY_NUM_THREADS = '1' to allow per-algorithm control

logger = logging.getLogger("Benchmark")
logger.setLevel(logging.INFO)
formatter = logging.Formatter('[%(asctime)s] %(message)s', datefmt='%Y-%m-%d %H:%M:%S')

ch = logging.StreamHandler()
ch.setFormatter(formatter)
logger.addHandler(ch)

fh = logging.FileHandler("benchmark_progress.log", mode='a')
fh.setFormatter(formatter)
logger.addHandler(fh)
years = [2021, 2022, 2023, 2024, 2025]
query_types = ["KNN_1", "KNN_10", "KNN_100", "Range_Small", "Range_Med", "Range_Large"]

def clean_df(df):
    df.columns = [c.strip() for c in df.columns]
    if 'Algo' in df.columns:
        df['Algo'] = df['Algo'].astype(str).str.strip()
    return df

def load_combined_data(year, res_dir, phase):
    df_list = []
    for algo in algos:
        if phase == "BaseQuery":
            file_path = os.path.join(res_dir, f"BaseQuery_{algo}.txt")
        else:
            file_path = os.path.join(res_dir, f"{year}_{phase}_{algo}.txt")
        if os.path.exists(file_path):
            try:
                df = pd.read_csv(file_path, sep="|", skipinitialspace=True)
                
                df_list.append(df)
            except: pass
    if df_list:
        return clean_df(pd.concat(df_list, ignore_index=True))
    return pd.DataFrame()

def generate_global_plots(bp, res_dir):
    out_dir = f"plots_and_tables/bp_{bp}"
    os.makedirs(out_dir, exist_ok=True)
    
    # Timeline sequences
    fwd_years = [2021, 2022, 2023, 2024, 2025]
    bwd_years = [2025, 2024, 2023, 2022, 2021]
    
    # 1. Memory Evolution
    plt.figure(figsize=(16, 6))
    
    # Calculate common X coordinates for tight year grouping
    global_x_coords = []
    global_boundaries = []
    curr_x = 0.0
    
    for y in fwd_years:
        df = load_combined_data(y, res_dir, "Forward_Batch")
        if not df.empty:
            vals = df[(df['Algo'] == algos[0]) & (df['Batch'].astype(str).str.strip() != 'SUMMARY')]['Mem_MB'].values
            num_batches = len(vals)
            year_x = []
            for _ in range(num_batches):
                year_x.append(curr_x)
                curr_x += 1.0
            global_x_coords.extend(year_x)
            last_x = year_x[-1] if year_x else curr_x
            global_boundaries.append(last_x + 0.2)
            curr_x = last_x + 0.4
            
    pivot = global_boundaries[-1] if global_boundaries else 0
    
    for y in bwd_years:
        df = load_combined_data(y, res_dir, "Backward_Batch")
        if not df.empty:
            vals = df[(df['Algo'] == algos[0]) & (df['Batch'].astype(str).str.strip() != 'SUMMARY')]['Mem_MB'].values
            num_batches = len(vals)
            year_x = []
            for _ in range(num_batches):
                year_x.append(curr_x)
                curr_x += 1.0
            global_x_coords.extend(year_x)
            last_x = year_x[-1] if year_x else curr_x
            global_boundaries.append(last_x + 0.2)
            curr_x = last_x + 0.4

    for algo in algos:
        mem_series = []
        for y in fwd_years:
            df = load_combined_data(y, res_dir, "Forward_Batch")
            if not df.empty and algo in df['Algo'].values:
                vals = df[(df['Algo'] == algo) & (df['Batch'].astype(str).str.strip() != 'SUMMARY')]['Mem_MB'].values
                mem_series.extend(vals)
                
        for y in bwd_years:
            df = load_combined_data(y, res_dir, "Backward_Batch")
            if not df.empty and algo in df['Algo'].values:
                vals = df[(df['Algo'] == algo) & (df['Batch'].astype(str).str.strip() != 'SUMMARY')]['Mem_MB'].values
                mem_series.extend(vals)
                
        if len(mem_series) > 0 and len(mem_series) == len(global_x_coords):
            linestyle = '--' if 'Rlog' in algo else '-'
            linewidth = 1.5 if 'Rlog' in algo else 2.5
            plt.plot(global_x_coords, mem_series, label=algo, linestyle=linestyle, linewidth=linewidth, marker='.')
            
    # Draw boundaries, custom X-ticks and backgrounds
    xticks_pos = []
    xtick_labels = []
    if len(global_boundaries) > 0:
        prev_b = -0.2
        plt.gca().set_xticks(global_x_coords, minor=True)
        
        year_colors = ['#e8f4f8', '#fff4e6', '#e6ffe6', '#f2e6ff', '#ffe6e6'] 
        
        for i, b in enumerate(global_boundaries):
            year_idx = i if i < 5 else 9 - i
            bg_color = year_colors[year_idx]
            
            if i < 5:
                label = f"{fwd_years[i]}\nInsert"
            else:
                label = f"{bwd_years[i-5]}\nDelete"
            
            plt.axvspan(prev_b, b, facecolor=bg_color, alpha=1.0)
            
            if b == pivot:
                plt.axvline(x=b, color='#d62728', linestyle='-', linewidth=2.5, label='Start Deletion (Pivot)', zorder=2)
            else:
                plt.axvline(x=b, color='#888888', linestyle='--', linewidth=1.5, alpha=0.5, zorder=1)
            
            xticks_pos.append((prev_b + b) / 2.0)
            xtick_labels.append(label)
            prev_b = b
            
    plt.xticks(xticks_pos, xtick_labels, rotation=0, ha='center', fontsize=9)
    plt.grid(False, which='major', axis='x')
    plt.grid(True, which='major', axis='y', alpha=0.3)

    max_x = global_boundaries[-1] if global_boundaries else 0
    plt.text(pivot / 2.0, 0.95, "Data Insertion Phase (Growing) ➔", transform=plt.gca().get_xaxis_transform(), ha='center', va='top', fontsize=12, fontweight='bold', color='#005b96', bbox=dict(facecolor='white', alpha=0.7, edgecolor='none', pad=3))
    plt.text(pivot + (max_x - pivot) / 2.0, 0.95, "Data Deletion Phase (Shrinking) ➔", transform=plt.gca().get_xaxis_transform(), ha='center', va='top', fontsize=12, fontweight='bold', color='#990000', bbox=dict(facecolor='white', alpha=0.7, edgecolor='none', pad=3))
    plt.text(pivot, 0.05, " ⬅ Start Deletion ", transform=plt.gca().get_xaxis_transform(), ha='left', va='bottom', fontsize=11, fontweight='bold', color='#d62728', bbox=dict(facecolor='white', alpha=0.8, edgecolor='red', boxstyle='round,pad=0.2'))
    
    plt.title(f"Global Logical Memory Evolution (BP={bp}%)")
    plt.xlabel("Timeline (Years Tightly Grouped)")
    plt.ylabel("Logical Memory Size (MB)")
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    plt.savefig(f"{out_dir}/Global_bp{bp}_Fig1_Memory.png", dpi=300)
    plt.close()

    # 2. Branching & Merging Time Evolution (Separated into 3 plots)
    
    # Pre-calculate data
    algo_fwd_times = {}
    algo_bwd_times = {}
    algo_cum_times = {}
    
    for algo in algos:
        fwd_times = []
        bwd_times = []
        for y in fwd_years:
            df = load_combined_data(y, res_dir, "Forward_Batch")
            if not df.empty and algo in df['Algo'].values:
                algo_df = df[df['Algo'] == algo]
                fwd_times.append(algo_df['Fork_ms'].sum() + algo_df['Commit_ms'].sum() + algo_df['Merge_ms'].sum())
            else: fwd_times.append(None)
            
        for y in bwd_years:
            df = load_combined_data(y, res_dir, "Backward_Batch")
            if not df.empty and algo in df['Algo'].values:
                algo_df = df[df['Algo'] == algo]
                bwd_times.append(algo_df['Fork_ms'].sum() + algo_df['Commit_ms'].sum() + algo_df['Merge_ms'].sum())
            else: bwd_times.append(None)
            
        algo_fwd_times[algo] = fwd_times
        algo_bwd_times[algo] = bwd_times
        
        cumulative_times = []
        acc = 0
        for t in fwd_times + bwd_times:
            if t is not None:
                acc += t
                cumulative_times.append(acc)
            else:
                cumulative_times.append(None)
        algo_cum_times[algo] = cumulative_times

    # 2A. Incremental Insertion: Branching Forward Phase Branching & Merging Time Merging Cost
    plt.figure(figsize=(10, 6))
    for algo in algos:
        if any(t is not None for t in algo_fwd_times[algo]):
            linestyle = '--' if 'Rlog' in algo else '-'
            linewidth = 1.5 if 'Rlog' in algo else 2.5
            plt.plot(fwd_years, algo_fwd_times[algo], marker='o', label=algo, linestyle=linestyle, linewidth=linewidth)
    
    plt.title(f"Incremental Insertion: Branching & Merging Cost (BP={bp}%)")
    plt.xlabel("Year")
    plt.ylabel("Total Time (ms) - Log Scale")
    plt.xticks(fwd_years, [str(y) for y in fwd_years])
    plt.grid(True, alpha=0.3)
    plt.yscale('log')
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    plt.savefig(f"{out_dir}/Global_bp{bp}_Fig2_Branching_Forward.png", dpi=300)
    plt.close()

    # 2B. Incremental Deletion: Branching Backward Phase Branching & Merging Time Merging Cost
    plt.figure(figsize=(10, 6))
    for algo in algos:
        if any(t is not None for t in algo_bwd_times[algo]):
            linestyle = '--' if 'Rlog' in algo else '-'
            linewidth = 1.5 if 'Rlog' in algo else 2.5
            plt.plot(bwd_years, algo_bwd_times[algo], marker='o', label=algo, linestyle=linestyle, linewidth=linewidth)
            
    plt.title(f"Incremental Deletion: Branching & Merging Cost (BP={bp}%)")
    plt.xlabel("Year")
    plt.ylabel("Total Time (ms) - Log Scale")
    plt.xticks(bwd_years, [str(y) for y in bwd_years])
    plt.grid(True, alpha=0.3)
    plt.yscale('log')
    plt.gca().invert_xaxis()
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    plt.savefig(f"{out_dir}/Global_bp{bp}_Fig2_Branching_Backward.png", dpi=300)
    plt.close()

    # 2C. Cumulative Branching & Merging Time
    plt.figure(figsize=(12, 6))
    x_labels_cum = [f"{y}\\nInsert" for y in fwd_years] + [f"{y}\\nDelete" for y in bwd_years]
    for algo in algos:
        if any(t is not None for t in algo_cum_times[algo]):
            linestyle = '--' if 'Rlog' in algo else '-'
            linewidth = 1.5 if 'Rlog' in algo else 2.5
            plt.plot(range(10), algo_cum_times[algo], marker='s', label=algo, linestyle=linestyle, linewidth=linewidth)

    plt.title(f"Cumulative Branching & Merging Time (BP={bp}%)")
    plt.xticks(range(10), x_labels_cum, rotation=0, ha='center', fontsize=9)
    plt.ylabel("Accumulated Time (ms) - Log")
    
    # Add backgrounds
    plt.axvspan(-0.5, 4.5, facecolor='#e4f1f8', alpha=0.6, zorder=0)
    plt.axvspan(4.5, 9.5, facecolor='#fcebe8', alpha=0.6, zorder=0)
    plt.axvline(x=4.5, color='#d62728', linestyle='-', linewidth=2, zorder=1)
    for x_line in range(1, 10):
        if x_line != 4.5:
            plt.axvline(x=x_line - 0.5, color='#aaaaaa', linestyle='--', linewidth=1, alpha=0.5, zorder=0)
            
    plt.grid(True, axis='y', alpha=0.3)
    plt.yscale('log')
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    
    # Annotations
    plt.text(2.0, 0.95, "Data Insertion ➔", transform=plt.gca().get_xaxis_transform(), ha='center', va='top', fontsize=12, fontweight='bold', color='#005b96', bbox=dict(facecolor='white', alpha=0.7, edgecolor='none'))
    plt.text(7.0, 0.95, "Data Deletion ➔", transform=plt.gca().get_xaxis_transform(), ha='center', va='top', fontsize=12, fontweight='bold', color='#990000', bbox=dict(facecolor='white', alpha=0.7, edgecolor='none'))
    plt.text(4.5, 0.05, " ⬅ Pivot ", transform=plt.gca().get_xaxis_transform(), ha='left', va='bottom', fontsize=11, fontweight='bold', color='#d62728', bbox=dict(facecolor='white', alpha=0.8, edgecolor='red', boxstyle='round,pad=0.2'))

    plt.tight_layout()
    plt.savefig(f"{out_dir}/Global_bp{bp}_Fig2_Branching_Cumulative.png", dpi=300)
    plt.close()

    # 3. Query Lifecycle (Unified single plots)
    phases_ordered = [("BaseQuery", 0)]
    for y in fwd_years: phases_ordered.append(("MidQuery", y))
    for y in bwd_years: phases_ordered.append(("EndQuery", y))
    
    # Map the 11 Query points to the end of their respective phases
    query_x_coords = [-1.0] + [b - 0.2 for b in global_boundaries]
    
    # xticks for labels (Aligned exactly with the dots)
    xticks_pos = query_x_coords[:]
    xtick_labels = ["Initial\nTree"]
    for i, b in enumerate(global_boundaries):
        if i < 5:
            xtick_labels.append(f"After\n{fwd_years[i]} Ins.")
        else:
            xtick_labels.append(f"After\n{bwd_years[i-5]} Del.")
        
    year_colors = ['#e8f4f8', '#fff4e6', '#e6ffe6', '#f2e6ff', '#ffe6e6'] 
    
    for qt in query_types:
        plt.figure(figsize=(16, 6))
        valid_plot = False
        for algo in algos:
            latencies = []
            for p_type, y in phases_ordered:
                df = load_combined_data(y, res_dir, p_type)
                if not df.empty and algo in df['Algo'].values:
                    algo_df = df[(df['Algo'] == algo) & (df['QType'].str.strip() == qt)]
                    if not algo_df.empty:
                        latencies.append(algo_df['Time_ms'].mean())
                    else: latencies.append(None)
                else: latencies.append(None)
            
            if any(l is not None for l in latencies):
                valid_plot = True
                linestyle = '--' if 'Rlog' in algo else '-'
                linewidth = 1.5 if 'Rlog' in algo else 2.5
                algo_idx = algos.index(algo)
                jitter = (algo_idx - len(algos)/2) * 0.015
                x_pos = [x + jitter for x in query_x_coords]
                markers = ['o', 's', '^', 'D', 'v', '<', '>', 'p']
                marker = markers[algo_idx % len(markers)]
                
                plt.plot(x_pos, latencies, marker=marker, markersize=8, alpha=0.8, label=algo, linestyle=linestyle, linewidth=linewidth)
                
        if valid_plot:
            plt.xticks(xticks_pos, xtick_labels, rotation=0, ha='center', fontsize=9)
            plt.axvspan(-1.8, -0.2, facecolor='#f0f0f0', alpha=1.0) # Initial block
            plt.axvline(x=-0.2, color='#888888', linestyle='--', linewidth=1.5, alpha=0.5)
            
            prev_b = -0.2
            for i, b in enumerate(global_boundaries):
                year_idx = i if i < 5 else 9 - i
                bg_color = year_colors[year_idx]
                plt.axvspan(prev_b, b, facecolor=bg_color, alpha=1.0)
                if b == pivot:
                    plt.axvline(x=b, color='#d62728', linestyle='-', linewidth=3, zorder=2, label='Start Deletion (Pivot)')
                else:
                    plt.axvline(x=b, color='#888888', linestyle='--', linewidth=1.5, alpha=0.5, zorder=1)
                prev_b = b
                
            plt.yscale('log')
            plt.text(pivot / 2.0, 0.95, "Data Insertion Phase ➔", transform=plt.gca().get_xaxis_transform(), ha='center', va='top', fontsize=12, fontweight='bold', color='#005b96', bbox=dict(facecolor='white', alpha=0.7, edgecolor='none'))
            max_x = global_boundaries[-1] if global_boundaries else 0
            plt.text(pivot + (max_x - pivot) / 2.0, 0.95, "Data Deletion Phase ➔", transform=plt.gca().get_xaxis_transform(), ha='center', va='top', fontsize=12, fontweight='bold', color='#990000', bbox=dict(facecolor='white', alpha=0.7, edgecolor='none'))
            plt.text(pivot, 0.05, " ⬅ Start Deletion ", transform=plt.gca().get_xaxis_transform(), ha='left', va='bottom', fontsize=11, fontweight='bold', color='#d62728', bbox=dict(facecolor='white', alpha=0.8, edgecolor='red', boxstyle='round,pad=0.2'))
            
            plt.title(f"Global {qt} Latency Lifecycle (BP={bp}%)")
            plt.ylabel("Average Time (ms) - Log Scale")
            plt.grid(True, axis='y', alpha=0.3)
            plt.grid(False, axis='x')
            plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
            plt.tight_layout()
            plt.savefig(f"{out_dir}/Global_bp{bp}_Fig3_Query_{qt}.png", dpi=300)
        plt.close()

    # 4. Query Nodes Lifecycle (Unified single plots)
    for qt in query_types:
        plt.figure(figsize=(16, 6))
        valid_plot = False
        for algo in algos:
            nodes = []
            for p_type, y in phases_ordered:
                df = load_combined_data(y, res_dir, p_type)
                if not df.empty and algo in df['Algo'].values and 'Nodes' in df.columns:
                    algo_df = df[(df['Algo'] == algo) & (df['QType'].str.strip() == qt)]
                    if not algo_df.empty:
                        nodes.append(algo_df['Nodes'].mean())
                    else: nodes.append(None)
                else: nodes.append(None)
            
            if any(n is not None and n > 0 for n in nodes):
                valid_plot = True
                linestyle = '--' if 'Rlog' in algo else '-'
                linewidth = 1.5 if 'Rlog' in algo else 2.5
                algo_idx = algos.index(algo)
                jitter = (algo_idx - len(algos)/2) * 0.015
                x_pos = [x + jitter for x in query_x_coords]
                markers = ['o', 's', '^', 'D', 'v', '<', '>', 'p']
                marker = markers[algo_idx % len(markers)]
                
                plt.plot(x_pos, nodes, marker=marker, markersize=8, alpha=0.8, label=algo, linestyle=linestyle, linewidth=linewidth)
                
        if valid_plot:
            plt.xticks(xticks_pos, xtick_labels, rotation=0, ha='center', fontsize=9)
            plt.axvspan(-1.8, -0.2, facecolor='#f0f0f0', alpha=1.0) # Initial block
            plt.axvline(x=-0.2, color='#888888', linestyle='--', linewidth=1.5, alpha=0.5)
            
            prev_b = -0.2
            for i, b in enumerate(global_boundaries):
                year_idx = i if i < 5 else 9 - i
                bg_color = year_colors[year_idx]
                plt.axvspan(prev_b, b, facecolor=bg_color, alpha=1.0)
                if b == pivot:
                    plt.axvline(x=b, color='#d62728', linestyle='-', linewidth=3, zorder=2, label='Start Deletion (Pivot)')
                else:
                    plt.axvline(x=b, color='#888888', linestyle='--', linewidth=1.5, alpha=0.5, zorder=1)
                prev_b = b
                
            plt.yscale('log')
            plt.text(pivot / 2.0, 0.95, "Data Insertion Phase ➔", transform=plt.gca().get_xaxis_transform(), ha='center', va='top', fontsize=12, fontweight='bold', color='#005b96', bbox=dict(facecolor='white', alpha=0.7, edgecolor='none'))
            max_x = global_boundaries[-1] if global_boundaries else 0
            plt.text(pivot + (max_x - pivot) / 2.0, 0.95, "Data Deletion Phase ➔", transform=plt.gca().get_xaxis_transform(), ha='center', va='top', fontsize=12, fontweight='bold', color='#990000', bbox=dict(facecolor='white', alpha=0.7, edgecolor='none'))
            plt.text(pivot, 0.05, " ⬅ Start Deletion ", transform=plt.gca().get_xaxis_transform(), ha='left', va='bottom', fontsize=11, fontweight='bold', color='#d62728', bbox=dict(facecolor='white', alpha=0.8, edgecolor='red', boxstyle='round,pad=0.2'))
            
            plt.title(f"Global {qt} Nodes Touched Lifecycle (BP={bp}%)")
            plt.ylabel("Average Nodes Touched - Log Scale")
            plt.grid(True, axis='y', alpha=0.3)
            plt.grid(False, axis='x')
            plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
            plt.tight_layout()
            plt.savefig(f"{out_dir}/Global_bp{bp}_Fig4_Nodes_{qt}.png", dpi=300)
        plt.close()

def print_interim_summary(bp, year, res_dir):
    logger.info(f"\n---> SUMMARY FOR {year} (BP={bp}%) <---")
    df_fwd = load_combined_data(year, res_dir, "Forward_Batch")
    df_bwd = load_combined_data(year, res_dir, "Backward_Batch")
    df_mid = load_combined_data(year, res_dir, "MidQuery")
    
    summary = []
    
    if not df_fwd.empty:
        agg_fwd = df_fwd.groupby('Algo')[['Fork_ms', 'Commit_ms', 'Merge_ms']].sum()
        agg_fwd['Inc_Ins_ms'] = agg_fwd['Fork_ms'] + agg_fwd['Commit_ms'] + agg_fwd['Merge_ms']
        
        agg_bwd = pd.DataFrame()
        if not df_bwd.empty:
            agg_bwd = df_bwd.groupby('Algo')[['Fork_ms', 'Commit_ms', 'Merge_ms']].sum()
            agg_bwd['Inc_Del_ms'] = agg_bwd['Fork_ms'] + agg_bwd['Commit_ms'] + agg_bwd['Merge_ms']
            
        mid_q = pd.DataFrame()
        if not df_mid.empty:
            df_mid['QType'] = df_mid['QType'].str.strip()
            mid_q = df_mid.pivot_table(index='Algo', columns='QType', values='Time_ms', aggfunc='mean')
            
        res = pd.concat([agg_fwd[['Inc_Ins_ms']], mid_q, agg_bwd[['Inc_Del_ms']] if not agg_bwd.empty else pd.DataFrame()], axis=1)
        res = res.reindex(algos).dropna(how='all').round(2)
        logger.info("\n" + res.to_string())
    logger.info("-" * 60 + "\n")

def run_benchmark():
    
    task_count = 0
    logger.info(f"Starting Real-World MVZD OSM Suite: {len(algos)} Algorithms.")
    
    dataset_name = os.path.basename(data_dir.strip('/')).replace('_workload', '')
    total_tasks = len(algos) * len(steps)
    task_count = 0
    
    for step in steps:
        step_dest = f"results_real_world/{dataset_name}/step_{step}"
        if not os.path.exists(step_dest): os.makedirs(step_dest)
        
        logger.info(f"=======================================================")
        logger.info(f"Phase: Event-Driven Workload (Step Size = {step})")
        logger.info(f"=======================================================")
        
        for algo in algos:
            task_count += 1
            logger.info(f"  [{task_count}/{total_tasks}] -> Running {algo} at step {step}...")
            
            start_t = time.time()
            if not os.path.exists("verification_results"): os.makedirs("verification_results")
            else:
                for f in os.listdir("verification_results"): os.remove(os.path.join("verification_results", f))
            
            algo_dest = os.path.join(step_dest, algo)
            if not os.path.exists(algo_dest): os.makedirs(algo_dest)
            
            real_algo = algo.replace("_Par", "")
            
            expected_end_file = os.path.join(algo_dest, f"2026_EndQuery_{real_algo}.txt")
            if os.path.exists(expected_end_file):
                logger.info(f"     [SKIP] {algo} at step {step} already completed.")
                continue
            
            env = os.environ.copy()
            if "_Par" in algo:
                if "PARLAY_NUM_THREADS" in env:
                    del env["PARLAY_NUM_THREADS"]
                commits_dir = "01_commits"
                if par_batch_size > 1:
                    commits_dir = f"03_commits_par_{par_batch_size}"
                    full_commits_path = os.path.join(data_dir, commits_dir)
                    if not os.path.exists(full_commits_path):
                        logger.info(f"Generating parallel batch dataset: {full_commits_path}")
                        subprocess.run(f"python3 scripts/generate_parallel_csv.py --dataset_dir {data_dir} --batch_size {par_batch_size}", shell=True, check=True)
            else:
                env["PARLAY_NUM_THREADS"] = "1"
                commits_dir = "01_commits"
            
            cmd = f"./verify_bench -algo {real_algo} -q_step {step} -start_year {start_year} -end_year {end_year} -dir {data_dir} -commits_dir {commits_dir}"
            process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, env=env)
            for line in iter(process.stdout.readline, ''):
                line = line.strip()
                if line and ("[" in line or "--- YEAR" in line or "Overlaps" in line):
                    logger.info(f"       {line}")
            
            process.wait()
            ret = process.returncode
            elapsed_s = time.time() - start_t
            
            if ret != 0:
                logger.warning(f"     [WARNING/OOM] {algo} failed (Time: {elapsed_s:.1f}s)")
            else:
                logger.info(f"     [OK] {algo} completed successfully in {elapsed_s:.1f}s")
            
            for f in os.listdir("verification_results"):
                if f.endswith(".txt") and real_algo in f:
                    shutil.move(os.path.join("verification_results", f), os.path.join(algo_dest, f))

def generate_hybrid_plots():
    import seaborn as sns
    import pandas as pd
    import matplotlib.pyplot as plt
    import glob
    import os
    import re
    
    dataset_name = os.path.basename(data_dir.strip('/')).replace('_workload', '')
    base_dir = f"results_real_world/{dataset_name}"
    
    colors = sns.color_palette("tab10", len(algos))
    style_map = {algo: {'color': colors[i], 'marker': 'o'} for i, algo in enumerate(algos)}
    
    print("Parsing physical step folders for Hybrid Workload...")
    
    hybrid_results = {algo: {'step': [], 'total_time_ms': []} for algo in algos}
    mem_results = {}
    
    for step in steps:
        step_dir = os.path.join(base_dir, f"step_{step}")
        fig_dir = os.path.join(step_dir, "figures")
        os.makedirs(fig_dir, exist_ok=True)
        
        for algo in algos:
            algo_dir = os.path.join(step_dir, algo)
            files = glob.glob(f"{algo_dir}/*_Commit_*.txt")
            
            total_commit_time = 0.0
            total_query_time = 0.0
            
            for f in files:
                try:
                    df = pd.read_csv(f, sep="|", skipinitialspace=True)
                    df.columns = [c.strip() for c in df.columns]
                    df['QType'] = df['QType'].str.strip()
                    
                    c_df = df[df['QType'] == 'Commit_Avg']
                    if not c_df.empty:
                        total_commit_time += c_df['Time_ms'].iloc[0] * step
                        # Also collect memory evolution (using step 500 as the atomic timeline)
                        if step == 500:
                            match_cnum = re.search(r"Commit_(\d+).txt", f)
                            if match_cnum:
                                cnum = int(match_cnum.group(1))
                                if algo not in mem_results:
                                    mem_results[algo] = {'commit': [], 'index_mem': [], 'raw_mem': []}
                                mem_results[algo]['commit'].append(cnum)
                                mem_results[algo]['index_mem'].append(c_df['Mem_MB'].iloc[0])
                                mem_results[algo]['raw_mem'].append(c_df['Nodes'].iloc[0]) # We repurposed Nodes column
                        
                    q_df = df[df['QType'] != 'Commit_Avg']
                    if not q_df.empty:
                        total_query_time += q_df['Time_ms'].sum()
                except Exception:
                    pass
            
            if total_commit_time > 0 or total_query_time > 0:
                hybrid_results[algo]['step'].append(step)
                hybrid_results[algo]['total_time_ms'].append(total_commit_time + total_query_time)
                
        # Generate plot for this step just to show something
        plt.figure(figsize=(10, 6))
        for algo in algos:
            if algo in hybrid_results and hybrid_results[algo]['step']:
                plt.bar(algo, hybrid_results[algo]['total_time_ms'][-1], color=style_map[algo]['color'])
        
        plt.xlabel('Algorithm', fontsize=12)
        plt.ylabel('Total Execution Time (ms)', fontsize=12)
        plt.title(f'Total Time at Step {step} ({dataset_name})', fontsize=14)
        plt.grid(True, linestyle='--', alpha=0.6, axis='y')
        plt.tight_layout()
        plt.savefig(os.path.join(fig_dir, f'Total_Time_Step_{step}.png'), dpi=300)
        plt.close()
        
    print("Generating Memory Overhead Evolution Plot...")
    plt.figure(figsize=(12, 6))
    plotted_raw = False
    for algo in algos:
        if algo in mem_results and mem_results[algo]['commit']:
            # Sort by commit
            sorted_idx = sorted(range(len(mem_results[algo]['commit'])), key=lambda k: mem_results[algo]['commit'][k])
            x = [mem_results[algo]['commit'][i] for i in sorted_idx]
            y_idx = [mem_results[algo]['index_mem'][i] for i in sorted_idx]
            y_raw = [mem_results[algo]['raw_mem'][i] for i in sorted_idx]
            
            plt.plot(x, y_idx, label=algo + " Index", color=style_map[algo]['color'], linewidth=2)
            
            if not plotted_raw:
                plt.plot(x, y_raw, label="Raw Data (Baseline)", color='black', linestyle='--', linewidth=2)
                plotted_raw = True

    plt.xlabel('Commit Sequence (Time)', fontsize=12)
    plt.ylabel('Memory Footprint (MB)', fontsize=12)
    plt.title(f'Multi-Version Index Memory Overhead vs Raw Data ({dataset_name})', fontsize=14)
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    plt.savefig(os.path.join(base_dir, "figures", "Memory_Overhead_Evolution.png"), dpi=300)
    plt.close()
    
    print(f"Hybrid simulation complete. Plot saved to {base_dir}/step_X/figures/")


import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run MVZD Benchmarks")
    parser.add_argument('--algos', nargs='+', default=["MVZD", "CPAMBB", "Rlog_1yr", "Rlog_2yr"], help="Algorithms to run")
    parser.add_argument('--step', type=int, default=100, help="Query step size")
    parser.add_argument('--start_year', type=int, default=2018)
    parser.add_argument('--end_year', type=int, default=2026)
    parser.add_argument('--par_batch_size', type=int, default=1000, help="Batch size for Parallel algorithms")
    parser.add_argument('--dir', type=str, default="dataset/bhutan_workload", help="Path to workload directory")
    args = parser.parse_args()

    # Override globals
    global algos, steps, start_year, end_year, data_dir, par_batch_size
    algos = args.algos
    steps = [args.step]
    start_year = args.start_year
    end_year = args.end_year
    data_dir = args.dir
    par_batch_size = args.par_batch_size

    run_benchmark()
    print("\n[All Done] Real-World Benchmark Completed! Generating Plots...")
    try:
        generate_hybrid_plots()
    except Exception as e:
        print(f"Error during plotting: {e}")

