import os
import shutil
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import logging
import time
import subprocess
from matplotlib.ticker import MaxNLocator

os.environ['PARLAY_NUM_THREADS'] = '1'

logger = logging.getLogger("Benchmark")
logger.setLevel(logging.INFO)
formatter = logging.Formatter('[%(asctime)s] %(message)s', datefmt='%Y-%m-%d %H:%M:%S')

ch = logging.StreamHandler()
ch.setFormatter(formatter)
logger.addHandler(ch)

fh = logging.FileHandler("benchmark_progress.log", mode='a')
fh.setFormatter(formatter)
logger.addHandler(fh)

bps = [100, 50, 25, 12.5, 6.25, 3.125]
algos = ["MVZD", "CPAMBB", "Rlog_1yr", "Rlog_2yr", "Rlog_3yr", "Rlog_4yr", "Rlog_5yr", "Rlog_NoSnap"]
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
    total_tasks = len(bps) * len(algos)
    task_count = 0
    logger.info(f"Starting Benchmark Suite: {len(bps)} BPs, {len(algos)} Algorithms. Total jobs = {total_tasks}")
    
    for bp in bps:
        dest = f"results_bp_{bp}"
        if not os.path.exists(dest): os.makedirs(dest)
        
        logger.info(f"=======================================================")
        logger.info(f"Phase: BP={bp}%")
        logger.info(f"=======================================================")
        
        for algo in algos:
            task_count += 1
            logger.info(f"  [{task_count}/{total_tasks}] -> Running {algo} (Full Lifecycle 2021 <-> 2025)...")
            
            start_t = time.time()
            if not os.path.exists("verification_results"): os.makedirs("verification_results")
            else:
                for f in os.listdir("verification_results"): os.remove(os.path.join("verification_results", f))
            
            expected_end_file = os.path.join(dest, f"2025_EndQuery_{algo}.txt")
            if os.path.exists(expected_end_file):
                logger.info(f"     [SKIP] {algo} already completed full lifecycle.")
                continue
            
            cmd = f"./verify_bench -bp {bp} -algo {algo}"
            process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
            for line in iter(process.stdout.readline, ''):
                line = line.strip()
                if line and ("[" in line or "--- YEAR" in line or "Batches:" in line):
                    logger.info(f"       {line}")
            
            process.wait()
            ret = process.returncode
            elapsed_s = time.time() - start_t
            
            if ret != 0:
                logger.warning(f"     [WARNING/OOM] {algo} failed for BP={bp} (Time: {elapsed_s:.1f}s)")
            else:
                logger.info(f"     [OK] {algo} completed successfully in {elapsed_s:.1f}s")
            
            for f in os.listdir("verification_results"):
                if f.endswith(".txt") and algo in f:
                    shutil.move(os.path.join("verification_results", f), os.path.join(dest, f))
        
        for year in years:
            print_interim_summary(bp, year, dest)
            
        try:
            generate_global_plots(bp, dest)
            logger.info(f"  [Plot] Successfully generated GLOBAL visual plots for BP={bp}% in 'plots_and_tables/bp_{bp}/'")
        except Exception as e:
            logger.error(f"  [Plot Error] Failed to generate global plots for BP={bp}%: {e}")


def parse_results():
    data_time = {algo: {"Method": algo, "Distribution": "bhutan"} for algo in algos}
    data_mem = {algo: {"Method": algo, "Distribution": "bhutan"} for algo in algos}
    data_nodes = {algo: {"Method": algo, "Distribution": "bhutan"} for algo in algos}
    
    for bp in bps:
        res_dir = f"results_bp_{bp}"
        if not os.path.exists(res_dir): continue
        
        for algo in algos:
            base_file = os.path.join(res_dir, f"BaseQuery_{algo}.txt")
            if os.path.exists(base_file):
                try:
                    df = clean_df(pd.read_csv(base_file, sep="|", skipinitialspace=True))
                    if len(df) > 0:
                        df['QType'] = df['QType'].str.strip()
                        for qtype in query_types:
                            qdf = df[df['QType'] == qtype]
                            if len(qdf) > 0:
                                avg_t = qdf['Time_ms'].mean()
                                avg_n = qdf['Nodes'].mean() if 'Nodes' in qdf.columns else 0
                                data_time[algo][f"Base Query {qtype}"] = f"{avg_t:.4f}"
                                data_nodes[algo][f"Base Query {qtype}"] = f"{avg_n:.1f}"
                except: pass

        for year in years:
            for algo in algos:
                fwd_file = os.path.join(res_dir, f"{year}_Forward_Batch_{algo}.txt")
                if os.path.exists(fwd_file):
                    try:
                        df = clean_df(pd.read_csv(fwd_file, sep="|", skipinitialspace=True))
                        if len(df) > 0:
                            total_ms = (df['Fork_ms'] + df['Commit_ms'] + df['Merge_ms']).sum()
                            data_time[algo][f"{year} Inc Ins {bp}%"] = f"{total_ms:.4f}"
                            data_mem[algo][f"{year} Inc Ins {bp}%"] = f"{df['Mem_MB'].iloc[-1]:.3f}"
                    except: pass
                
                bwd_file = os.path.join(res_dir, f"{year}_Backward_Batch_{algo}.txt")
                if os.path.exists(bwd_file):
                    try:
                        df = clean_df(pd.read_csv(bwd_file, sep="|", skipinitialspace=True))
                        if len(df) > 0:
                            total_ms = (df['Fork_ms'] + df['Commit_ms'] + df['Merge_ms']).sum()
                            data_time[algo][f"{year} Inc Del {bp}%"] = f"{total_ms:.4f}"
                            data_mem[algo][f"{year} Inc Del {bp}%"] = f"{df['Mem_MB'].iloc[-1]:.3f}"
                    except: pass

                mid_file = os.path.join(res_dir, f"{year}_MidQuery_{algo}.txt")
                if os.path.exists(mid_file):
                    try:
                        df = clean_df(pd.read_csv(mid_file, sep="|", skipinitialspace=True))
                        if len(df) > 0:
                            df['QType'] = df['QType'].str.strip()
                            for qtype in query_types:
                                qdf = df[df['QType'] == qtype]
                                if len(qdf) > 0:
                                    avg_t = qdf['Time_ms'].mean()
                                    avg_n = qdf['Nodes'].mean() if 'Nodes' in qdf.columns else 0
                                    data_time[algo][f"{year} Query after inc ins {qtype}"] = f"{avg_t:.4f}"
                                    data_nodes[algo][f"{year} Query after inc ins {qtype}"] = f"{avg_n:.1f}"
                    except: pass

                end_file = os.path.join(res_dir, f"{year}_EndQuery_{algo}.txt")
                if os.path.exists(end_file):
                    try:
                        df = clean_df(pd.read_csv(end_file, sep="|", skipinitialspace=True))
                        if len(df) > 0:
                            df['QType'] = df['QType'].str.strip()
                            for qtype in query_types:
                                qdf = df[df['QType'] == qtype]
                                if len(qdf) > 0:
                                    avg_t = qdf['Time_ms'].mean()
                                    avg_n = qdf['Nodes'].mean() if 'Nodes' in qdf.columns else 0
                                    data_time[algo][f"{year} Query after inc del {qtype}"] = f"{avg_t:.4f}"
                                    data_nodes[algo][f"{year} Query after inc del {qtype}"] = f"{avg_n:.1f}"
                    except: pass
                            
    cols = ["Distribution", "Method"]
    for qt in query_types: cols.append(f"Base Query {qt}")
    for year in years:
        for bp in bps: cols.append(f"{year} Inc Ins {bp}%")
        for qt in query_types: cols.append(f"{year} Query after inc ins {qt}")
        for bp in bps: cols.append(f"{year} Inc Del {bp}%")
        for qt in query_types: cols.append(f"{year} Query after inc del {qt}")
        
    df_time = pd.DataFrame(list(data_time.values()))
    df_mem = pd.DataFrame(list(data_mem.values()))
    df_nodes = pd.DataFrame(list(data_nodes.values()))
    
    cols_time = [c for c in cols if c in df_time.columns]
    cols_mem = [c for c in cols if c in df_mem.columns]
    cols_nodes = [c for c in cols if c in df_nodes.columns]
    
    df_time = df_time[cols_time]
    df_mem = df_mem[cols_mem]
    df_nodes = df_nodes[cols_nodes]
    
    df_time.to_csv("final_report_time.csv", index=False)
    df_mem.to_csv("final_report_memory.csv", index=False)
    df_nodes.to_csv("final_report_nodes.csv", index=False)
    logger.info("\n[SUCCESS] Generated final_report_time.csv, final_report_memory.csv, and final_report_nodes.csv!")
    generate_publication_tables()



def format_ms(v):
    try:
        v = float(v)
        return f"{v:.3f}"
    except:
        return v

def generate_publication_tables():
    import matplotlib.pyplot as plt
    import matplotlib.patches as patches
    import pandas as pd
    import numpy as np
    import re
    import os

    logger.info("  [Table] Generating publication-quality tables...")

    if not os.path.exists("final_report_time.csv"):
        logger.warning("  [Table] final_report_time.csv not found. Skipping table generation.")
        return

    # 1. Create the intermediate Table_Year_*.csv
    df = pd.read_csv("final_report_time.csv")
    build_times = {"MVZD": 38.736, "CPAMBB": 112.695}
    for algo in ["Rlog_1yr", "Rlog_2yr", "Rlog_3yr", "Rlog_4yr", "Rlog_5yr", "Rlog_NoSnap"]:
        build_times[algo] = 113.361

    our_bps = ["100%", "50%", "25%", "12.5%", "6.25%", "3.125%"]
    queries = ["KNN_1", "KNN_10", "KNN_100", "Range_Small", "Range_Med", "Range_Large"]

    for year in range(2021, 2026):
        out_file = f"Table_Year_{year}.csv"
        with open(out_file, "w") as f:
            f.write("Trees,Build")
            f.write(",Incremental Insert" + "," * (len(our_bps) - 1))
            f.write(",Query after Inc. Ins." + "," * (len(queries) - 1))
            f.write(",Incremental Delete" + "," * (len(our_bps) - 1))
            f.write(",Query after Inc. Del." + "," * (len(queries) - 1))
            f.write("\n")
            f.write(",")
            for bp in our_bps: f.write(f",{bp}")
            for q in queries: f.write(f",{q}")
            for bp in our_bps: f.write(f",{bp}")
            for q in queries: f.write(f",{q}")
            f.write("\n")
            
            algos_ordered = ["MVZD", "CPAMBB", "Rlog_1yr", "Rlog_2yr", "Rlog_3yr", "Rlog_4yr", "Rlog_5yr", "Rlog_NoSnap"]
            for algo in algos_ordered:
                row_df = df[df['Method'] == algo]
                if row_df.empty: continue
                row_data = row_df.iloc[0]
                f.write(f"{algo},{build_times[algo]:.4f}")
                for bp in our_bps:
                    col_name = f"{year} Inc Ins {bp}"
                    val = row_data[col_name] if col_name in row_data else "N/A"
                    f.write(f",{val}")
                for q in queries:
                    col_name = f"{year} Query after inc ins {q}"
                    val = row_data[col_name] if col_name in row_data else "N/A"
                    f.write(f",{val}")
                for bp in our_bps:
                    col_name = f"{year} Inc Del {bp}"
                    val = row_data[col_name] if col_name in row_data else "N/A"
                    f.write(f",{val}")
                for q in queries:
                    col_name = f"{year} Query after inc del {q}"
                    val = row_data[col_name] if col_name in row_data else "N/A"
                    f.write(f",{val}")
                f.write("\n")

    # 2. Draw Individual Yearly Tables
    for year in range(2021, 2026):
        try:
            df_year = pd.read_csv(f"Table_Year_{year}.csv", header=None)
            data = df_year.values.tolist()
            has_build = (year == 2021)
            if not has_build:
                for r in range(len(data)): data[r].pop(1)

            cols = len(data[0])
            rows = len(data)

            if has_build:
                widths = [3.0, 1.3] + [1.8]*24
                ins_start, q_ins_start, del_start, q_del_start = 2, 8, 14, 20
            else:
                widths = [3.0] + [1.8]*24
                ins_start, q_ins_start, del_start, q_del_start = 1, 7, 13, 19

            total_w = sum(widths)
            cum_w = [0] + list(np.cumsum(widths))
            row_h = 0.8
            total_h = rows * row_h

            fig, ax = plt.subplots(figsize=(total_w * 0.55, (total_h + 2.0) * 0.55))
            ax.set_xlim(0, total_w)
            ax.set_ylim(-2.0, total_h)
            ax.axis('off')

            min_vals = {}
            for c in range(1, cols):
                col_vals = []
                for r in range(2, rows):
                    try:
                        col_vals.append(float(data[r][c]))
                    except: pass
                if col_vals: min_vals[c] = min(col_vals)

            for r in range(rows):
                y = total_h - (r + 1) * row_h
                algo_name = str(data[r][0]) if r >= 2 else ""
                compact_ins = False
                compact_del = False
                if r >= 2:
                    m = re.match(r"Rlog_(\d+)yr", algo_name)
                    if m:
                        N = int(m.group(1))
                        if (year - 2020) % N == 0: compact_ins = True
                        if (5 + (2026 - year)) % N == 0: compact_del = True

                for c in range(cols):
                    x = cum_w[c]
                    w = widths[c]
                    raw_val = str(data[r][c])
                    if pd.isna(data[r][c]) or raw_val == "nan": raw_val = ""
                    
                    is_header = r < 2
                    bg_color = 'white'
                    font_weight = 'normal'
                    text_color = 'black'
                    is_fastest = False
                    is_compact_cell = False
                    display_val = raw_val
                    
                    if is_header:
                        bg_color = '#eaeaea'
                        font_weight = 'bold'
                        if r == 1:
                            display_val = display_val.replace("Range_", "R_").replace("KNN_", "K_").replace("%", "")
                    elif c == 0:
                        bg_color = '#f9f9f9'
                        font_weight = 'bold'
                        if compact_ins and compact_del: display_val = f"{algo_name} *[Both]"
                        elif compact_ins: display_val = f"{algo_name} *[Ins]"
                        elif compact_del: display_val = f"{algo_name} *[Del]"
                    else:
                        if raw_val != "N/A" and raw_val != "":
                            try:
                                fval = float(raw_val)
                                display_val = format_ms(fval)
                                min_v = min_vals.get(c, None)
                                if min_v is not None:
                                    if min_v == 0:
                                        is_fastest = (fval == 0)
                                        ratio = 1.0 if fval == 0 else 999.0
                                    else:
                                        is_fastest = (fval <= min_v * 1.001)
                                        ratio = fval / min_v
                                    if ratio < 3.0: bg_color = '#78d978' 
                                    elif ratio < 10.0: bg_color = '#a3e6a3' 
                                    elif ratio < 35.0: bg_color = '#d9f2d9' 
                                    elif ratio < 100.0: bg_color = '#e0e0e0' 
                                    else: bg_color = '#b0b0b0' 
                                    if is_fastest: font_weight = 'bold'
                            except: pass
                            
                    bbox_props = None
                    if r >= 2 and raw_val != "N/A" and raw_val != "":
                        if compact_ins and ins_start <= c < ins_start + 6: is_compact_cell = True
                        if compact_del and del_start <= c < del_start + 6: is_compact_cell = True
                        if is_compact_cell:
                            bbox_props = dict(boxstyle='square,pad=0.25', facecolor='none', edgecolor='black', linewidth=1.2)

                    if r == 0:
                        if c == ins_start: w = sum(widths[ins_start:ins_start+6])
                        elif c == q_ins_start: w = sum(widths[q_ins_start:q_ins_start+6])
                        elif c == del_start: w = sum(widths[del_start:del_start+6])
                        elif c == q_del_start: w = sum(widths[q_del_start:q_del_start+6])
                        elif c not in [0, 1] if has_build else c != 0: continue
                            
                    rect = patches.Rectangle((x, y), w, row_h, linewidth=0.5, edgecolor='#dddddd', facecolor=bg_color)
                    ax.add_patch(rect)
                    
                    if bbox_props:
                        ax.text(x + w/2, y + row_h/2, display_val, ha='center', va='center', fontsize=12 if not is_header else 13, fontweight=font_weight, color=text_color, bbox=bbox_props)
                    else:
                        ax.text(x + w/2, y + row_h/2, display_val, ha='center', va='center', fontsize=12 if not is_header else 13, fontweight=font_weight, color=text_color)

                    if is_fastest and not is_header:
                        text_len = len(display_val) * 0.17 
                        ax.plot([x + w/2 - text_len/2, x + w/2 + text_len/2], [y + row_h/2 - 0.20, y + row_h/2 - 0.20], color='black', linewidth=1.5)

            if has_build:
                group_boundaries = [0, cum_w[1], cum_w[2], cum_w[8], cum_w[14], cum_w[20], total_w]
            else:
                group_boundaries = [0, cum_w[1], cum_w[7], cum_w[13], cum_w[19], total_w]

            for gb in group_boundaries:
                ax.plot([gb, gb], [0, total_h], color='black', linewidth=2.0 if gb in [0, total_w] else 1.2)
                
            ax.plot([0, total_w], [total_h, total_h], color='black', linewidth=2) 
            ax.plot([0, total_w], [total_h - row_h, total_h - row_h], color='black', linewidth=1) 
            ax.plot([0, total_w], [total_h - 2*row_h, total_h - 2*row_h], color='black', linewidth=1.5) 
            ax.plot([0, total_w], [0, 0], color='black', linewidth=2) 

            legend_y = -0.8
            legend_x = 0.5
            ax.text(legend_x, legend_y, "The fastest time is in", fontsize=12, va='center')
            legend_x += 4.5
            ax.text(legend_x, legend_y, "bold and underlined", fontsize=12, fontweight='bold', va='center')
            ax.plot([legend_x - 0.1, legend_x + 3.8], [legend_y - 0.2, legend_y - 0.2], color='black', linewidth=1.5)
            legend_x += 4.2
            boxes = [('#78d978', '< 3x the fastest'), ('#a3e6a3', '< 10x the fastest'), ('#d9f2d9', '< 35x the fastest'), ('#e0e0e0', '< 100x the fastest'), ('#b0b0b0', '> 100x the fastest')]
            for color, label in boxes:
                rect = patches.Rectangle((legend_x, legend_y - 0.2), 0.8, 0.4, linewidth=0.5, edgecolor='#cccccc', facecolor=color)
                ax.add_patch(rect)
                legend_x += 1.0
                ax.text(legend_x, legend_y, f": {label}", fontsize=12, va='center')
                legend_x += 3.8
                
            ax.text(0.5, -1.6, "* Indicates the algorithm triggered a full tree compaction. The specific phase is tagged as [Ins], [Del], or [Both].", ha='left', va='center', fontsize=12, style='italic', color='#555555')
            ax.text(total_w/2, total_h + 0.3, f"Running time (in milliseconds) on Bhutan dataset (Year {year}). Lower is better.", ha='center', va='bottom', fontsize=14, fontweight='bold', color='#333333')
            plt.tight_layout()
            plt.savefig(f"Paper_Table_Year_{year}_Beautiful.png", dpi=300, bbox_inches='tight')
            plt.close()
        except Exception as e:
            logger.error(f"Failed to generate table for {year}: {e}")

    # 3. Draw Super Combined Table
    try:
        all_data = []
        df2021 = pd.read_csv("Table_Year_2021.csv", header=None)
        headers = df2021.values.tolist()[:2]
        headers[0].insert(0, "Year")
        headers[1].insert(0, "")
        all_data.extend(headers)

        year_rows_map = {} 
        min_vals_per_year = {}
        current_row_idx = 2

        for year in range(2021, 2026):
            df_year = pd.read_csv(f"Table_Year_{year}.csv", header=None)
            data = df_year.values.tolist()[2:]
            
            min_vals = {}
            cols_orig = len(data[0])
            for c in range(1, cols_orig):
                col_vals = []
                for r in range(len(data)):
                    try: col_vals.append(float(data[r][c]))
                    except: pass
                if col_vals: min_vals[c] = min(col_vals)
            min_vals_per_year[year] = min_vals

            for r in range(len(data)):
                if year != 2021: data[r][1] = "-"
                data[r].insert(0, str(year))
                all_data.append(data[r])
                year_rows_map[current_row_idx] = year
                current_row_idx += 1

        cols = len(all_data[0])
        rows = len(all_data)
        widths = [1.2, 2.8, 1.3] + [1.8]*24
        total_w = sum(widths)
        cum_w = [0] + list(np.cumsum(widths))
        row_h = 0.8
        total_h = rows * row_h

        fig, ax = plt.subplots(figsize=(total_w * 0.55, (total_h + 2.0) * 0.55))
        ax.set_xlim(0, total_w)
        ax.set_ylim(-2.0, total_h)
        ax.axis('off')

        for r in range(rows):
            y = total_h - (r + 1) * row_h
            year = year_rows_map.get(r, None)
            algo_name = str(all_data[r][1]) if r >= 2 else ""
            
            compact_ins = False
            compact_del = False
            if r >= 2:
                m = re.match(r"Rlog_(\d+)yr", algo_name)
                if m:
                    N = int(m.group(1))
                    if (year - 2020) % N == 0: compact_ins = True
                    if (5 + (2026 - year)) % N == 0: compact_del = True

            for c in range(cols):
                x = cum_w[c]
                w = widths[c]
                raw_val = str(all_data[r][c])
                if pd.isna(all_data[r][c]) or raw_val == "nan": raw_val = ""
                
                is_header = r < 2
                bg_color = 'white'
                font_weight = 'normal'
                text_color = 'black'
                is_fastest = False
                is_compact_cell = False
                display_val = raw_val
                
                if is_header:
                    bg_color = '#eaeaea'
                    font_weight = 'bold'
                    if r == 1: display_val = display_val.replace("Range_", "R_").replace("KNN_", "K_").replace("%", "")
                elif c == 0:
                    bg_color = '#eaeaea' 
                    font_weight = 'bold'
                elif c == 1:
                    bg_color = '#f9f9f9'
                    font_weight = 'bold'
                    if compact_ins and compact_del: display_val = f"{algo_name} *[Both]"
                    elif compact_ins: display_val = f"{algo_name} *[Ins]"
                    elif compact_del: display_val = f"{algo_name} *[Del]"
                else:
                    if raw_val != "N/A" and raw_val != "-" and raw_val != "":
                        try:
                            fval = float(raw_val)
                            display_val = format_ms(fval)
                            min_v = min_vals_per_year[year].get(c-1, None)
                            if min_v is not None:
                                if min_v == 0:
                                    is_fastest = (fval == 0)
                                    ratio = 1.0 if fval == 0 else 999.0
                                else:
                                    is_fastest = (fval <= min_v * 1.001)
                                    ratio = fval / min_v
                                if ratio < 3.0: bg_color = '#78d978' 
                                elif ratio < 10.0: bg_color = '#a3e6a3' 
                                elif ratio < 35.0: bg_color = '#d9f2d9' 
                                elif ratio < 100.0: bg_color = '#e0e0e0' 
                                else: bg_color = '#b0b0b0' 
                                if is_fastest: font_weight = 'bold'
                        except: pass
                        
                bbox_props = None
                if r >= 2 and raw_val != "N/A" and raw_val != "-" and raw_val != "":
                    if compact_ins and 3 <= c < 9: is_compact_cell = True
                    if compact_del and 15 <= c < 21: is_compact_cell = True
                    if is_compact_cell:
                        bbox_props = dict(boxstyle='square,pad=0.25', facecolor='none', edgecolor='black', linewidth=1.2)

                if r == 0:
                    if c == 3: w = sum(widths[3:9])
                    elif c == 9: w = sum(widths[9:15])
                    elif c == 15: w = sum(widths[15:21])
                    elif c == 21: w = sum(widths[21:27])
                    elif c not in [0, 1, 2]: continue
                
                if c == 0 and r >= 2:
                    if (r - 2) % 8 == 0:
                        rect = patches.Rectangle((x, y - 7*row_h), w, 8*row_h, linewidth=0.5, edgecolor='#dddddd', facecolor=bg_color)
                        ax.add_patch(rect)
                        ax.text(x + w/2, y - 3*row_h, str(year), ha='center', va='center', rotation=90, fontsize=20, fontweight='bold')
                else:
                    rect = patches.Rectangle((x, y), w, row_h, linewidth=0.5, edgecolor='#dddddd', facecolor=bg_color)
                    ax.add_patch(rect)
                    if bbox_props:
                        ax.text(x + w/2, y + row_h/2, display_val, ha='center', va='center', fontsize=12 if not is_header else 13, fontweight=font_weight, color=text_color, bbox=bbox_props)
                    else:
                        ax.text(x + w/2, y + row_h/2, display_val, ha='center', va='center', fontsize=12 if not is_header else 13, fontweight=font_weight, color=text_color)

                if is_fastest and not is_header:
                    text_len = len(display_val) * 0.17 
                    ax.plot([x + w/2 - text_len/2, x + w/2 + text_len/2], [y + row_h/2 - 0.20, y + row_h/2 - 0.20], color='black', linewidth=1.5)

        group_boundaries = [0, cum_w[1], cum_w[2], cum_w[3], cum_w[9], cum_w[15], cum_w[21], total_w]
        for gb in group_boundaries:
            ax.plot([gb, gb], [0, total_h], color='black', linewidth=2.0 if gb in [0, total_w] else 1.2)
            
        for i in range(6):
            y_line = total_h - 2*row_h - i*8*row_h
            ax.plot([0, total_w], [y_line, y_line], color='black', linewidth=2.0 if i in [0, 5] else 1.5)

        ax.plot([0, total_w], [total_h, total_h], color='black', linewidth=2) 
        ax.plot([0, total_w], [total_h - row_h, total_h - row_h], color='black', linewidth=1) 

        legend_y = -0.8
        legend_x = 0.5
        ax.text(legend_x, legend_y, "The fastest time for each year dataset is in", fontsize=12, va='center')
        legend_x += 7.3
        ax.text(legend_x, legend_y, "bold and underlined", fontsize=12, fontweight='bold', va='center')
        ax.plot([legend_x - 0.1, legend_x + 3.8], [legend_y - 0.2, legend_y - 0.2], color='black', linewidth=1.5)
        legend_x += 4.2

        boxes = [('#78d978', '< 3x the fastest'), ('#a3e6a3', '< 10x the fastest'), ('#d9f2d9', '< 35x the fastest'), ('#e0e0e0', '< 100x the fastest'), ('#b0b0b0', '> 100x the fastest')]
        for color, label in boxes:
            rect = patches.Rectangle((legend_x, legend_y - 0.2), 0.8, 0.4, linewidth=0.5, edgecolor='#cccccc', facecolor=color)
            ax.add_patch(rect)
            legend_x += 1.0
            ax.text(legend_x, legend_y, f": {label}", fontsize=12, va='center')
            legend_x += 3.8
            
        ax.text(0.5, -1.6, "* Indicates the algorithm triggered a full tree compaction. The specific phase is tagged as [Ins], [Del], or [Both].", ha='left', va='center', fontsize=12, style='italic', color='#555555')
        ax.text(total_w/2, total_h + 0.5, "Running time (in milliseconds) on Bhutan dataset across 5 years of evolution. Lower is better.", ha='center', va='bottom', fontsize=16, fontweight='bold', color='#333333')

        plt.tight_layout()
        plt.savefig("Paper_Table_All_Years_Combined.png", dpi=300, bbox_inches='tight')
        plt.close()
        logger.info("  [Table] Successfully generated all publication tables!")
    except Exception as e:
        logger.error(f"Failed to generate super combined table: {e}")


if __name__ == "__main__":
    run_benchmark()
    parse_results()
