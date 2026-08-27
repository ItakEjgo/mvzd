import re

with open("run_experiments.py", "r") as f:
    text = f.read()

# Replace generate_plots
start_plot = text.find("def generate_plots(bp, year, res_dir):")
end_plot = text.find("def print_interim_summary(bp, year, res_dir):")

new_global_plots = """def generate_global_plots(bp, res_dir):
    out_dir = f"plots_and_tables/bp_{bp}"
    os.makedirs(out_dir, exist_ok=True)
    
    # Timeline sequences
    fwd_years = [2021, 2022, 2023, 2024, 2025]
    bwd_years = [2025, 2024, 2023, 2022, 2021]
    
    # 1. Memory Evolution
    plt.figure(figsize=(16, 6))
    for algo in algos:
        mem_series = []
        boundaries = []
        current_len = 0
        
        # Forward
        for y in fwd_years:
            df = load_combined_data(y, res_dir, "Forward_Batch")
            if not df.empty and algo in df['Algo'].values:
                vals = df[(df['Algo'] == algo) & (df['Batch'].astype(str).str.strip() != 'SUMMARY')]['Mem_MB'].values
                mem_series.extend(vals)
                current_len += len(vals)
                boundaries.append(current_len)
        
        pivot = current_len
        
        # Backward
        for y in bwd_years:
            df = load_combined_data(y, res_dir, "Backward_Batch")
            if not df.empty and algo in df['Algo'].values:
                vals = df[(df['Algo'] == algo) & (df['Batch'].astype(str).str.strip() != 'SUMMARY')]['Mem_MB'].values
                mem_series.extend(vals)
                current_len += len(vals)
                boundaries.append(current_len)
                
        if len(mem_series) > 0:
            linestyle = '--' if 'Rlog' in algo else '-'
            linewidth = 1.5 if 'Rlog' in algo else 2.5
            plt.plot(mem_series, label=algo, linestyle=linestyle, linewidth=linewidth)
            
    # Draw boundaries
    if len(boundaries) > 0:
        for b in boundaries[:-1]:
            if b == pivot:
                plt.axvline(x=b - 1, color='red', linestyle='-', linewidth=2, label='Start Deletion')
            else:
                plt.axvline(x=b - 1, color='grey', linestyle=':', alpha=0.5)
                
    plt.gca().xaxis.set_major_locator(MaxNLocator(integer=True))
    plt.title(f"Global Logical Memory Evolution (BP={bp}%)")
    plt.xlabel("Continuous Batch Sequence (2021 -> 2025 -> 2021)")
    plt.ylabel("Logical Memory Size (MB)")
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(f"{out_dir}/Global_bp{bp}_Fig1_Memory.png", dpi=300)
    plt.close()

    # 2. Update Time Evolution
    fig, axes = plt.subplots(1, 2, figsize=(18, 6))
    
    for algo in algos:
        fwd_times = []
        bwd_times = []
        for y in fwd_years:
            df = load_combined_data(y, res_dir, "Forward_Batch")
            if not df.empty and algo in df['Algo'].values:
                algo_df = df[df['Algo'] == algo]
                fwd_times.append(algo_df['Fork_ms'].sum() + algo_df['Commit_ms'].sum() + algo_df['Merge_ms'].sum())
            else: fwd_times.append(0)
            
        for y in bwd_years:
            df = load_combined_data(y, res_dir, "Backward_Batch")
            if not df.empty and algo in df['Algo'].values:
                algo_df = df[df['Algo'] == algo]
                bwd_times.append(algo_df['Fork_ms'].sum() + algo_df['Commit_ms'].sum() + algo_df['Merge_ms'].sum())
            else: bwd_times.append(0)
            
        linestyle = '--' if 'Rlog' in algo else '-'
        linewidth = 1.5 if 'Rlog' in algo else 2.5
        axes[0].plot(fwd_years, fwd_times, marker='o', label=algo, linestyle=linestyle, linewidth=linewidth)
        axes[1].plot(bwd_years, bwd_times, marker='o', label=algo, linestyle=linestyle, linewidth=linewidth)
        
    axes[0].set_title("Forward Phase Update Time")
    axes[0].set_xlabel("Year")
    axes[0].set_ylabel("Total Time (ms)")
    axes[0].grid(True, alpha=0.3)
    axes[0].set_yscale('log')
    
    axes[1].set_title("Backward Phase Update Time")
    axes[1].set_xlabel("Year")
    axes[1].grid(True, alpha=0.3)
    axes[1].set_yscale('log')
    axes[1].invert_xaxis()
    axes[1].legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    
    plt.suptitle(f"Global Update Time Evolution (BP={bp}%)")
    plt.tight_layout()
    plt.savefig(f"{out_dir}/Global_bp{bp}_Fig2_Update.png", dpi=300)
    plt.close()

    # 3. Query Lifecycle
    phases_ordered = [("BaseQuery", 0)]
    for y in fwd_years: phases_ordered.append(("MidQuery", y))
    for y in bwd_years: phases_ordered.append(("EndQuery", y))
    
    x_labels = ["Base"] + [f"{y} Mid" for y in fwd_years] + [f"{y} End" for y in bwd_years]
    
    for qt in query_types:
        plt.figure(figsize=(16, 6))
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
            
            linestyle = '--' if 'Rlog' in algo else '-'
            linewidth = 1.5 if 'Rlog' in algo else 2.5
            plt.plot(range(len(x_labels)), latencies, marker='o', label=algo, linestyle=linestyle, linewidth=linewidth)
            
        plt.xticks(range(len(x_labels)), x_labels, rotation=45)
        plt.yscale('log')
        plt.title(f"Global {qt} Latency Lifecycle (BP={bp}%)")
        plt.ylabel("Average Time (ms) - Log Scale")
        plt.grid(True, alpha=0.3)
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.tight_layout()
        plt.savefig(f"{out_dir}/Global_bp{bp}_Fig3_Query_{qt}.png", dpi=300)
        plt.close()

"""
text = text[:start_plot] + new_global_plots + text[end_plot:]

# Replace the loop call
old_call = """        for year in [2021, 2022, 2023, 2024, 2025]:
            print_interim_summary(bp, year, dest)
            generate_plots(bp, year, dest)
            logger.info(f"  [Plot] Successfully generated visual plots for Year {year} (BP={bp}%) in 'plots_and_tables/'")"""

new_call = """        for year in [2021, 2022, 2023, 2024, 2025]:
            print_interim_summary(bp, year, dest)
        
        generate_global_plots(bp, dest)
        logger.info(f"  [Plot] Successfully generated GLOBAL visual plots for BP={bp}% in 'plots_and_tables/'")"""

text = text.replace(old_call, new_call)

with open("run_experiments.py", "w") as f:
    f.write(text)

