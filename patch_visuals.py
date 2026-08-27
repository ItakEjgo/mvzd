import re

with open("run_experiments.py", "r") as f:
    text = f.read()

# Replace the Memory plot section
old_mem_plot = """    # Draw boundaries
    if len(boundaries) > 0:
        for b in boundaries[:-1]:
            if b == pivot:
                plt.axvline(x=b - 1, color='red', linestyle='-', linewidth=2, label='Start Deletion')
            else:
                plt.axvline(x=b - 1, color='grey', linestyle=':', alpha=0.5)
                
    plt.gca().xaxis.set_major_locator(MaxNLocator(integer=True))
    plt.title(f"Global Logical Memory Evolution (BP={bp}%)")
    plt.xlabel("Continuous Batch Sequence (2021 -> 2025 -> 2021)")"""

new_mem_plot = """    # Draw boundaries and custom X-ticks
    xticks_pos = []
    xtick_labels = []
    if len(boundaries) > 0:
        prev_b = 0
        for i, b in enumerate(boundaries):
            if b == pivot:
                plt.axvline(x=b - 1, color='red', linestyle='-', linewidth=2, label='Start Deletion')
            else:
                plt.axvline(x=b - 1, color='grey', linestyle=':', alpha=0.5)
            
            # Midpoint for the label
            mid = (prev_b + b - 1) / 2.0
            xticks_pos.append(mid)
            
            # Construct label
            batch_count = b - prev_b
            if i < 5: 
                label = f"{fwd_years[i]} Fwd\\n({batch_count} batches)"
            else: 
                label = f"{bwd_years[i-5]} Bwd\\n({batch_count} batches)"
            xtick_labels.append(label)
            
            prev_b = b
            
    plt.xticks(xticks_pos, xtick_labels, rotation=45, fontsize=9)
    plt.title(f"Global Logical Memory Evolution (BP={bp}%)")
    plt.xlabel("Timeline (Years & Batches)")"""

text = text.replace(old_mem_plot, new_mem_plot)

# Replace Update plot to include Cumulative
old_update = """    # 2. Update Time Evolution
    fig, axes = plt.subplots(1, 2, figsize=(18, 6))"""

new_update = """    # 2. Update Time Evolution (Includes Cumulative)
    fig, axes = plt.subplots(1, 3, figsize=(24, 6))"""

text = text.replace(old_update, new_update)

old_update_loop = """        axes[0].plot(fwd_years, fwd_times, marker='o', label=algo, linestyle=linestyle, linewidth=linewidth)
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
    plt.close()"""

new_update_loop = """        axes[0].plot(fwd_years, fwd_times, marker='o', label=algo, linestyle=linestyle, linewidth=linewidth)
        axes[1].plot(bwd_years, bwd_times, marker='o', label=algo, linestyle=linestyle, linewidth=linewidth)
        
        # Cumulative total
        cumulative_times = []
        acc = 0
        for t in fwd_times + bwd_times:
            acc += t
            cumulative_times.append(acc)
        axes[2].plot(range(10), cumulative_times, marker='s', label=algo, linestyle=linestyle, linewidth=linewidth)
        
    axes[0].set_title("Forward Phase Update Time (Spiky)")
    axes[0].set_xlabel("Year")
    axes[0].set_ylabel("Time (ms) - Log")
    axes[0].grid(True, alpha=0.3)
    axes[0].set_yscale('log')
    
    axes[1].set_title("Backward Phase Update Time (Spiky)")
    axes[1].set_xlabel("Year")
    axes[1].grid(True, alpha=0.3)
    axes[1].set_yscale('log')
    axes[1].invert_xaxis()
    
    x_labels_cum = [f"{y}F" for y in fwd_years] + [f"{y}B" for y in bwd_years]
    axes[2].set_title("Cumulative Update Time (Smooth Overhead)")
    axes[2].set_xticks(range(10))
    axes[2].set_xticklabels(x_labels_cum, rotation=45)
    axes[2].set_ylabel("Accumulated Time (ms) - Log")
    axes[2].grid(True, alpha=0.3)
    axes[2].set_yscale('log')
    axes[2].legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    
    plt.suptitle(f"Global Update Time Evolution (BP={bp}%)")
    plt.tight_layout()
    plt.savefig(f"{out_dir}/Global_bp{bp}_Fig2_Update.png", dpi=300)
    plt.close()"""

text = text.replace(old_update_loop, new_update_loop)

with open("run_experiments.py", "w") as f:
    f.write(text)
