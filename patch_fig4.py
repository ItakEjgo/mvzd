import re

with open("run_experiments.py", "r") as f:
    text = f.read()

# Find the end of Fig3 logic
target = """        plt.tight_layout()
        plt.savefig(f"{out_dir}/Global_bp{bp}_Fig3_Query_{qt}.png", dpi=300)
        plt.close()"""

new_fig4 = """        plt.tight_layout()
        plt.savefig(f"{out_dir}/Global_bp{bp}_Fig3_Query_{qt}.png", dpi=300)
        plt.close()

    # 4. Query Nodes Lifecycle
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
                        valid_plot = True
                    else: nodes.append(None)
                else: nodes.append(None)
            
            linestyle = '--' if 'Rlog' in algo else '-'
            linewidth = 1.5 if 'Rlog' in algo else 2.5
            if any(n is not None and n > 0 for n in nodes):
                plt.plot(range(len(x_labels)), nodes, marker='o', label=algo, linestyle=linestyle, linewidth=linewidth)
            
        if valid_plot:
            plt.xticks(range(len(x_labels)), x_labels, rotation=45)
            plt.yscale('log')
            plt.title(f"Global {qt} Nodes Touched Lifecycle (BP={bp}%)")
            plt.ylabel("Average Nodes Touched - Log Scale")
            plt.grid(True, alpha=0.3)
            plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
            plt.tight_layout()
            plt.savefig(f"{out_dir}/Global_bp{bp}_Fig4_Nodes_{qt}.png", dpi=300)
        plt.close()"""

text = text.replace(target, new_fig4)

with open("run_experiments.py", "w") as f:
    f.write(text)

