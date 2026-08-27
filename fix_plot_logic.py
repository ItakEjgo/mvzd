import re

with open("run_experiments.py", "r") as f:
    text = f.read()

# 1. Fix load_combined_data
old_load = """def load_combined_data(year, res_dir, phase):
    df_list = []
    for algo in algos:
        file_path = os.path.join(res_dir, f"{year}_{phase}_{algo}.txt")"""

new_load = """def load_combined_data(year, res_dir, phase):
    df_list = []
    for algo in algos:
        if phase == "BaseQuery":
            file_path = os.path.join(res_dir, f"BaseQuery_{algo}.txt")
        else:
            file_path = os.path.join(res_dir, f"{year}_{phase}_{algo}.txt")"""

text = text.replace(old_load, new_load)

# 2. Fix generate_plots directory
old_out = """def generate_plots(bp, year, res_dir):
    out_dir = "plots_and_tables"
    os.makedirs(out_dir, exist_ok=True)"""

new_out = """def generate_plots(bp, year, res_dir):
    out_dir = f"plots_and_tables/bp_{bp}"
    os.makedirs(out_dir, exist_ok=True)"""

text = text.replace(old_out, new_out)

with open("run_experiments.py", "w") as f:
    f.write(text)

