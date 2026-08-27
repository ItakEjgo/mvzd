import re

with open("run_experiments.py", "r") as f:
    text = f.read()

old_code = """    cols = [c for c in cols if c in df_time.columns]
    df_time = df_time[cols]
    df_mem = df_mem[cols]
    df_nodes = df_nodes[[c for c in cols if c in df_nodes.columns]]"""

new_code = """    cols_time = [c for c in cols if c in df_time.columns]
    cols_mem = [c for c in cols if c in df_mem.columns]
    cols_nodes = [c for c in cols if c in df_nodes.columns]
    
    df_time = df_time[cols_time]
    df_mem = df_mem[cols_mem]
    df_nodes = df_nodes[cols_nodes]"""

text = text.replace(old_code, new_code)

with open("run_experiments.py", "w") as f:
    f.write(text)

