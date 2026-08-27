import re

with open("run_experiments.py", "r") as f:
    text = f.read()

# Remove the try/except generate_plots block
text = re.sub(r'            try:\n                generate_plots\(bp, year, dest\).*?logger\.error.*?\{e\}\"\)', '', text, flags=re.DOTALL)

# Insert the global plot call after the year loop
old_loop_end = """        for year in years:
            print_interim_summary(bp, year, dest)"""

new_loop_end = """        for year in years:
            print_interim_summary(bp, year, dest)
            
        try:
            generate_global_plots(bp, dest)
            logger.info(f"  [Plot] Successfully generated GLOBAL visual plots for BP={bp}% in 'plots_and_tables/bp_{bp}/'")
        except Exception as e:
            logger.error(f"  [Plot Error] Failed to generate global plots for BP={bp}%: {e}")"""

text = text.replace(old_loop_end, new_loop_end)

with open("run_experiments.py", "w") as f:
    f.write(text)
