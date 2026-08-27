import re

with open("run_experiments.py", "r") as f:
    text = f.read()

text = re.sub(r'bps\s*=\s*\[.*?\]', 'bps = [100, 50, 25, 12.5, 6.25, 3.125]', text)

with open("run_experiments.py", "w") as f:
    f.write(text)

