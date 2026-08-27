import re

with open("run_experiments.py", "r") as f:
    text = f.read()

text = re.sub(r'bps\s*=\s*.*', 'bps = [0.01, 0.1, 1, 10]', text)

with open("run_experiments.py", "w") as f:
    f.write(text)

