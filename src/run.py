import subprocess

# Define your commands
commands = [
    # ["python3.11", "read_pbf.py", '14']
    ["bash", "run_latency_insert.sh", '10000_2'],
    ["bash", "run_latency_insert.sh", '100000_2'],
    ["bash", "run_latency_insert.sh", '1000000_2'],
    ["bash", "run_latency_delete.sh", '1000000_2'],
    ["bash", "run_latency_delete.sh", '100000_2'],
    ["bash", "run_latency_delete.sh", '10000_2']
]

# Start each command in parallel and capture stdout and stderr
processes = []
for cmd in commands:
    processes.append(subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE))

# Wait for all processes to complete and capture the output
for p in processes:
    stdout, stderr = p.communicate()  # Capture output and errors
    print(f"Output:\n{stdout.decode()}")  # Decode and print stdout
    if stderr:
        print(f"Error: {stderr.decode()}")  # Decode and print stderr if any

print("All commands finished.")
