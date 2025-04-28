import subprocess

file = "10000000_2"

# Define your commands
commands = [
    # ["python3.11", "read_pbf.py", '14']
    ["bash", "run_latency_insert_parallel.sh", file, '1'],
    ["bash", "run_latency_insert_parallel.sh", file, '2'],
    ["bash", "run_latency_insert_parallel.sh", file, '3'],
    ["bash", "run_latency_insert_parallel.sh", file, '4'],
    ["bash", "run_latency_insert_parallel.sh", file, '5'],
    ["bash", "run_latency_insert_parallel.sh", file, '6'],
    ["bash", "run_latency_insert_parallel.sh", file, '7'],
    ["bash", "run_latency_insert_parallel.sh", file, '8'],
    ["bash", "run_latency_insert_parallel.sh", file, '9'],
    ["bash", "run_latency_insert_parallel.sh", file, '10'],
    ["bash", "run_latency_delete_parallel.sh", file, '1'],
    ["bash", "run_latency_delete_parallel.sh", file, '2'],
    ["bash", "run_latency_delete_parallel.sh", file, '3'],
    ["bash", "run_latency_delete_parallel.sh", file, '4'],
    ["bash", "run_latency_delete_parallel.sh", file, '5'],
    ["bash", "run_latency_delete_parallel.sh", file, '6'],
    ["bash", "run_latency_delete_parallel.sh", file, '7'],
    ["bash", "run_latency_delete_parallel.sh", file, '8'],
    ["bash", "run_latency_delete_parallel.sh", file, '9'],
    ["bash", "run_latency_delete_parallel.sh", file, '10']
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
