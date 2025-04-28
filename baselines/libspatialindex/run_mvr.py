#!/usr/bin/env python3
import argparse
import subprocess

synthetic_dir = '/data/bhuan102/mvzd-data-processing/processed_data/synthetic_data'
real_dir = '/data/bhuan102/mvzd-data-processing/processed_data/real_data'

range_count_query = '/data/bhuan102/mvzd-data-processing/range_count_query'

data_type= 'varden'
data_size = '10000_2'
data_name = '1.in'

def comb(dir, data_type, data_size, file_name):
    return dir + '/' + data_type + '/' + data_size + '/' + data_name

def main():
    parser = argparse.ArgumentParser(
        description='Run the binary "main" with various command line options.'
    )
    parser.add_argument('-i', metavar='<Path-to-Input>', type=str,
                        help='Path to the input file', default=comb(synthetic_dir, data_type, data_size, data_name))
    parser.add_argument('-t', metavar='<Task-Name>', type=str,
                        help='Name of the task', default='build')
    parser.add_argument('-b', metavar='<Path-to-Batch-File>', type=str,
                        help='Path to the batch file')
    parser.add_argument('-bf', metavar='<batch-fraction>', type=float,
                        help='Batch fraction', default='10')
    parser.add_argument('-r', metavar='<Path-to-Range-Query>', type=str,
                        help='Path to the range query', default=comb(range_count_query, data_type, data_size, data_name + '-2.qry'))
    parser.add_argument('-real', metavar='<Is-Real-Dataset?>', type=str,
                        help='Is it a real dataset? (e.g., yes/no)', default='0')

    args = parser.parse_args()

    # Build the command list to execute
    cmd = ['./main']
    if args.i:
        cmd.extend(['-i', args.i])
    if args.t:
        cmd.extend(['-t', args.t])
    if args.b:
        cmd.extend(['-b', args.b])
    if args.bf is not None:
        cmd.extend(['-bf', str(args.bf)])
    if args.r:
        cmd.extend(['-r', args.r])
    if args.real:
        cmd.extend(['-real', args.real])

    print("Running command:", " ".join(cmd))
    
    try:
        # Use stdout and stderr instead of capture_output for Python 3.6
        result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                universal_newlines=True, check=True)
        print("Standard Output:\n", result.stdout)
        if result.stderr:
            print("Standard Error:\n", result.stderr)
    except subprocess.CalledProcessError as e:
        print("An error occurred while running the command:")
        print(e)

if __name__ == '__main__':
    main()
