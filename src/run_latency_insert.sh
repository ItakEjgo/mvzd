#!/bin/bash

U_100M="/data/bhuan102/mvzd-data-processing/processed_data/synthetic_data/uniform/100000000_2/1.in"
V_100M="/data/bhuan102/mvzd-data-processing/processed_data/synthetic_data/varden/100000000_2/1.in"

MVR_dir="/data/bhuan102/cpam-quadtree/CPAM/baselines/libspatialindex"
    
make

# numactl -i all $MVR_dir/main -i $U_100M -a mvrtree -t batch-insert > debug/mvr-100M-U-insert.log
# numactl -i all $MVR_dir/main -i $U_100M -a mvrtree -t batch-delete > debug/mvr-100M-U-delete.log
# numactl -i all $MVR_dir/main -i $U_100M -a mv3rtree -t batch-insert > debug/mv3r-100M-U-insert.log
# numactl -i all $MVR_dir/main -i $U_100M -a mv3rtree -t batch-delete > debug/mv3r-100M-U-delete.log

numactl -i all $MVR_dir/main -i $V_100M -a mvrtree -t batch-insert > debug/mvr-100M-V-insert.log
numactl -i all $MVR_dir/main -i $V_100M -a mvrtree -t batch-delete > debug/mvr-100M-V-delete.log
numactl -i all $MVR_dir/main -i $V_100M -a mv3rtree -t batch-insert > debug/mv3r-100M-V-insert.log
numactl -i all $MVR_dir/main -i $V_100M -a mv3rtree -t batch-delete > debug/mv3r-100M-V-delete.log