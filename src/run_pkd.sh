#!/bin/bash


data_dir_varden="/data/bhuan102/mvzd-data-processing/processed_data/synthetic_data/varden"
data_dir_uniform="/data/bhuan102/mvzd-data-processing/processed_data/synthetic_data/uniform"

file="1.in"
dir=$1

pkd_new_dir="../baselines/SpaceTreeLib/build"

export PARLAY_NUM_THREADS


for i in 1 2 4 8 16 32 64 128; do
    export PARLAY_NUM_THREADS=$i
    echo "---------------EXP_SPLIT--------------------------"
    echo "using $i threads"
    numactl -i all $pkd_new_dir/test -p $data_dir_varden/$dir/$file -d 2 -t 0 -q 4 -r 3 -i 1
done