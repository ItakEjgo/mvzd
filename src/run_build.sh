#!/bin/bash

data_dir_varden="/data/bhuan102/mvzd-data-processing/processed_data/synthetic_data/varden"
data_dir_uniform="/data/bhuan102/mvzd-data-processing/processed_data/synthetic_data/uniform"

file="1.in"
dir=$1

pkd_new_dir="/data/bhuan102/SpaceTreeLib/build"

range_count_dir_varden="/data/bhuan102/mvzd-data-processing/range_count_query/varden"
range_count_dir_uniform="/data/bhuan102/mvzd-data-processing/range_count_query/uniform"
    
make

file="1.in"
dir=$1
# for i in 1 2 4 8 16 32 64; do
# for i in 1 2 4 8 16 24 48 96; do
for i in 1 2 4 8 16 32 64 128; do
    export PARLAY_NUM_THREADS=$i
    echo "---------------EXP_SPLIT--------------------------"
    echo "using $i threads" 
    numactl -i all ./main -i $data_dir_varden/$dir/$file -t build -a mvzd
done