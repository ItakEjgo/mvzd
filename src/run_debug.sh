#!/bin/bash

data_dir_varden="/data/bhuan102/mvzd-data-processing/processed_data/synthetic_data/varden"
data_dir_uniform="/data/bhuan102/mvzd-data-processing/processed_data/synthetic_data/uniform"
data_dir_real="/data/bhuan102/mvzd-data-processing/processed_data/real_data/median"

range_count_dir_varden="/data/bhuan102/mvzd-data-processing/range_count_query/varden"
range_count_dir_uniform="/data/bhuan102/mvzd-data-processing/range_count_query/uniform"

MVR_dir="/data/bhuan102/libspatialindex"
pkd_new_dir="/data/bhuan102/SpaceTreeLib/build"
    
make

export PARLAY_NUM_THREADS
# export PARLAY_NUM_THREADS=1
file="1.in"
dir=$1

numactl -i all ./main -i $data_dir_varden/$dir/$file -t spatial-diff -a combined -r $range_count_dir_varden/$dir/1.in-0.qry
# numactl -i all ./main -i $data_dir_varden/$dir/$file -t diff -a combined -r $range_count_dir_varden/$dir/1.in-0.qry

# for i in 1 2 4 8 16 24 48 96; do
#     export PARLAY_NUM_THREADS=$i
#     echo "---------------EXP_SPLIT--------------------------"
#     echo "using $i threads"
#     numactl -i all ./main -i $data_dir_varden/$dir/$file -t diff -a combined
# done

# echo "processing varden-$dir"
# LD_PRELOAD=/usr/local/lib/libjemalloc.so.2 numactl -i all ./main -i $data_dir_varden/$dir/$file -t debug -a mvzd
# echo "processing uniform-$dir"
# LD_PRELOAD=/usr/local/lib/libjemalloc.so.2 numactl -i all ./main -i $data_dir_uniform/$dir/$file -t debug -a mvzd
# numactl -i all ./main -i $data_dir_uniform/$dir/$file -t diff -a mvzd

# for file in `ls $data_dir_real`; do
#     echo "dealing with $file:"
#     LD_PRELOAD=/usr/local/lib/libjemalloc.so.2 numactl -i all ./main -i $data_dir_real/$file -t debug -a mvzd -real 1
# done

# for i in 1 2 4 8 16 24 48 96; do
# only use 1 core
    # i=96
    # export PARLAY_NUM_THREADS=$i 
    # echo "[INFO] Running with Threads: $i"
# export PARLAY_NUM_THREADS
    # file="1.in"
    # dir=$1
    # Tree construction
    # echo "Build Latency Test:"
    # echo "---------------EXP_SPLIT--------------------------"
    # echo "Varden Results"
    # LD_PRELOAD=/usr/local/lib/libjemalloc.so.2 numactl -i all ./main -i $data_dir_varden/$dir/$file -t debug -a mvzd

    # echo "---------------EXP_SPLIT--------------------------"
    # echo "Uniform Results"

    # LD_PRELOAD=/usr/local/lib/libjemalloc.so.2 numactl -i all ./main -i $data_dir_uniform/$dir/$file -t debug -a mvzd
    
# done
