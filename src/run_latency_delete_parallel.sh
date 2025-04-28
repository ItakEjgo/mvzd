#!/bin/bash

data_dir_varden="/data/bhuan102/mvzd-data-processing/SpaceTreeLib/build/generated_data/ss_varden_bigint"
data_dir_uniform="/data/bhuan102/mvzd-data-processing/SpaceTreeLib/build/generated_data/uniform_bigint"
MVR_dir="/data/bhuan102/libspatialindex"
    
make

# only use 1 core
export PARLAY_NUM_THREADS=1 
# export PARLAY_NUM_THREADS
file="1.in"

dir=$1
percent=$2
# Tree construction
# echo "Build Latency Test:"
# for dir in `ls $data_dir_varden`; do
#     # build test
#     echo "---------------EXP_SPLIT--------------------------"
#     echo "build latency (sequential) for varden-$dir:"
#     numactl -i all ./main -i $data_dir_varden/$dir/$file -t build -a combined
#     $MVR_dir/main -i $data_dir_varden/$dir/$file -t build
# done

# for dir in `ls $data_dir_uniform`; do
#     # build test
#     echo "---------------EXP_SPLIT--------------------------"
#     echo "build latency (sequential) for uniform-$dir:"
#     numactl -i all ./main -i $data_dir_uniform/$dir/$file -t build -a combined
#     $MVR_dir/main -i $data_dir_uniform/$dir/$file -t build
# done

# Batch insertion
# echo "Batch-insert Latency Test:"
# # for dir in `ls $data_dir_varden`; do
#     # batch-insert test
#     echo "---------------EXP_SPLIT--------------------------"
#     echo "insert latency (sequential) for varden-$dir-$percent%:"
#     numactl -i all ./main -i $data_dir_varden/$dir/$file -t batch-insert -a combined -bf $percent
#     $MVR_dir/main -i $data_dir_varden/$dir/$file -t batch-insert -bf $percent
# # done

# # for dir in `ls $data_dir_uniform`; do
#     # batch-insert test
#     echo "---------------EXP_SPLIT--------------------------"
#     echo "insert latency (sequential) for uniform-$dir-$percent%:"
#     numactl -i all ./main -i $data_dir_uniform/$dir/$file -t batch-insert -a combined -bf $percent
#     $MVR_dir/main -i $data_dir_uniform/$dir/$file -t batch-insert -bf $percent
# # done

# Batch deletion
# echo "Batch-delete Latency Test:"
# for dir in `ls $data_dir_varden`; do
    # batch-delete test
    echo "---------------EXP_SPLIT--------------------------"
    echo "delete latency (sequential) for varden-$dir-$percent%:"
    # for percent in 1 2 3 4 5 6 7 8 9 10; do
    numactl -i all ./main -i $data_dir_varden/$dir/$file -t batch-delete -a combined -bf $percent
    $MVR_dir/main -i $data_dir_varden/$dir/$file -t batch-delete -bf $percent
    # done
# done

# for dir in `ls $data_dir_uniform`; do
    # batch-delete test
    echo "---------------EXP_SPLIT--------------------------"
    echo "delete latency (sequential) for uniform-$dir-$percent%:"
    numactl -i all ./main -i $data_dir_uniform/$dir/$file -t batch-delete -a combined -bf $percent
    $MVR_dir/main -i $data_dir_uniform/$dir/$file -t batch-delete -bf $percent
# done
