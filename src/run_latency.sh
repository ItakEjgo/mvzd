#!/bin/bash

data_dir_varden="/data/bhuan102/mvzd-data-processing/SpaceTreeLib/build/generated_data/ss_varden_bigint"
data_dir_uniform="/data/bhuan102/mvzd-data-processing/SpaceTreeLib/build/generated_data/uniform_bigint"
data_dir_real="/data/bhuan102/mvzd-data-processing/processed_data/real_data"
MVR_dir="/data/bhuan102/libspatialindex"
    
make

# only use 1 core
export PARLAY_NUM_THREADS=1 
# export PARLAY_NUM_THREADS
# file="1.in"
file="mongolia-250101.osm.txt"

dir=$1
for file in `ls $data_dir_real/$dir`; do
    echo "[INFO] Dealing with $file"

    echo "Build Latency Test:"
    echo "---------------EXP_SPLIT--------------------------"
    echo "build latency (sequential) for real-$dir:"
    numactl -i all ./main -i $data_dir_real/$dir/$file -t build -a combined -real 1
    $MVR_dir/main -i $data_dir_real/$dir/$file -t build -real 1

    percent=10

    echo "Batch-insert Latency Test:"
    echo "---------------EXP_SPLIT--------------------------"
    echo "insert latency (sequential) for real-$dir-$percent:"
    numactl -i all ./main -i $data_dir_real/$dir/$file -t batch-insert -a combined -bf $percent -real 1
    $MVR_dir/main -i $data_dir_real/$dir/$file -t batch-insert -bf $percent -real 1

    echo "Batch-delete Latency Test:"
    echo "---------------EXP_SPLIT--------------------------"
    echo "delete latency (sequential) for real-$dir-$percent%:"
    numactl -i all ./main -i $data_dir_real/$dir/$file -t batch-delete -a combined -bf $percent -real 1
    $MVR_dir/main -i $data_dir_real/$dir/$file -t batch-delete -bf $percent -real 1
done

    # Tree construction
    # for dir in `ls $data_dir_varden`; do
        # build test
        # echo "Build Latency Test:"
        # echo "---------------EXP_SPLIT--------------------------"
        # echo "build latency (sequential) for varden-$dir:"
        # numactl -i all ./main -i $data_dir_varden/$dir/$file -t build -a combined -real 1
        # $MVR_dir/main -i $data_dir_varden/$dir/$file -t build -real 1

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
    # for dir in `ls $data_dir_varden`; do
    # batch-insert test
    # percent=10
    # echo "---------------EXP_SPLIT--------------------------"
    # echo "insert latency (sequential) for varden-$dir-$percent:"
    # for percent in 1 2 3 4 5 6 7 8 9 10; do
        # echo "batch $percent%:"
        # numactl -i all ./main -i $data_dir_varden/$dir/$file -t batch-insert -a combined -bf $percent -real 1
        # $MVR_dir/main -i $data_dir_varden/$dir/$file -t batch-insert -bf $percent -real 1
    # done
# done

# for dir in `ls $data_dir_uniform`; do
    # batch-insert test
    # echo "---------------EXP_SPLIT--------------------------"
    # echo "insert latency (sequential) for uniform-$dir:"
    # for percent in 1 2 3 4 5 6 7 8 9 10; do
    #     echo "batch $percent%:"
    #     numactl -i all ./main -i $data_dir_uniform/$dir/$file -t batch-insert -a combined -bf $percent
    #     $MVR_dir/main -i $data_dir_uniform/$dir/$file -t batch-insert -bf $percent
    # done
# done

# Batch deletion
    # echo "Batch-delete Latency Test:"
# for dir in `ls $data_dir_varden`; do
    # batch-delete test
    # echo "---------------EXP_SPLIT--------------------------"
    # echo "delete latency (sequential) for varden-$dir-$percent%:"
    # for percent in 1 2 3 4 5 6 7 8 9 10; do
    # echo "batch $percent%:"
    # numactl -i all ./main -i $data_dir_varden/$dir/$file -t batch-delete -a combined -bf $percent -real 1
    # $MVR_dir/main -i $data_dir_varden/$dir/$file -t batch-delete -bf $percent -real 1
    # done
# done

# # for dir in `ls $data_dir_uniform`; do
#     # batch-delete test
#     echo "---------------EXP_SPLIT--------------------------"
#     echo "delete latency (sequential) for uniform-$dir:"
#     for percent in 1 2 3 4 5 6 7 8 9 10; do
#         echo "batch $percent%:"
#         numactl -i all ./main -i $data_dir_uniform/$dir/$file -t batch-delete -a combined -bf $percent
#         $MVR_dir/main -i $data_dir_uniform/$dir/$file -t batch-delete -bf $percent
#     done
# # done
