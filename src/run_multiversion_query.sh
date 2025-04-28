#!/bin/bash

data_dir_varden="/data/bhuan102/mvzd-data-processing/processed_data/synthetic_data/varden"
data_dir_uniform="/data/bhuan102/mvzd-data-processing/processed_data/synthetic_data//uniform"
MVR_dir="/data/bhuan102/libspatialindex"

pkd_new_dir="/data/bhuan102/SpaceTreeLib/build"

range_count_dir_varden="/data/bhuan102/mvzd-data-processing/range_count_query/varden"
range_count_dir_uniform="/data/bhuan102/mvzd-data-processing/range_count_query/uniform"
    
make

export PARLAY_NUM_THREADS
file="1.in"
dir=$1
echo "---------------EXP_SPLIT--------------------------"
# numactl -i all $pkd_new_dir/test -p $data_dir_varden/$dir/$file -d 2 -t 0 -q 4 -r 3 -i 1 -s 0
# echo "Varden Results"
# # LD_PRELOAD=/usr/local/lib/libjemalloc.so.2 numactl -i all ./main -i $data_dir_varden/$dir/$file -t build -a mvzd
# numactl -i all ./main -i $data_dir_varden/$dir/$file -t range-count -a mvzd -r $range_count_dir_varden/$dir/1.in-0.qry
numactl -i all ./main -i $data_dir_varden/$dir/$file -t multi-version-query-test -a mvzd -mv $range_count_dir_varden/$dir
diff mvzd_range_count-2-on-0.txt mvzd_range_count-2-on-6.txt
# diff mvzd_range_count-2-on-0.txt mvzd_range_count-2-on-4.txt
# diff mvzd_range_count-2-on-0.txt mvzd_range_count-2-on-6.txt