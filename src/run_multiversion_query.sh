#!/bin/bash

#bhutan dir
data_dir="/data/bhuan102/mvzd-data-processing/raw_data/real_data/coordinate_data/bhutan-180101.txt"
diff_dir="/data/bhuan102/mvzd-data-processing/processed_data/real_data/bhutan/"

#japan dir
# data_dir="/data/bhuan102/mvzd-data-processing/raw_data/real_japan/coordinate_data/japan-140101.txt"
# diff_dir="/data/bhuan102/mvzd-data-processing/processed_data/real_data/multi-version/"

MVR_dir="/data/bhuan102/cpam-quadtree/CPAM/baselines/libspatialindex"
Rtree_dir="/data/bhuan102/cpam-quadtree/CPAM/baselines/boostRtree"
    
make

export PARLAY_NUM_THREADS
# echo "[CPAM]:"
# numactl -i all ./main -i $data_dir -a combined -mv $diff_dir -real 1 -t multi-version-test

# numactl -i all ./main -i $data_dir -a cpambb -mv $diff_dir -real 1 -t multi-version-test > debug/bhutan-cpambb-case-study.log
# numactl -i all ./main -i $data_dir -a cpamz -mv $diff_dir -real 1 -t multi-version-test > debug/bhutan-cpamz-case-study.log
# echo "[ZDTree]:"
numactl -i all ./main -i $data_dir -a mvzd -mv $diff_dir -real 1 -t multi-version-test > debug/bhutan-zdtree-case-study-new.log
# echo "[boostRTree]:"
# numactl -i all $Rtree_dir/main -i $data_dir -a mvzd -mv $diff_dir -real 1 -t multi-version-test > debug/bhutan-rtree-case-study-japan.log
# echo "[MVRTree]:"
# numactl -i all $MVR_dir/main -i $data_dir -a mvrtree -mv $diff_dir -real 1 -t multi-version-test > debug/bhutan-mvrtree-case-study-japan.log
# echo "[MV3RTree]:"
# numactl -i all $MVR_dir/main -i $data_dir -a mv3rtree -mv $diff_dir -real 1 -t multi-version-test > debug/bhutan-mv3rtree-case-study-japan.log

# export PARLAY_NUM_THREADS=1
# echo "[CPAM]:"
# numactl -i all ./main -i $data_dir -a cpambb -mv $diff_dir -real 1 -t multi-version-test > debug/bhutan-cpambb-case-study-sequential.log
# numactl -i all ./main -i $data_dir -a cpamz -mv $diff_dir -real 1 -t multi-version-test > debug/bhutan-cpamz-case-study-sequential.log
# echo "[ZDTree Sequential]:"
# numactl -i all ./main -i $data_dir -a mvzd -mv $diff_dir -real 1 -t multi-version-test > debug/bhutan-zdtree-case-study-sequential-japan.log