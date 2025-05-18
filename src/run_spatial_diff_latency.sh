#!/bin/bash

make

Unifile="/data/bhuan102/mvzd-data-processing/processed_data/synthetic_data/uniform/100000000_2/1.in"
Uniqry0="/data/bhuan102/mvzd-data-processing/range_count_query/uniform/100000000_2/1.in-0.qry"
Uniqry1="/data/bhuan102/mvzd-data-processing/range_count_query/uniform/100000000_2/1.in-1.qry"
Uniqry2="/data/bhuan102/mvzd-data-processing/range_count_query/uniform/100000000_2/combined_selected.qry"

Varfile="/data/bhuan102/mvzd-data-processing/processed_data/synthetic_data/varden/100000000_2/1.in"
Varqry0="/data/bhuan102/mvzd-data-processing/range_count_query/varden/100000000_2/1.in-0.qry"
Varqry1="/data/bhuan102/mvzd-data-processing/range_count_query/varden/100000000_2/1.in-1.qry"
Varqry2="/data/bhuan102/mvzd-data-processing/range_count_query/varden/100000000_2/combined_selected.qry"
rtree_binary="../"

# for i in 1 2 4 8 16 24 48 96; do
#     export PARLAY_NUM_THREADS=$i 
#     echo "[INFO] Running with Threads: $i"
#     echo "[INFO] Dealing Small Region Query:"
#     numactl -i all ./main -i $file -a combined -t range-report -r $qry0
#     echo "[INFO] Dealing Median Region Query:"
#     numactl -i all ./main -i $file -a combined -t range-report -r $qry1 
# done

# exp on CPU0
PHYSICAL_CPUS=(0 4 8 12 16 20 24 28 32 36 40 44 48 52 56 60 64 68 72 76 80 84 88 92 96 100 104 108)
PHYSICAL_CORES=$(seq -s, 0 95)

i=96
export PARLAY_NUM_THREADS=$i
echo "[INFO] Running with Threads: $i"
numactl --physcpubind=$PHYSICAL_CORES -i all ./main -i $Unifile -a combined -t spatial-commit-merge -r $Uniqry2 > output/100M-U-Commit-Merge.log
numactl --physcpubind=$PHYSICAL_CORES -i all ./main -i $Varfile -a combined -t spatial-commit-merge -r $Varqry2 > output/100M-V-Commit-Merge.log



# PHYSICAL_CPUS=(0 2 4 6 8 10 12 14 16 18 20 22 24 26 28 30)

# for i in $(seq 1 28); do
# i=1
# export PARLAY_NUM_THREADS=$i 
# # CPU_LIST=$(IFS=, ; echo "${PHYSICAL_CPUS[*]:0:$i}")

# echo "[INFO] Running with Threads: $i on CPU0"

# echo "[INFO] Dealing Small Region Query:"
# numactl --membind=0 --physcpubind=0 ./test -i $Unifile -a combined -t spatial-diff -r $Uniqry2 > output/100M-U-comb-CPAMBB-Plain.log
# numactl --membind=0 --physcpubind=0 ../baselines/boostRtree/main -i $Unifile -a combined -t spatial-diff -r $Uniqry2 > output/100M-U-rtree-comb.log
# numactl --membind=0 --physcpubind=0 ./test -i $Varfile -a combined -t spatial-diff -r $Varqry2 > output/100M-V-comb-CPAMBB-Plain.log
# numactl --membind=0 --physcpubind=0 ../baselines/boostRtree/main -i $Varfile -a combined -t spatial-diff -r $Varqry2 > output/100M-V-rtree-comb.log

# echo "[INFO] Dealing Median Region Query:"
# numactl --membind=0 --physcpubind=0 ./test -i $Unifile -a combined -t spatial-diff -r $Uniqry1 > output/100M-U-1-comb.log
# numactl --membind=0 --physcpubind=0 ../baselines/boostRtree/main -i $Unifile -a combined -t spatial-diff -r $Uniqry1 > output/100M-U-1-rtree-comb.log
# numactl --membind=0 --physcpubind=0 ./test -i $Varfile -a combined -t spatial-diff -r $Varqry1 > output/100M-V-1-comb.log
# numactl --membind=0 --physcpubind=0 ../baselines/boostRtree/main -i $Varfile -a combined -t spatial-diff -r $Varqry1 > output/100M-V-1-rtree-comb.log
# done
