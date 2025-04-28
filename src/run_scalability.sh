#!/bin/bash

make

file="/data/bhuan102/mvzd-data-processing/processed_data/synthetic_data/uniform/100000000_2/1.in"
qry0="/data/bhuan102/mvzd-data-processing/range_count_query/uniform/100000000_2/1.in-0.qry"
qry1="/data/bhuan102/mvzd-data-processing/range_count_query/uniform/100000000_2/1.in-1.qry"

for i in 1 2 4 8 16 24 48 96; do
    export PARLAY_NUM_THREADS=$i 
    echo "[INFO] Running with Threads: $i"
    echo "[INFO] Dealing Small Region Query:"
    numactl -i all ./main -i $file -a combined -t range-report -r $qry0
    echo "[INFO] Dealing Median Region Query:"
    numactl -i all ./main -i $file -a combined -t range-report -r $qry1 
done

# exp on CPU0
PHYSICAL_CPUS=(0 4 8 12 16 20 24 28 32 36 40 44 48 52 56 60 64 68 72 76 80 84 88 92)
# PHYSICAL_CPUS=(0 2 4 6 8 10 12 14 16 18 20 22 24 26 28 30)

# for i in $(seq 1 28); do
#     export PARLAY_NUM_THREADS=$i 
#     CPU_LIST=$(IFS=, ; echo "${PHYSICAL_CPUS[*]:0:$i}")

#     echo "[INFO] Running with Threads: $i on CPU0"
#     echo "[INFO] Dealing Small Region Query:"
#     numactl --membind=0 --cpunodebind=0 taskset -c $CPU_LIST ./main -i $file -a combined -t range-report -r $qry0
#     echo "[INFO] Dealing Median Region Query:"
#     numactl --membind=0 --cpunodebind=0 taskset -c $CPU_LIST ./main -i $file -a combined -t range-report -r $qry1
# done
