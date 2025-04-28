make

real_japan="/data/bhuan102/mvzd-data-processing/processed_data/real_data/multi-version/japan-20140101.in"
mv_dir="/data/bhuan102/mvzd-data-processing/processed_data/real_data/multi-version/"

# Use all cores
export PARLAY_NUM_THREADS

numactl -i all ./main -i $real_japan -real 1 -a mvzd -t multi-version-test -mv $mv_dir