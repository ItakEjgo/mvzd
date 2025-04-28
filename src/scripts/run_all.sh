#!/bin/bash

data_dir="/colddata/bhuan102/quadtree"
data_type="varden"
data_folder="10000000_2"
# data_folder="100000000_2"
# data_folder="1000000000_2"
pkd_old_dir="/data/bhuan102/KDtree/build"
pkd_new_dir="/data/bhuan102/SpaceTreeLib/build"

mkdir -p ../exp_res/$data_type/$data_folder
data_folder=$data_type/$data_folder

make

export PARLAY_NUM_THREADS

batch_percent=10
# mkdir -p ../exp_res/$data_folder/bach_$batch_percent

rangeSize=0
for rangeSize in 0 1 2; do
    for file in `ls $data_dir/$data_folder`; do
        outdir=../exp_res/$data_folder/region_${rangeSize}/
        mkdir -p $outdir
        #### range count experiment
        #   pkd-tree
        numactl -i all $pkd_new_dir/test -p $data_dir/$data_folder/$file -d 2 -t 0 -q 4 -r 3 -i 1 -s $rangeSize > $outdir/${file}_pkd_new_range_count.time
        mv pkd_range_count.txt $outdir/${file}_pkd_new_range_count.res

        numactl -i all $pkd_old_dir/test -p $data_dir/$data_folder/$file -d 2 -t 0 -q 4 -r 3 -i 1 -s $rangeSize > $outdir/${file}_pkd_old_range_count.time
        mv pkd_range_count.txt $outdir/${file}_pkd_old_range_count.res

        diff $outdir/${file}_pkd_old_range_count.res $outdir/${file}_pkd_new_range_count.res

        #### range report experiment
        #   pkd-tree
        numactl -i all $pkd_new_dir/test -p $data_dir/$data_folder/$file -d 2 -t 0 -q 8 -r 3 -i 1 -s $rangeSize > $outdir/${file}_pkd_new_range_report.time
        mv pkd_range_report.txt $outdir/${file}_pkd_new_range_report.res

        numactl -i all $pkd_old_dir/test -p $data_dir/$data_folder/$file -d 2 -t 0 -q 8 -r 3 -i 1 -s $rangeSize > $outdir/${file}_pkd_old_range_report.time
        mv pkd_range_report.txt $outdir/${file}_pkd_old_range_report.res

        diff $outdir/${file}_pkd_old_range_report.res $outdir/${file}_pkd_new_range_report.res
    done
done