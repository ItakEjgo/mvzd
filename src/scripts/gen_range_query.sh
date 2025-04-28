#!/bin/bash

# data_dir="/colddata/bhuan102/quadtree"
# data_dir="/data/bhuan102/quadtree"

data_dir="/data/bhuan102/mvzd-data-processing/SpaceTreeLib/build/generated_data"

# data_type="varden"
# data_type="ss_varden_bigint"
data_type="uniform_bigint"

# data_folder="10000_2"
# data_folder="100000_2"
# data_folder="1000000_2"
data_folder="10000000_2"
# data_folder="10000000_2"
# data_folder="1000000000_2"

# pkd_dir="/data/bhuan102/SpaceTreeLib/build"
pkd_dir="/data/bhuan102/mvzd-data-processing/SpaceTreeLib/build"


export PARLAY_NUM_THREADS
query_dir="/data/bhuan102/mvzd-data-processing"

range_count_dir=$query_dir/range_count_query/$data_type
range_report_dir=$query_dir/range_report_query/$data_type

echo $range_count_dir
echo $range_report_dir

mkdir -p $range_count_dir/$data_folder
mkdir -p $range_report_dir/$data_folder

for rangeSize in 0 1 2; do
    file="1.in"
# for rangeSize in 0; do
    # for file in `ls $data_dir/$data_type/$data_folder`; do
        #### range count experiment
        #   pkd-tree
        numactl -i all $pkd_dir/test -p $data_dir/$data_type/$data_folder/$file -d 2 -t 0 -x 1 -q 2 -xz $rangeSize >> log.txt
        cp range_count2.qry $range_report_dir/$data_folder/$file-$rangeSize.qry
        mv range_count2.qry $range_count_dir/$data_folder/$file-$rangeSize.qry

        # rm pkd_range_count.txt
        # echo $range_count_dir/$data_folder/$file-$rangeSize.qry
        # mv pkd_range_count.txt $outdir/${file}_pkd_range_count.res

        # numactl -i all $pkd_dir/test -p $data_dir/$data_type/$data_folder/$file -d 2 -t 0 -x 1 -q 4 -xz $rangeSize >> log.txt       
         # mv range_report2.qry range_report.qry
        # mv range_report2.qry $range_report_dir/$data_folder/$file-$rangeSize-01M.qry
        # rm pkd_range_report.txt
        # mv pkd_range_report.txt $outdir/${file}_pkd_range_report.res
        #   zd-tree
        # numactl -i all ./main $data_dir/$data_folder/$file 4 $batch_percent 0 > $outdir/${file}_zd_range_count.time
        # mv zd_range_count.txt $outdir/${file}_zd_range_count.res
        # diff $outdir/${file}_pkd_range_count.res $outdir/${file}_zd_range_count.res 
    # done
done