#!/bin/bash

data_dir="/data/bhuan102/quadtree"
data_type="real"
data_folder="japan2014"
pkd_dir="/data/bhuan102/KDtree/build"

make

mkdir -p ../exp_res/$data_type/$data_folder
data_folder=$data_type/$data_folder

export PARLAY_NUM_THREADS

mkdir -p ../exp_res/$data_folder/

# Real Test
for file in `ls $data_dir/$data_folder`
do
    outdir=../exp_res/$data_folder/${file}
    mkdir -p $outdir
    numactl -i all ./main $data_dir/$data_folder/$file 4 10 1 > ${outdir}/${file}_zd_size.time
    # mv *.res $outdir
    echo "size: $file finished"
done

# range report experiment
# file="1.in"
# for file in `ls $data_dir/$data_folder`
# do
#     $pkd_dir/test -p $data_dir/$data_folder/$file -d 2 -t 0 -q 8 -r 3 -i 1 > ../exp_res/$data_folder/${file}_pkd_range_report.time
#     mv pkd_range_report.txt ../exp_res/$data_folder/${file}_pkd_range_report.res

#     ./main $data_dir/$data_folder/$file 8 > ../exp_res/$data_folder/${file}_zd_range_report.time
#     mv zd_range_report.txt ../exp_res/$data_folder/${file}_zd_range_report.res

#     echo "range report: $file"
#     diff ../exp_res/$data_folder/${file}_pkd_range_report.res ../exp_res/$data_folder/${file}_zd_range_report.res
# done

# knn report test
# for file in `ls $data_dir/$data_folder`
# do
#     $pkd_dir/test -p $data_dir/$data_folder/$file -d 2 -t 0 -q 1 -r 3 -i 1 > ../exp_res/$data_folder/${file}_pkd_knn_report.time
#     # mv pkd_knn_report.txt ../exp_res/$data_folder/${file}_pkd_knn_report.res

#     ./main $data_dir/$data_folder/$file 1 > ../exp_res/$data_folder/${file}_zd_knn_report.time
#     mv zd_knn_report.txt ../exp_res/$data_folder/${file}_zd_knn_report.res

#     echo "knn: $file"
# done