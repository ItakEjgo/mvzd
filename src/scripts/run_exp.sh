#!/bin/bash

data_dir="/colddata/bhuan102/quadtree"
# data_dir="/data/bhuan102/quadtree"
# data_type="uniform"
data_type="varden"
# data_type="syn"
# data_type="real"
# data_folder="japan2014"
data_folder="10000000_2"
# data_folder="100000000_2"
# data_folder="1000000000_2"
# data_folder="heavy_duplicate"
pkd_dir="/data/bhuan102/SpaceTreeLib/build"

mkdir -p ../exp_res/$data_type/$data_folder
data_folder=$data_type/$data_folder

make clean
make

export PARLAY_NUM_THREADS

batch_percent=10
# mkdir -p ../exp_res/$data_folder/bach_$batch_percent

# for rangeSize in 0 1 2; do
# for rangeSize in 0 1; do
    rangeSize=0
    for file in `ls $data_dir/$data_folder`; do
        outdir=../exp_res/$data_folder/region_${rangeSize}/
        mkdir -p $outdir
        #### range count experiment
        #   pkd-tree
        numactl -i all $pkd_dir/test -p $data_dir/$data_folder/$file -d 2 -t 0 -q 4 -r 3 -i 1 -s $rangeSize > $outdir/${file}_pkd_range_count.time
        mv range_count2.qry range_count.qry
        mv pkd_range_count.txt $outdir/${file}_pkd_range_count.res

        # numactl -i all $pkd_dir/test -p $data_dir/$data_folder/$file -d 2 -t 0 -q 8 -r 3 -i 1 -s $rangeSize > $outdir/${file}_pkd_range_report.time
        # mv range_report2.qry range_report.qry
        # mv pkd_range_report.txt $outdir/${file}_pkd_range_report.res
        #   zd-tree
        numactl -i all ./main $data_dir/$data_folder/$file 4 $batch_percent 0 > $outdir/${file}_zd_range_count.time
        mv zd_range_count.txt $outdir/${file}_zd_range_count.res
        diff $outdir/${file}_pkd_range_count.res $outdir/${file}_zd_range_count.res 
        # mv zd_range_report.txt ../exp_res/$data_folder/bach_${batch_percent}/${file}_zd_range_report.res
        # numactl -i all ./main $data_dir/$data_folder/$file 4 $batch_percent 0 > $outdir/${file}_cpamaug_range.time
        # numactl -i all ./main $data_dir/$data_folder/$file 4 $batch_percent 0 > $outdir/${file}_zd_range_count.time
        # mv *.txt $outdir
        # mv br_range_count.txt ../exp_res/$data_folder/bach_${batch_percent}/${file}_br_range_count_morton.res
        # mv br_range_report.txt ../exp_res/$data_folder/bach_${batch_percent}/${file}_br_range_report.res
        # diff $outdir/${file}_pkd_range_count.res $outdir/cpamaug_hilbert_range_count.res
        # diff $outdir/${file}_pkd_range_count.res $outdir/cpamaug_zorder_range_report.res
        # diff $outdir/${file}_pkd_range_count.res $outdir/br_hilbert_range_count.res
        # diff $outdir/${file}_pkd_range_count.res $outdir/br_zorder_range_count.res
        # diff $outdir/${file}_pkd_range_report.res $outdir/br_hilbert_range_report.res
        # diff $outdir/${file}_pkd_range_report.res $outdir/br_zorder_range_report.res

        # diff ../exp_res/$data_folder/bach_${batch_percent}/${file}_pkd_range_count.res ../exp_res/$data_folder/bach_${batch_percent}/${file}_zd_range_count.res
        # diff ../exp_res/$data_folder/bach_${batch_percent}/${file}_pkd_range_count.res ../exp_res/$data_folder/bach_${batch_percent}/${file}_br_range_count_morton.res
        # diff ../exp_res/$data_folder/bach_${batch_percent}/${file}_pkd_range_report.res ../exp_res/$data_folder/bach_${batch_percent}/${file}_zd_range_report.res
        # diff ../exp_res/$data_folder/bach_${batch_percent}/${file}_pkd_range_report.res ../exp_res/$data_folder/bach_${batch_percent}/${file}_br_range_report.res
        # diff $outdir/${file}_pkd_range_count.res $outdir/cpamaug_zorder_range_count.res
        # diff $outdir/${file}_pkd_range_report.res $outdir/cpamaug_zorder_range_report.res

        echo "range count: $file finished"

        #---------------------------------------

        #### range report experiment
        #   pkd-tree
        # numactl -i all $pkd_dir/test -p $data_dir/$data_folder/$file -d 2 -t 0 -q 8 -r 3 -i 1 > ../exp_res/$data_folder/bach_${batch_percent}/${file}_pkd_range_report.time
        # mv pkd_range_report.txt ../exp_res/$data_folder/bach_${batch_percent}/${file}_pkd_range_report.res
        #   zd-tree
        # numactl -i all ./main $data_dir/$data_folder/$file 8 $batch_percent > ../exp_res/$data_folder/bach_${batch_percent}/${file}_zd_range_report.time
        # mv zd_range_report.txt ../exp_res/$data_folder/bach_${batch_percent}/${file}_zd_range_report.res
        # diff ../exp_res/$data_folder/bach_${batch_percent}/${file}_pkd_range_report.res ../exp_res/$data_folder/bach_${batch_percent}/${file}_zd_range_report.res
        # echo "range report: $file finished."
        #---------------------------------------
        #### knn report experiment
        # numactl -i all ./main $data_dir/$data_folder/$file 1 $batch_percent > ../exp_res/$data_folder/${file}_zd_knn_report.time
    done
# done

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