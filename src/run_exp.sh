#!/bin/bash

# 10M input & query
# data_dir1="/data/bhuan102/quadtree/varden/10000000_2"
data_dir1="/data/bhuan102/mvzd-data-processing/range_count_query/ss_varden_bigint/10000000_2"
# count_dir1="/data/bhuan102/mvzd-data-processing/range_count_query/varden/10000000_2"
count_dir1="/data/bhuan102/mvzd-data-processing/range_count_query/ss_varden_bigint/10000000_2"
# report_dir1="/data/bhuan102/mvzd-data-processing/range_report_query/varden/10000000_2"
report_dir1="/data/bhuan102/mvzd-data-processing/range_report_query/ss_varden_bigint/10000000_2"

# 1B input & query
data_dir2="/data/bhuan102/quadtree/varden/1000000000_2"
count_dir2="/data/bhuan102/mvzd-data-processing/range_count_query/varden/1000000000_2"
report_dir2="/data/bhuan102/mvzd-data-processing/range_report_query/varden/1000000000_2"

make clean
make

echo "10M 1.in range count (1M querys)"
echo "range 0:"
numactl -i all ./main -i $data_dir1/1.in -t range-count -a combined -r $count_dir1/1.in-0.qry
echo "range 1:"
numactl -i all ./main -i $data_dir1/1.in -t range-count -a combined -r $count_dir1/1.in-1.qry
echo "range 2:"
numactl -i all ./main -i $data_dir1/1.in -t range-count -a combined -r $count_dir1/1.in-2.qry

echo "10M 1.in range report (1M querys)"
echo "range 0:"
numactl -i all ./main -i $data_dir1/1.in -t range-report -a combined -r $report_dir1/1.in-0.qry
echo "range 1:"
numactl -i all ./main -i $data_dir1/1.in -t range-report -a combined -r $report_dir1/1.in-1.qry
echo "range 2:"
numactl -i all ./main -i $data_dir1/1.in -t range-report -a combined -r $report_dir1/1.in-2.qry

# echo "10M 1.in build"
# numactl -i all ./main -i $data_dir1/1.in -t build -a combined 

# echo "10M 2.in build"
# numactl -i all ./main -i $data_dir1/2.in -t build -a combined 

# echo "1B 1.in build"
# numactl -i all ./main -i $data_dir2/1.in -t build -a combined 

# echo "1B 2.in build"
# numactl -i all ./main -i $data_dir2/2.in -t build -a combined 

# echo "10M 1.in batch-insert"
# echo "10 percent:"
# numactl -i all ./main -i $data_dir1/1.in -t batch-insert -a combined -bf 10
# echo "30 percent:"
# numactl -i all ./main -i $data_dir1/1.in -t batch-insert -a combined -bf 30
# echo "50 percent:"
# numactl -i all ./main -i $data_dir1/1.in -t batch-insert -a combined -bf 50

# echo "10M 2.in batch-insert"
# echo "10 percent:"
# numactl -i all ./main -i $data_dir1/2.in -t batch-insert -a combined -bf 10
# echo "30 percent:"
# numactl -i all ./main -i $data_dir1/2.in -t batch-insert -a combined -bf 30
# echo "50 percent:"
# numactl -i all ./main -i $data_dir1/2.in -t batch-insert -a combined -bf 50

# echo "1B 1.in batch-insert"
# echo "10 percent:"
# numactl -i all ./main -i $data_dir2/1.in -t batch-insert -a combined -bf 10
# echo "30 percent:"
# numactl -i all ./main -i $data_dir2/1.in -t batch-insert -a combined -bf 30
# echo "50 percent:"
# numactl -i all ./main -i $data_dir2/1.in -t batch-insert -a combined -bf 50

# echo "1B 2.in batch-insert"
# echo "10 percent:"
# numactl -i all ./main -i $data_dir2/2.in -t batch-insert -a combined -bf 10
# echo "30 percent:"
# numactl -i all ./main -i $data_dir2/2.in -t batch-insert -a combined -bf 30
# echo "50 percent:"
# numactl -i all ./main -i $data_dir2/2.in -t batch-insert -a combined -bf 50

# echo "10M 1.in batch-delete"
# echo "10 percent:"
# numactl -i all ./main -i $data_dir1/1.in -t batch-delete -a combined -bf 10
# echo "30 percent:"
# numactl -i all ./main -i $data_dir1/1.in -t batch-delete -a combined -bf 30
# echo "50 percent:"
# numactl -i all ./main -i $data_dir1/1.in -t batch-delete -a combined -bf 50

# echo "10M 2.in batch-delete"
# echo "10 percent:"
# numactl -i all ./main -i $data_dir1/2.in -t batch-delete -a combined -bf 10
# echo "30 percent:"
# numactl -i all ./main -i $data_dir1/2.in -t batch-delete -a combined -bf 30
# echo "50 percent:"
# numactl -i all ./main -i $data_dir1/2.in -t batch-delete -a combined -bf 50

# echo "1B 1.in batch-delete"
# echo "10 percent:"
# numactl -i all ./main -i $data_dir2/1.in -t batch-delete -a combined -bf 10
# echo "30 percent:"
# numactl -i all ./main -i $data_dir2/1.in -t batch-delete -a combined -bf 30
# echo "50 percent:"
# numactl -i all ./main -i $data_dir2/1.in -t batch-delete -a combined -bf 50

# echo "1B 2.in batch-delete"
# echo "10 percent:"
# numactl -i all ./main -i $data_dir2/2.in -t batch-delete -a combined -bf 10
# echo "30 percent:"
# numactl -i all ./main -i $data_dir2/2.in -t batch-delete -a combined -bf 30
# echo "50 percent:"
# numactl -i all ./main -i $data_dir2/2.in -t batch-delete -a combined -bf 50


# echo "10M 1.in range count (1M querys)"
# echo "range 0:"
# numactl -i all ./main -i $data_dir1/1.in -t range-count -a combined -r $count_dir1/1.in-0.qry
# echo "range 1:"
# numactl -i all ./main -i $data_dir1/1.in -t range-count -a combined -r $count_dir1/1.in-1.qry
# echo "range 2:"
# numactl -i all ./main -i $data_dir1/1.in -t range-count -a combined -r $count_dir1/1.in-2.qry

# echo "10M 2.in range count (1M querys)"
# echo "range 0:"
# numactl -i all ./main -i $data_dir1/2.in -t range-count -a combined -r $count_dir1/2.in-0.qry
# echo "range 1:"
# numactl -i all ./main -i $data_dir1/2.in -t range-count -a combined -r $count_dir1/2.in-1.qry
# echo "range 2:"
# numactl -i all ./main -i $data_dir1/2.in -t range-count -a combined -r $count_dir1/2.in-2.qry

# echo "10M 1.in range report (1M querys)"
# echo "range 0:"
# numactl -i all ./main -i $data_dir1/1.in -t range-report -a combined -r $report_dir1/1.in-0.qry
# echo "range 1:"
# numactl -i all ./main -i $data_dir1/1.in -t range-report -a combined -r $report_dir1/1.in-1.qry
# echo "range 2:"
# numactl -i all ./main -i $data_dir1/1.in -t range-report -a combined -r $report_dir1/1.in-2.qry

# echo "10M 2.in range report (1M querys)"
# echo "range 0:"
# numactl -i all ./main -i $data_dir1/2.in -t range-report -a combined -r $report_dir1/2.in-0.qry
# echo "range 1:"
# numactl -i all ./main -i $data_dir1/2.in -t range-report -a combined -r $report_dir1/2.in-1.qry
# echo "range 2:"
# numactl -i all ./main -i $data_dir1/2.in -t range-report -a combined -r $report_dir1/2.in-2.qry

# echo "1B 1.in range count (10000 querys)"
# numactl -i all ./main -i $data_dir2/1.in -t range-count -a combined -r $count_dir2/1.in-0.qry

# echo "1B 2.in range count (10000 querys)"
# numactl -i all ./main -i $data_dir2/2.in -t range-count -a combined -r $count_dir2/2.in-0.qry

# echo "1B 1.in range report (10000 querys)"
# numactl -i all ./main -i $data_dir2/1.in -t range-report -a combined -r $report_dir2/1.in-0.qry

# echo "1B 2.in range report (10000 querys)"
# numactl -i all ./main -i $data_dir2/2.in -t range-report -a combined -r $report_dir2/2.in-0.qry


