#include <bits/stdc++.h>

#include <cpam/cpam.h>
#include <parlay/primitives.h>
#include "geobase.h"
#include <parlay/internal/get_time.h>
#include <parlay/hash_table.h>
#include <cpam/parse_command_line.h>
#include "zdtree.hpp"
#include "morton.hpp"
#include "helper/time_loop.h"

#include "zdtest.hpp"
// #include "seq_zdtree.hpp"

using namespace std;

// #define DEBUG

#ifdef SMALL_TEST
size_t leaf_size = 32;
#else
size_t leaf_size = 32;
#endif
size_t maxSize = 100;
geobase::Bounding_Box largest_mbr;
// Breakdown time evaluation
break_down zd_build_break_down;
break_down cpam_build_break_down;

void line_splitter(){
	cout << "-------------------------------------------------------" << endl;
}


// void test(string input, int q_type, int batch_percent, int use_real = 0){
// 	ifstream fin(input);
// 	parlay::sequence<geobase::Point> P;
// 	largest_mbr = geobase::read_pts(P, fin, use_real);
// 	if (!use_real){
// 		string count_qry_file = "range_count.qry";
// 		// string count_qry_file = "test.qry";
// 		string report_qry_file = "range_report.qry";
// 		// string report_qry_file = "test.qry";
// 		auto range_count_querys = geobase::read_range_query(count_qry_file, 4, maxSize);
// 		auto range_report_querys = geobase::read_range_query(report_qry_file, 8, maxSize);

// 		/* Big Test */
// 		//CPAMBB
// 		// CPAMBB::build_test(P, 0);
// 		// CPAMBB::insert_test(P, batch_percent, 0);
// 		// CPAMBB::delete_test(P, batch_percent, 0);
// 		// CPAMBB::range_count_test(P, range_count_querys);
// 		// CPAMBB::range_count_test(P, range_report_querys);
// 		// CPAMBB::range_report_test(P, range_report_querys);
// 		// line_splitter();
	
// 		// // BRTree, static binary rtree
// 		// BRTest::build_test(P, 0);	//	Test Zorder
// 		// BRTest::range_count_test(P, range_count_querys, 0);	
// 		// BRTest::range_report_test(P, range_report_querys, 0);
// 		// // cout << "-------------------------------------------------------" << endl;
// 		// line_splitter();
// 		// // BRTest::build_test(P, 1);	//	Test Hilbert
// 		// // BRTest::range_count_test(P, range_count_querys, 1);
// 		// // BRTest::range_report_test(P, range_report_querys, 1);
// 		// // cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << endl;

// 		// CPAMZ::build_test(P, 0);
// 		// CPAMZ::insert_test(P, batch_percent, 0);
// 		// CPAMZ::delete_test(P, batch_percent, 0);
// 		// // CPAMZ::range_count_test(P, range_count_querys, 0);
// 		// CPAMZ::range_report_test(P, range_report_querys, 0);
// 		// // CPAMZ::range_report_test(P, range_report_querys, 1);
// 		// line_splitter();
// 		// // cout << "-------------------------------------------------------" << endl;
// 		// // CPAMZ::build_test(P, 1);
// 		// // CPAMZ::insert_test(P, batch_percent, 1);
// 		// // CPAMZ::delete_test(P, batch_percent, 1);
// 		// // CPAMZ::range_report_test(P, range_report_querys, 1);
// 		// // // CPAMZ::range_report_test(P, range_count_querys, 1);
// 		ZDTest::build_test(P);
// 		// ZDTest::batch_insert_test(P, batch_percent);
// 		ZDTest::range_count_test(P, range_count_querys);
// 		// ZDTest::range_report_test(P, range_report_querys);
// 		// // ZDTest::knn_test(P);
// 	}
// 	else{	//	real-test

// 		parlay::sequence<geobase::Point> insert_pts, delete_pts;
// 		auto insert_dir = "real_insert.in";
// 		auto delete_dir = "real_delete.in";
// 		ifstream fin_insert(insert_dir), fin_delete(delete_dir);
// 		read_pts(insert_pts, fin_insert, 1);
// 		read_pts(delete_pts, fin_delete, 1);
// 		cout << "original data: " << P.size() << endl;
// 		cout << "inserted data: " << insert_pts.size() << endl;
// 		cout << "deleted data: " << delete_pts.size() << endl;

// 		// CPAMZ::build_test(P, 0);
// 		// CPAMZ::insert_test(P, batch_percent, 0);
// 		// CPAMZ::delete_test(P, batch_percent, 0);
// 		// CPAMZ::insert_test(P, insert_pts, 0);
// 		// CPAMZ::delete_test(P, delete_pts, 0);
// 		// cout << "-------------------------------------------------------" << endl;
// 		// CPAMZ::build_test(P, 1);
// 		// CPAMZ::insert_test(P, batch_percent, 1);
// 		// CPAMZ::delete_test(P, batch_percent, 1);
// 		// CPAMZ::insert_test(P, insert_pts, 1);
// 		// CPAMZ::delete_test(P, delete_pts, 1);
// 		// // CPAMZ::range_report_test(P, range_count_querys, 1);
// 		// cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << endl;
// 		ZDTest::build_test(P);
// 		// ZDTest::batch_insert_test(P, batch_percent);
// 		ZDTest::batch_insert_test(P, insert_pts, delete_pts);
// 		// cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << endl;
// 	}
// 	// ZDTest::CPAMZ_build_test(P);
// 	// ZDTest::CPAMZ_insert_test(P, batch_percent);
// 	// ZDTest::CPAMZ_delete_test(P, batch_percent);
// 	// ZDTest::zdtree_test(P, range_querys, q_type);
// 	// for (auto i = 0; i < 10; i++){
// 	// 	auto gen_size = 10000 + 10000 * i;
// 	// 	ZDTest::mutation_test(P, gen_size);
// 	// }
// 	// ZDTest::mutation_test(P, 1000000);
// 	// ZDTest::merge_test(P);
// 	// for (auto i = 1; i <= 100000; i *= 10){
// 	// for (auto i = 1; i <= 50; i++){
// 	// 	ZDTest::node_count_test(P, i);
// 	// }
// }

// available tasks:
// build, batch-insert, batch-delete, range-count, range-report, knn
void run(int argc, char** argv){
	cpam::commandLine cmd(argc, argv, "[-i <Path-to-Input>] [-o <Path-to-Output>] [-t <Task-Name>] [-a <Algorithm-Name>] "
									  "[-b <Path-to-Batch-file>] [-bf <batch-fraction>] "
									  "[-r <Path-to-Range-Query>] [-real <Is-Real-Dataset?>] "
									  "[-mv <Dir-to-Multi-Version>] [-s <Single-Point-Query>]"
									  );
	if (!cmd.getOption("-t")){
		cout << "[ERROR]: <Task-Name> is not specified." << endl;
		return;
	}

	if (!cmd.getOption("-i")){
		cout << "[ERROR]: <Path-to-Input> is not specified." << endl;
		return;
	}

	if (!cmd.getOption("-a")){
		cout << "[ERROR]: <Algorithm-Name> is not specified." << endl;
		return;
	}

	string task = cmd.getOptionValue("-t");
	string algo = cmd.getOptionValue("-a");
	string input_file = cmd.getOptionValue("-i");
	int is_real = cmd.getOptionIntValue("-real", 0);

	/* read input file */
	ifstream fin(input_file);
	parlay::sequence<geobase::Point> P;
	largest_mbr = geobase::read_pts(P, fin, is_real);	//	change to true if id is contained.

	if (task == "debug"){
		cout << "total points: " << P.size() << endl;
		get_sorted_points(P);
		count_duplicate_zvalues(P);
		return;
	}

	/* build test */
	if (task == "build"){
		if (algo == "mvzd"){
			ZDTest::build_test(P);
		}
		else if (algo == "cpambb"){
			CPAMBB::build_test(P);
		}
		else if (algo == "cpamz"){
			CPAMZ::build_test(P);
		}
		else if (algo == "combined"){
			ZDTest::build_test(P);
			line_splitter();
			CPAMZ::build_test(P);
			line_splitter();
			CPAMBB::build_test(P);
		}
		return;
	}

	/* batch-insert */
	if (task == "batch-insert"){
		// tested batch-size

		parlay::sequence<size_t> batch_sizes = {
			100000,
			400000,
			700000,
			1000000,
			4000000,
			7000000,
			10000000,
			40000000,
			70000000,
			100000000,
			400000000,
			700000000,
			1000000000
		};
		// parlay::sequence<size_t> batch_sizes = {
		// 	10000,
		// 	40000,
		// 	70000,
		// 	100000,
		// 	400000,
		// 	700000,
		// 	1000000,
		// 	4000000,
		// 	7000000,
		// 	10000000,
		// 	40000000,
		// 	70000000,
		// 	100000000
		// };

		if (algo == "combined"){
			ZDTest::batch_insert_test(P, batch_sizes);
			CPAMZ::batch_insert_test(P, batch_sizes);
			CPAMBB::batch_insert_test(P, batch_sizes);
		}
	}

	/* batch-delete */
	if (task == "batch-delete"){
		// tested batch-size
		parlay::sequence<size_t> batch_sizes = {
			100000,
			400000,
			700000,
			1000000,
			4000000,
			7000000,
			10000000,
			40000000,
			70000000,
			100000000,
			400000000,
			700000000,
			1000000000
		};
		// parlay::sequence<size_t> batch_sizes = {
		// 	10000,
		// 	40000,
		// 	70000,
		// 	100000,
		// 	400000,
		// 	700000,
		// 	1000000,
		// 	4000000,
		// 	7000000,
		// 	10000000,
		// 	40000000,
		// 	70000000,
		// 	100000000
		// };

		if (algo == "combined"){
			ZDTest::batch_delete_test(P, batch_sizes);
			CPAMZ::batch_delete_test(P, batch_sizes);
			CPAMBB::batch_delete_test(P, batch_sizes);
		}
	}

	/* range-count*/
	if (task == "range-count"){
		if (!cmd.getOption("-r")){
			cout << "[Error]: <Path-to-Range-Query> is not specified." << endl;
		}
		else{
			string query_file = cmd.getOptionValue("-r");
			auto s = cmd.getOptionIntValue("-s", 0);
			
			auto [cnt, range_count_querys] = geobase::read_range_query(query_file, 8, maxSize);
			if (s == 1){
			// test single point
				for (size_t i = 0; i < range_count_querys.size(); i++){
					range_count_querys[i] = geobase::Bounding_Box(P[i], P[i]);
				}
			}
			// range_count_querys = range_count_querys.substr(0, 1000);

			// line_splitter();
			// cout << "[zdtree]: " << endl;
			// ZDTest::range_count_test(P, range_count_querys, cnt);
			line_splitter();
			cout << "[cpambb]: " << endl;
			CPAMBB::range_count_test(P, range_count_querys, cnt);
			line_splitter();
			cout << "[cpamz]: " << endl;
			CPAMZ::range_count_test(P, range_count_querys, cnt);
		}
		return;
	}

	/* range-report */
	if (task == "range-report"){
		if (!cmd.getOption("-r")){
			cout << "[Error]: <Path-to-Range-Query> is not specified." << endl;
		}
		else{
			string query_file = cmd.getOptionValue("-r");
			auto s = cmd.getOptionIntValue("-s", 0);

			auto [cnt, range_report_querys] = geobase::read_range_query(query_file, 8, maxSize);
			if (s == 1){
				// test single point
				for (size_t i = 0; i < range_report_querys.size(); i++){
					range_report_querys[i] = geobase::Bounding_Box(P[i], P[i]);
				}
			}
			// range_report_querys = range_report_querys.substr(0, 1000);
			// line_splitter();
			// cout << "[zdtree]: " << endl;
			// ZDTest::range_report_test(P, range_report_querys, cnt);
			line_splitter();
			cout << "[cpambb]: " << endl;
			CPAMBB::range_report_test(P, range_report_querys, cnt);
			// line_splitter();
			// cout << "[cpamz]: " << endl;
			// CPAMZ::range_report_test(P, range_report_querys, cnt);
		}
		return;
	}

	if (task == "knn"){
		if (algo == "mvzd"){
			ZDTest::knn_test(P);
		}
		else if (algo == "cpambb"){
			CPAMBB::knn_test(P);
		}
		else if (algo == "combined"){
			parlay::sequence<size_t> k_vals = {
				1,
				2,
				5,
				10,
				20,
				50,
				70,
				100
			};
			for (size_t k: k_vals){
				cout << "[INFO] k = " << k << endl;
				ZDTest::knn_test(P, k, 50000);
				CPAMBB::knn_test(P, k, 50000);
			}
		}
		else{
			cout << "unsupported" << endl;
		}
	}

	/* multi-version test */
	if (task == "multi-version-test"){
		if (!cmd.getOption("-mv")){
			cout << "[ERROR]: <Dir-to-Multi-Version> is not specified." << endl;
		}
		else{
			if (algo == "mvzd"){
				string mv_dir = cmd.getOptionValue("-mv");
				ZDTest::multi_version_test(P, mv_dir);
			}
			else{
				cout << "unsupported" << endl;
			}
		}
		return;
	}

	/* multi-version with mixed query */
	if (task == "multi-version-query-test"){
		if (!cmd.getOption("-mv")){
			cout << "[ERROR]: <Dir-to-Multi-Version> is not specified." << endl;
		}
		else{
			if (algo == "mvzd"){
				string query_dir = cmd.getOptionValue("-mv");	// for multi-version-query, the -mv path is the query input
				ZDTest::multi_version_query_test(P, query_dir);
				// CPAMBB::multi_version_query_test(P, query_dir);
			}
			else{
				cout << "unsupported" << endl;
			}
		}
	}

	/* diff test */
	if (task == "diff"){
		if (algo == "mvzd"){
			ZDTest::diff_test(P);
		}
		else if (algo == "cpambb"){	//	CPAM-BB
			CPAMBB::diff_test(P);
		}
		else if (algo == "cpamz"){	//	CPAMZ
			CPAMZ::diff_test(P);
		}
		else if (algo == "combined"){
			ZDTest::diff_test(P);
			CPAMBB::diff_test(P);
			CPAMZ::diff_test(P);
		}
	}

	if (task == "spatial-diff"){
		string query_file = cmd.getOptionValue("-r");
		auto [cnt, queries] = geobase::read_range_query(query_file, 8, maxSize);
		// size_t batch_size = 10000;
        size_t ratio = 5;
		// queries = queries.substr(0, 50000);
		// queries = queries.substr(15, 1);
		// queries = queries.substr(109, 1);
		// queries = queries.substr(49999, 1);

		parlay::sequence<size_t> batch_sizes = {
            10000,
            20000,
            50000,
            100000,
            200000,
            500000,
            1000000,
            2000000,
            5000000,
            10000000,
            20000000, 
            50000000,
            100000000
        };
        parlay::sequence<size_t> ratios = {1, 2, 3, 4, 5, 6, 7, 8, 9};
		// queries = queries.substr(660, 1);
        
        cout << "[INFO]: Exp for Batch Size" << endl;
		ZDTest::spatial_diff_test_latency(P, queries, batch_sizes, ratio);
		// CPAMZ::spatial_diff_test_latency(P, queries, batch_sizes, ratio);
		CPAMBB::spatial_diff_test_latency(P, queries, batch_sizes, ratio);
		
		line_splitter();

		ZDTest::plain_spatial_diff_test_latency(P, queries, batch_sizes, ratio);
		ZDTest::plain_spatial_diff_test_latency(P, queries, batch_sizes, ratio, true);
		// CPAMZ::plain_spatial_diff_test_latency(P, queries, batch_sizes, ratio);
		CPAMBB::plain_spatial_diff_test_latency(P, queries, batch_sizes, ratio);

        // for (auto &batch_size: batch_sizes){
        //     if (batch_size > P.size()) break;
        //     cout << "[INFO] Dealing with Batch Size: " << batch_size << endl;
		// 	cout << "[INFO] Proposed Solutions:" << endl;
		// 	ZDTest::spatial_diff_test_latency(P, queries, batch_size, ratio);
        // 	CPAMZ::spatial_diff_test_latency(P, queries, batch_size, ratio);
			// CPAMBB::spatial_diff_test_latency(P, queries, batch_size, ratio);
			
		// 	cout << "[INFO] Plain Solutions:" << endl;

		// 	ZDTest::plain_spatial_diff_test_latency(P, queries, batch_size, ratio);
		// 	ZDTest::plain_spatial_diff_test_latency(P, queries, batch_size, ratio, true);
		// 	CPAMZ::plain_spatial_diff_test_latency(P, queries, batch_size, ratio);
		// 	CPAMBB::plain_spatial_diff_test_latency(P, queries, batch_size, ratio);
        // }

		// line_splitter();

        // batch_size = P.size() / 10 < 10000000 ? P.size() / 10 : 10000000;
        // cout << "[INFO]: Exp for Insert/Delete Ratio" << endl;
        // for (auto &ratio: ratios){
        //     cout << "[INFO] Dealing with Ratio: " << ratio << endl;
		// 	ZDTest::spatial_diff_test_latency(P, queries, batch_size, ratio);
        // 	CPAMZ::spatial_diff_test_latency(P, queries, batch_size, ratio);
		// 	CPAMBB::spatial_diff_test_latency(P, queries, batch_size, ratio);

		// 	cout << "[INFO] Plain Solutions:" << endl;

		// 	ZDTest::plain_spatial_diff_test_latency(P, queries, batch_size, ratio);
		// 	ZDTest::plain_spatial_diff_test_latency(P, queries, batch_size, ratio, true);
		// 	CPAMZ::plain_spatial_diff_test_latency(P, queries, batch_size, ratio);
		// 	CPAMBB::plain_spatial_diff_test_latency(P, queries, batch_size, ratio);
        // }

	}

	if (task == "spatial-diff-old"){
		// cout << "dealing with spatial diff" << endl;
		if (!cmd.getOption("-r")){
			cout << "[Error]: <Path-to-Range-Query> is not specified." << endl;
		}
		else{
			parlay::sequence<size_t> batch_sizes = {
				10000,
				20000,
				50000,
				100000,
				200000,
				500000,
				1000000,
				2000000,
				5000000,
				10000000,
				20000000, 
				50000000,
				100000000
				// 400000000,
				// 700000000,
				// 100000000,
			};

			string query_file = cmd.getOptionValue("-r");
			auto [cnt, range_report_querys] = geobase::read_range_query(query_file, 8, maxSize);
			// bool early_end = false;
			
			// if (query_file[query_file.length() - 5] != '0'){
			// 	early_end = true;
			// }

			// range_report_querys = range_report_querys.substr(0, 1);
			line_splitter();
			cout << "[zdtree fix ratio]: " << endl;
			// ZDTest::spatial_diff_test(P, range_report_querys, batch_sizes, early_end);
			ZDTest::spatial_diff_test_fix_ratio(P, range_report_querys, batch_sizes);
			line_splitter();
			cout << "[cpambb fix ratio]: " << endl;
			// CPAMBB::spatial_diff_test(P, range_report_querys, batch_sizes, early_end);
			CPAMBB::spatial_diff_test_fix_ratio(P, range_report_querys, batch_sizes);
			line_splitter();
			cout << "[cpamz fix ratio]: " << endl;
			CPAMZ::spatial_diff_test_fix_ratio(P, range_report_querys, batch_sizes);
			cout << endl;
			line_splitter();
			cout << endl;


			parlay::sequence<size_t> ratios = {1, 2, 3, 4, 5, 6, 7, 8, 9};
			/* Fix 10% modifications */
			line_splitter();
			cout << "[zdtree fix ratio]: " << endl;
			// ZDTest::spatial_diff_test(P, range_report_querys, batch_sizes, early_end);
			ZDTest::spatial_diff_test_fix_size(P, range_report_querys, ratios);
			line_splitter();
			cout << "[cpambb fix ratio]: " << endl;
			// CPAMBB::spatial_diff_test(P, range_report_querys, batch_sizes, early_end);
			CPAMBB::spatial_diff_test_fix_size(P, range_report_querys, ratios);
			line_splitter();
			cout << "[cpamz fix ratio]: " << endl;
			CPAMZ::spatial_diff_test_fix_size(P, range_report_querys, ratios);


			// cout << "Q: "; print_mbr(range_report_querys[0]);
			// for (size_t i = 0; i < ret1.size(); i++){
			// 	if (ret1[i].id != ret2[i].id){
			// 		// cout << ret1[i].id << " " << ret2[i].id << " " << ret3[i].id << endl;
			// 		cout << ret1[i].id << " " << ret2[i].id << endl;
			// 		cout << ret1[i].x << " " << ret1[i].y << endl;
			// 		// cout << ret2[i].x << " " << ret2[i].y << endl;
			// 		// cout << ret2[i].x << " " << ret2[i].y << endl;
			// 		int x; cin >> x;
			// 	}
			// }
		}
	}


	
	
}

int main(int argc, char **argv) {
	ios::sync_with_stdio(0); cin.tie(0);
	run(argc, argv);
	// geobase::kBoundedQueue<geobase::Point, geobase::nn_pair> kbq;
	// auto test_p = geobase::Point(0, 2032, 66199);
	// cout << bitset<64>(test_p.morton_id) << endl;
	// gc_test();
	// BRTest::hilbert_test();
	// test(argv[1], atoi(argv[2]), atoi(argv[3]), atoi(argv[4]));
	return 0;
}