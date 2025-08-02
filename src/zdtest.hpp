
#include <bits/stdc++.h>

#include <cpam/cpam.h>
#include <parlay/primitives.h>
#include "geobase.h"
#include <parlay/internal/get_time.h>
#include <parlay/hash_table.h>
#include "zdtree.hpp"
// #include "seq_zdtree.hpp"
#include "morton.hpp"
#include "helper/time_loop.h"

#include "hilbert.h"
#include "binary_rtree.hpp"

#define TEST	//	print for correctness check

using namespace std;
using namespace geobase;

extern size_t leaf_size;
extern Bounding_Box largest_mbr;
extern size_t maxSize;

double zd_leaf_copy_time, zd_inte_copy_time;
double leaf_time, inte_time;
size_t visited_leaf, visited_inte;


namespace CPAMBB{


	template<typename PT>
	void multi_version_test(PT P, string dir, int start_year = 14, int version_num = 5){
		auto cur_year = start_year;
		
		parlay::sequence<geobase::Point> P_delete[version_num], P_insert[version_num], P_update[version_num], P_updove[version_num];

		for (auto i = 0; i != version_num; i++){
			// cout << "[INFO] Year: " << cur_year << "-" << cur_year + 1 << " status:" << endl;
			auto delete_file_name = dir + to_string(cur_year) + "-" + to_string(cur_year + 1) + "-delete.txt";
			auto insert_file_name = dir + to_string(cur_year) + "-" + to_string(cur_year + 1) + "-insert.txt";
			auto update_file_name = dir + to_string(cur_year) + "-" + to_string(cur_year + 1) + "-update.txt";
			auto updove_file_name = dir + to_string(cur_year) + "-" + to_string(cur_year + 1) + "-update_remove.txt";
			// cout << delete_file_name << endl << insert_file_name << endl;
			ifstream fin_delete(delete_file_name);
			ifstream fin_insert(insert_file_name);
			ifstream fin_update(update_file_name);
			ifstream fin_updove(updove_file_name);

			auto delete_mbr = read_pts(P_delete[i], fin_delete, 1);
			auto insert_mbr = read_pts(P_insert[i], fin_insert, 1);
			auto update_mbr = read_pts(P_update[i], fin_update, 1);
			auto updove_mbr = read_pts(P_updove[i], fin_updove, 1);

			if (P_update[i].size() != P_updove[i].size()){
				cout << "[ERROR]: inconsistent # of update pts!" << endl;
			}

			P_delete[i].append(P_updove[i]);
			P_insert[i].append(P_update[i]);

			cur_year += 1;
			delete_mbr = insert_mbr; // useless, just remove warning
			update_mbr = updove_mbr;
		}

		vector<CPAMBB::zmap> all_versions;
		CPAMBB::zmap tree;

		/* Build initial version */
		auto build_avg = time_loop(
            3, 1.0, [&]() {
				tree.clear();
			},
            [&]() {
				tree = CPAMBB::map_init(P);	// initi
            },
   	    	[&](){
			});

		cout << "[cpambb init build time]: " << fixed << setprecision(6) << build_avg << " Seconds" << endl;

		all_versions.emplace_back(tree);

		auto f_noop = [&](const auto &et){	return 0; };

		std::unordered_map<size_t, bool> mmp = {}, num_mmp = {};
		double cur_mem = 0.0, prev_mem = 0.0; 
		size_t cur_inte_num = 0, cur_leaf_num = 0, cur_leaf_sz = 0;
		size_t pre_inte_num = 0, pre_leaf_num = 0, pre_leaf_sz = 0;
		cur_leaf_sz += pre_leaf_sz;
		
		prev_mem = 1.0 * tree.size_in_bytes(f_noop, mmp);
		tie(pre_inte_num, pre_leaf_num, pre_leaf_sz) = tree.node_stats(num_mmp);

		cout << "[init-version memory]: " << prev_mem / 1024.0 / 1024.0 << " MB" << endl;
		cout << "[init-version node nums]: " << pre_inte_num << " interior nodes, " << pre_leaf_num << " leaf nodes" << endl;
		
		vector<zmap> new_ver(version_num);

		for (auto i = 0; i < version_num; i++){
			cout << "dealing with version " << i + 1 << ":" << endl;

			auto commit_avg = time_loop(
				3, 1.0, 
				[&]() {
					new_ver[i].clear();
				},
				[&]() {
					new_ver[i] = CPAMBB::map_commit(all_versions[i], P_insert[i], P_delete[i]);
				},
				[&](){});

			all_versions.emplace_back(new_ver[i]);
			cur_mem = 0;
			mmp.clear();
			num_mmp.clear();
			cur_inte_num = cur_leaf_num = 0;

			for (size_t j = 0; j < all_versions.size(); j++){
				cur_mem += 1.0 * all_versions[j].size_in_bytes(f_noop, mmp); 	// accumulate all version memories, shared pointers only count once
				auto [tmp_inte_num, tmp_leaf_num, tmp_leaf_sz] = all_versions[j].node_stats(num_mmp);
				cur_inte_num += tmp_inte_num, cur_leaf_num += tmp_leaf_num;
			}
			cout << "[new ver commit time]: " << fixed << setprecision(6) << commit_avg << " Seconds" << endl;
			cout << "[cpambb memory usage]: " << (cur_mem - prev_mem) / 1024.0 / 1024.0  << " MB" << endl;
			cout << "[cpambb node nums]: " << cur_inte_num - pre_inte_num << " interior nodes, " << cur_leaf_num - pre_leaf_num << " leaf nodes" << endl;
			prev_mem = cur_mem;
			pre_inte_num = cur_inte_num;
			pre_leaf_num = cur_leaf_num;
		}
	}

	template<typename PT, typename RQ>
	void plain_spatial_diff_test_latency(PT &P, RQ &range_queries, parlay::sequence<size_t> &batch_sizes, size_t &insert_ratio){
		/*  build tree */
		auto cpambb0 = CPAMBB::map_init(P);	//	initial version
		auto max_batch_size = batch_sizes[batch_sizes.size() - 1];
		/* get insert, delete points */
		auto P_test = geobase::shuffle_point(P, max_batch_size);
		auto [P_insert_set, P_delete_set] = geobase::split_insert_delete(P_test, insert_ratio, P.size());

		for (auto &batch_size: batch_sizes){
			cout << "[INFO] Batch Size: " << batch_size << endl;
			auto insert_num = batch_size / 10 * insert_ratio;
			auto delete_num = batch_size / 10 * (10 - insert_ratio);

			auto P_insert = P_insert_set.substr(0, insert_num);
			auto P_delete = P_delete_set.substr(0, delete_num);

			auto P_newver = geobase::collect_newver_point(P, P_insert, P_delete);
			auto cpambb1 = CPAMBB::map_init(P_newver);

			parlay::sequence<size_t> addCnt(range_queries.size());
			parlay::sequence<size_t> removeCnt(range_queries.size());

			auto l_pts = parlay::sequence<Point>::uninitialized(2 * maxSize);
			auto r_pts = parlay::sequence<Point>::uninitialized(2 * maxSize);

			for (size_t i = 0; i < range_queries.size(); i++){
				auto avg_time = time_loop(
					3, 1.0,
					[&](){},
					[&](){
						diff_type ret_diff(maxSize, maxSize);
						CPAMBB::plain_map_spatial_diff(cpambb0, cpambb1, range_queries[i], ret_diff, l_pts, r_pts);
						ret_diff.compact();
						addCnt[i] = ret_diff.add.size();
						removeCnt[i] = ret_diff.remove.size();
					},
					[&]{}
				);
				cout << fixed << setprecision(6) << i << " " << avg_time << endl;
			}

			// auto avg_time = time_loop(
			// 	3, 1.0,
			// 	[&](){},
			// 	[&](){
			// 		for (size_t i = 0; i < range_queries.size(); i++){
			// 			diff_type ret_diff(maxSize, maxSize);
			// 			auto l_pts = parlay::sequence<Point>::uninitialized(2 * maxSize);
			// 			auto r_pts = parlay::sequence<Point>::uninitialized(2 * maxSize);
			// 			CPAMBB::plain_map_spatial_diff(cpambb0, cpambb1, range_queries[i], ret_diff, l_pts, r_pts);
			// 			ret_diff.compact();
			// 			addCnt[i] = ret_diff.add.size();
			// 			removeCnt[i] = ret_diff.remove.size();
			// 		}
			// 	},
			// 	[&]{}
			// );
			// cout << fixed << setprecision(6) << "[cpambb-plain] spatial-diff time (avg): " << avg_time << endl;
			#ifdef TEST
				string file_name = "output/cpambb_spatial_diff_plain-" + to_string(batch_size); 
				ofstream spatialDiffOut(file_name);
				for (size_t i = 0; i < range_queries.size(); i++){
					spatialDiffOut<< addCnt[i] << " " << removeCnt[i] << endl;
				}
			#endif
		}
	}

	template<typename PT, typename RQ>
    void spatial_diff_test_latency(PT &P, RQ &range_queries, parlay::sequence<size_t> &batch_sizes, size_t &insert_ratio){
        /*  build tree */
		auto cpambb0 = CPAMBB::map_init(P);	//	initial version
		auto max_batch_size = batch_sizes[batch_sizes.size() - 1];

        /* get insert, delete points */
        auto P_test = geobase::shuffle_point(P, max_batch_size);
        auto [P_insert_set, P_delete_set] = geobase::split_insert_delete(P_test, insert_ratio, P.size());

		for (auto &batch_size: batch_sizes){
			cout << "[INFO] Batch Size: " << batch_size << endl;
			auto insert_num = batch_size / 10 * insert_ratio;
			auto delete_num = batch_size / 10 * (10 - insert_ratio);

			auto P_insert = P_insert_set.substr(0, insert_num);
			auto P_delete = P_delete_set.substr(0, delete_num);

			auto cpambb1 = CPAMBB::map_delete(P_delete, cpambb0); 
			auto cpambb2 = CPAMBB::map_insert(P_insert, cpambb1);	//	new	version
        
        	parlay::sequence<size_t> addCnt(range_queries.size());
        	parlay::sequence<size_t> removeCnt(range_queries.size());

			auto l_pts = parlay::sequence<Point>::uninitialized(2 * maxSize);
			auto r_pts = parlay::sequence<Point>::uninitialized(2 * maxSize);
			for (size_t i = 0; i < range_queries.size(); i++){
				auto avg_time = time_loop(
					3, 1.0,
					[&](){},
					[&](){
						// for (size_t i = 0; i < range_queries.size(); i++){
							diff_type ret_diff(maxSize, maxSize);
							CPAMBB::plain_map_spatial_diff(cpambb0, cpambb2, range_queries[i], ret_diff, l_pts, r_pts);
							ret_diff.compact();
							addCnt[i] = ret_diff.add.size();
							removeCnt[i] = ret_diff.remove.size();
						// }
					},
					[&]{}
				);
				cout << fixed << setprecision(6) << i << " " << avg_time << endl;
			}
        	// auto avg_time = time_loop(
            // 	3, 1.0,
            // 	[&](){},
            // 	[&](){
            //     	for (size_t i = 0; i < range_queries.size(); i++){
			// 			diff_type ret_diff(maxSize, maxSize);
			// 			auto l_pts = parlay::sequence<Point>::uninitialized(2 * maxSize);
			// 			auto r_pts = parlay::sequence<Point>::uninitialized(2 * maxSize);
			// 			CPAMBB::plain_map_spatial_diff(cpambb0, cpambb2, range_queries[i], ret_diff, l_pts, r_pts);
			// 			ret_diff.compact();
			// 			addCnt[i] = ret_diff.add.size();
			// 			removeCnt[i] = ret_diff.remove.size();
            //     	}
            // 	},
            // 	[&]{}
        	// );
        	// cout << fixed << setprecision(6) << "[cpambb] spatial-diff time (avg): " << avg_time << endl;
			
        
			#ifdef TEST
        		string file_name = "output/cpambb_spatial_diff-" + to_string(batch_size); 
        		ofstream spatialDiffOut(file_name);
        		for (size_t i = 0; i < range_queries.size(); i++){
            		spatialDiffOut<< addCnt[i] << " " << removeCnt[i] << endl;
        		}
			#endif
		}
    }

	/*	50% insertion, 50% deletion	*/
	template<typename PT, typename RQ>
	auto spatial_diff_test_fix_size(PT P,  RQ &range_queries, parlay::sequence<size_t> &batch_sizes, bool use_hilbert = false){
		auto cpambb0 = CPAMBB::map_init(P);	//	initial version

		for (auto &batch_size: batch_sizes){
			if (batch_size > P.size()) break;
			cout << "[ratio (insert/delete)]: " << batch_size << "/" << 10 - batch_size << endl;	

			auto insert_num = P.size() / 100 * batch_size;
			auto delete_num = P.size() / 100 * (10 - batch_size);

			auto P_insert = P.substr(0, insert_num);
			auto P_delete = P.substr(P.size() - delete_num, delete_num);
	
			parlay::parallel_for(0, P_insert.size(), [&](size_t j){
				P_insert[j].id += P.size();
			});

			cout << "# of insertion/deletion: " << P_insert.size() << ", " << P_delete.size();
	
			auto cpambb1 = CPAMBB::map_delete(P_delete, cpambb0); 
			auto cpambb2 = CPAMBB::map_insert(P_insert, cpambb1);	//	new	version
	
			parlay::sequence<size_t> addCnt(range_queries.size());
			parlay::sequence<size_t> removeCnt(range_queries.size());

			auto diff_avg = time_loop(
				3, 1.0, [&]() {},
				[&]() {
					parlay::parallel_for(0, range_queries.size(), [&](int i){
						diff_type ret_diff(maxSize, maxSize);
						CPAMBB::map_spatial_diff(cpambb0, cpambb2, range_queries[i], ret_diff);
						ret_diff.compact();
						addCnt[i] = ret_diff.add.size();
						removeCnt[i] = ret_diff.remove.size();
					});
				},
			[&](){} );
			cout << fixed << setprecision(6) << "[CPAMBB] spatial-diff time (avg): " << diff_avg << endl;
		}
	}
	

	/*	50% insertion, 50% deletion	*/
	template<typename PT, typename RQ>
	auto spatial_diff_test_fix_ratio(PT P,  RQ &range_queries, parlay::sequence<size_t> &batch_sizes, bool use_hilbert = false){
		auto cpambb0 = CPAMBB::map_init(P);	//	initial version

		for (auto &batch_size: batch_sizes){
			if (batch_size > P.size()) break;
			cout << "[batch-size]: " << batch_size << endl;
	
			auto P_insert = P.substr(0, batch_size / 2);
			auto P_delete = P.substr(P.size() - batch_size / 2, batch_size / 2);
	
			parlay::parallel_for(0, P_insert.size(), [&](size_t j){
				P_insert[j].id += P.size();
			});

			cout << "# of insertion/deletion: " << P_insert.size() << ", " << P_delete.size();
	
			auto cpambb1 = CPAMBB::map_delete(P_delete, cpambb0); 
			auto cpambb2 = CPAMBB::map_insert(P_insert, cpambb1);	//	new	version
	
			parlay::sequence<size_t> addCnt(range_queries.size());
			parlay::sequence<size_t> removeCnt(range_queries.size());

			auto diff_avg = time_loop(
				3, 1.0, [&]() {},
				[&]() {
					parlay::parallel_for(0, range_queries.size(), [&](int i){
						diff_type ret_diff(maxSize, maxSize);
						CPAMBB::map_spatial_diff(cpambb0, cpambb2, range_queries[i], ret_diff);
						ret_diff.compact();
						addCnt[i] = ret_diff.add.size();
						removeCnt[i] = ret_diff.remove.size();
					});
				},
			[&](){} );
			cout << fixed << setprecision(6) << "[CPAMBB] spatial-diff time (avg): " << diff_avg << endl;
		}
	}

	template<typename PT, typename RQ>
	auto spatial_diff_test(PT P, RQ &range_queries, parlay::sequence<size_t> &batch_sizes, bool &early_end, bool use_hilbert = false){
		auto cpambb0 = CPAMBB::map_init(P);	//	initial version

		for (auto &batch_size: batch_sizes){
			if (batch_size > P.size()) batch_size = P.size();
			cout << "[batch-size]: " << batch_size << endl;

			auto P_insert = P.substr(0, batch_size);
			auto P_delete = P.substr(P.size() - batch_size, batch_size);

			parlay::parallel_for(0, P_insert.size(), [&](size_t j){
				P_insert[j].id += P.size();
			});

			auto cpambb1 = CPAMBB::map_delete(P_delete, cpambb0); 
			auto cpambb2 = CPAMBB::map_insert(P_insert, cpambb0);	//	new	version

			parlay::sequence<size_t> addCnt(range_queries.size());
			parlay::sequence<size_t> removeCnt(range_queries.size());
			// auto addCnt = 0;
			// auto removeCnt = 0;
			// parlay::sequence<Point> ret;

			// for (size_t i = 0; i < range_queries.size(); i++){
			// 	auto diff_avg = time_loop(
		    // 	3, 1.0, [&]() {},
		    // 	[&]() {
			// 		auto [add, remove] = CPAMBB::map_spatial_diff(cpambb1, cpambb2, range_queries[i]);
			// 		addCnt[i] = add.size();
			// 		removeCnt[i] = remove.size();
			// 			// ret = remove;
			// 	},
			// 	[&](){} );
			// 	cout << fixed << setprecision(6) << diff_avg << endl;
			// }

	    	// auto diff_avg = time_loop(
		    // 	3, 1.0, [&]() {},
		    // 	[&]() {
			// 		parlay::parallel_for(0, range_queries.size(), [&](int i){
			// 		// for (size_t i = 0; i < range_queries.size(); i++){
			// 			// auto [add, remove] = CPAMBB::map_diff(cpambb0, cpambb1);
			// 			auto [add, remove] = CPAMBB::map_spatial_diff(cpambb1, cpambb2, range_queries[i]);
			// 			addCnt[i] = add.size();
			// 			removeCnt[i] = remove.size();
			// 			// ret = remove;
			// 		// }
			// 		});
		    // 	},
	    	// [&](){} );
			// // cout << addCnt << ", " << removeCnt << endl;
			// cout << fixed << setprecision(6) << "[CPAMBB] spatial-diff time (avg): " << diff_avg << endl;

			if (!early_end){
				decltype(cpambb0) commit_ver;
	    		auto commit_avg = time_loop(
		    		3, 1.0, [&]() {
						commit_ver.clear();
					},
		    		[&]() {
						commit_ver = CPAMBB::map_commit(cpambb0, P_insert, P_delete);
		    		},
	    			[&](){} 
				);

				parlay::sequence<Point> conflict_insert, conflict_update, conflict_delete;
				decltype(cpambb0) merge_ver;
				auto merge_avg = time_loop(
		    		3, 1.0, [&]() {
						merge_ver.clear();
						conflict_insert.clear();
						conflict_update.clear();
						conflict_delete.clear();
					},
		    		[&]() {
						tie(merge_ver, conflict_insert, conflict_update, conflict_delete) = CPAMBB::map_merge(cpambb0, cpambb1, cpambb2);
		    		},
	    			[&](){} 
				);
				cout << "[INFO] commit, merge size: " << commit_ver.size() << ", " << merge_ver.size() << endl; 
				cout << fixed << setprecision(6) << "[CPAMBB]: spatial-commit time (avg): " << commit_avg << endl;
				cout << fixed << setprecision(6) << "[CPAMBB]: spatial-merge time (avg): " << merge_avg << endl;
			}

			// string file_name = "output/cpambb_spatial_diff-" + to_string(batch_size); 
			// ofstream spatialDiffOut(file_name);
			// for (size_t i = 0; i < range_queries.size(); i++){
			// 	spatialDiffOut<< addCnt[i] << " " << removeCnt[i] << endl;
			// }
		}

		// return ret;
	}


	template<typename PT>
	void multi_version_query_test(PT P, string query_dir, int batch_percent = 10, int version_num = 6){
		// build zdtree initial version
		auto CPAMZ = CPAMBB::map_init(P);
		cout << "build finished" << endl;

		auto num_insert_version = version_num / 2;
		auto num_delete_version = version_num - num_insert_version;
		auto batch_size = P.size() * batch_percent / 100;
		

		parlay::sequence<decltype(CPAMZ)> all_versions(7);
		all_versions[0] = CPAMZ;
		// insert 3 versions
		for (auto i = 0; i != num_insert_version; i++){
			auto P_insert = P.substr(i * batch_size, batch_size);

			parlay::parallel_for(0, P_insert.size(), [&](size_t j){
				P_insert[j].id += (i + 1) * P.size();
			});

			all_versions[i + 1] = CPAMBB::map_insert(P_insert, all_versions[i]); 
		}
		cout << "insert finished" << endl;	
		//	delete 3 versions
		for (auto i = 0; i != num_delete_version; i++){
			auto P_delete = P.substr(i * batch_size, batch_size);

			parlay::parallel_for(0, P_delete.size(), [&](size_t j){
				P_delete[j].id += (i + 1) * P.size();
			});
			
			all_versions[i + 4] = CPAMBB::map_delete(P_delete, all_versions[i + 3]); 
			// all_versions.push_back(new_version);
		}
		cout << "delete finished" << endl;	
		
		// range count query test
		for (auto i = 0; i != 1; i++){	//	0, 1, 2 represent small, median, and large regions
			auto cur_query_dir = query_dir + "/1.in-" + to_string(i) + ".qry";
			auto [cnt, range_count_querys] = read_range_query(cur_query_dir, 8, maxSize);

			for (size_t j = 0; j < all_versions.size(); j++){

				parlay::sequence<size_t> rangeCnt(range_count_querys.size());
				auto rangeCnt_avg = time_loop(
					3, 1.0, [&]() {},
					[&]() {
						parlay::parallel_for(
							0, range_count_querys.size(),
							[&]( size_t k ) {
								rangeCnt[k] = CPAMBB::range_count(all_versions[j], range_count_querys[k]); 
						});
					},
				[&](){} );

				cout << fixed << setprecision(6) << "CPAM-BB range count time (avg) for region " << i << " on version " << j << ": " << rangeCnt_avg << endl;
				auto output_name = "CPAM-BB_range_count-" + to_string(i) + "-on-" + to_string(j) + ".txt";
				ofstream regionCntOut(output_name);
				for (size_t k = 0; k < range_count_querys.size(); k++){
					regionCntOut << rangeCnt[k] << endl;
				}
			}
		}

		// range report query test
		for (auto i = 0; i != 1; i++){	//	0, 1, 2 represent small, median, and large regions
			auto cur_query_dir = query_dir + "/1.in-" + to_string(i) + ".qry";
			auto [cnt, range_report_querys] = read_range_query(cur_query_dir, 8, maxSize);

			for (size_t j = 0; j < all_versions.size(); j++){

				parlay::sequence<size_t> rangeCnt(range_report_querys.size());
				auto rangeReport_avg = time_loop(
					3, 1.0, [&]() {},
					[&]() {
						parlay::parallel_for(
							0, range_report_querys.size(),
							[&]( size_t k ) {
		    					rangeCnt[k] = CPAMBB::range_report(all_versions[j], range_report_querys[k]).size();
						});
					},
				[&](){} );

				cout << fixed << setprecision(6) << "CPAM-BB range report time (avg) for region " << i << " on version " << j << ": " << rangeReport_avg << endl;
				auto output_name = "CPAM-BB_range_report-" + to_string(i) + "-on-" + to_string(j) + ".txt";
				ofstream regionCntOut(output_name);
				for (size_t k = 0; k < range_report_querys.size(); k++){
					regionCntOut << rangeCnt[k] << endl;
				}
			}
		}
	}

	template<typename PT>
	void diff_test(PT P, int batch_percent = 10, bool use_hilbert = false){
		auto cpambb0 = CPAMBB::map_init(P);	//	initial version

		auto batch_size = P.size() * batch_percent / 100;	//	insertion 10%
		auto P_insert = P.substr(0, batch_size);
		parlay::parallel_for(0, P_insert.size(), [&](size_t j){
			P_insert[j].id += P.size();
		});

		auto P_delete = P.substr(0, 2 * batch_size);

		auto cpambb1 = CPAMBB::map_insert(P_insert, cpambb0);	//	new	version
		cpambb1 = CPAMBB::map_delete(P_delete, cpambb1); 

		auto add_sz = 0, remove_sz = 0;
	    auto cpam_diff_avg = time_loop(
		    3, 1.0, [&]() {},
		    [&]() {
				auto [add, remove] = CPAMBB::map_diff(cpambb0, cpambb1);
				add_sz = add.size();
				remove_sz = remove.size();
		    },
	    [&](){} );
		

		cout << "add size: " << add_sz << endl;
		cout << "remove size: " << remove_sz << endl;
		if (use_hilbert) cout << "[Hilbert-CPAMBB]: ";
		else cout << "[Zorder-CPAMBB]: ";
		cout << fixed << setprecision(6) << "CPAMBB diff time (avg): " << cpam_diff_avg << endl;
	}


	template<typename PT>
	void build_test(PT P, bool use_hilbert = false){
		CPAMBB::zmap tree;

		auto cpam_build_avg = time_loop(
			3, 1.0, [&](){
				tree.clear();
			},
			[&](){
				tree = CPAMBB::map_init(P);
			},
		[&](){} );
		// auto [mem_inte_nodes, mem_leaf_nodes] = CPAMBB::size_in_bytes();
		auto [num_inte_nodes, num_leaf_nodes, leaf_size] = tree.node_stats();
		// cout << "leaf sz = " << 1.0 * leaf_size / 1024.0 / 1024.0 << " MB" << endl;

		auto f_noop = [&](const auto &et){
			return 0;
		};

		cout << "[cpambb memory usage]: " << endl <<
			"[# of inte nodes]: " << num_inte_nodes << endl << 
			"[# of leaf nodes]: " << num_leaf_nodes << endl <<
			"[tree size]: " << 1.0 * tree.size_in_bytes(f_noop) / 1024.0 / 1024.0 << " MB" << endl;
			// "[memory usage for inte nodes]: " << 1.0 * mem_inte_nodes / 1024.0 / 1024.0 << " MB" << endl <<
			// "[memory usage for leaf nodes]: " << 1.0 * mem_leaf_nodes / 1024.0 / 1024.0 << " MB"  << endl;

		// cout << "cpambb print stats: " << endl;
		// CPAMBB::print_stats();

		if (use_hilbert) cout << "[Hilbert-CPAMBB]: ";
		else cout << "[Zorder-CPAMBB]: ";
		cout << fixed << setprecision(6) << "build time (avg): " << cpam_build_avg << endl;
	}

	template<class PT, class RQ>
    void range_count_test(PT P, RQ querys, parlay::sequence<size_t> &cnt, bool use_hilbert = false){
	    auto tree = CPAMBB::map_init(P, use_hilbert);

		parlay::sequence<size_t> rangeCnt(querys.size());

		// auto rangeReport_avg = time_loop(
		// 	1, 0.0, [&]() {},
		// 	[&]() {
		// 		parlay::parallel_for(
		// 			0, querys.size(),
		// 			[&]( size_t i ) {
		//     			rangeCnt[i] = CPAMBB::range_count(CPAMZ, querys[i], use_hilbert);
		// 		});
		// 	},
		// [&](){} );
		// if (use_hilbert) cout << "[Hilbert-CPAMBB]: ";
		// else cout << "[Zorder-CPAMBB]: ";
		// cout << fixed << setprecision(6) << "range count time (avg): " << rangeReport_avg << endl;
		// bool ok = true;
		// parlay::parallel_for(0, querys.size(), [&](size_t i){
		// 	if (cnt[i] != rangeCnt[i]){
		// 		ok = false;
		// 	}
		// });
		// if (!ok){
		// 	cout << "[ERROR] incorrect range count result !!!";
		// }

		/* Latency Test */
		for (size_t i = 0; i < querys.size(); i++){
			auto avg_time = time_loop(
				3, 1.0, 
				[&]() {},
				[&]() {					
					rangeCnt[i] = CPAMBB::range_count(tree, querys[i], use_hilbert);
				},
				[&](){} );
			if (rangeCnt[i] != cnt[i]){
				cout << "[ERROR] Incorrect" << endl;
			}
			else{
				cout << fixed << setprecision(6) << rangeCnt[i] << " " << avg_time << endl;
			}
		}
	
		// auto output_file = use_hilbert ? "cpamaug_hilbert_range_count.res" : "cpamaug_zorder_range_count.res";
		// ofstream regionCntOut(output_file);
		// for (size_t i = 0; i < querys.size(); i++){
		// 	// cout << rangeCnt[i] << endl;
		// 	regionCntOut << rangeCnt[i] << endl;
		// 	if (cnt[i] != rangeCnt[i]){
		// 		cout << "[ERROR] incorrect range count result " << rangeCnt[i] << "-" << cnt[i] << endl;
		// 	}
		// }
    }


	template<class PT, class RQ>
    void range_report_test(PT P, RQ querys, parlay::sequence<size_t> &cnt, bool use_hilbert = false, size_t par_for_granularity = 100){
	    auto tree = CPAMBB::map_init(P, use_hilbert);

		// parlay::sequence<Bounding_Box> q2(querys.size());
		parlay::sequence<size_t> rangeCnt(querys.size());
		parlay::sequence<parlay::sequence<Point> > rangeReport(querys.size());
		for (size_t i = 0; i < querys.size(); i++){
			rangeReport[i].resize(cnt[i]);
		}

		// size_t tot_inte = 0, tot_leaf = 0; double tot_time = 0;
		// for (size_t i = 0; i < querys.size(); i++){
		auto avg_time = time_loop(
			3, 1.0, 
			[&]() {},
			[&]() {					
				// for (size_t i = 0; i < querys.size(); i++){
				parlay::parallel_for(0, querys.size(), [&](size_t i){
					rangeCnt[i] = CPAMBB::range_report(tree, querys[i], rangeReport[i], use_hilbert);
					// rangeCnt[i] = CPAMBB::range_report(tree, querys[(i + 233) % querys.size()], rangeReport[(i + 233) % querys.size()], use_hilbert);
					// rangeCnt[i] = CPAMBB::range_report(tree, querys[(i + 666) % querys.size()], rangeReport[(i + 666) % querys.size()], use_hilbert);
					// rangeCnt[i] = CPAMBB::range_report(tree, querys[(i + 123) % querys.size()], rangeReport[(i + 123) % querys.size()], use_hilbert);
					// rangeCnt[i] = CPAMBB::range_report(tree, querys[(i + 321) % querys.size()], rangeReport[(i + 321) % querys.size()], use_hilbert);
					// rangeCnt[i] = CPAMBB::range_report(tree, querys[(i + 777) % querys.size()], rangeReport[(i + 777) % querys.size()], use_hilbert);
					// rangeCnt[i] = CPAMBB::range_report(tree, querys[(i + 555) % querys.size()], rangeReport[(i + 555) % querys.size()], use_hilbert);
					// rangeCnt[i] = CPAMBB::range_report(tree, querys[(i + 12) % querys.size()], rangeReport[(i + 12) % querys.size()], use_hilbert);
					// rangeCnt[i] = CPAMBB::range_report(tree, querys[(i + 6) % querys.size()], rangeReport[(i + 6) % querys.size()], use_hilbert);
					// rangeCnt[i] = CPAMBB::range_report(tree, querys[(i + 88) % querys.size()], rangeReport[(i + 88) % querys.size()], use_hilbert);
					// rangeCnt[i] = CPAMBB::range_report(tree, querys[(i + 100) % querys.size()], rangeReport[(i + 100) % querys.size()], use_hilbert);
					// rangeCnt[i] = CPAMBB::range_report(tree, querys[i], rangeReport[i], use_hilbert);
					// tot_inte += visited_inte;
					// tot_leaf += visited_leaf;
				});
				// }
			},
			[&](){} 
		);

		if (use_hilbert) cout << "[Hilbert-CPAMBB]: ";
		else cout << "[CPAMBB]: ";
		cout << fixed << setprecision(6) << "range report time (avg): " << avg_time << endl;

			// parlay::parallel_for(0, querys.size(), [&](size_t i){
			// 	q2[i] = Bounding_Box(rangeReport[i][0], rangeReport[i][0]);
			// });
			// tot_inte += visited_inte, tot_leaf += visited_leaf;
			// print_mbr(querys[i]);
			// cout << "[INFO] visited inte: " << visited_inte << endl;
			// cout << "[INFO] visited leaf: " << visited_leaf << endl;
			// cout << "[INFO] cnt = " << rangeCnt[i] << endl;
			// cout << "[INFO] query " << i << " time: " << fixed << setprecision(6) << avg_time << endl;
			// cout << "[INFO] totoal visited inte: " << tot_inte << endl;
			// cout << "[INFO] totoal visited leaf: " << tot_leaf << endl;
			// cout << "[INFO] tot query time: " << fixed << setprecision(6) << avg_time << endl;
			// tot_time += avg_time;

		/* range report sample test. */
		// for (size_t i = 0; i < querys.size(); i++){
		// 	auto avg_time = time_loop(
		// 		3, 1.0, 
		// 		[&]() {
		// 		},
		// 		[&]() {					
		// 			rangeCnt[i] = CPAMBB::range_report(tree, querys[i], rangeReport[i], use_hilbert);
		// 		},
		// 		[&](){} 
		// 	);

		// 	if (rangeCnt[i] != cnt[i]){
		// 		cout << "[ERROR] Incorrect" << rangeCnt[i] << " " << cnt[i] << endl;
		// 	}
		// 	else{
		// 		cout << fixed << setprecision(6) << rangeCnt[i] << " " << avg_time << endl;
		// 	}
		// }

		// }
		// cout << "[INFO] totoal visited inte: " << tot_inte << endl;
		// cout << "[INFO] totoal visited leaf: " << tot_leaf << endl;

		// size_t tot_inte, tot_leaf;

		// auto rangeReport_avg = time_loop(
		// 	3, 1.0, [&]() {
		// 		tot_inte = 0;
		// 		tot_leaf = 0;
		// 	},
		// 	[&]() {
		// 		parlay::parallel_for(
		// 			0, querys.size(),
		// 			[&]( size_t i ) {
		// 				visited_inte = 0;
		// 				visited_leaf = 0;
		//     			// rangeCnt[i] = CPAMBB::range_report(CPAMZ, querys[i], use_hilbert).size();
		// 				cout << "[start] " << i << endl;
		//     			rangeCnt[i] = CPAMBB::range_report(CPAMZ, querys[i], rangeReport[i], use_hilbert);
		// 				cout << i << ": " << rangeCnt[i] << endl;
		// 				cout << "[visted inte]: " << visited_inte << endl;
		// 				cout << "[visted leaf]: " << visited_leaf << endl;
		// 				tot_inte += visited_inte;
		// 				tot_leaf += visited_leaf;
		// 				print_mbr(querys[i]);
		// 				for (size_t j = 0; j < rangeCnt[i]; j++){
		// 					auto pt = rangeReport[i][j];
		// 					cout << pt.x << " " << pt.y << endl;
		// 				}
		// 		});
		// 	},
		// [&](){
		// } );

		// cout << "[tot visted inte]: " << tot_inte << endl;
		// cout << "[tot visted leaf]: " << tot_leaf << endl;

		// if (use_hilbert) cout << "[Hilbert-CPAMBB]: ";
		// else cout << "[Zorder-CPAMBB]: ";
		// cout << fixed << setprecision(6) << "range report time (avg): " << rangeReport_avg << endl;

		// /* Correctness Check */
		// bool ok = true;
		// parlay::parallel_for(0, querys.size(), [&](size_t i){
		// 	if (cnt[i] != rangeCnt[i]){
		// 		ok = false;
		// 	}
		// });
		// if (!ok){
		// 	cout << "[ERROR] incorrect range count result !!!" << endl;
		// }

		// auto output_file = use_hilbert ? "cpamaug_hilbert_range_report.res" : "cpamaug_zorder_range_report.res";
		// ofstream regionCntOut(output_file);
		// for (size_t i = 0; i < querys.size(); i++){
		// 	// cout << rangeCnt[i] << endl;
		// 	regionCntOut << rangeCnt[i] << endl;
		// }
    }

	template<class PT>
    void knn_test(PT P, size_t k = 10, size_t q_num = 50000, bool use_hilbert = false){
	    auto tree = CPAMBB::map_init(P, use_hilbert);
		
		auto knn_sqrdis = parlay::sequence<size_t>::uninitialized(q_num);
		
		auto avg_time = time_loop(
			3, 1.0, 
			[&]() {},
			[&]() {					
				for (size_t i = 0; i < q_num; i++){
					knn_sqrdis[i] = CPAMBB::knn(tree, P[i], k).top().second;
				}
			},
			[&](){} );
		cout << fixed << setprecision(6) << "[CPAMBB] KNN Latency: " << avg_time << endl;


		// auto rangeReport_avg = time_loop(
		// 	3, 1.0, [&]() {},
		// 	[&]() {
		// 		parlay::parallel_for(
		// 			0, q_num,
		// 			[&]( size_t i ) {
		//     			knn_sqrdis[i] = CPAMBB::knn(tree, P[i], k).top().second;
		// 		});`
		// 	},
		// [&](){} );

		// if (use_hilbert) cout << "[Hilbert-CPAMBB]: ";
		// else cout << "[Zorder-CPAMBB]: ";
		// cout << fixed << setprecision(6) << "knn report time (avg): " << rangeReport_avg << endl;

		// /* Correctness Check */

		// auto output_file = "cpambb-knn.res";
		// ofstream RES(output_file);
		// for (size_t i = 0; i < q_num; i++){
		// 	RES << knn_sqrdis[i] << endl;
		// }
    }

	

	template<typename PT>
	void batch_insert_test(PT P, parlay::sequence<size_t> &batch_sizes, bool use_hilbert = false){
		auto n = P.size();
		auto m1 = CPAMBB::map_init(P, use_hilbert);	//	build original tree
		decltype(m1) m2;

		auto rand_p = shuffle_point(P);

		for (auto &num_processed: batch_sizes){
			if (num_processed > P.size()) num_processed = P.size();
			auto P2 = rand_p.substr(0, num_processed);

	    	parlay::parallel_for (0, P2.size(), [&](int i){
		    	P2[i].id = n + i;
	    	});
			bool print_flag = true;
	    	auto cpam_insert_avg = time_loop(
		    	3, 1.0, [&]() {
					m2.clear();
				},
		    	[&]() {
					m2 = CPAMBB::map_insert(P2, m1, use_hilbert);
		    	},
	    	[&](){
				if (print_flag){
					cout << "# of points: " << m2.size() << endl;
					print_flag = false;
				}			
			});

			cout << "[batch_size]: " << num_processed << endl;
			if (use_hilbert) cout << "[Hilbert-CPAMBB]: ";
			else cout << "[Zorder-CPAMBB]: ";
			cout << fixed << setprecision(6) << "batch insert time (avg): " << cpam_insert_avg << endl;
		}

	}

	template<typename PT>
	void batch_delete_test(PT P, parlay::sequence<size_t> &batch_sizes, bool use_hilbert = false){
		auto m1 = CPAMBB::map_init(P, use_hilbert);	//	build original tree
		decltype(m1) m2;
		// auto rand_p = shuffle_point(P);

		for (auto &num_processed: batch_sizes){
			if (num_processed > P.size()) num_processed = P.size();
			auto P2 = P.substr(0, num_processed);
			bool print_flag = true;

	    	auto cpam_insert_avg = time_loop(
		    	3, 1.0, [&]() {
					m2.clear();
				},
		    	[&]() {
					m2 = CPAMBB::map_delete(P2, m1, use_hilbert);
		    	},
	    	[&](){
				if (print_flag){
					cout << "# of points: " << m2.size() << endl;
					print_flag = false;
				}
			} );

			cout << "[batch_size]: " << num_processed << endl;
			if (use_hilbert) cout << "[Hilbert-CPAMBB]: ";
			else cout << "[Zorder-CPAMBB]: ";
			cout << fixed << setprecision(6) << "batch delete time (avg): " << cpam_insert_avg << endl;
		}
	}


}

//	use CPAM (without bounding box) to build
namespace CPAMZ{

	template<typename PT>
	void multi_version_test(PT P, string dir, int start_year = 14, int version_num = 5){
		auto cur_year = start_year;
		
		parlay::sequence<geobase::Point> P_delete[version_num], P_insert[version_num], P_update[version_num], P_updove[version_num];

		for (auto i = 0; i != version_num; i++){
			// cout << "[INFO] Year: " << cur_year << "-" << cur_year + 1 << " status:" << endl;
			auto delete_file_name = dir + to_string(cur_year) + "-" + to_string(cur_year + 1) + "-delete.txt";
			auto insert_file_name = dir + to_string(cur_year) + "-" + to_string(cur_year + 1) + "-insert.txt";
			auto update_file_name = dir + to_string(cur_year) + "-" + to_string(cur_year + 1) + "-update.txt";
			auto updove_file_name = dir + to_string(cur_year) + "-" + to_string(cur_year + 1) + "-update_remove.txt";
			// cout << delete_file_name << endl << insert_file_name << endl;
			ifstream fin_delete(delete_file_name);
			ifstream fin_insert(insert_file_name);
			ifstream fin_update(update_file_name);
			ifstream fin_updove(updove_file_name);

			auto delete_mbr = read_pts(P_delete[i], fin_delete, 1);
			auto insert_mbr = read_pts(P_insert[i], fin_insert, 1);
			auto update_mbr = read_pts(P_update[i], fin_update, 1);
			auto updove_mbr = read_pts(P_updove[i], fin_updove, 1);

			if (P_update[i].size() != P_updove[i].size()){
				cout << "[ERROR]: inconsistent # of update pts!" << endl;
			}

			P_delete[i].append(P_updove[i]);
			P_insert[i].append(P_update[i]);

			cur_year += 1;
			delete_mbr = insert_mbr; // useless, just remove warning
			update_mbr = updove_mbr;
		}

		vector<Morton::zmap> all_versions;
		Morton::zmap tree;

		/* Build initial version */
		auto build_avg = time_loop(
            3, 1.0, [&]() {
				tree.clear();
			},
            [&]() {
				tree = Morton::CPAMZ_init(P);	// initi
            },
   	    	[&](){
			});

		cout << "[cpamz init build time]: " << fixed << setprecision(6) << build_avg << " Seconds" << endl;

		all_versions.emplace_back(tree);

		auto f_noop = [&](const auto &et){	return 0; };

		std::unordered_map<size_t, bool> mmp = {};
		double cur_mem = 0.0, prev_mem = 0.0; 
		prev_mem = 1.0 * tree.size_in_bytes(f_noop, mmp);

		cout << "[init-version memory]: " << prev_mem / 1024.0 / 1024.0 << " MB" << endl;
		
		vector<Morton::zmap> new_ver(version_num);

		for (auto i = 0; i < version_num; i++){
			cout << "dealing with version " << i + 1 << ":" << endl;

			auto commit_avg = time_loop(
				3, 1.0, 
				[&]() {
					new_ver[i].clear();
				},
				[&]() {
					new_ver[i] = Morton::CPAMZ_commit(all_versions[i], P_insert[i], P_delete[i]);
				},
				[&](){});

			all_versions.emplace_back(new_ver[i]);
			cur_mem = 0;
			mmp.clear();
			for (size_t j = 0; j < all_versions.size(); j++){
				cur_mem += 1.0 * all_versions[j].size_in_bytes(f_noop, mmp); 	// accumulate all version memories, shared pointers only count once
			}
			cout << "[new ver commit time]: " << fixed << setprecision(6) << commit_avg << " Seconds" << endl;
			cout << "[cpamz memory usage]: " << (cur_mem - prev_mem) / 1024.0 / 1024.0  << " MB" << endl;
			prev_mem = cur_mem;
		}
	}


	template<typename PT, typename RQ>
    void plain_spatial_diff_test_latency(PT &P, RQ &range_queries, parlay::sequence<size_t> &batch_sizes, size_t &insert_ratio){
        /*  build tree */
		auto cpamz0 = Morton::CPAMZ_init(P);	//	initial version
		auto max_batch_size = batch_sizes[batch_sizes.size() - 1];
		      
        /* get insert, delete points */
        auto P_test = geobase::shuffle_point(P, max_batch_size);
        auto [P_insert_set, P_delete_set] = geobase::split_insert_delete(P_test, insert_ratio, P.size());
		
		for (auto &batch_size: batch_sizes){
			cout << "[INFO] Batch Size: " << batch_size << endl;
			auto insert_num = batch_size / 10 * insert_ratio;
			auto delete_num = batch_size / 10 * (10 - insert_ratio);

			auto P_insert = P_insert_set.substr(0, insert_num);
			auto P_delete = P_delete_set.substr(0, delete_num);

			auto P_newver = geobase::collect_newver_point(P, P_insert, P_delete);

			auto cpamz1 = Morton::CPAMZ_init(P_newver);
        
			parlay::sequence<size_t> addCnt(range_queries.size());
			parlay::sequence<size_t> removeCnt(range_queries.size());
	
			auto avg_time = time_loop(
				3, 1.0,
				[&](){},
				[&](){
					for (size_t i = 0; i < range_queries.size(); i++){
						auto [add, remove] = Morton::plain_map_spatial_diff(cpamz0, cpamz1, range_queries[i]);
						addCnt[i] = add.size();
						removeCnt[i] = remove.size();
					}
				},
				[&]{}
			);
			cout << fixed << setprecision(6) << "[cpamz-plain] spatial-diff time (avg): " << avg_time << endl;
			
			#ifdef TEST
				string file_name = "output/cpamz_spatial_diff_plain-" + to_string(batch_size); 
				ofstream spatialDiffOut(file_name);
				for (size_t i = 0; i < range_queries.size(); i++){
					spatialDiffOut<< addCnt[i] << " " << removeCnt[i] << endl;
				}
			#endif
		}
    }

	template<typename PT, typename RQ>
    void spatial_diff_test_latency(PT &P, RQ &range_queries, parlay::sequence<size_t> &batch_sizes, size_t &insert_ratio){
        /*  build tree */
		auto cpamz0 = Morton::CPAMZ_init(P);	//	initial version
		auto max_batch_size = batch_sizes[batch_sizes.size() - 1];

        /* get insert, delete points */
        auto P_test = geobase::shuffle_point(P, max_batch_size);
        auto [P_insert_set, P_delete_set] = geobase::split_insert_delete(P_test, insert_ratio, P.size());

		for (auto &batch_size: batch_sizes){
			cout << "[INFO] Batch Size: " << batch_size << endl;
			
			auto insert_num = batch_size / 10 * insert_ratio;
			auto delete_num = batch_size / 10 * (10 - insert_ratio);

			auto P_insert = P_insert_set.substr(0, insert_num);
			auto P_delete = P_delete_set.substr(0, delete_num);

			auto cpamz1 = Morton::CPAMZ_delete(P_delete, cpamz0);	//	new	version
			auto cpamz2 = Morton::CPAMZ_insert(P_insert, cpamz1); 
			
			parlay::sequence<size_t> addCnt(range_queries.size());
			parlay::sequence<size_t> removeCnt(range_queries.size());

			auto avg_time = time_loop(
				3, 1.0,
				[&](){},
				[&](){
					for (size_t i = 0; i < range_queries.size(); i++){
						auto [add, remove] = Morton::map_spatial_diff(cpamz0, cpamz2, range_queries[i]);
						addCnt[i] = add.size();
						removeCnt[i] = remove.size();
					}
				},
				[&]{}
			);
			cout << fixed << setprecision(6) << "[cpamz] spatial-diff time (avg): " << avg_time << endl;
			
			#ifdef TEST
				string file_name = "output/cpamz_spatial_diff-" + to_string(batch_size); 
				ofstream spatialDiffOut(file_name);
				for (size_t i = 0; i < range_queries.size(); i++){
					spatialDiffOut<< addCnt[i] << " " << removeCnt[i] << endl;
				}
			#endif
		}
    }

	/*	different ratio: 1-9, 2-8, 3-7, ..., 8-2, 9-1 */
	template<typename PT, typename RQ>
	auto spatial_diff_test_fix_size(PT P,  RQ &range_queries, parlay::sequence<size_t> &batch_sizes, bool use_hilbert = false){
		auto cpamz0 = Morton::CPAMZ_init(P);	//	initial version
	
		for (auto &batch_size: batch_sizes){
			if (batch_size > P.size()) break;
			cout << "[ratio (insert/delete)]: " << batch_size << "/" << 10 - batch_size << endl;	

			auto insert_num = P.size() / 100 * batch_size;
			auto delete_num = P.size() / 100 * (10 - batch_size);

			auto P_insert = P.substr(0, insert_num);
			auto P_delete = P.substr(P.size() - delete_num, delete_num);
			parlay::parallel_for(0, P_insert.size(), [&](size_t j){
				P_insert[j].id += P.size();
			});

			cout << "# of insertion/deletion: " << P_insert.size() << ", " << P_delete.size();

			auto cpamz1 = Morton::CPAMZ_delete(P_delete, cpamz0);	//	new	version
			auto cpamz2 = Morton::CPAMZ_insert(P_insert, cpamz1); 

			parlay::sequence<size_t> addCnt(range_queries.size());
			parlay::sequence<size_t> removeCnt(range_queries.size());

				auto diff_avg = time_loop(
				3, 1.0, [&]() {},
				[&]() {
					parlay::parallel_for(0, range_queries.size(), [&](int i){
						auto [add, remove] = Morton::map_spatial_diff(cpamz0, cpamz2, range_queries[i]);
						addCnt[i] = add.size();
						removeCnt[i] = remove.size();
					});
				},
			[&](){} );
			cout << fixed << setprecision(6) << "[CPAMZ]: spatial-diff time (avg): " << diff_avg << endl;
		}
	}

	/*	50% insertion, 50% deletion	*/
	template<typename PT, typename RQ>
	auto spatial_diff_test_fix_ratio(PT P,  RQ &range_queries, parlay::sequence<size_t> &batch_sizes, bool use_hilbert = false){
		auto cpamz0 = Morton::CPAMZ_init(P);	//	initial version
	
		for (auto &batch_size: batch_sizes){
			if (batch_size > P.size()) break;
			cout << "[batch-size]: " << batch_size << endl;
			
			auto P_insert = P.substr(0, batch_size / 2);
			auto P_delete = P.substr(P.size() - batch_size / 2, batch_size / 2);
			parlay::parallel_for(0, P_insert.size(), [&](size_t j){
				P_insert[j].id += P.size();
			});

			cout << "# of insertion/deletion: " << P_insert.size() << ", " << P_delete.size();

			auto cpamz1 = Morton::CPAMZ_delete(P_delete, cpamz0);	//	new	version
			auto cpamz2 = Morton::CPAMZ_insert(P_insert, cpamz1); 

			parlay::sequence<size_t> addCnt(range_queries.size());
			parlay::sequence<size_t> removeCnt(range_queries.size());

		 	auto diff_avg = time_loop(
		    	3, 1.0, [&]() {},
		    	[&]() {
					parlay::parallel_for(0, range_queries.size(), [&](int i){
						auto [add, remove] = Morton::map_spatial_diff(cpamz0, cpamz2, range_queries[i]);
						addCnt[i] = add.size();
						removeCnt[i] = remove.size();
					});
		    	},
	    	[&](){} );
			cout << fixed << setprecision(6) << "[CPAMZ]: spatial-diff time (avg): " << diff_avg << endl;
		}
	}
	


	template<typename PT, typename RQ>
	auto spatial_diff_test(PT P, RQ &range_queries, parlay::sequence<size_t> &batch_sizes, bool &early_end, bool use_hilbert = false){
		auto cpamz0 = Morton::CPAMZ_init(P);	//	initial version
	
		// parlay::sequence<Point> ret;

		for (auto &batch_size: batch_sizes){
			if (batch_size > P.size()) batch_size = P.size();
			cout << "[batch-size]: " << batch_size << endl;
			
			auto P_insert = P.substr(0, batch_size);
			auto P_delete = P.substr(P.size() - batch_size, batch_size);
			parlay::parallel_for(0, P_insert.size(), [&](size_t j){
				P_insert[j].id += P.size();
			});

			auto cpamz1 = Morton::CPAMZ_delete(P_delete, cpamz0);	//	new	version
			auto cpamz2 = Morton::CPAMZ_insert(P_insert, cpamz0); 

			parlay::sequence<size_t> addCnt(range_queries.size());
			parlay::sequence<size_t> removeCnt(range_queries.size());

			// for (size_t i = 0; i < range_queries.size(); i++){
			// 	auto diff_avg = time_loop(
		    // 	3, 1.0, [&]() {},
		    // 	[&]() {
			// 		auto [add, remove] = Morton::map_spatial_diff(cpamz1, cpamz2, range_queries[i]);
			// 		addCnt[i] = add.size();
			// 		removeCnt[i] = remove.size();
			// 			// ret = remove;
			// 	},
			// 	[&](){} );
			// 	cout << fixed << setprecision(6) << diff_avg << endl;
			// }

	    	// auto diff_avg = time_loop(
		    // 	3, 1.0, [&]() {},
		    // 	[&]() {
			// 		parlay::parallel_for(0, range_queries.size(), [&](int i){
			// 		// for (size_t i = 0; i < range_queries.size(); i++){
			// 			auto [add, remove] = Morton::map_spatial_diff(cpamz1, cpamz2, range_queries[i]);
			// 			addCnt[i] = add.size();
			// 			removeCnt[i] = remove.size();
			// 			// ret = add;
			// 			// for (size_t j = 0; j < 10; j++) cout << remove[j].id << " "; cout << endl;
			// 		// }
			// 		});
		    // 	},
	    	// [&](){} );
			// // cout << addCnt[0] << ", " << removeCnt[0] << endl;
			// cout << fixed << setprecision(6) << "[CPAMZ]: spatial-diff time (avg): " << diff_avg << endl;

			if (!early_end){
				decltype(cpamz0) commit_ver;
				auto commit_avg = time_loop(
					3, 1.0, [&]() {
						commit_ver.clear();
					},
					[&]() {
						commit_ver = Morton::CPAMZ_commit(cpamz0, P_insert, P_delete);
					},
					[&](){} 
				);

				decltype(cpamz0) merge_ver;
				parlay::sequence<Point> conflict_insert, conflict_update, conflict_delete;
				// cout << cpamz0.size() << ", " << cpamz1.size() << ", " << cpamz2.size() << endl;
				auto merge_avg = time_loop(
					3, 1.0, [&]() {
						merge_ver.clear();
						conflict_insert.clear();
						conflict_update.clear();
						conflict_delete.clear();
					},
					[&]() {
						tie(merge_ver, conflict_insert, conflict_update, conflict_delete)  = Morton::CPAMZ_merge(cpamz0, cpamz1, cpamz2);
					},
					[&](){} 
				);
				cout << "[INFO] commit, merge size: " << commit_ver.size() << ", " << merge_ver.size() << endl; 
				cout << fixed << setprecision(6) << "[CPAMZ]: spatial-commit time (avg): " << commit_avg << endl;
				cout << fixed << setprecision(6) << "[CPAMZ]: spatial-merge time (avg): " << merge_avg << endl;
			}
			
			// if (use_hilbert) cout << "[Hilbert-CPAMZ]: ";
			// else cout << "[Zorder-CPAMZ]: ";

			// string file_name = "output/cpamz_spatial_diff-" + to_string(batch_size); 
			// ofstream spatialDiffOut(file_name);
			// for (size_t i = 0; i < range_queries.size(); i++){
			// 	spatialDiffOut<< addCnt[i] << " " << removeCnt[i] << endl;
			// }
		}
		// return ret;
	}


	template<typename PT>
	void diff_test(PT P, int batch_percent = 10, bool use_hilbert = false){
		auto tree = Morton::CPAMZ_init(P);	//	initial version

		auto batch_size = P.size() * batch_percent / 100;	//	insertion 10%
		auto P_insert = P.substr(0, batch_size);
		parlay::parallel_for(0, P_insert.size(), [&](size_t j){
			P_insert[j].id += P.size();
		});

		auto P_delete = P.substr(0, 2 * batch_size);

		auto new_ver = Morton::CPAMZ_insert(P_insert, tree);	//	new	version
		new_ver = Morton::CPAMZ_delete(P_delete, new_ver); 

		auto add_sz = 0, remove_sz = 0;
	    auto cpam_diff_avg = time_loop(
		    3, 1.0, [&]() {},
		    [&]() {
				auto [add, remove] = Morton::map_diff(tree, new_ver);
				add_sz = add.size();
				remove_sz = remove.size();
		    },
	    [&](){} );
		

		cout << "add size: " << add_sz << endl;
		cout << "remove size: " << remove_sz << endl;

		if (use_hilbert) cout << "[Hilbert-CPAMZ]: ";
		else cout << "[Zorder-CPAMZ]: ";
		cout << fixed << setprecision(6) << "spatial diff time (avg): " << cpam_diff_avg << endl;
	}

	template<typename PT>
	void build_test(PT P, bool use_hilbert = false){
		Morton::zmap tree;
		auto cpam_build_avg = time_loop(
			3, 1.0, [&](){
				tree.clear();
			},
			[&](){
				tree = Morton::CPAMZ_init(P);
			},
		[&](){} );

		// auto [mem_inte_nodes, mem_leaf_nodes] = Morton::size_in_bytes();
		auto [num_inte_nodes, num_leaf_nodes, leaf_size] = tree.node_stats();

		auto f_noop = [&](const auto &et){
			return 0;
		};

		cout << "[cpamz memory usage]: " << endl <<
			"[# of inte nodes]: " << num_inte_nodes << endl << 
			"[# of leaf nodes]: " << num_leaf_nodes << endl <<
			"[tree size]: " << 1.0 * tree.size_in_bytes(f_noop) / 1024.0 / 1024.0 << " MB" << endl;
			// "[memory usage for inte nodes]: " << 1.0 * mem_inte_nodes / 1024.0 / 1024.0 << " MB" << endl <<
			// "[memory usage for leaf nodes]: " << 1.0 * mem_leaf_nodes / 1024.0 / 1024.0 << " MB"  << endl;

		if (use_hilbert) cout << "[Hilbert-CPAMZ]: ";
		else cout << "[Zorder-CPAMZ]: ";
		cout << fixed << setprecision(6) << "build time (avg): " << cpam_build_avg << endl;
	}
	
	template<typename PT>
	void batch_insert_test(PT P, parlay::sequence<size_t> &batch_sizes, bool use_hilbert = false){
		auto n = P.size();
		auto m1 = Morton::CPAMZ_init(P, use_hilbert);	//	build original tree
		decltype(m1) m2;

		auto rand_p = shuffle_point(P);

		for (auto &num_processed: batch_sizes){
			if (num_processed > P.size()) num_processed = P.size();

			auto P2 = rand_p.substr(0, num_processed);
	    	parlay::parallel_for (0, P2.size(), [&](int i){
		    	P2[i].id = n + i;
	    	});

			bool print_flag = true;
	    	auto cpam_insert_avg = time_loop(
		    	3, 1.0, [&]() {
					m2.clear();
				},
		    	[&]() {
					m2 = Morton::CPAMZ_insert(P2, m1, use_hilbert);
		    	},
	    	[&](){
				if (print_flag){
					cout << "# of points: " << m2.size() << endl;
					print_flag = false;
				}
			} );

			// auto [mem_inte_nodes, mem_leaf_nodes] = Morton::size_in_bytes();
			// auto [num_inte_nodes, num_leaf_nodes, leaf_size] = tree.node_stats();

			// cout << "[cpamz memory usage]: " << endl <<
			// 	"[# of inte nodes]: " << num_inte_nodes << endl << 
			// 	"[# of leaf nodes]: " << num_leaf_nodes << endl <<
			// 	"[memory usage for inte nodes]: " << 1.0 * mem_inte_nodes / 1024.0 / 1024.0 << " MB" << endl <<
			// 	"[memory usage for leaf nodes]: " << 1.0 * mem_leaf_nodes / 1024.0 / 1024.0 << " MB"  << endl;

			cout << "[batch_size]: " << num_processed << endl;
			if (use_hilbert) cout << "[Hilbert-CPAMZ]: ";
			else cout << "[Zorder-CPAMZ]: ";
			cout << fixed << setprecision(6) << "batch insert time (avg): " << cpam_insert_avg << endl;
		}
	}

	template<typename PT>
	void batch_delete_test(PT P, parlay::sequence<size_t> &batch_sizes, bool use_hilbert = false){
		auto m1 = Morton::CPAMZ_init(P, use_hilbert);	//	build original tree
		decltype(m1) m2;
		// auto rand_p = shuffle_point(P);

		for (auto &num_processed: batch_sizes){
			if (num_processed > P.size()) num_processed = P.size();
			auto P2 = P.substr(0, num_processed);

			bool print_flag = true;
	    	auto cpam_insert_avg = time_loop(
		    	3, 1.0, [&]() {
					m2.clear();
				},
		    	[&]() {
					m2 = Morton::CPAMZ_delete(P2, m1, use_hilbert);
		    	},
	    	[&](){
				if (print_flag){
					cout << "# of points: " << m2.size() << endl;
					print_flag = false;
				}			
			} );
			
			cout << "[batch_size]: " << num_processed << endl;
			if (use_hilbert) cout << "[Hilbert-CPAMZ]: ";
			else cout << "[Zorder-CPAMZ]: ";
			cout << fixed << setprecision(6) << "batch delete time (avg): " << cpam_insert_avg << endl;
		}
	}

	template<class PT, class RQ>
    void range_report_test(PT &P, RQ querys, parlay::sequence<size_t> &cnt, bool use_hilbert = false){
	    auto CPAMZ = Morton::CPAMZ_init(P, use_hilbert);
		// for (auto pt: P){
		// 	cout << pt.morton_id << endl;
		// }

		parlay::sequence<size_t> rangeCnt(querys.size());

		for (size_t i = 0; i < querys.size(); i++){
			auto avg_time = time_loop(
				3, 1.0, 
				[&]() {},
				[&]() {					
					rangeCnt[i] = Morton::range_report(CPAMZ, querys[i], use_hilbert).size();
				},
				[&](){} );
			if (rangeCnt[i] != cnt[i]){
				cout << "[ERROR] Incorrect" << endl;
			}
			else{
				cout << fixed << setprecision(6) << rangeCnt[i] << " " << avg_time << endl;
			}
		}

		// for (size_t i = 0; i < querys.size(); i++){
		// 	auto ret = Morton::range_report(CPAMZ, querys[i], use_hilbert);
		// 	for (auto pt: ret){
		// 		cout << pt.x << " " << pt.y << endl;
		// 	}
		// 	rangeCnt[i] = ret.size();
		// }

		// auto rangeReport_avg = time_loop(
		// 	3, 1.0, [&]() {},
		// 	[&]() {
		// 		parlay::parallel_for(
		// 			0, querys.size(),
		// 			[&]( size_t i ) {
		//     			rangeCnt[i] = Morton::range_report(CPAMZ, querys[i], use_hilbert).size();
		// 		});
		// 	},
		// [&](){} );

		// if (use_hilbert) cout << "[Hilbert-CPAMZ]: ";
		// else cout << "[Zorder-CPAMZ]: ";
		// cout << fixed << setprecision(6) << "range report time (avg): " << rangeReport_avg << endl;

		// bool ok = true;
		// parlay::parallel_for(0, querys.size(), [&](size_t i){
		// 	if (cnt[i] != rangeCnt[i]){
		// 		ok = false;
		// 	}
		// });
		// if (!ok){
		// 	cout << "[ERROR] incorrect range count result !!!";
		// }	

		// auto output_file = use_hilbert ? "cpam_hilbert_range_report.res" : "cpam_zorder_range_report.res";
		// ofstream regionCntOut(output_file);
		// for (size_t i = 0; i < querys.size(); i++){
		// 	regionCntOut << rangeCnt[i] << endl;
		// 	if (cnt[i] != rangeCnt[i]){
		// 		cout << "[ERROR] incorrect range count result " << rangeCnt[i] << "-" << cnt[i] << endl;
		// 	}
		// }
    }

	template<class PT, class RQ>
    void range_count_test(PT &P, RQ querys, parlay::sequence<size_t> &cnt, bool use_hilbert = false){
	    range_report_test(P, querys, cnt);
    }

}

namespace ZDTest{

	template<typename PT>
	void spatial_join_test(PT &P1, PT &P2, FT &point_dis){
		auto P_set = get_sorted_points(P1);
		ZDTree::Tree zdtree(leaf_size);
		zdtree.build(P_set);	// initial version for P1

		P_set = get_sorted_points(P2);
		ZDTree::Tree zdtree2(leaf_size);
		zdtree2.build(P_set);

		parlay::sequence<pair<geobase::Point, geobase::Point> > join_res = {};

		auto avg_time = time_loop(
			3, 1.0,
			[&](){
				join_res.clear();
			},
			[&](){
				zdtree.two_version_spatial_join(zdtree.root, zdtree2.root, point_dis, join_res, largest_mbr, largest_mbr);
			},
			[&]{}
		);
		cout << fixed << setprecision(6) << "[zdtree] spatial-join time (avg): " << avg_time << endl;

		#ifdef TEST
			auto cmp = [](const auto &lhs, const auto &rhs){
				if (lhs.first.id != rhs.first.id) return lhs.first.id < rhs.first.id;
				return lhs.second.id < rhs.second.id;
			};
			join_res = parlay::sort(join_res, cmp);
			parlay::sequence<pair<geobase::Point, geobase::Point> > bf_res = {};
			
			for (auto &pt1: P1){
				for (auto &pt2: P2){
					if (dcmp(geobase::point_point_sqrdis(pt1, pt2) - point_dis * point_dis) <= 0){
						bf_res.emplace_back(pt1, pt2);
					}
				}
			}

			bf_res = parlay::sort(bf_res, cmp);
			if (join_res != bf_res){
				cout << "[ERROR] Incorrect Join Results" << endl;
				cout << "Join Res: " << join_res.size() << endl;
				for (auto &par: join_res){
					cout << fixed << setprecision(6) << sqrt(point_point_sqrdis(par.first, par.second)) << " | " << par.first << ", " << par.second << endl;
				}
				cout << "BF Res: " << bf_res.size() << endl;
				for (auto &par: bf_res){
					cout << fixed << setprecision(6) << sqrt(point_point_sqrdis(par.first, par.second)) << " | " << par.first << ", " << par.second << endl;
				}
			}
			else{
				cout << "[INFO] Correct Join Restuls" << endl;
				cout << "join/bf size: " << join_res.size() << ", " << bf_res.size() << endl;
			}
		#endif
	}

	template<typename PT, typename RQ>
    void plain_spatial_diff_test_latency(PT &P, RQ &range_queries, parlay::sequence<size_t> &batch_sizes, size_t &insert_ratio, bool dual_traverse = false){
        /*  build tree */
		auto P_set = get_sorted_points(P);
		ZDTree::Tree zdtree(leaf_size);
		zdtree.build(P_set);	// initial version

		auto max_batch_size = batch_sizes[batch_sizes.size() - 1];
        /* get insert, delete points */
        auto P_test = geobase::shuffle_point(P, max_batch_size);
        auto [P_insert_set, P_delete_set] = geobase::split_insert_delete(P_test, insert_ratio, P.size());

        /* get insert, delete points */
		for (auto &batch_size: batch_sizes){
			cout << "[INFO] Batch Size: " << batch_size << endl;
			auto insert_num = batch_size / 10 * insert_ratio;
			auto delete_num = batch_size / 10 * (10 - insert_ratio);

			auto P_insert = P_insert_set.substr(0, insert_num);
			auto P_delete = P_delete_set.substr(0, delete_num);

			auto P_newver = geobase::collect_newver_point(P, P_insert, P_delete);
	
			P_set = get_sorted_points(P_newver);
			ZDTree::Tree newtree(leaf_size);
			newtree.build(P_set);
	
			parlay::sequence<size_t> addCnt(range_queries.size());
			parlay::sequence<size_t> removeCnt(range_queries.size());
			// parlay::sequence<size_t> pts1Cnt(range_queries.size());
			// parlay::sequence<size_t> pts2Cnt(range_queries.size());
			// map<size_t, size_t> pts_map, diff_map;

			auto pts1 = parlay::sequence<Point>::uninitialized(2 * maxSize); 
			auto pts2 = parlay::sequence<Point>::uninitialized(2 * maxSize); 
			for (size_t i = 0; i < range_queries.size(); i++){
				auto avg_time = time_loop(
					3, 1.0,
					[&](){},
					[&](){
						if (!dual_traverse){
							size_t cnt1 = 0, cnt2 = 0;
							zdtree.range_report(range_queries[i], largest_mbr, cnt1, pts1);
							pts1.resize(cnt1);
							newtree.range_report(range_queries[i], largest_mbr, cnt2, pts2);
							pts2.resize(cnt2);
							diff_type ret_diff(maxSize, maxSize);
							merge_pts(pts1, pts2, ret_diff);
							ret_diff.compact();
							addCnt[i] = ret_diff.add.size();
							removeCnt[i] = ret_diff.remove.size();
						}
						else{
							diff_type ret_diff(maxSize, maxSize);
							zdtree.spatial_two_version_diff(zdtree.root, newtree.root, range_queries[i], largest_mbr, ret_diff);
							ret_diff.compact();
							addCnt[i] = ret_diff.add.size();
							removeCnt[i] = ret_diff.remove.size();
						}
					},
					[&]{}
				);
				cout << fixed << setprecision(6) << i << " " << avg_time << endl;
			}
			
	
			// auto avg_time = time_loop(
			// 	3, 1.0,
			// 	[&](){
			// 	},
			// 	[&](){
			// 		for (size_t i = 0; i < range_queries.size(); i++){
			// 			auto pts1 = parlay::sequence<Point>::uninitialized(2 * maxSize); 
			// 			auto pts2 = parlay::sequence<Point>::uninitialized(2 * maxSize); 
			// 			if (!dual_traverse){
			// 				size_t cnt1 = 0, cnt2 = 0;
			// 				zdtree.range_report(range_queries[i], largest_mbr, cnt1, pts1);
			// 				pts1.resize(cnt1);
			// 				newtree.range_report(range_queries[i], largest_mbr, cnt2, pts2);
			// 				pts2.resize(cnt2);
			// 				diff_type ret_diff(maxSize, maxSize);
			// 				merge_pts(pts1, pts2, ret_diff);
			// 				ret_diff.compact();
			// 				addCnt[i] = ret_diff.add.size();
			// 				removeCnt[i] = ret_diff.remove.size();
			// 			}
			// 			else{
			// 				diff_type ret_diff(maxSize, maxSize);
			// 				zdtree.spatial_two_version_diff(zdtree.root, newtree.root, range_queries[i], largest_mbr, ret_diff);
			// 				ret_diff.compact();
			// 				addCnt[i] = ret_diff.add.size();
			// 				removeCnt[i] = ret_diff.remove.size();
			// 			}
			// 		}
			// 	},
			// 	[&]{}
			// );
			// if (!dual_traverse){
			// 	cout << fixed << setprecision(6) << "[zdtree-plain] spatial-diff time (avg): " << avg_time << endl;
			// }
			// else{
			// 	cout << fixed << setprecision(6) << "[zdtree-plain-dual] spatial-diff time (avg): " << avg_time << endl;
			// }
			
			#ifdef TEST
				string file_name = "output/zd_spatial_diff_plain-" + to_string(batch_size); 
				if (dual_traverse){
					file_name = file_name + "-dual";
				}
				ofstream spatialDiffOut(file_name);
				for (size_t i = 0; i < range_queries.size(); i++){
					spatialDiffOut<< addCnt[i] << " " << removeCnt[i] << endl;
				}
			#endif
		}
    }


	template<typename PT, typename RQ>
    void spatial_diff_test_latency(PT &P, RQ &range_queries, parlay::sequence<size_t> &batch_sizes, size_t &insert_ratio){
        /*  build tree */
		auto P_set = get_sorted_points(P);
		ZDTree::Tree zdtree(leaf_size);
		zdtree.build(P_set);	// initial version

		auto max_batch_size = batch_sizes[batch_sizes.size() - 1];
        /* get insert, delete points */
        auto P_test = geobase::shuffle_point(P, max_batch_size);
        auto [P_insert_set, P_delete_set] = geobase::split_insert_delete(P_test, insert_ratio, P.size());
		
		for (auto &batch_size: batch_sizes){
			if (batch_size > P.size()) batch_size = P.size();

			cout << "[INFO] Batch Size: " << batch_size << endl;
			auto insert_num = batch_size / 10 * insert_ratio;
			auto delete_num = batch_size / 10 * (10 - insert_ratio);

			auto P_insert = P_insert_set.substr(0, insert_num);
			auto P_delete = P_delete_set.substr(0, delete_num);
			
			// get version 1 by deletion
			auto P_delete_sorted = get_sorted_points(P_delete);
			auto new_ver = zdtree.multi_version_batch_delete_sorted(P_delete_sorted, zdtree.root);	//	new version
			
			// get version 2 by insertion
			auto P_insert_sorted = get_sorted_points(P_insert);
			auto new_ver2 = zdtree.multi_version_batch_insert_sorted(P_insert_sorted, new_ver);
			
			parlay::sequence<size_t> addCnt(range_queries.size());
			parlay::sequence<size_t> removeCnt(range_queries.size());

			diff_type ret_diff(maxSize, maxSize);

			for (size_t i = 0; i < range_queries.size(); i++){
				auto avg_time = time_loop(
					3, 1.0,
					[&](){},
					[&](){
						ret_diff.reset(maxSize, maxSize);
						zdtree.spatial_two_version_diff(zdtree.root, new_ver2, range_queries[i], largest_mbr, ret_diff);
						ret_diff.add.resize(ret_diff.add_cnt);
						ret_diff.remove.resize(ret_diff.remove_cnt);
						addCnt[i] = ret_diff.add.size();
						removeCnt[i] = ret_diff.remove.size();
					},
					[&]{}
				);
				cout << fixed << setprecision(6) << i << " " << avg_time << endl;
			}

			// auto avg_time = time_loop(
			// 	3, 1.0,
			// 	[&](){},
			// 	[&](){
			// 		for (size_t i = 0; i < range_queries.size(); i++){
			// 			// cout << "processing: " << i << endl;
			// 			// print_mbr(range_queries[i]);
			// 			// auto[add, remove] = zdtree.spatial_two_version_diff(zdtree.root, new_ver2, range_queries[i], largest_mbr);
			// 			ret_diff.reset(maxSize, maxSize);
			// 			zdtree.spatial_two_version_diff(zdtree.root, new_ver2, range_queries[i], largest_mbr, ret_diff);
			// 			ret_diff.add.resize(ret_diff.add_cnt);
			// 			ret_diff.remove.resize(ret_diff.remove_cnt);
			// 			// print_Pset_info(ret_diff.add, "add");
			// 			// print_Pset_info(ret_diff.remove, "remove");
			// 			addCnt[i] = ret_diff.add.size();
			// 			removeCnt[i] = ret_diff.remove.size();
			// 		}
			// 	},
			// 	[&]{}
			// );
			// cout << fixed << setprecision(6) << "[zdtree] spatial-diff time (avg): " << avg_time << endl;
			
			#ifdef TEST
				string file_name = "output/zd_spatial_diff-" + to_string(batch_size); 
				ofstream spatialDiffOut(file_name);
				for (size_t i = 0; i < range_queries.size(); i++){
					spatialDiffOut<< addCnt[i] << " " << removeCnt[i] << endl;
				}
			#endif

		}
    }

	/*	50% insertion, 50% deletion	*/
	template<typename PT, typename RQ>
	auto spatial_diff_test_fix_size(PT P,  RQ &range_queries, parlay::sequence<size_t> &batch_sizes, bool use_hilbert = false){
		auto P_set = get_sorted_points(P);
		ZDTree::Tree zdtree(leaf_size);
		zdtree.build(P_set);	// initial version

		parlay::sequence<size_t> addCnt(range_queries.size());
		parlay::sequence<size_t> removeCnt(range_queries.size());

		for (auto &batch_size: batch_sizes){
			if (batch_size > P.size()) break;
			cout << "[ratio (insert/delete)]: " << batch_size << "/" << 10 - batch_size << endl;

			auto insert_num = P.size() / 100 * batch_size;	//	0.5 insert
			auto delete_num = P.size() / 100 * (10 - batch_size);	//	0.5 delete

			auto P_insert = P.substr(0, insert_num);
			parlay::parallel_for(0, P_insert.size(), [&](size_t j){
				P_insert[j].id += P.size();
			});

			auto P_delete = P.substr(P.size() - delete_num, delete_num);

			cout << "# of insertion/deletion: " << P_insert.size() << ", " << P_delete.size();

			// get version 1 by deletion
			auto P_delete_sorted = get_sorted_points(P_delete);
			auto new_ver = zdtree.multi_version_batch_delete_sorted(P_delete_sorted, zdtree.root);	//	new version

			// get version 2 by insertion
			auto P_insert_sorted = get_sorted_points(P_insert);
			auto new_ver2 = zdtree.multi_version_batch_insert_sorted(P_insert_sorted, new_ver);

	    	auto diff_avg = time_loop(
		    	3, 1.0, [&]() {},
		    	[&]() {
					parlay::parallel_for(0, range_queries.size(), [&](int i){
						diff_type ret_diff(maxSize, maxSize);
						zdtree.spatial_two_version_diff(zdtree.root, new_ver2, range_queries[i], largest_mbr, ret_diff);
						addCnt[i] = ret_diff.add.size();
						removeCnt[i] = ret_diff.remove.size();
					});
		    	},
	    		[&](){} 
			);
			cout << fixed << setprecision(6) << "[zdtree]: spatial diff time (avg): " << diff_avg << endl;
		}
	}

	/*	50% insertion, 50% deletion	*/
	template<typename PT, typename RQ>
	auto spatial_diff_test_fix_ratio(PT P,  RQ &range_queries, parlay::sequence<size_t> &batch_sizes, bool use_hilbert = false){
		auto P_set = get_sorted_points(P);
		ZDTree::Tree zdtree(leaf_size);
		zdtree.build(P_set);	// initial version

		parlay::sequence<size_t> addCnt(range_queries.size());
		parlay::sequence<size_t> removeCnt(range_queries.size());

		for (auto &batch_size: batch_sizes){
			if (batch_size > P.size()) break;
			cout << "[batch-size]: " << batch_size << endl;

			auto insert_num = batch_size / 2;	//	0.5 insert
			auto delete_num = batch_size - insert_num;	//	0.5 delete

			auto P_insert = P.substr(0, insert_num);
			parlay::parallel_for(0, P_insert.size(), [&](size_t j){
				P_insert[j].id += P.size();
			});

			auto P_delete = P.substr(P.size() - delete_num, delete_num);

			cout << "# of insertion/deletion: " << P_insert.size() << ", " << P_delete.size();

			// get version 1 by deletion
			auto P_delete_sorted = get_sorted_points(P_delete);
			auto new_ver = zdtree.multi_version_batch_delete_sorted(P_delete_sorted, zdtree.root);	//	new version

			// get version 2 by insertion
			auto P_insert_sorted = get_sorted_points(P_insert);
			auto new_ver2 = zdtree.multi_version_batch_insert_sorted(P_insert_sorted, new_ver);

	    	auto diff_avg = time_loop(
		    	3, 1.0, [&]() {},
		    	[&]() {
					parlay::parallel_for(0, range_queries.size(), [&](int i){
						diff_type ret_diff(maxSize, maxSize);
						zdtree.spatial_two_version_diff(zdtree.root, new_ver2, range_queries[i], largest_mbr, ret_diff);
						addCnt[i] = ret_diff.add.size();
						removeCnt[i] = ret_diff.remove.size();
					});
		    	},
	    		[&](){} 
			);
			cout << fixed << setprecision(6) << "[zdtree]: spatial diff time (avg): " << diff_avg << endl;
		}
	}



	template<typename PT, typename RQ>
	auto spatial_diff_test(PT P,  RQ &range_queries, parlay::sequence<size_t> &batch_sizes, bool &early_end, bool use_hilbert = false){
		auto P_set = get_sorted_points(P);
		ZDTree::Tree zdtree(leaf_size);
		zdtree.build(P_set);	// initial version

		parlay::sequence<size_t> addCnt(range_queries.size());
		parlay::sequence<size_t> removeCnt(range_queries.size());

		for (auto &batch_size: batch_sizes){
			if (batch_size > P.size()) batch_size = P.size();
			cout << "[batch-size]: " << batch_size << endl;

			auto P_insert = P.substr(0, batch_size);
			parlay::parallel_for(0, P_insert.size(), [&](size_t j){
				P_insert[j].id += P.size();
			});

			auto P_delete = P.substr(P.size() - batch_size, batch_size);

			// get version 1 by deletion
			auto P_delete_sorted = get_sorted_points(P_delete);
			auto new_ver = zdtree.multi_version_batch_delete_sorted(P_delete_sorted, zdtree.root);	//	new version

			// get version 2 by insertion
			auto P_insert_sorted = get_sorted_points(P_insert);
			auto new_ver2 = zdtree.multi_version_batch_insert_sorted(P_insert_sorted, zdtree.root);

			// cout << P_insert.size() << ", " << P_delete.size() << endl;
			// cout << zdtree.collect_records(new_ver).size() << ", " << zdtree.collect_records(new_ver2).size() << endl;

			// for (size_t i = 0; i < range_queries.size(); i++){
			// 	auto diff_avg = time_loop(
		    // 	3, 1.0, [&]() {},
		    // 	[&]() {
			// 			auto[add, remove] = zdtree.spatial_two_version_diff(new_ver, new_ver2, range_queries[i], largest_mbr);
			// 			addCnt[i] = add.size();
			// 			removeCnt[i] = remove.size();
			// 			// ret = add;
		    // 	},
	    	// 	[&](){} );
			// 	cout << fixed << setprecision(6) << diff_avg << endl;
			// }

	    	// auto diff_avg = time_loop(
		    // 	3, 1.0, [&]() {},
		    // 	[&]() {
			// 		parlay::parallel_for(0, range_queries.size(), [&](int i){
			// 			auto[add, remove] = zdtree.spatial_two_version_diff(new_ver, new_ver2, range_queries[i], largest_mbr);
			// 			addCnt[i] = add.size();
			// 			removeCnt[i] = remove.size();
			// 			// ret = add;
			// 		});
		    // 	},
	    	// 	[&](){} 
			// );
			// cout << fixed << setprecision(6) << "[zdtree]: spatial diff time (avg): " << diff_avg << endl;

			// cout << "[INFO] add set size, P_insert size: " << add.size() << ", " << P_insert.size() << endl;
			// cout << "[INFO] remove set size, P_delete size: " << remove.size() << ", " << P_delete.size() << endl;
			if (!early_end){
				decltype(zdtree.root) commit_ver;
	    		auto commit_avg = time_loop(
		    		3, 1.0, [&]() {
						commit_ver.reset();
					},
		    		[&]() {
						commit_ver = zdtree.commit(zdtree.root, P_insert, P_delete);
		    		},
	    			[&](){} 
				);

				// cout << "init size = " << zdtree.collect_records(zdtree.root).size() << endl;
				// cout << "commit size = " << zdtree.collect_records(commit_ver).size() << endl;
				// cout << "commit finished." << endl;

				decltype(zdtree.root) merge_ver;
				parlay::sequence<Point> conflict_insert, conflict_update, conflict_delete;

				auto merge_avg = time_loop(
		    		3, 1.0, [&]() {
						merge_ver.reset();
						conflict_insert.clear();
						conflict_update.clear();
						conflict_delete.clear();
					},
		    		[&]() {
						tie(merge_ver, conflict_insert, conflict_update, conflict_delete) = zdtree.merge(zdtree.root, new_ver, new_ver2);
		    		},
	    			[&](){} 
				);

				cout << "[INFO] commit, merge size: " << zdtree.collect_records(commit_ver).size() << ", " << zdtree.collect_records(merge_ver).size() << endl;
				cout << fixed << setprecision(6) << "[zdtree]: spatial commit time (avg): " << commit_avg << endl;
				cout << fixed << setprecision(6) << "[zdtree]: spatial merge time (avg): " << merge_avg << endl;
			}

			// string file_name = "output/zd_spatial_diff-" + to_string(batch_size);
			// ofstream spatialDiffOut(file_name);
			// for (size_t i = 0; i < range_queries.size(); i++){
			// 	spatialDiffOut<< addCnt[i] << " " << removeCnt[i] << endl;
			// }
		}

		// return ret;
	}

	template<typename PT>
	void diff_test(PT P, int batch_percent = 10, bool use_hilbert = false){
		auto P_set = get_sorted_points(P);
		ZDTree::Tree zdtree(leaf_size);
		zdtree.build(P_set);	// initial version

		auto batch_size = P.size() * batch_percent / 100;	//	insertion 10%

		auto P_insert = P.substr(0, batch_size);
		parlay::parallel_for(0, P_insert.size(), [&](size_t j){
			P_insert[j].id += P.size();
		});
		
		auto P_insert_sorted = get_sorted_points(P_insert);
		auto new_ver = zdtree.multi_version_batch_insert_sorted(P_insert_sorted, zdtree.root);

		auto P_delete = P.substr(0, 2 * batch_size);
		// auto P_delete = P.substr(0, 10);
		auto P_delete_sorted = get_sorted_points(P_delete);
		new_ver = zdtree.multi_version_batch_delete_sorted(P_delete_sorted, new_ver);	//	new version

		auto add_sz = 0, remove_sz = 0;
	    auto zd_diff_avg = time_loop(
		    3, 1.0, [&]() {},
		    [&]() {
				auto [add, remove] = zdtree.diff(zdtree.root, new_ver, 64);
				add_sz = add.size();
				remove_sz = remove.size();
		    },
	    [&](){} );

		cout << "add size: " << add_sz << endl;
		cout << "remove size: " << remove_sz << endl;
		cout << fixed << setprecision(6) << "[zdtree] diff time (avg): " << zd_diff_avg << endl;
	}
	
	

	template<typename PT>
	void multi_version_query_test(PT P, string query_dir, int batch_percent = 10, int version_num = 6){
		// build zdtree initial version
		auto P_set = get_sorted_points(P);
		ZDTree::Tree zdtree(leaf_size);
		zdtree.build(P_set);
		cout << "[INFO] Tree build finished." << endl;
		zdtree.multi_version_roots.emplace_back(zdtree.root);

		auto num_insert_version = version_num / 2;
		auto num_delete_version = version_num - num_insert_version;
		auto batch_size = P.size() * batch_percent / 100;
		
		shared_ptr<ZDTree::BaseNode> new_ver = zdtree.root;
		// insert 3 versions
		for (auto i = 0; i != num_insert_version; i++){
			auto P_insert = P.substr(i * batch_size, batch_size);

			parlay::parallel_for(0, P_insert.size(), [&](size_t j){
				P_insert[j].id += (i + 1) * P.size();
			});

			auto P_insert_sorted = get_sorted_points(P_insert);
			new_ver = zdtree.multi_version_batch_insert_sorted(P_insert_sorted, new_ver);
			zdtree.multi_version_roots.emplace_back(new_ver);
		}
		cout << "[INFO] Tree insertion finished." << endl;

		//	delete 3 versions
		for (auto i = 0; i != num_delete_version; i++){
			auto P_delete = P.substr(i * batch_size, batch_size);

			parlay::parallel_for(0, P_delete.size(), [&](size_t j){
				P_delete[j].id += (i + 1) * P.size();
			});

			auto P_delete_sorted = get_sorted_points(P_delete);
			new_ver = zdtree.multi_version_batch_delete_sorted(P_delete_sorted, new_ver);
			zdtree.multi_version_roots.emplace_back(new_ver);
		}

		cout << "[INFO] Tree deletion finished." << endl;

		// for (size_t i = 0; i < zdtree.multi_version_roots.size(); i++){
		// 	auto cur_hash = zdtree.tree_hash(zdtree.multi_version_roots[i]);
		// 	auto cur_hash2 = zdtree.tree_hash_non_associative(zdtree.multi_version_roots[i]);
		// 	cout << "version " << i << " hash is: " << cur_hash << ", " << cur_hash2 << endl;
		// }

		// zdtree.print_leaf(zdtree.multi_version_roots[0]);
		// zdtree.print_leaf(zdtree.multi_version_roots[2]);
		// return;

		// range count query test
		for (auto i = 0; i != 3; i++){	//	0, 1, 2 represent small, median, and large regions
			auto cur_query_dir = query_dir + "/1.in-" + to_string(i) + ".qry";
			auto [cnt, range_count_querys] = read_range_query(cur_query_dir, 8, maxSize);
			// range_count_querys = range_count_querys.substr(0, 1000);
			// cout << "test i" << endl;

			// cout << zdtree.multi_version_roots.size() << endl;
			for (size_t j = 0; j < zdtree.multi_version_roots.size(); j++){
				parlay::sequence<size_t> rangeCnt(range_count_querys.size());
				auto rangeCnt_avg = time_loop(
					3, 1, [&]() {},
					// 1, 0, [&]() {},
					[&]() {
						parlay::parallel_for(
							0, range_count_querys.size(),
							[&]( size_t k ) {
								rangeCnt[k] = zdtree.range_count(zdtree.multi_version_roots[j], range_count_querys[k], largest_mbr);
						});
					},
				[&](){} );

				cout << fixed << setprecision(6) << "zdtree range count time (avg) for region " << i << " on version " << j << ": " << rangeCnt_avg << endl;
				auto output_name = "mvzd_range_count-" + to_string(i) + "-on-" + to_string(j) + ".txt";
				ofstream regionCntOut(output_name);
				for (size_t k = 0; k < range_count_querys.size(); k++){
					regionCntOut << rangeCnt[k] << endl;
				}
			}
		}

		// range report test
		for (auto i = 0; i != 2; i++){
			auto cur_query_dir = query_dir + "/1.in-" + to_string(i) + ".qry";
			auto [cnt, range_report_querys] = read_range_query(cur_query_dir, 8, maxSize);

			for (size_t j = 0; j < zdtree.multi_version_roots.size(); j++){
				parlay::sequence<parlay::sequence<Point> > rangeReport(range_report_querys.size());
				parlay::sequence<size_t> rangeReportCnt(range_report_querys.size(), 0);
				for (size_t k = 0; k < range_report_querys.size(); k++){
					rangeReport[k].resize(maxSize);	
				}
				auto rangeReport_avg = time_loop(
					3, 1.0, [&]() {},
					[&]() {
						parlay::parallel_for(
							0, range_report_querys.size(),
							[&]( size_t k ) {
								rangeReportCnt[k] = 0;
								zdtree.range_report(zdtree.multi_version_roots[j], range_report_querys[k], largest_mbr, rangeReportCnt[k], rangeReport[k]);
						});
					},
				[&](){} );
				cout << fixed << setprecision(6) << "zdtree range report time (avg) for region " << i << " on version " << j << ": " << rangeReport_avg << endl;
				auto output_name = "mvzd_range_report-" + to_string(i) + "-on-" + to_string(j) + ".txt";
				ofstream regionCntOut(output_name);
				for (size_t k = 0; k < range_report_querys.size(); k++){
					regionCntOut << rangeReportCnt[k] << endl;
				}
			}
		}
		
		
	}

	template<typename PT>
	void multi_version_test(PT P, string dir, int start_year = 14, int version_num = 5){
		auto cur_year = start_year;
		
		parlay::sequence<geobase::Point> P_delete[version_num], P_insert[version_num], P_update[version_num], P_updove[version_num];

		for (auto i = 0; i != version_num; i++){
			// cout << "[INFO] Year: " << cur_year << "-" << cur_year + 1 << " status:" << endl;
			auto delete_file_name = dir + to_string(cur_year) + "-" + to_string(cur_year + 1) + "-delete.txt";
			auto insert_file_name = dir + to_string(cur_year) + "-" + to_string(cur_year + 1) + "-insert.txt";
			auto update_file_name = dir + to_string(cur_year) + "-" + to_string(cur_year + 1) + "-update.txt";
			auto updove_file_name = dir + to_string(cur_year) + "-" + to_string(cur_year + 1) + "-update_remove.txt";
			// cout << delete_file_name << endl << insert_file_name << endl;
			ifstream fin_delete(delete_file_name);
			ifstream fin_insert(insert_file_name);
			ifstream fin_update(update_file_name);
			ifstream fin_updove(updove_file_name);

			auto delete_mbr = read_pts(P_delete[i], fin_delete, 1);
			auto insert_mbr = read_pts(P_insert[i], fin_insert, 1);
			auto update_mbr = read_pts(P_update[i], fin_update, 1);
			auto updove_mbr = read_pts(P_updove[i], fin_updove, 1);

			// cout << "# of delete pts: " << P_delete[i].size() << endl;
			// cout << "# of insert pts: " << P_insert[i].size() << endl;

			if (P_update[i].size() != P_updove[i].size()){
				cout << "[ERROR]: inconsistent # of update pts!" << endl;
			}
			// cout << "# of update pts: " << P_update[i].size() << endl;

			P_delete[i].append(P_updove[i]);
			P_insert[i].append(P_update[i]);

			cur_year += 1;
			delete_mbr = insert_mbr; // useless, just remove warning
			update_mbr = updove_mbr;
		}

		ZDTree::Tree zdtree(leaf_size);

		/* Build initial version */
		auto build_avg = time_loop(
            3, 1.0, [&]() {
				zdtree.clear();
			},
            [&]() {
                auto P_set = get_sorted_points(P);
                zdtree.build(P_set);
            },
   	    	[&](){
			});

		zdtree.multi_version_roots.emplace_back(zdtree.root);	//	store the initial version 0
		auto tot = zdtree.num_of_nodes();
		auto stat = zdtree.get_tree_statistics();

		cout << "[zdtree init build time]: " << fixed << setprecision(6) << build_avg << " Seconds" << endl;
		cout << "[zdtree num of total tree nodes]: " << tot << endl;
		cout << "[zdtree node nums]: " << stat.num_inte_nodes << " interior nodes, " << stat.num_leaf_nodes << " leaf nodes" << endl;
		cout << "[zdtree memory usage]: " << stat.get_total() << " MB" << endl;
		// cout << "[memory usage for inte nodes]: " << 1.0 * stat.mem_inte_nodes / 1024.0 / 1024.0 << " MB" << endl 
		//      << "[memory usage for leaf nodes]: " << 1.0 * stat.mem_leaf_nodes / 1024.0 / 1024.0 << " MB"  << endl;
		
		shared_ptr<ZDTree::BaseNode> new_ver = zdtree.root;

		for (auto i = 0; i < version_num; i++){
			cout << "dealing with version " << i + 1 << ":" << endl;

			auto commit_avg = time_loop(
				3, 1.0, 
				[&]() {
					new_ver.reset();
				},
				[&]() {
					new_ver = zdtree.commit(zdtree.multi_version_roots[i], P_insert[i], P_delete[i]);
				},
				[&](){});

			zdtree.multi_version_roots.emplace_back(new_ver);
			cout << "[new ver commit time]: " << fixed << setprecision(6) << commit_avg << " Seconds" << endl;

			tot = zdtree.num_of_nodes();
			auto cur_stat = zdtree.get_tree_statistics();

			cout << "[num of total tree nodes]: " << tot << endl;
			cout << "[zdtree memory usage]: " << cur_stat.get_total() - stat.get_total() << " MB" << endl;
			cout << "[zdtree tot inte memory usage]: " << 1.0 * cur_stat.mem_inte_nodes / 1024.0 / 1024.0 << " MB" << endl;
			cout << "[zdtree tot leaf memory usage]: " << 1.0 * cur_stat.mem_leaf_nodes / 1024.0 / 1024.0 << " MB" << endl;
			cout << "[zdtree increased nodes]: " << cur_stat.num_inte_nodes - stat.num_inte_nodes << " interior nodes, " <<
													cur_stat.num_leaf_nodes - stat.num_leaf_nodes << " leaf nodes" << endl;
			stat = cur_stat;

			// cout << "[memory usage for inte nodes]: " << 1.0 * stat.mem_inte_nodes / 1024.0 / 1024.0 << " MB" << endl 
			//  	<< "[memory usage for leaf nodes]: " << 1.0 * stat.mem_leaf_nodes / 1024.0 / 1024.0 << " MB"  << endl;	
		}
	}

    template<class PT>
    void build_test(PT P){
		// auto node_num = 0;
		// unordered_map<int, int> leaf_map = {};
		ZDTree::Tree zdtree(leaf_size);

        auto zdtree_build_avg = time_loop(
            3, 1.0, [&]() {
				zdtree.clear();
			},
            [&]() {
                auto P_set = get_sorted_points(P);
                zdtree.build(P_set);
				zdtree.multi_version_roots.emplace_back(zdtree.root);
				// node_num = zdtree.num_of_nodes();
				// leaf_map = zdtree.leaf_size_status();
            },
        [&](){
		});
		auto stat = zdtree.get_tree_statistics();
		// auto node_num = zdtree.num_of_nodes();

		// cout << "total # of nodes: " << node_num << endl;
		cout << "[zdtree stats]: " << endl <<
			"[# of inte nodes]: " << stat.num_inte_nodes << endl <<
			"[# of leaf nodes]: " << stat.num_leaf_nodes << endl <<
			"[memory usage for inte nodes]: " << 1.0 * stat.mem_inte_nodes / 1024.0 / 1024.0 << " MB" << endl <<
			"[memory usage for leaf nodes]: " << 1.0 * stat.mem_leaf_nodes / 1024.0 / 1024.0 << " MB"  << endl;
        cout << fixed << setprecision(6) << "[zdtree]: build time (avg): " << zdtree_build_avg << endl;
    }

	// template<typename PT>
	// void CPAMZ_build_test(PT P){
	// 	auto cpam_build_avg = time_loop(
	// 		3, 1.0, [&](){},
	// 		[&](){
	// 			auto P_set = get_unsorted_address(P);
	// 			auto CPAMZ = Morton::CPAMZ_init(P_set);
	// 		},
	// 	[&](){} );
	// 	cout << fixed << setprecision(6) << "CPAM build time (avg): " << cpam_build_avg << endl;
	// }

    // template<typename PT>
    // void mutation_test(PT P, int gen_size = 10000){
    //     // auto n = P.size();
    //     auto P_set = get_sorted_address(P);
    //     ZDTree::Tree zdtree(leaf_size);
    //     zdtree.build(P_set);
    //     srand((int)time(0));

    //     set<int> delete_id {}, update_id = {}, insert_id = {};
    //     // for (auto i = 0; i < gen_size; i++){
    //     //     auto id = rand() % n;
    //     //     delete_id.insert(id); 
    //     // }
    //     // for (auto i = 0; i < gen_size; i++){
    //     //     auto id = rand() % n;
    //     //     if (delete_id.find(id) == delete_id.end()) update_id.insert(id);
    //     // }
    //     parlay::sequence<Point> P_insert = {}, P_delete = {}, P_update = {};
	// 	// biased test
	// 	for (auto &pt: P){
	// 		if (pt.x >= 45400 && pt.x < 45600  && pt.y >= 46200 && pt.y < 46400) {
	// 			update_id.insert(pt.id);
	// 		}
	// 	}
	// 	cout << "selected size: " << update_id.size() << endl;

    //     // P_insert.resize(gen_size);
    //     // for (auto i = 0; i < gen_size; i++){
    //     //     insert_id.insert(n + i);
    //     //     P_insert[i] = Point(n + i, rand() % (int)largest_mbr.second.x, rand() % (int)largest_mbr.second.y);
    //     // }

    //     // P_delete.resize(delete_id.size());
    //     auto cur = 0;
    //     // for (auto &id: delete_id){
    //     //     P_delete[cur++] = P[id];
    //     // }

    //     P_update.resize(update_id.size());
    //     cur = 0;
    //     for (auto &id: update_id){
    //         P_update[cur++] = P[id];
    //     }
        
    //     shared_ptr<ZDTree::BaseNode> new_ver = nullptr;

    //     P_set = get_sorted_address(P_delete);   //  delete P_delete
    //     new_ver = zdtree.multi_version_batch_delete_sorted(P_set, zdtree.root);

    //     P_set = get_sorted_address(P_update);   //  delete P_update first 
    //     new_ver = zdtree.multi_version_batch_delete_sorted(P_set, new_ver);   
    //     parlay::parallel_for(0, P_update.size(), [&](int i){    //  update point coordinates
    //         P_update[i].x += rand() % 10 + 1;
    //         P_update[i].y += rand() % 10 + 1;
    //     });
    //     P_set = get_sorted_address(P_update);
    //     new_ver = zdtree.multi_version_batch_insert_sorted(P_set, new_ver); //  insert modified points

    //     P_set = get_sorted_address(P_insert);
    //     new_ver = zdtree.multi_version_batch_insert_sorted(P_set, new_ver); //  insert P_insert
    //     // auto [add, remove] = zdtree.diff(zdtree.root, new_ver, 64);
	// 	auto time_s = chrono::high_resolution_clock::now();
    //     auto [add, remove] = zdtree.two_version_diff(zdtree.root, new_ver);
	// 	auto time_t = chrono::high_resolution_clock::now();
	// 	auto duration = chrono::duration_cast<chrono::microseconds>(time_t - time_s);
	// 	auto tot_time = duration.count();
	// 	cout << "total time (micro-seconds): " << tot_time << endl;
	// 	cout << fixed << setprecision(6) << "breakdown time (case0-4): " << zdtree.case0_time << "(" << 100.0 * zdtree.case0_time / tot_time << "%), " 
	// 																	 << zdtree.case1_time << "(" << 100.0 * zdtree.case1_time / tot_time << "%), " 
	// 									 								 << zdtree.case2_time << "(" << 100.0 * zdtree.case2_time / tot_time << "%), " 
	// 																	 << zdtree.case3_time << "(" << 100.0 * zdtree.case3_time / tot_time << "%), "
	// 																	 << zdtree.case4_time << "(" << 100.0 * zdtree.case4_time / tot_time << "%)" << endl;

	// 	cout << "gen size: " << gen_size << " [case0-4], "; 
    //     cout << zdtree.case0 << ", " << zdtree.case1 << ", " << zdtree.case2 << ", " << zdtree.case3 << ", " << zdtree.case4 << endl;
	// 	auto [t_insert, t_delete, t_update] = zdtree.filter_diff_results(add, remove);
	// 	set<int> ret_insert_id = {}, ret_delete_id = {}, ret_update_id = {};
	// 	for (auto &pt: t_insert) ret_insert_id.insert(pt->id);
	// 	for (auto &pt: t_delete) ret_delete_id.insert(pt->id);
	// 	for (auto &pt: t_update) ret_update_id.insert(pt->id);
    //     cout << "insert check: " << ret_insert_id.size() << ", " << insert_id.size() << " " << (ret_insert_id == insert_id) << endl;
    //     cout << "delete check: " << ret_delete_id.size() << ", " << delete_id.size() << " " << (ret_delete_id == delete_id) << endl;
    //     cout << "update check: " << ret_update_id.size() << ", " << update_id.size() << " " << (ret_update_id == update_id) << endl;
    // }

	// template<typename PT>
    // void merge_test(PT P, int gen_size = 10000){
    //     int n = P.size();
    //     auto P_set = get_sorted_address(P);
    //     ZDTree::Tree zdtree(leaf_size);
    //     zdtree.build(P_set);

    //     set<int> delete_id1 = {0, 1, 2, 3, 4, 5}, update_id1 = {6, 7, 8, 9, 10}, insert_id1 = {n, n + 1, n + 2, n + 3, n + 4, n + 5};
    //     set<int> delete_id2 = {0, 1, 2, 6, 7, 8}, update_id2 = {3, 4, 5, 9, 11}, insert_id2 = {n, n + 1, n + 2, n + 6, n + 7, n + 8};

    //     parlay::sequence<Point> P_insert1 = {}, P_delete1 = {}, P_update1 = {};
    //     parlay::sequence<Point> P_insert2 = {}, P_delete2 = {}, P_update2 = {};

	// 	for (auto &id: insert_id1) {
	// 		auto new_p = Point(id, P[id - n].x + 1, P[id - n].x + 1);
	// 		P_insert1.emplace_back(new_p);
	// 	}
	// 	for (auto &id: insert_id2) {
	// 		auto new_p = Point(id, P[id - n].x + 2, P[id - n].x + 2);
	// 		P_insert2.emplace_back(new_p);
	// 	}

	// 	for (auto &id: delete_id1) P_delete1.emplace_back(P[id]);
	// 	for (auto &id: delete_id2) P_delete2.emplace_back(P[id]);

	// 	for (auto &id: update_id1) P_update1.emplace_back(P[id]);
	// 	for (auto &id: update_id2) P_update2.emplace_back(P[id]);

    //     shared_ptr<ZDTree::BaseNode> new_ver = nullptr;
    //     shared_ptr<ZDTree::BaseNode> ver1 = nullptr;
    //     shared_ptr<ZDTree::BaseNode> ver2 = nullptr;

    //     P_set = get_sorted_address(P_delete1);   //  delete P_delete
    //     new_ver = zdtree.multi_version_batch_delete_sorted(P_set, zdtree.root);
    //     P_set = get_sorted_address(P_update1);   //  delete P_update first 
    //     new_ver = zdtree.multi_version_batch_delete_sorted(P_set, new_ver);   
    //     parlay::parallel_for(0, P_update1.size(), [&](int i){    //  update point coordinates
    //         P_update1[i].x += 1;
    //         P_update1[i].y += 1;
    //     });
    //     P_set = get_sorted_address(P_update1);
    //     new_ver = zdtree.multi_version_batch_insert_sorted(P_set, new_ver); //  insert modified points
    //     P_set = get_sorted_address(P_insert1);
    //     ver1 = zdtree.multi_version_batch_insert_sorted(P_set, new_ver); //  insert P_insert, finish modifying branch 1

    //     P_set = get_sorted_address(P_delete2);   //  delete P_delete
    //     new_ver = zdtree.multi_version_batch_delete_sorted(P_set, zdtree.root);
    //     P_set = get_sorted_address(P_update2);   //  delete P_update first 
    //     new_ver = zdtree.multi_version_batch_delete_sorted(P_set, new_ver);   
    //     parlay::parallel_for(0, P_update2.size(), [&](int i){    //  update point coordinates
    //         P_update2[i].x += 2;
    //         P_update2[i].y += 2;
    //     });
    //     P_set = get_sorted_address(P_update2);
    //     new_ver = zdtree.multi_version_batch_insert_sorted(P_set, new_ver); //  insert modified points
    //     P_set = get_sorted_address(P_insert2);
    //     ver2 = zdtree.multi_version_batch_insert_sorted(P_set, new_ver); //  insert P_insert, finish modifying branch 1

	// 	for (auto &pt: P_update1) {cout << pt.x << ", " << pt.y << " ";} cout << endl;
	// 	for (auto &pt: P_update2) {cout << pt.x << ", " << pt.y << " ";} cout << endl;

	// 	auto [conflict_insert, conflict_update, conflict_delate] = zdtree.merge(zdtree.root, ver1, ver2);

	// 	cout << "insert conflict: "; for (auto &pt: conflict_insert) {cout << pt->id << " ";} cout << endl;
	// 	cout << "update conflict: "; for (auto &pt: conflict_update) {cout << pt->id << " ";} cout << endl;
	// 	cout << "delate conflict: "; for (auto &pt: conflict_delate) {cout << pt->id << " ";} cout << endl;
		
    // }


	// template<typename PT>
	// void CPAMZ_insert_test(PT P, size_t batch_percent){
	// 	auto n = P.size();
	// 	auto num_processed = batch_percent * n / 100;
	// 	auto P2 = P.substr(0, num_processed);
	//     parlay::parallel_for (0, P2.size(), [&](int i){
	// 	    P2[i].id = n + i;
	//     });
	// 	auto P_set = get_unsorted_address(P);
	// 	auto m1 = Morton::CPAMZ_init(P_set);	//	build original tree
	// 	decltype(m1) m2;
	//     auto cpam_insert_avg = time_loop(
	// 	    3, 1.0, [&]() {},
	// 	    [&]() {
	// 		    P_set = get_sorted_address(P2);
	// 			m2 = Morton::CPAMZ_insert(P_set, m1);
	// 	    },
	//     [&](){} );
	// 	cout << "cpam original size: " << m1.size() << endl;
	// 	cout << "cpam new size: " << m2.size() << endl;
	// 	cout << fixed << setprecision(6) << "CPAM insert time (avg): " << cpam_insert_avg << endl;
	// }

	// template<typename PT>
	// void CPAMZ_delete_test(PT P, size_t batch_percent){
	// 	auto n = P.size();
	// 	auto num_processed = batch_percent * n / 100;
	// 	auto P2 = P.substr(0, num_processed);
	//     parlay::parallel_for (0, P2.size(), [&](int i){
	// 	    P2[i].id = n + i;
	//     });
	//     auto P3 = P;	// 1+x% 
	//     P3.resize(n + num_processed);
	//     parlay::parallel_for (0, n + num_processed, [&](size_t i){
	// 	    if (i < n) P3[i] = P[i];
	// 	    else P3[i] = P2[i - n];
	//     });
	// 	auto P_set = get_unsorted_address(P3);
	// 	auto m1 = Morton::CPAMZ_init(P_set);	//	build original tree
	// 	decltype(m1) m2;
	//     auto cpam_insert_avg = time_loop(
	// 	    3, 1.0, [&]() {},
	// 	    [&]() {
	// 		    P_set = get_sorted_address(P2);
	// 			m2 = Morton::CPAMZ_delete(P_set, m1);
	// 	    },
	//     [&](){} );
	// 	cout << "cpam original size: " << m1.size() << endl;
	// 	cout << "cpam new size: " << m2.size() << endl;
	// 	cout << fixed << setprecision(6) << "CPAM delete time (avg): " << cpam_insert_avg << endl;
	// }

	// template<typename PT>
	// void node_count_test(PT P, size_t batch_percent){
	// 	auto n = P.size();
	//     auto num_processed = batch_percent * n / 100;
	//     // auto num_processed = batch_num;
	// 	cout << "num processed: " << batch_percent << "%, " << num_processed << endl;
	//     // auto num_processed = 1000;
	// 	auto P2 = P.substr(0, num_processed);	//	first x%
	// 	auto P3 = P.substr(num_processed, n - num_processed); 	//	1-x%
		
	// 	ZDTree::Tree zdtree(leaf_size);
	// 	SEQ_ZDTree::Tree seq_zdtree(leaf_size);

	// 	auto P_set = get_sorted_address(P);
	// 	zdtree.build(P_set);	//	original tree;
	// 	seq_zdtree.build(P_set);
	// 	cout << "original hash: " << zdtree.tree_hash(zdtree.root) << ", " << seq_zdtree.tree_hash(seq_zdtree.root) << endl;

	// 	shared_ptr<ZDTree::BaseNode> new_ver = nullptr;
	// 	shared_ptr<SEQ_ZDTree::BaseNode> seq_new_ver = nullptr;

	// 	auto zdtree_time_cost = time_loop(
	// 		3, 1.0, [&](){
	// 			zdtree.clear();
	// 			P_set = get_sorted_address(P3);
	// 			zdtree.build(P_set);
	// 			zdtree.multi_version_roots.clear();
	// 			zdtree.multi_version_roots.emplace_back(zdtree.root);
	// 		}, 
	// 		[&](){
	// 			P_set = get_sorted_address(P2);
	// 			new_ver = zdtree.multi_version_batch_insert_sorted(P_set, zdtree.root);
	// 		},
	// 	[&](){} ); 
	// 	zdtree.multi_version_roots.emplace_back(new_ver);

	// 	auto seq_zdtree_time_cost = time_loop(
	// 		3, 1.0, [&](){
	// 			seq_zdtree.clear();
	// 			P_set = get_sorted_address(P3);
	// 			seq_zdtree.build(P_set);
	// 			seq_new_ver = seq_zdtree.root;
	// 			seq_zdtree.multi_version_roots.clear();
	// 			seq_zdtree.multi_version_roots.emplace_back(seq_zdtree.root);
	// 		}, 
	// 		[&](){
	// 			P_set = get_sorted_address(P2);
	// 			// seq_new_ver = seq_zdtree.multi_version_batch_insert_sorted(P_set, seq_zdtree.root);
	// 			for (size_t i = 0; i < P_set.size(); i++){
	// 				parlay::sequence<Point*> P_cur{P_set[i]};
	// 				seq_new_ver = seq_zdtree.multi_version_batch_insert_sorted(P_cur, seq_new_ver);
	// 				seq_zdtree.multi_version_roots.emplace_back(seq_new_ver);
	// 			}
	// 		},
	// 	[&](){} ); 
	// 	seq_zdtree.multi_version_roots.emplace_back(seq_new_ver);

	// 	cout << "after hash (par v.s. seq): " << zdtree.tree_hash(new_ver) << " " << seq_zdtree.tree_hash(seq_new_ver) << endl;
	// 	cout << fixed << setprecision(6) << "insert time (par v.s. seq): " << zdtree_time_cost << " " << seq_zdtree_time_cost << endl;

	// 	// auto record_check = [&](auto &lhs, auto &rhs){
	// 	// 	if (lhs.size() != rhs.size()) {
	// 	// 		cout << "size does not match!" << endl;
	// 	// 		return false;
	// 	// 	}
	// 	// 	for (size_t i = 0; i < lhs.size(); i++){
	// 	// 		if (lhs[i]->x != rhs[i]->x || lhs[i]->y != rhs[i]->y){
	// 	// 			cout << "order does not match!" << endl;
	// 	// 			return false;
	// 	// 		}
	// 	// 	}
	// 	// 	return true;
	// 	// };
	// 	// auto l_records = zdtree.collect_records(zdtree.root);
	// 	// auto r_records = seq_zdtree.collect_records(seq_zdtree.root);
	// 	// cout << record_check(l_records, r_records) << endl;
	// 	// l_records = zdtree.collect_records(new_ver);
	// 	// r_records = seq_zdtree.collect_records(seq_new_ver);
	// 	// cout << record_check(l_records, r_records) << endl;
	// 	cout << "zdtree tot node count (par v.s. seq): " << zdtree.node_count() << " " << seq_zdtree.node_count() << endl;
	// 	zdtree.multi_version_roots.clear();
	// 	zdtree.multi_version_roots.emplace_back(zdtree.root);
	// 	zdtree.multi_version_roots.emplace_back(new_ver);

	// 	seq_zdtree.multi_version_roots.clear();
	// 	seq_zdtree.multi_version_roots.emplace_back(seq_zdtree.root);
	// 	seq_zdtree.multi_version_roots.emplace_back(seq_new_ver);
	// 	cout << "zdtree tot node count (two-version par v.s. seq): " << zdtree.node_count() << " " << seq_zdtree.node_count() << endl;
		
	// 	cout << endl;
	// 	// P_set = get_sorted_address(P);

	// 	// auto P_set2 = parlay::sort(r_records,  [&](auto lhs, auto rhs){
	// 	// 	auto msd = 0;
	// 	// 	if (geobase::less_msb(static_cast<unsigned int>(lhs->x) ^ static_cast<unsigned int>(rhs->x), static_cast<unsigned int>(lhs->y) ^ static_cast<unsigned int>(rhs->y))) 
	// 	// 		msd = 1;
	// 	// 	return !msd ? lhs->x < rhs->x : lhs->y < rhs->y;
	// 	// });

	// 	// for (auto i = 0; i < n; i++){
	// 	// 	// if (!(*l_records[i] == *P_set[i])){
	// 	// 	// 	cout << "zdtree order not correct!" << endl;
	// 	// 	// 	break;
	// 	// 	// }
	// 	// 	if (!(*r_records[i] == *P_set[i])){
	// 	// 	// if (!(*P_set2[i] == *P_set[i])){
	// 	// 		// cout << P_set[i]->x << ", " << P_set[i]->y << "; " << P_set2[i]->x << ", " << P_set2[i]->y << endl;
	// 	// 		cout << "seq_zdtree order not correct!" << endl;
	// 	// 		break;
	// 	// 	}
	// 	// }
	// }

    template<typename PT>
    void batch_insert_test(PT P, parlay::sequence<size_t> &batch_sizes){
	    auto n = P.size();

	    ZDTree::Tree zdtree(leaf_size);
	    auto P_set = get_sorted_points(P); // build original tree

	    zdtree.build(P_set);

		auto rand_p = shuffle_point(P);

		shared_ptr<ZDTree::BaseNode> new_ver = nullptr;
		
		for (auto &num_processed: batch_sizes){
			if (num_processed > P.size()) num_processed = P.size();
	    	auto P2 = rand_p.substr(0, num_processed);	// first x%
	    	parlay::parallel_for (0, P2.size(), [&](int i){
		    	P2[i].id = n + i;
	    	});

			bool print_flag = true;
	    	auto zdtree_insert_avg = time_loop(
		    	3, 1.0, [&]() {
					new_ver.reset();
				},
		    	[&]() {
					// parlay::internal::timer t("debug", true);
			    	P_set = get_sorted_points(P2);
					// t.next("sort time");
			    	new_ver = zdtree.multi_version_batch_insert_sorted(P_set, zdtree.root);
					// t.next("insert time");
		    	},
	    	[&](){
				if (print_flag){
					cout << "# of points: " << zdtree.collect_records(new_ver).size() << endl;
					print_flag = false;
				}
			} );

			cout << "[batch_size]: " << num_processed << endl;		
	    	cout << fixed << setprecision(6) << "[zdtree]: batch insert time (avg): " << zdtree_insert_avg << endl;
		}
    }

    template<typename PT>
    void batch_delete_test(PT P, parlay::sequence<size_t> &batch_sizes){
	    ZDTree::Tree zdtree(leaf_size);
	    auto P_set = get_sorted_points(P); // build original tree
	    zdtree.build(P_set);
	    shared_ptr<ZDTree::BaseNode> new_ver = nullptr;

		// auto rand_p = shuffle_point(P);
		

		for (auto &num_processed: batch_sizes){
			if (num_processed > P.size()) num_processed = P.size();
	    	auto P2 = P.substr(0, num_processed);	// first x%
			// bool print_flag = true;
	    	auto zdtree_delete_avg = time_loop(
		    	3, 1.0, [&]() {
					new_ver.reset();
					// zdtree.visited_leaf = 0;
					// zdtree.visited_inte = 0;
					// zd_leaf_copy_time = 0.0;
					// zd_inte_copy_time = 0.0;
					// cout << "before size: " << zdtree.collect_records(zdtree.root).size() << endl;
				},
		    	[&]() {
					// parlay::internal::timer t("debug", true);
			    	// P_set = get_sorted_address(P2);
			    	P_set = get_sorted_points(P2);
					// t.next("sort time");
			    	new_ver = zdtree.multi_version_batch_delete_sorted(P_set, zdtree.root);
					// t.next("delete time");
		    	},
	    	[&](){
				// cout << "after size: " << zdtree.collect_records(zdtree.root).size() << endl;
				// if (print_flag){
				// 	cout << "[touched inte]: " << zdtree.visited_inte << endl;
				// 	cout << "[touched leaf]: " << zdtree.visited_leaf << endl;
				// 	print_flag = false;
				// }
			} );
			cout << "# of points: " << zdtree.collect_records(new_ver).size() << endl;
			cout << "[batch_size]: " << num_processed << endl;
	    	// cout << fixed << setprecision(6) << "[leaf copy time]: " << zd_leaf_copy_time << endl;
	    	// cout << fixed << setprecision(6) << "[inte copy time]: " << zd_inte_copy_time << endl;
	    	cout << fixed << setprecision(6) << "[zdtree]: batch delete time (avg): " << zdtree_delete_avg << endl;
		}
    }


    // template<typename PT>
    // void batch_insert_test(PT P, size_t batch_percent){
	//     auto n = P.size();
	//     auto num_processed = batch_percent * n / 100;
	//     auto P2 = P.substr(0, num_processed);	// first x%
	//     parlay::parallel_for (0, P2.size(), [&](int i){
	// 	    P2[i].id = n + i;
	//     });
	//     auto P3 = P;	// 1+x% 
	//     P3.resize(n + num_processed);
	//     parlay::parallel_for (0, n + num_processed, [&](size_t i){
	// 	    if (i < n) P3[i] = P[i];
	// 	    else P3[i] = P2[i - n];
	//     });
	//     auto P4 = P.substr(num_processed, n - num_processed);
	
	//     ZDTree::Tree zdtree(leaf_size);
	//     // auto P_set = get_sorted_address(P); // build original tree
	//     auto P_set = get_sorted_points(P); // build original tree

	//     zdtree.build(P_set);
	//     // auto tree_hash = zdtree.tree_hash(zdtree.root);
	//     // cout << "tree hash: " << tree_hash << ", node cnt: " << zdtree.node_cnt << ", leaf cnt: " << zdtree.leaf_cnt << ", record cnt: " << zdtree.record_cnt << endl;
	//     // zdtree.clear();
	//     // P_set = get_sorted_address(P4);	
	//     // P_set = get_sorted_points(P4);	
	//     // zdtree.build(P_set);	// build for (1-x)%

	//     shared_ptr<ZDTree::BaseNode> new_ver = nullptr;
	//     auto zdtree_insert_avg = time_loop(
	// 	    3, 1.0, [&]() {},
	// 	    [&]() {
	// 		    // P_set = get_sorted_address(P2);
	// 		    P_set = get_sorted_points(P2);
	// 		    new_ver = zdtree.multi_version_batch_insert_sorted(P_set, zdtree.root);
	// 	    },
	//     [&](){} );

	//     // zdtree.node_cnt = zdtree.leaf_cnt = zdtree.record_cnt = 0;
	//     // auto old_ver_hash = zdtree.tree_hash(zdtree.root);
	//     // cout << "before insertion: " << old_ver_hash << ", node cnt: " << zdtree.node_cnt << ", leaf cnt: " << zdtree.leaf_cnt << ", record cnt: " << zdtree.record_cnt << endl;
	//     // zdtree.node_cnt = zdtree.leaf_cnt = zdtree.record_cnt = 0;
	//     // auto new_ver_hash = zdtree.tree_hash(new_ver);
	//     // cout << "after insertion " << new_ver_hash << ", node cnt: " << zdtree.node_cnt << ", leaf cnt: " << zdtree.leaf_cnt << ", record cnt: " << zdtree.record_cnt << endl;
	//     // auto [add, remove] = zdtree.diff(zdtree.root, new_ver, 64);
	//     // cout << "add, remove size = " << add.size() << ", " << remove.size() << endl;
	// 	// if (tree_hash != new_ver_hash || old_ver_hash == new_ver_hash){
	// 	// 	cout << "[ERROR] Insertion Failed (inconsistent). " << endl;
	// 	// }
	//     cout << fixed << setprecision(6) << "[zdtree]: insert time (avg): " << zdtree_insert_avg << endl;
	//     zdtree.clear();

	//     // P_set = get_sorted_address(P3);
	//     P_set = get_sorted_points(P3);
	//     zdtree.build(P_set); // build for 1+x%
	//     parlay::hashtable<parlay::hash_numeric<int> > table(P2.size(), parlay::hash_numeric<int>{});

	//     auto zdtree_delete_avg = time_loop(
	// 	    3, 1.0, [&]() { },
	// 	    [&]() {
	// 		    // P_set = get_sorted_address(P2);
	// 		    P_set = get_sorted_points(P2);
	// 		    new_ver = zdtree.multi_version_batch_delete_sorted(P_set, zdtree.root);
	// 	    },
	//     [&](){} );

	//     // zdtree.node_cnt = zdtree.leaf_cnt = zdtree.record_cnt = 0;
	//     // old_ver_hash = zdtree.tree_hash(zdtree.root);
	//     // cout << "before deletion: " << old_ver_hash << ", node cnt: " << zdtree.node_cnt << ", leaf cnt: " << zdtree.leaf_cnt << ", record cnt: " << zdtree.record_cnt << endl;
	//     // zdtree.node_cnt = zdtree.leaf_cnt = zdtree.record_cnt = 0;
	//     // new_ver_hash = zdtree.tree_hash(new_ver);
	//     // cout << "after deletion " << new_ver_hash << ", node cnt: " << zdtree.node_cnt << ", leaf cnt: " << zdtree.leaf_cnt << ", record cnt: " << zdtree.record_cnt << endl;

	//     // tie(add, remove) = zdtree.diff(zdtree.root, new_ver, 64);
	//     // cout << "add, remove size = " << add.size() << ", " << remove.size() << endl;
	// 	// if (tree_hash != new_ver_hash || old_ver_hash == new_ver_hash){
	// 	// 	cout << "[ERROR] Deletion Failed (inconsistent). " << endl;
	// 	// }

	//     cout << fixed << setprecision(6) << "[zdtree]: delete time (avg): " << zdtree_delete_avg << endl;
    // }

    template<typename PT>
    void batch_insert_test(PT P, PT P_insert, PT P_delete){
	    ZDTree::Tree zdtree(leaf_size);
	    // auto P_set = get_sorted_address(P); // build original tree
	    auto P_set = get_sorted_points(P); // build original tree
	    zdtree.build(P_set);
		zdtree.multi_version_roots.emplace_back(zdtree.root);

		auto tot_node_count = zdtree.num_of_nodes();
		cout << "[initial node count]: " << tot_node_count << endl;

	    shared_ptr<ZDTree::BaseNode> new_ver = nullptr;
	    auto zdtree_insert_avg = time_loop(
		    3, 1.0, [&]() {},
		    [&]() {
			    P_set = get_sorted_points(P_insert);
			    new_ver = zdtree.multi_version_batch_insert_sorted(P_set, zdtree.root);
		    },
	    [&](){} );

		zdtree.multi_version_roots.emplace_back(new_ver);

	    // zdtree.node_cnt = zdtree.leaf_cnt = zdtree.record_cnt = 0;
	    // auto old_ver_hash = zdtree.tree_hash(zdtree.root);
	    // cout << "before insertion: " << old_ver_hash << ", node cnt: " << zdtree.node_cnt << ", leaf cnt: " << zdtree.leaf_cnt << ", record cnt: " << zdtree.record_cnt << endl;
	    // zdtree.node_cnt = zdtree.leaf_cnt = zdtree.record_cnt = 0;
	    // auto new_ver_hash = zdtree.tree_hash(new_ver);
	    // cout << "after insertion: " << new_ver_hash << ", node cnt: " << zdtree.node_cnt << ", leaf cnt: " << zdtree.leaf_cnt << ", record cnt: " << zdtree.record_cnt << endl;
	    cout << fixed << setprecision(6) << "[zdtree]: insert time (avg): " << zdtree_insert_avg << endl;

		tot_node_count = zdtree.num_of_nodes();
		cout << "[after_insertion]: " << tot_node_count << endl;
		
	    // zdtree.clear();
		// P_set = get_sorted_points(P);
	    // zdtree.build(P_set); 

	    parlay::hashtable<parlay::hash_numeric<int> > table(P_delete.size(), parlay::hash_numeric<int>{});

		auto inserted_ver = new_ver;
	    auto zdtree_delete_avg = time_loop(
		    3, 1.0, [&]() { },
		    [&]() {
			    P_set = get_sorted_points(P_delete);
			    // new_ver = zdtree.multi_version_batch_delete_sorted(P_set, zdtree.root);
			    new_ver = zdtree.multi_version_batch_delete_sorted(P_set, inserted_ver);
		    },
	    [&](){} );

		zdtree.multi_version_roots.emplace_back(new_ver);		

	    // zdtree.node_cnt = zdtree.leaf_cnt = zdtree.record_cnt = 0;
	    // old_ver_hash = zdtree.tree_hash(zdtree.root);
	    // cout << "before deletion: " << old_ver_hash << ", node cnt: " << zdtree.node_cnt << ", leaf cnt: " << zdtree.leaf_cnt << ", record cnt: " << zdtree.record_cnt << endl;
	    // zdtree.node_cnt = zdtree.leaf_cnt = zdtree.record_cnt = 0;
	    // new_ver_hash = zdtree.tree_hash(new_ver);
	    // cout << "after deletion " << new_ver_hash << ", node cnt: " << zdtree.node_cnt << ", leaf cnt: " << zdtree.leaf_cnt << ", record cnt: " << zdtree.record_cnt << endl;

	    cout << fixed << setprecision(6) << "[zdtree]: delete time (avg): " << zdtree_delete_avg << endl;

		tot_node_count = zdtree.num_of_nodes();
		cout << "[after_deletion]: " << tot_node_count << endl;

		zdtree.multi_version_roots[1].reset();
		tot_node_count = zdtree.num_of_nodes();
		cout << "[after_remove_intermediate]: " << tot_node_count << endl;

    }


	template<class PT, class RQ>
	void range_count_test(PT P, RQ querys, parlay::sequence<size_t> &cnt){
	    size_t n = P.size();
	    auto P2 = P.substr(0, P.size() / 2);   // insert a half
	    for (size_t i = 0; i < P2.size(); i++){
		    P2[i].id = n + i;
	    }
	    ZDTree::Tree zdtree(leaf_size);

	    auto P_set = get_sorted_points(P);
	    zdtree.build(P_set);	// 	build for P

		parlay::sequence<size_t> rangeCnt(querys.size());
		for (size_t i = 0; i < querys.size(); i++){
			auto avg_time = time_loop(
				3, 1.0, 
				[&]() {},
				[&]() {					
					rangeCnt[i] = zdtree.range_count(querys[i], largest_mbr);
				},
				[&](){} );
			if (rangeCnt[i] != cnt[i]){
				cout << "[ERROR] Incorrect" << endl;
			}
			else{
				cout << fixed << setprecision(6) << rangeCnt[i] << " " << avg_time << endl;
			}
		}
		
		// auto rangeCnt_avg = time_loop(
		// 	3, 1.0, [&]() {},
		// 	[&]() {
		// 		parlay::parallel_for(
		// 			0, querys.size(),
		// 			[&]( size_t i ) {
		// 				// rangeCnt[i] = zdtree.range_count(querys[i]);
		// 				rangeCnt[i] = zdtree.range_count(querys[i], largest_mbr);
		// 		});
		// 	},
		// [&](){} );
		// cout << fixed << setprecision(6) << "[zdtree] range count time (avg): " << rangeCnt_avg << endl;

		// bool ok = true;
		// parlay::parallel_for(0, querys.size(), [&](size_t i){
		// 	if (cnt[i] != rangeCnt[i]){
		// 		ok = false;
		// 	}
		// });
		// if (!ok){
		// 	cout << "[ERROR] incorrect range count result !!!" << endl;
		// }

		// ofstream regionCntOut("zd_range_count.txt");
		
		// for (size_t i = 0; i < querys.size(); i++){
		// 	regionCntOut << rangeCnt[i] << endl;
		// 	if (cnt[i] != rangeCnt[i]){
		// 		cout << "[ERROR] incorrect range count result " << rangeCnt[i] << "-" << cnt[i] << endl;
		// 	}
		// }
		// cout << "range count finished." << endl;
	}

	template<class PT, class RQ>
	void range_report_test(PT P, RQ querys, parlay::sequence<size_t> &cnt){
	    ZDTree::Tree zdtree(leaf_size);
	    auto P_set = get_sorted_points(P);
	    zdtree.build(P_set);	// 	build for P

		parlay::sequence<parlay::sequence<Point> > rangeReport(querys.size());
		parlay::sequence<size_t> rangeReportCnt(querys.size(), 0);

		for (size_t i = 0; i < querys.size(); i++){
			rangeReport[i].resize(cnt[i]);	
		}

		for (size_t i = 0; i < querys.size(); i++){
			// print_mbr(querys[i]);
			auto avg_time = time_loop(
				3, 1.0, 
				[&]() {
					rangeReportCnt[i] = 0;
				},
				[&]() {					
					zdtree.range_report(querys[i], largest_mbr, rangeReportCnt[i], rangeReport[i]);
				},
				[&](){} );
			if (rangeReportCnt[i] != cnt[i]){
				cout << "[ERROR] Incorrect" << endl;
				cout << rangeReportCnt[i] << " " << cnt[i] << endl;
			}
			else{
				cout << fixed << setprecision(6) << rangeReportCnt[i] << " " << avg_time << endl;
			}
		}

		// auto rangeReport_avg = time_loop(
		// 	3, 1.0, [&]() {},
		// 	[&]() {
		// 		parlay::parallel_for(
		// 			0, querys.size(),
		// 			[&]( size_t i ) {
		// 				rangeReportCnt[i] = 0;
		// 				// zdtree.range_report(querys[i], rangeReportCnt[i], rangeReport[i]);
		// 				zdtree.range_report(querys[i], largest_mbr, rangeReportCnt[i], rangeReport[i]);
		// 		});
		// 	},
		// [&](){} );

		// cout << fixed << setprecision(6) << "[zdtree] range report time (avg): " << rangeReport_avg << endl;

		// bool ok = true;
		// parlay::parallel_for(0, querys.size(), [&](size_t i){
		// 	if (cnt[i] != rangeReportCnt[i]){
		// 		ok = false;
		// 	}
		// });
		// if (!ok){
		// 	cout << "[ERROR] incorrect range count result !!!" << endl;
		// }

		// ofstream regionReportOut("zd_range_report.txt");
		// for (size_t i = 0; i < querys.size(); i++){
		// 	regionReportOut << rangeReportCnt[i] << endl;
		// }
		// cout << "range report finshied." << endl;
	}


	template<class PT>
	void knn_test(PT P, size_t k = 10, size_t q_num = 50000){
	    ZDTree::Tree zdtree(leaf_size);
	    // auto P_set = get_sorted_address(P);
	    auto P_set = get_sorted_points(P);
	    zdtree.build(P_set);	// 	build for P

		parlay::sequence<FT> knn_sqrdis(q_num);

		auto knnReport_avg = time_loop(
			3, 1.0, [&]() {},
			[&]() {
				for (size_t i = 0; i < q_num; i++){
					knn_sqrdis[i] = zdtree.knn_report(k, P[i], largest_mbr).top().second;
				}
				// parlay::parallel_for(
				// 	0, q_num,
				// 	[&]( size_t i ) {
				// 		knn_sqrdis[i] = zdtree.knn_report(k, P[i], largest_mbr).top().second;
				// });
			},
		[&](){} );
		cout << fixed << setprecision(6) << "[zdtree] knn report time (avg): " << knnReport_avg << endl;

		ofstream knnReportOut("zd-knn.res");
		for (size_t i = 0; i < q_num; i++){
			knnReportOut << knn_sqrdis[i] << endl;
			// size_t id = i < 5 ? i : n / 2 + i;
			// knnReportOut << knn_sqrdis[i] << ", brute force: " << geobase::knn_bf(k, *P_set[i], P) << endl;
			// knnReportOut << knn_sqrdis[id] << ", brute force: " << knn_bf(k, P_set[id], P_set) << endl;
		}
		// cout << "knn finished." << endl;
	}


    template<class PT, class RQ>
    void zdtree_test(PT P, RQ querys, int q_type){
	    size_t n = P.size();
	    auto P2 = P.substr(0, P.size() / 2);   // insert a half
	    for (size_t i = 0; i < P2.size(); i++){
		    P2[i].id = n + i;
	    }
	    ZDTree::Tree zdtree(leaf_size);
	    auto P_set = get_sorted_address(P);
	    zdtree.build(P_set);	// 	build for P

	    // auto phash = zdtree.tree_hash(zdtree.root);	// hash value for construct P directly
	    // cout << "build zdtree finished." << endl;

	    /* insertion test */
	    // if (q_type == 8){
		    // P_set = get_sorted_address(P2);
		    // zdtree.multi_version_batch_insert_sorted(P_set, zdtree.root); // insert P2 to P
		    // cout << "insertion finished." << endl;
		    // // auto p2hash = zdtree.tree_hash(zdtree.root); // hash value for insert P2
		    // zdtree.multi_version_batch_delete_sorted(P_set, zdtree.root);
		    // cout << "deletion finished." << endl;
		    // auto p4hash = zdtree.tree_hash(zdtree.root); // hash value for delete P2
		    // cout << "phash = " << phash << ", p4hash = " << p4hash << endl;	// these two should be the same
		    // cout << "p3hash = " << p3hash << ", p2hash = " << p2hash << endl;	// these two should be the same
	    // }

	    /* in order traverse */
	    // zdtree.in_order_traverse(zdtree.root, 0, 0, 32, true);
	    // zdtree.in_order_traverse(zdtree.root);
	    // cout << "traverse finished." << endl;

	    /* range count test */
	    // if (q_type == 4){
		    parlay::sequence<size_t> rangeCnt(querys.size());
		    auto rangeCnt_avg = time_loop(
			    3, 1.0, [&]() {},
        	    [&]() {
				    parlay::parallel_for(
    	        	    0, querys.size(),
					    [&]( size_t i ) {
						    // rangeCnt[i] = zdtree.range_count(querys[i]);
						    rangeCnt[i] = zdtree.range_count(querys[i], largest_mbr);
            	    });
			    },
		    [&](){} );
		    cout << fixed << setprecision(6) << "zdtree range count time (avg): " << rangeCnt_avg << endl;
		    ofstream regionCntOut("zd_range_count.txt");
		    for (size_t i = 0; i < querys.size(); i++){
			    regionCntOut << rangeCnt[i] << endl;
		    }
			// cout << "range count finished." << endl;
	    // }

	    /* range report test */
	    // if (q_type == 8){
		    parlay::sequence<parlay::sequence<Point*> > rangeReport(querys.size());
		    parlay::sequence<size_t> rangeReportCnt(querys.size(), 0);
		    for (size_t i = 0; i < querys.size(); i++){
			    rangeReport[i].resize(maxSize);	
		    }
		    auto rangeReport_avg = time_loop(
			    3, 1.0, [&]() {},
        	    [&]() {
				    parlay::parallel_for(
    	        	    0, querys.size(),
					    [&]( size_t i ) {
						    rangeReportCnt[i] = 0;
						    // zdtree.range_report(querys[i], rangeReportCnt[i], rangeReport[i]);
						    zdtree.range_report(querys[i], largest_mbr, rangeReportCnt[i], rangeReport[i]);
            	    });
			    },
		    [&](){} );
		    cout << fixed << setprecision(6) << "zdtree range report time (avg): " << rangeReport_avg << endl;
		    ofstream regionReportOut("zd_range_report.txt");
		    for (size_t i = 0; i < querys.size(); i++){
			    regionReportOut << rangeReportCnt[i] << endl;
		    }
			// cout << "range report finshied." << endl;
	    // }

	    /* knn report test */
	    // if (q_type == 1){
		    size_t k = 10;
		    P_set.resize(n);
		    parlay::parallel_for(0, n, [&](size_t i){
			    P_set[i] = i < P.size() ? &P[i] : &P2[i - P.size()];
		    });
		    size_t knn_test_size = n;
		    parlay::sequence<FT> knn_sqrdis(knn_test_size);
		    auto knnReport_avg = time_loop(
			    3, 1.0, [&]() {},
        	    [&]() {
				    parlay::parallel_for(
    	        	    0, knn_test_size,
					    [&]( size_t i ) {
						    // knn_sqrdis[i] = zdtree.knn_report(k, &P[i]).top().second;
						    // knn_sqrdis[i] = zdtree.knn_report(k, &P[i], largest_mbr).top().second;
						    knn_sqrdis[i] = zdtree.knn_report(k, P_set[i], largest_mbr).top().second;
            	    });
			    },
		    [&](){} );
		    cout << fixed << setprecision(6) << "zdtree knn report time (avg): " << knnReport_avg << endl;
		    ofstream knnReportOut("zd_knn_report.txt");
		    for (size_t i = 0; i < 10; i++){
			    size_t id = i < 5 ? i : n / 2 + i;
			    // knnReportOut << knn_sqrdis[i] << ", brute force: " << geobase::knn_bf(k, *P_set[i], P) << endl;
			    knnReportOut << knn_sqrdis[id] << ", brute force: " << knn_bf(k, *P_set[id], P_set) << endl;
		    }
			cout << "knn finished." << endl;
	    // }
    }

    template<class PT, class RQ>
    void CPAMZ_test(PT &P, RQ querys){
	    auto CPAMZ = Morton::CPAMZ_init(P);
    #ifdef DEBUG
	    cout << "CPAMZ size = " << CPAMZ.size() << endl;
    #endif
	    parlay::internal::timer t;
	    parlay::sequence<Point> range_query_res;
	    for (auto query_mbr: querys){
		    range_query_res = Morton::range_report(CPAMZ, query_mbr);
	    }
	    cout << "Z-CPAM query time (avg): " << t.next_time() << endl;
    #ifdef DEBUG
	    // auto range_query_res = Morton::range_report(CPAMZ[0], query_mbr);
	    cout << "[Z-CPAM RESULT] range query results: " << range_query_res.size() << "| ";
	    for (auto p: range_query_res){
		    cout << "(" << p.x << ", " << p.y << ")" << " ";
	    }
	    cout << endl;	
    #endif
    }

    template<class PT, class RQ>
    void zMAP_test(PT &P, RQ querys){
	    auto zMAP = Morton::zMAP_init(P);
    #ifdef DEBUG
	    cout << "zMAP size = " << zMAP.size() << endl;
    #endif
	    parlay::internal::timer t;
	    parlay::sequence<Point> range_query_res;
	    for (auto query_mbr: querys){
		    range_query_res = Morton::zMAP_range_report(zMAP, query_mbr);
	    }
	    cout << "zMAP query time (avg): " << t.next_time() << endl;
	    // auto range_query_res = Morton::range_report(CPAMZ[0], query_mbr);
    #ifdef DEBUG
	    cout << "[Z-Map RESULT] range query results: " << range_query_res.size() << "| ";
	    for (auto p: range_query_res){
		    cout << "(" << p.x << ", " << p.y << ")" << " ";
	    }
	    cout << endl;	
    #endif
    }
}