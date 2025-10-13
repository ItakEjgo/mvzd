#pragma once 

#include <bits/stdc++.h>

#include <cpam/cpam.h>
// #include <pam/pam.h>
#include <parlay/primitives.h>
#include "geobase.h"
#include <parlay/internal/get_time.h>
#include "pam/utils.h"

#define BR_MBR
#define SEQ


namespace CPAMBB_BOX{
	using namespace std;
	using namespace geobase;
	using parlay::sequence;
	using parlay::par_do;
	using parlay::par_do_if;

	using key_type = pair<unsigned long long, unsigned long long>;	// morton_id, id
	using val_type = Rectangle;
	using aug_type = pair<Bounding_Box, size_t>;

	//	CPAM entry
	struct entry {
		using key_t = key_type;
		using val_t = val_type;
		using aug_t = aug_type;

		static inline bool comp(key_t a, key_t b) { return a < b; }
		static aug_t get_empty() { return make_pair(Bounding_Box{Point(1e60, 1e60), Point(-1, -1)}, 0); }
		static aug_t from_entry(key_t k, val_t v) { return make_pair(Bounding_Box(Point(v.x_low, v.y_low), Point(v.x_high, v.y_high)), 1); }
		static aug_t combine(aug_t a, aug_t b) { return make_pair(merge_mbr(a.first, b.first), a.second + b.second); }
	};

	// using zmap = cpam::aug_map<entry, 32>;
	using zmap = cpam::aug_map<entry, 16>;
	using par = std::tuple<entry::key_t, entry::val_t>;

	// template<class T, class MBR>
	// auto filter_range(T &tree, MBR query_mbr, bool use_hilbert = false){
	// 	auto BL = query_mbr.first;
	// 	auto UR = query_mbr.second;

	// 	auto fpt = [&](par cur){ 
	// 		if (BL.x <= get<1>(cur).x && get<1>(cur).x <= UR.x &&
	// 			BL.y <= get<1>(cur).y && get<1>(cur).y <= UR.y){
	// 			return true;
	// 		}
	// 		return false;
	// 	};

	// 	auto fbb = [&](auto cur){ 
	// 		return !mbr_exclude_mbr(query_mbr, cur.first);
	// 	};

	// 	// return zmap::aug_filter(tree, fbb);
	// 	// return  zmap::filter(zmap::aug_filter(tree, fbb), fpt);
	// 	return zmap::aug_filter2(tree, fpt, fbb);
	// }

	// template<typename M, typename DIFF>
	// auto map_spatial_diff(M &lhs, M &rhs, Bounding_Box &query, DIFF &ret_diff){
	// 	auto filtered_lhs = filter_range(lhs, query);
	// 	auto filtered_rhs = filter_range(rhs, query);
	// 	ret_diff.add = zmap::values(zmap::map_difference(filtered_rhs, filtered_lhs));
	// 	ret_diff.add_cnt = ret_diff.add.size();
	// 	ret_diff.remove = zmap::values(zmap::map_difference(filtered_lhs, filtered_rhs));
	// 	ret_diff.remove_cnt = ret_diff.remove.size();
	// 	return;
	// }

	// template<typename M>
	// auto map_diff(M &lhs, M &rhs){
	// 	auto add = zmap::values(zmap::map_difference(rhs, lhs));
	// 	auto remove = zmap::values(zmap::map_difference(lhs, rhs));
	// 	return make_tuple(add, remove);
	// }

	template<class BOX>
	auto map_init(BOX &P, bool use_hilbert = false){
		size_t n = P.size();
		parlay::sequence<unsigned long long> Z_vals(n);

		parlay::parallel_for(0, n, [&](int i){
			Z_vals[i] = P[i].get_centroid().interleave_bits();
		});
		parlay::sequence<par> entries(n);
		parlay::parallel_for(0, n, [&](int i){
			entries[i] = {{Z_vals[i], P[i].id}, P[i]};
		});
		zmap m1(entries);
		return m1;
	}	

	template<typename PT, typename M>
	auto map_insert(PT &P, M &mmp, bool use_hilbert = false){
		size_t n = P.size();
		parlay::sequence<unsigned long long> Z_vals(n);

		parlay::internal::timer t("debug", false);
		parlay::parallel_for(0, n, [&](int i){
			Z_vals[i] = P[i].get_centroid().interleave_bits();
		});
		auto insert_pts = parlay::sequence<par>::uninitialized(n);
		parlay::parallel_for(0, n, [&](int i){
			insert_pts[i] = {{Z_vals[i], P[i].id}, P[i]};
			// insert_pts[i] = {P[i]->id, P[i]};
		});
		t.next("init time");
		auto m2 = zmap::multi_insert(mmp, insert_pts);
		t.next("insert time");

		return m2;
	}

	template<typename PT, typename M>
	auto map_delete(PT &P, M &mmp, bool use_hilbert = false){
		size_t n = P.size();
		parlay::sequence<unsigned long long> Z_vals(n);

		parlay::parallel_for(0, n, [&](int i){
			Z_vals[i] = P[i].get_centroid().interleave_bits();
		});
		parlay::sequence<par> delete_pts(n);
		parlay::parallel_for(0, n, [&](int i){
			delete_pts[i] = {{Z_vals[i], P[i].id}, P[i]};
		});
		auto m2 = zmap::multi_delete(mmp, delete_pts);

		return m2;
	}

	template<class T, class MBR>
	auto intersects_with(T &tree, MBR &query_mbr, parlay::sequence<Rectangle> &out, bool use_hilbert = false){
		// auto ret = zmap::values(filter_range(tree, query_mbr, use_hilbert));
		auto f = [&](auto &cur){ 
			return mbr_mbr_relation(cur, query_mbr);
		};

		int64_t ret = 0;
		// for (auto i = 0; i < 20; i++){
			// ret = 0;
			// zmap::range_report_filter(tree, f, ret, out);
		zmap::intersects_filter(tree, f, ret, out);
		// }
		return ret;
	}

	// template<class T, class MBR>
	// auto range_count(T &zCPAM, MBR &query_mbr, bool use_hilbert = false){
	// 	auto f = [&](auto cur){ 
	// 		return mbr_mbr_relation(cur, query_mbr);
	// 	};

	// 	auto f2 = [&](auto cur){ 
	// 		return point_in_mbr(cur, query_mbr);
	// 	};

	// 	// auto res = zmap::range_count_filter(zCPAM, f, f2);
	// 	auto res = zmap::range_count_filter2(zCPAM, f, f2);
	// 	return res;
	// }

	// template<class T, class MBR>
	// auto range_report(T &tree, MBR &query_mbr, parlay::sequence<Point> &out, bool use_hilbert = false){
	// 	// auto ret = zmap::values(filter_range(tree, query_mbr, use_hilbert));
	// 	auto f = [&](auto &cur){ 
	// 		return mbr_mbr_relation(cur, query_mbr);
	// 	};

	// 	int64_t ret = 0;
	// 	// for (auto i = 0; i < 20; i++){
	// 		// ret = 0;
	// 		// zmap::range_report_filter(tree, f, ret, out);
	// 		zmap::range_report_filter2(tree, f, ret, out);
	// 	// }
	// 	return ret;
	// }

	// template<typename M, typename DIFF>
	// auto plain_map_spatial_diff(M &lhs, M &rhs, Bounding_Box &query, DIFF &ret_diff, parlay::sequence<Point> &l_pts, parlay::sequence<Point> &r_pts){
	// 	// auto l_pts = zmap::values(filter_range(lhs, query));
	// 	// auto r_pts = zmap::values(filter_range(rhs, query));

	// 	auto ret1 = range_report(lhs, query, l_pts);
	// 	auto ret2 = range_report(rhs, query, r_pts);
	// 	l_pts.resize(ret1);
	// 	r_pts.resize(ret2);
	// 	// print_Pset_info(l_pts, "lpts");
	// 	// for (auto &pt: l_pts){
	// 	// 	cout << pt.morton_id << endl;
	// 	// }
	// 	// print_Pset_info(r_pts, "rpts");
	// 	// for (auto &pt: r_pts){
	// 	// 	cout << pt.morton_id << endl;
	// 	// }
	// 	// cout << "debug: " << l_pts.size() << ", " << r_pts.size() << endl;
	// 	// auto [add, remove] = merge_pts(l_pts, r_pts);
	// 	// cout << "l size = " << l_pts.size() << endl;
	// 	// cout << "r size = " << r_pts.size() << endl;
		
	// 	merge_pts(l_pts, r_pts, ret_diff);
		
	// 	// geobase::print_Pset_info(add, "add");
	// 	// geobase::print_Pset_info(remove, "remove");
	// 	return;
	// }

	// template<typename T>
	// auto knn(T &tree, geobase::Point &query_point, size_t &k){
	// 	auto f = [&](auto cur_pt){ return point_point_sqrdis(cur_pt, query_point); };

	// 	auto f2 = [&](auto cur_mbr){ return point_mbr_sqrdis(query_point, cur_mbr); };

	// 	priority_queue<nn_pair, vector<nn_pair>, nn_pair_cmp> nn_res;
	// 	zmap::knn_filter(tree, f, f2, k, nn_res);
	// 	return nn_res;
	// }

	// template<typename PT, typename M>
	// auto map_commit(M &mmp, PT &P_insert, PT &P_delete){
	// 	// parlay::internal::timer t("CPAM-BB breakdown", true);
	// 	auto new_ver = map_delete(P_delete, mmp);	//	new	version
	// 	// t.next("delete time");
	// 	new_ver = map_insert(P_insert, new_ver); 
	// 	// t.next("insert time");
	// 	return new_ver;
	// }

	// template<typename M>
	// auto map_merge(M &base, M &v1, M &v2){
	// 	parlay::internal::timer t("merge", false);
	// 	auto [add1, remove1] = map_diff(base, v1);
	// 	auto [add2, remove2] = map_diff(base, v2); 
	// 	t.next("diff time");

	// 	auto [insert1, delete1, update1] = geobase::filter_diff_results(add1, remove1);
	// 	auto [insert2, delete2, update2] = geobase::filter_diff_results(add2, remove2);

	// 	auto [no_conflict_insert, conflict_insert] = geobase::merge_by_id_with_conflict(insert1, insert2);	// insert conflict
	// 	auto [no_conflict_update, conflict_update] = geobase::merge_by_id_with_conflict(update1, update2);	// update conflict
	// 	//	delate = delete + update
	// 	auto [no_conflict_delete, conflict_delete] = geobase::merge_by_id_with_conflict(delete1, update2);	//	delete and update conflict
	// 	auto [no_conflict_delete1, conflict_delete1] = geobase::merge_by_id_with_conflict(delete2, update1);	//	delete and update conflict, another case

	// 	no_conflict_delete.append(no_conflict_delete1);
	// 	conflict_delete.append(conflict_delete1);

	// 	t.next("conflict time");

	// 	auto new_ver = map_commit(base, no_conflict_insert, no_conflict_delete);
	// 	t.next("commit time");

		
	// 	return make_tuple(new_ver, conflict_insert, conflict_update, conflict_delete);
	// }

	// //	return size of interior nodes and sizeof leaf nodes size, respectively
	// auto size_in_bytes(){
	// 	size_t inte_used = zmap::GC::used_node();
	// 	size_t internal_nodes_space = sizeof(typename zmap::GC::regular_node) * inte_used;
	// 	auto [used, unused] = parlay::internal::get_default_allocator().stats();
	// 	return make_tuple(internal_nodes_space, used);
	// }

	auto print_stats(){
		zmap::GC::print_stats();
	}

}