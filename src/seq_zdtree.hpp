#pragma once 

#include <bits/stdc++.h>

#include <cpam/cpam.h>
#include <parlay/primitives.h>
#include "geobase.h"
#include <parlay/internal/get_time.h>
#include "pam/utils.h"

// #define USE_MBR
// #define SEQ

namespace SEQ_ZDTree{

	using namespace std;
	using namespace geobase;
	using parlay::sequence;
	using parlay::par_do;
	using parlay::par_do_if;

	struct BaseNode{
		BaseNode(){}

		virtual~BaseNode() = default;
		virtual bool is_leaf(){ return false; }
		virtual size_t get_num_points(){ return 0; }
	};

	struct InteNode: BaseNode{
		shared_ptr<BaseNode> l_son, r_son;
		size_t num_pts;

		InteNode(): l_son(nullptr), r_son(nullptr), num_pts(0){}
		// InteNode(InteNode &x): l_son(x->l_son), r_son(x->r_son), num_pts(x->num_pts){}
		
		virtual bool is_leaf(){ return false; }
		virtual size_t get_num_points(){ return num_pts; }
	};

	struct LeafNode: BaseNode{
		sequence<Point*> records = {};
		// parlay::slice<Point*, Point*> records {nullptr, nullptr};

		template<typename Records>
		LeafNode(Records &r){
			records.resize(r.size());
			for (size_t i = 0; i < r.size(); i++){
				records[i] = r[i];
			}
		}
		
		template<typename Records>
		LeafNode(Records &p, size_t l, size_t r){
			records.resize(r - l);
			for (size_t i = l; i < r; i++){
				records[i - l] = p[i];
			}
		}

		// LeafNode(LeafNode &x): records(x->records){}

		virtual bool is_leaf(){ return true; }
		virtual size_t get_num_points(){ return records.size(); }

		void print_records(){
			cout << records.size() << endl;
			for (size_t i = 0; i < records.size(); i++){
				cout << "(" << records[i]->x << ", " << records[i]->y << ")" << endl; 
			}
		}
	};

	class Tree{
	private:
		size_t visited_leaf = 0;
		size_t visited_inte = 0;
		size_t granularity_cutoff = 1000;

	public:
		size_t inte_node_cnt = 0;
		size_t leaf_node_cnt = 0;
		size_t leaf_size = 1;
		size_t node_cnt = 0;
		size_t leaf_cnt = 0;
		size_t record_cnt = 0;
		size_t case0 = 0;
		size_t case1 = 0;
		size_t case2 = 0;
		size_t case3 = 0;
		size_t case4 = 0;
		FT case0_time = 0;
		FT case1_time = 0;
		FT case2_time = 0;
		FT case3_time = 0;
		FT case4_time = 0;
		shared_ptr<BaseNode> root;
		vector<shared_ptr<BaseNode> > multi_version_roots = {root};

		Tree(size_t _leaf_sz);
		// ~Tree();

		shared_ptr<BaseNode> node_copy(shared_ptr<BaseNode> &x);

		// entrance of building zdtree & tree construction
		shared_ptr<BaseNode> build(sequence<Point*> &P, size_t l, size_t r, size_t b);
		void build(sequence<Point*> &P);
		void clear();

		void in_order_traverse(shared_ptr<BaseNode> &x);
		void in_order_traverse(shared_ptr<BaseNode> &x, unsigned x_prefix, unsigned y_prefix, size_t b, bool x_splitter);
		int  tree_hash(shared_ptr<BaseNode> &x);
		void tree_compress(shared_ptr<BaseNode> &x);

		void tree_node_count(shared_ptr<BaseNode> &x, set<uintptr_t> &s);	//	store the nodes of a root in a set s
		size_t node_count();	//	count the total number of nodes among all versions

		void merge_nodes(shared_ptr<BaseNode> &lhs, shared_ptr<BaseNode> &rhs, shared_ptr<InteNode> &cur); 	// merge two sons to current node.
		void delete_merge_nodes(shared_ptr<BaseNode> &L, shared_ptr<BaseNode> &R, InteNode* cur_node);

		shared_ptr<InteNode> create_internal(shared_ptr<BaseNode> &L, shared_ptr<BaseNode> &R);	//	create an internal node, do not store pointers to original records.
		shared_ptr<LeafNode> create_leaf(sequence<Point*> &P, size_t l, size_t r, size_t b); // create a leaf, store all (pointers of) records.

		// size_t num_covered_points(shared_ptr<BaseNode> x);

		// 	insertion
		void batch_insert_sorted(sequence<Point*> &P);
		void batch_insert_sorted_node(shared_ptr<BaseNode> &x, sequence<Point*> &P, size_t l, size_t r, size_t b);
		// 	persistent insertion
		shared_ptr<BaseNode> multi_version_batch_insert_sorted(sequence<Point*> &P, shared_ptr<BaseNode> &old_root);
		void multi_version_batch_insert_sorted_node(shared_ptr<BaseNode> &x, sequence<Point*> &P, size_t l, size_t r, size_t b);


		//	deletion
		void batch_delete_sorted(sequence<Point*> &P);
		void batch_delete_sorted_node(shared_ptr<BaseNode> &x, sequence<Point*> &P, size_t l, size_t r, size_t b);
		template <class HashTable>
		void batch_delete_sorted(sequence<Point*> &P, HashTable& ht);
		template <class HashTable>
		void batch_delete_sorted_node(shared_ptr<BaseNode> &x, sequence<Point*> &P, size_t l, size_t r, size_t b, HashTable &ht);
		//	multi-version deletion
		shared_ptr<BaseNode> multi_version_batch_delete_sorted(sequence<Point*> &P, shared_ptr<BaseNode> &old_root);
		void multi_version_batch_delete_sorted_node(shared_ptr<BaseNode> &x, sequence<Point*> &P, size_t l, size_t r, size_t b);

		auto filter_diff_results(sequence<Point*> &add, sequence<Point*> &remove);
		
		auto collect_records(shared_ptr<BaseNode> &x);
		// auto diff(shared_ptr<BaseNode> &tree1, shared_ptr<BaseNode> &tree2);
		auto diff(shared_ptr<BaseNode> &node1, shared_ptr<BaseNode> &node2, size_t b);
		auto leaf_leaf_diff(shared_ptr<BaseNode> &leaf1, shared_ptr<BaseNode> &leaf2);
		auto leaf_inte_diff(sequence<Point*> &P, size_t l, size_t r, shared_ptr<BaseNode> &inte, size_t b);
		auto two_version_diff(shared_ptr<BaseNode> &node1, shared_ptr<BaseNode> &node2);

		//	commit, based on op do multi-version insertion, deletion, and update
		// auto commit(size_t op, shared_ptr<BaseNode> &old_version, sequence<Point*> &P_old, sequence<Point*> &P_new);
		// merge, works like git merge, return a conflict set
		auto merge(shared_ptr<BaseNode> &base, shared_ptr<BaseNode> &node1, shared_ptr<BaseNode> &node2);
		
		
		// intervals are left close right open: i.e., [L, R)
		// 2D range report
		#ifdef USE_MBR
		template <class Out> 
		void range_report_node(shared_ptr<BaseNode> &x, Bounding_Box &query_mbr, size_t &cnt, Out &out);
		template <class Out>
		void range_report(Bounding_Box &query_mbr, size_t &cnt, Out &out);
		#endif

		template <class Out> 
		void range_report_node(shared_ptr<BaseNode> &x, Bounding_Box &query_mbr, Bounding_Box &cur_mbr, FT x_prefix, FT y_prefix, size_t b, bool x_splitter, size_t &cnt, Out &out);
		template <class Out>
		void range_report(Bounding_Box &query_mbr, Bounding_Box &cur_mbr, size_t &cnt, Out &out);

		// 2D range count
		#ifdef USE_MBR
		size_t range_count_node(shared_ptr<BaseNode> &x, Bounding_Box &query_mbr);
		size_t range_count(Bounding_Box &query_mbr);
		size_t range_count_node(shared_ptr<BaseNode> &x, Bounding_Box &query_mbr, Bounding_Box &cur_mbr, bool x_splitter);
		#endif

		size_t range_count_node(shared_ptr<BaseNode> &x, Bounding_Box &query_mbr, Bounding_Box &cur_mbr, FT x_prefix, FT y_prefix, size_t b, bool x_splitter);
		size_t range_count(Bounding_Box &query_mbr, Bounding_Box &cur_mbr);

		// 2D k nearest neighbor report
		#ifdef USE_MBR
		template<class T> 
		void knn_report_node(shared_ptr<BaseNode> &x, size_t &k, Point* query_point, T &nn_res);
		auto knn_report(size_t &k, Point* query_point);
		#endif

		template<class T> 
		void knn_report_node(shared_ptr<BaseNode> &x, size_t &k, Point* query_point, Bounding_Box &cur_mbr, FT x_prefix, FT y_prefix, size_t b, bool x_splitter, T &nn_res);
		auto knn_report(size_t &k, Point* query_point, Bounding_Box &cur_mbr);
	};
		
	//	filter inserted points, deleted points, and updated points from add and remove sets. The updated points contain the new coordinates
	auto Tree::filter_diff_results(sequence<Point*> &add, sequence<Point*> &remove){
		auto id_cmp = [&](auto lhs, auto rhs){ return lhs->id < rhs->id; };
		auto sorted_add = parlay::sort(add, id_cmp);
		auto sorted_remove = parlay::sort(remove, id_cmp);

		sequence<Point*> insert_points = {}, delete_points = {}, update_points = {};
		size_t i = 0, j = 0;
		while (i < sorted_add.size() && j < sorted_remove.size()){
			if (sorted_add[i]->id < sorted_remove[j]->id){	//	we should insert sorted_add[i]
				insert_points.emplace_back(sorted_add[i++]);
			}
			else if (sorted_add[i]->id == sorted_remove[j]->id){	//	exist in both add and remove, indicate it is an update point
				update_points.emplace_back(sorted_add[i++]);
				j++;
			}
			else{	//	we should remove sorted_remove[j]
				delete_points.emplace_back(sorted_remove[j++]);
			}
		}
		while (i < sorted_add.size()) insert_points.emplace_back(sorted_add[i++]);
		while (j < sorted_remove.size()) delete_points.emplace_back(sorted_remove[j++]);
		return make_tuple(insert_points, delete_points, update_points);
	}
	
	auto Tree::collect_records(shared_ptr<BaseNode> &x){
		sequence<Point*> ret = {};
		if (!x){
			return ret;
		}
		if (x->is_leaf()){
			auto cur_leaf = static_cast<LeafNode*>(x.get());
			ret = cur_leaf->records;
			return ret;
		}
		auto cur_inte = static_cast<InteNode*>(x.get());
		sequence<Point*> R;
		auto collect_left = [&](){ret = collect_records(cur_inte->l_son); }; 
		auto collect_right = [&](){R = collect_records(cur_inte->r_son); };
		par_do_if(x->get_num_points() >= granularity_cutoff,
			collect_left,
			collect_right
		);
		ret.append(R);
		return ret;
	}

	auto Tree::leaf_leaf_diff(shared_ptr<BaseNode> &leaf1, shared_ptr<BaseNode> &leaf2){
		sequence<Point*> add = {}, remove = {};
		auto cur_leaf1 = static_cast<LeafNode*>(leaf1.get());
		auto cur_leaf2 = static_cast<LeafNode*>(leaf2.get());

		auto not_in_leaf2 = [&](auto pt){
			// return ht.find(x->id) == static_cast<int>(x->id);	// id is unsigned, but hash table is signed (needs empty -1).
			for (size_t i = 0; i < cur_leaf2->records.size(); i++){
				if (pt->id == cur_leaf2->records[i]->id && *pt == *cur_leaf2->records[i] ) return false;
			}
			return true;
		};
		remove = parlay::filter(cur_leaf1->records, not_in_leaf2);

		auto not_in_leaf1 = [&](auto pt){
			// return ht.find(x->id) == static_cast<int>(x->id);	// id is unsigned, but hash table is signed (needs empty -1).
			for (size_t i = 0; i < cur_leaf1->records.size(); i++){
				if (pt->id == cur_leaf1->records[i]->id && *pt == *cur_leaf1->records[i]) return false;
			}
			return true;
		};
		add = parlay::filter(cur_leaf2->records, not_in_leaf1);
		return make_tuple(add, remove);	
	}

	//	note, this returns {add, remove} w.r.t. leaf, pay attention to the order (especially when swapped) 
	//	leaf must not be nullptr
	auto Tree::leaf_inte_diff(sequence<Point*> &P, size_t l, size_t r, shared_ptr<BaseNode> &x, size_t b){
		sequence<Point*> add = {}, remove = {};
		if (!x){	//	inte could be nullptr if only one branch exists
			remove = P.substr(l, r - l);
			return make_tuple(add, remove);
		}
		if (x->is_leaf()){	//	Base case when we meet two leaf nodes
			auto cur_leaf = static_cast<LeafNode*>(x.get());
			for (size_t i = l; i < r; i++){	// in P[l..r] but not in the leaf
				bool flag = true; 
				for (size_t j = 0; j < cur_leaf->records.size(); j++){
					if (P[i]->id == cur_leaf->records[j]->id && *P[i] == *cur_leaf->records[j]){
						flag = false;
						break;
					}
				}
				if (flag) remove.emplace_back(P[i]);
			}

			for (size_t j = 0; j < cur_leaf->records.size(); j++){	// in the leaf but not in P[l..r]
				bool flag = true; 
				for (size_t i = l; i < r; i++){	
					if (P[i]->id == cur_leaf->records[j]->id && *P[i] == *cur_leaf->records[j]){
						flag = false;
						break;
					}
				}
				if (flag) add.emplace_back(cur_leaf->records[j]);
			}
			return make_tuple(add, remove);
		}
		//	divide into two subcases
		auto splitter = split_by_bit(P, l, r, b);
		auto cur_inte = static_cast<InteNode*>(x.get());
		sequence<Point*> L_add, L_remove, R_add, R_remove;

		auto diff_left = [&](){ tie(L_add, L_remove) = leaf_inte_diff(P, l, splitter, cur_inte->l_son, b - 1); };
		auto diff_right = [&](){ tie(R_add, R_remove) = leaf_inte_diff(P, splitter, r, cur_inte->r_son, b - 1); };
		par_do_if(r - l + x->get_num_points() >= granularity_cutoff,
			diff_left,
			diff_right);
		L_add.append(R_add);
		L_remove.append(R_remove);
		return make_tuple(L_add, L_remove);
	}

	auto Tree::diff(shared_ptr<BaseNode> &node1, shared_ptr<BaseNode> &node2, size_t b){
		if (node1 == node2){	//	case0, no difference, all empty or point to the same (shared) node
			utils::fetch_and_add(&case0, 1);
			sequence<Point*> add = {}, remove = {};
			return make_tuple(add, remove);
		}
		if (!node1){	//	case1, node1 is empty, we need to add all points in node2
			utils::fetch_and_add(&case1, 1);
			sequence<Point*> add = {}, remove = {};
			auto time_s = chrono::high_resolution_clock::now();
			add = collect_records(node2);
			auto time_t = chrono::high_resolution_clock::now();
			auto duration = chrono::duration_cast<chrono::microseconds>(time_t - time_s);
			case1_time += duration.count();
			return make_tuple(add, remove);
		}
		if (!node2){	//	case1, node2 is empty, we need to remove all points in node1
			utils::fetch_and_add(&case1, 1);
			sequence<Point*> add = {}, remove = {};
			auto time_s = chrono::high_resolution_clock::now();
			remove = collect_records(node1); 
			auto time_t = chrono::high_resolution_clock::now();
			auto duration = chrono::duration_cast<chrono::microseconds>(time_t - time_s);
			case1_time += duration.count();
			return make_tuple(add, remove);
		}
		// both are not empty, but they correspond to the same region
		if (node1->is_leaf() && node2->is_leaf()){	// case2, two leafs
			utils::fetch_and_add(&case2, 1);
			auto time_s = chrono::high_resolution_clock::now();
			auto [add, remove] = leaf_leaf_diff(node1, node2);
			auto time_t = chrono::high_resolution_clock::now();
			auto duration = chrono::duration_cast<chrono::microseconds>(time_t - time_s);
			case2_time += duration.count();
			return make_tuple(add, remove);
		}
		if (node1->is_leaf() && !node2->is_leaf()){	//	case3, a leaf and an inte
			utils::fetch_and_add(&case3, 1);
			auto P = static_cast<LeafNode*>(node1.get())->records;
			auto time_s = chrono::high_resolution_clock::now();
			auto [add, remove] = leaf_inte_diff(P, 0, P.size(), node2, b);	
			auto time_t = chrono::high_resolution_clock::now();
			auto duration = chrono::duration_cast<chrono::microseconds>(time_t - time_s);
			case3_time += duration.count();
			return make_tuple(add, remove);
		}
		if (!node1->is_leaf() && node2->is_leaf()){	//	case3, an inte and a leaf
			utils::fetch_and_add(&case3, 1);
			auto P = static_cast<LeafNode*>(node2.get())->records;
			auto time_s = chrono::high_resolution_clock::now();
			auto [remove, add] = leaf_inte_diff(P, 0, P.size(), node1, b);	//	inversed due to the second one is leaf 
			auto time_t = chrono::high_resolution_clock::now();
			auto duration = chrono::duration_cast<chrono::microseconds>(time_t - time_s);
			case3_time += duration.count();
			return make_tuple(add, remove);
		}
		//	case4, two inte
		utils::fetch_and_add(&case4, 1);
		sequence<Point*> L_add, L_remove, R_add, R_remove;
		auto cur_inte1 = static_cast<InteNode*>(node1.get());
		auto cur_inte2 = static_cast<InteNode*>(node2.get());
		auto diff_left = [&](){tie(L_add, L_remove) = diff(cur_inte1->l_son, cur_inte2->l_son, b - 1); };	// both go left
		auto diff_right = [&](){tie(R_add, R_remove) = diff(cur_inte1->r_son, cur_inte2->r_son, b - 1);	};	// both go right
		// SEQ
		diff_left();
		diff_right();
		// par_do_if(node1->get_num_points() + node2->get_num_points() >= granularity_cutoff,
		// 	diff_left,
		// 	diff_right);

		L_add.append(R_add);
		L_remove.append(R_remove);
		return make_tuple(L_add, L_remove);
	}

	auto Tree::two_version_diff(shared_ptr<BaseNode> &node1, shared_ptr<BaseNode> &node2){
		case0 = 0, case1 = 0, case2 = 0, case3 = 0, case4 = 0;
		case0_time = 0, case1_time = 0, case2_time = 0, case3_time = 0, case4_time = 0;
		auto [add, remove] = diff(node1, node2, 64);
		return make_tuple(add, remove);
	}

	//	current return conflict only
	auto Tree::merge(shared_ptr<BaseNode> &base, shared_ptr<BaseNode> &node1, shared_ptr<BaseNode> &node2){
		auto [add1, remove1] = two_version_diff(base, node1);
		auto [add2, remove2] = two_version_diff(base, node2); 

		for (auto &pt: add1) {cout << pt->id << ":" << pt->x << "," << pt->y << " ";} cout << endl;
		for (auto &pt: remove1) {cout << pt->id << ":" << pt->x << "," << pt->y << " ";} cout << endl;

		auto [insert1, delete1, update1] = filter_diff_results(add1, remove1);
		cout << insert1.size() << " " << delete1.size() << " " << update1.size() << endl;
		for (auto &pt: update1) {cout << pt->id << ":" << pt->x << "," << pt->y << " ";} cout << endl;
		auto [insert2, delete2, update2] = filter_diff_results(add2, remove2);
		cout << insert2.size() << " " << delete2.size() << " " << update2.size() << endl;
		for (auto &pt: update2) {cout << pt->id << ":" << pt->x << "," << pt->y << " ";} cout << endl;

		auto [no_conflict_insert, conflict_insert] = merge_by_id_with_conflict(insert1, insert2);	// insert conflict
		auto [no_conflict_update, conflict_update] = merge_by_id_with_conflict(update1, update2);	// update conflict
		//	delate = delete + update
		auto [no_conflict_delate, conflict_delate] = merge_by_id_with_conflict(delete1, update2);	//	delete and update conflict
		auto [no_conflict_delate1, conflict_delate1] = merge_by_id_with_conflict(delete2, update1);	//	delete and update conflict, another case
		no_conflict_delate.append(no_conflict_delate1);
		conflict_delate.append(conflict_delate1);

		return make_tuple(conflict_insert, conflict_update, conflict_delate);
	}

	Tree::Tree(size_t _leaf_sz){
		leaf_size = _leaf_sz;
	}

	void Tree::clear(){
		root.reset();
	}

	void Tree::delete_merge_nodes(shared_ptr<BaseNode> &L, shared_ptr<BaseNode> &R, InteNode* cur_node){
		// deal with MBR, covered points of parent
		if (L == nullptr) {
			cur_node->num_pts = R->get_num_points();
		}
		else if (R == nullptr){
			cur_node->num_pts = L->get_num_points();
		} 
		else {
			cur_node->num_pts = L->get_num_points() + 
								R->get_num_points();
		}
		cur_node->l_son = move(L);
		cur_node->r_son = move(R);
	}

	void Tree::merge_nodes(shared_ptr<BaseNode> &L, shared_ptr<BaseNode> &R, shared_ptr<InteNode> &cur_node){
		// deal with MBR, covered points of parent
		if (L == nullptr) {
			cur_node->num_pts = R->get_num_points();
		}
		else if (R == nullptr){
			cur_node->num_pts = L->get_num_points();
		} 
		else {
			cur_node->num_pts = L->get_num_points() + 
								R->get_num_points();
		}
		cur_node->l_son = move(L);
		cur_node->r_son = move(R);
	}

	shared_ptr<InteNode> Tree::create_internal(shared_ptr<BaseNode> &L, shared_ptr<BaseNode> &R){
		shared_ptr<InteNode> cur_node(new InteNode());
		#ifdef CNT_NODE
		utils::fetch_and_add(&inte_node_cnt, 1);
		#endif
		// augmented changes happen here
		merge_nodes(L, R, cur_node);

		return cur_node;		
	}

	shared_ptr<LeafNode> Tree::create_leaf(sequence<Point*> &P, size_t l, size_t r, size_t b){
		// auto cur_records = parlay::make_slice(&P[l], &P[r]);
		// shared_ptr<LeafNode> cur_node(new LeafNode(cur_records));
		#ifdef CNT_NODE
		utils::fetch_and_add(&leaf_node_cnt, 1);
		#endif
		shared_ptr<LeafNode> cur_node(new LeafNode(P, l, r));
		// auto cur_records = parlay::make_slice(&P[l], &P[r]);
		// shared_ptr<LeafNode> cur_node(new LeafNode(cur_records));
		return cur_node;
	}

	shared_ptr<BaseNode> Tree::build(sequence<Point*> &P, size_t l, size_t r, size_t b){
		if (!b || (r - l <= leaf_size)){
			return create_leaf(P, l, r, b);
		}
		auto splitter = split_by_bit(P, l, r, b);
		shared_ptr<BaseNode> L = nullptr;
		shared_ptr<BaseNode> R = nullptr;
		if (l < splitter) L = build(P, l, splitter, b - 1);
		if (splitter < r) R = build(P, splitter, r, b - 1);
		return create_internal(L, R);
	}

	void Tree::build(sequence<Point*> &P){
		if (!P.size()) return;
		root = build(P, 0, P.size(), 64);
	}

	void Tree::batch_insert_sorted_node(shared_ptr<BaseNode> &x, sequence<Point*> &P, size_t l, size_t r, size_t b){
		if (x == nullptr){
			x = build(P, l, r, b);
			return;
		}
		if (x->is_leaf()){
			auto cur_leaf = static_cast<LeafNode*>(x.get());
			auto cur_records = P.substr(l, r - l);
			if (!b || cur_leaf->records.size() + cur_records.size() <= leaf_size){	// current leaf is not full
				cur_leaf->records = morton_merge(cur_leaf->records, cur_records);
				return;
			}
			else{
				auto new_points = morton_merge(cur_leaf->records, cur_records);
				x = build(new_points, 0, new_points.size(), b);
				return;
			}
		}
		auto splitter = split_by_bit(P, l, r, b);
		auto cur_inte = static_cast<InteNode*>(x.get());
		if (l < splitter){
			batch_insert_sorted_node(cur_inte->l_son, P, l, splitter, b - 1);
		}
		if (splitter < r){
			batch_insert_sorted_node(cur_inte->r_son, P, splitter, r, b - 1);
		}
		
		delete_merge_nodes(cur_inte->l_son, cur_inte->r_son, cur_inte);
	}

	void Tree::batch_insert_sorted(sequence<Point*> &P){
		if (!P.size()) return;
		if (root == nullptr) build(P);
		else batch_insert_sorted_node(root, P, 0, P.size(), 64);
	}

	shared_ptr<BaseNode> Tree::node_copy(shared_ptr<BaseNode> &x){
		if (!x) return nullptr;
		if (x->is_leaf()){
			// cout << "copying leaf node" << endl;
			#ifdef CNT_NODE
			utils::fetch_and_add(&leaf_node_cnt, 1);	
			#endif
			auto cur_leaf = static_cast<LeafNode*>(x.get());
			shared_ptr<LeafNode> new_leaf(new LeafNode(cur_leaf->records));
			return new_leaf;
		}
		else{
			#ifdef CNT_NODE
			utils::fetch_and_add(&inte_node_cnt, 1);	
			#endif
			// cout << "copying inte node" << endl;
			auto cur_inte = static_cast<InteNode*>(x.get());
			shared_ptr<InteNode> new_inte(new InteNode());
			new_inte->num_pts = cur_inte->num_pts;
			new_inte->l_son = cur_inte->l_son;
			new_inte->r_son = cur_inte->r_son;
			return new_inte;
		}
	}

	void Tree::multi_version_batch_insert_sorted_node(shared_ptr<BaseNode> &x, sequence<Point*> &P, size_t l, size_t r, size_t b){
		// cout << "current is leaf? " << x->is_leaf() << ", ";
		// cout << "insert " << r - l << " points, l = " << l << ", r = " << r << ", b = " << b << endl;
		if (x == nullptr){
			x = build(P, l, r, b);
			return;
		}
		if (x->is_leaf()){
			auto cur_leaf = static_cast<LeafNode*>(x.get());
			auto cur_records = P.substr(l, r - l);
			if (!b || cur_leaf->records.size() + cur_records.size() <= leaf_size){	// current leaf is not full
				cur_leaf->records = morton_merge(cur_leaf->records, cur_records);
				return;
			}
			else{
				// auto new_points = parlay::merge(cur_leaf->records, parlay::make_slice(&P[l], &P[r]), less);
				// auto new_points = parlay::merge(cur_leaf->records, cur_records);
				auto new_points = morton_merge(cur_leaf->records, cur_records);
				x = build(new_points, 0, new_points.size(), b);
				return;
			}
		}
		auto splitter = split_by_bit(P, l, r, b);
		auto cur_inte = static_cast<InteNode*>(x.get());
		auto new_l_son(cur_inte->l_son);
		auto new_r_son(cur_inte->r_son);

		if (l < splitter){ 
			new_l_son = node_copy(cur_inte->l_son);
			multi_version_batch_insert_sorted_node(new_l_son, P, l, splitter, b - 1);
		};
		if (splitter < r){ 
			new_r_son = node_copy(cur_inte->r_son);
			multi_version_batch_insert_sorted_node(new_r_son, P, splitter, r, b - 1); 
		};
			
		delete_merge_nodes(new_l_son, new_r_son, cur_inte);
	}

	shared_ptr<BaseNode> Tree::multi_version_batch_insert_sorted(sequence<Point*> &P, shared_ptr<BaseNode> &old_root){
		auto new_root = node_copy(old_root);
		// multi_version_roots.emplace_back(new_root);
		multi_version_batch_insert_sorted_node(new_root, P, 0, P.size(), 64);
		return new_root;
	}

	// TODO: deletion has bugs for hash table. Deletion tree is not identical as build/insertion/ 
	template<class HashTable>
	void Tree::batch_delete_sorted_node(shared_ptr<BaseNode> &x, sequence<Point*> &P, size_t l, size_t r, size_t b, HashTable &ht){
		if (x->is_leaf()){
			auto cur_leaf = static_cast<LeafNode*>(x.get());
			parlay::sequence<bool> to_remove(cur_leaf->records.size(), true);
			auto not_in = [&](auto x){
				// return ht.find(x->id) == static_cast<int>(x->id);	// id is unsigned, but hash table is signed (needs empty -1).
				for (size_t i = l; i < r; i++){
					if (x->id == P[i]->id) return false;
				}
				return true;
			};
			cur_leaf->records = parlay::filter(cur_leaf->records, not_in);
			if (!cur_leaf->records.size()) x.reset();
			return;
		}
		auto splitter = split_by_bit(P, l, r, b);
		auto cur_inte = static_cast<InteNode*>(x.get());
		#ifdef SEQ
		if (l < splitter){
			batch_delete_sorted_node(cur_inte->l_son, P, l, splitter, b - 1, ht);
		}
		if (splitter < r){
			batch_delete_sorted_node(cur_inte->r_son, P, splitter, r, b - 1, ht);
		}
		#else
		auto delete_left = [&](){ if (l < splitter){ batch_delete_sorted_node(cur_inte->l_son, P, l, splitter, b - 1, ht); }; };
		auto delete_right = [&](){ if (splitter < r){ batch_delete_sorted_node(cur_inte->r_son, P, splitter, r, b - 1, ht); }; };
		par_do_if(r - l >= granularity_cutoff,
			delete_left,
			delete_right); 
		#endif

		auto less = [&](auto lhs, auto rhs){
			auto msd = 0;
			if (geobase::less_msb(static_cast<unsigned int>(lhs->x) ^ static_cast<unsigned int>(rhs->x), static_cast<unsigned int>(lhs->y) ^ static_cast<unsigned int>(rhs->y))) 
				msd = 1;
			return !msd ? lhs->x < rhs->x : lhs->y < rhs->y;
		};

		if (cur_inte->l_son == nullptr && cur_inte->r_son == nullptr) x.reset();
		else{
			// delete_merge_nodes(cur_inte->l_son, cur_inte->r_son, cur_inte);
			// if (cur_inte->l_son == nullptr) x = move(cur_inte->r_son);
			// else if (cur_inte->r_son == nullptr) x = move(cur_inte->l_son);
			if (cur_inte->l_son == nullptr) {
				if (cur_inte->r_son->get_num_points() <= leaf_size) x = move(cur_inte->r_son);
			}
			else if (cur_inte->r_son == nullptr) {
				if (cur_inte->l_son->get_num_points() <= leaf_size) x = move(cur_inte->l_son);
			}
			else if (cur_inte->l_son->get_num_points() + cur_inte->r_son->get_num_points() <= leaf_size){
			// if (cur_inte->l_son->get_num_points() + cur_inte->r_son->get_num_points() <= leaf_size){
				auto L = static_cast<LeafNode*>(cur_inte->l_son.get());
				auto R = static_cast<LeafNode*>(cur_inte->r_son.get());
				auto cur_records = parlay::merge(L->records, R->records, less);
				x = create_leaf(cur_records, 0, cur_records.size(), 0);
			}
			else{
				// cur_inte->num_pts = new_l_son->get_num_points() + new_r_son->get_num_points();
				// cur_inte->l_son = move(new_l_son);
				// cur_inte->r_son = move(new_r_son);
				delete_merge_nodes(cur_inte->l_son, cur_inte->r_son, cur_inte);
			}
		}
	}

	void Tree::batch_delete_sorted_node(shared_ptr<BaseNode> &x, sequence<Point*> &P, size_t l, size_t r, size_t b){
		if (x == nullptr){
			return;
		}
		auto not_in = [&](auto x){
			for (size_t i = l; i < r; i++){
				if (x->id == P[i]->id) return false;
			}
			return true;
		};
		if (x->is_leaf()){
			auto cur_leaf = static_cast<LeafNode*>(x.get());
			cur_leaf->records = parlay::filter(cur_leaf->records, not_in);
			return;
		}
		auto splitter = split_by_bit(P, l, r, b);
		auto cur_inte = static_cast<InteNode*>(x.get());
		#ifdef SEQ
		if (l < splitter){
			batch_delete_sorted_node(cur_inte->l_son, P, l, splitter, b - 1);
		}
		if (splitter < r){
			batch_delete_sorted_node(cur_inte->r_son, P, splitter, r, b - 1);
		}
		#else
		auto delete_left = [&](){ if (l < splitter){ batch_delete_sorted_node(cur_inte->l_son, P, l, splitter, b - 1); }; };
		auto delete_right = [&](){ if (splitter < r){ batch_delete_sorted_node(cur_inte->r_son, P, splitter, r, b - 1); }; };
		par_do_if(r - l >= granularity_cutoff,
			delete_left,
			delete_right); 
		#endif

		delete_merge_nodes(cur_inte->l_son, cur_inte->r_son, cur_inte);
	}

	void Tree::multi_version_batch_delete_sorted_node(shared_ptr<BaseNode> &x, sequence<Point*> &P, size_t l, size_t r, size_t b){
		if (x->is_leaf()){
			auto cur_leaf = static_cast<LeafNode*>(x.get());
			parlay::sequence<bool> to_remove(cur_leaf->records.size(), true);
			auto not_in = [&](auto x){
				// return ht.find(x->id) == static_cast<int>(x->id);	// id is unsigned, but hash table is signed (needs empty -1).
				for (size_t i = l; i < r; i++){
					if (x->id == P[i]->id) return false;
				}
				return true;
			};
			cur_leaf->records = parlay::filter(cur_leaf->records, not_in);
			if (!cur_leaf->records.size()) x.reset();
			return;
		}
		auto splitter = split_by_bit(P, l, r, b);
		auto cur_inte = static_cast<InteNode*>(x.get());
		auto new_l_son(cur_inte->l_son);
		auto new_r_son(cur_inte->r_son);

		auto delete_left = [&](){ 
			if (l < splitter){ 
				new_l_son = node_copy(cur_inte->l_son);
				multi_version_batch_delete_sorted_node(new_l_son, P, l, splitter, b - 1); 
			}; 
		};
		auto delete_right = [&](){
			if (splitter < r){ 
				new_r_son = node_copy(cur_inte->r_son);
				multi_version_batch_delete_sorted_node(new_r_son, P, splitter, r, b - 1); 
			}; 
		};

		par_do_if(r - l >= granularity_cutoff,
			delete_left,
			delete_right); 

		auto less = [&](auto lhs, auto rhs){
			auto msd = 0;
			if (geobase::less_msb(static_cast<unsigned int>(lhs->x) ^ static_cast<unsigned int>(rhs->x), static_cast<unsigned int>(lhs->y) ^ static_cast<unsigned int>(rhs->y))) 
				msd = 1;
			return !msd ? lhs->x < rhs->x : lhs->y < rhs->y;
		};

		cur_inte->l_son = new_l_son;
		cur_inte->r_son = new_r_son;
		if (!new_l_son && !new_r_son) x.reset();
		else{
			if (!new_l_son) {
				if (new_r_son->get_num_points() <= leaf_size) x = move(new_r_son);
			}
			else if (!new_r_son) {
				if (new_l_son->get_num_points() <= leaf_size) x = move(new_l_son);
			}
			else if (new_l_son->get_num_points() + new_r_son->get_num_points() <= leaf_size){
			// if (cur_inte->l_son->get_num_points() + cur_inte->r_son->get_num_points() <= leaf_size){
				auto L = static_cast<LeafNode*>(new_l_son.get());
				auto R = static_cast<LeafNode*>(new_r_son.get());
				auto cur_records = parlay::merge(L->records, R->records, less);
				x = create_leaf(cur_records, 0, cur_records.size(), 0);
			}
			else{
				// cur_inte->num_pts = new_l_son->get_num_points() + new_r_son->get_num_points();
				// cur_inte->l_son = move(new_l_son);
				// cur_inte->r_son = move(new_r_son);
				delete_merge_nodes(new_l_son, new_r_son, cur_inte);
			}
		}
	}

	shared_ptr<BaseNode> Tree::multi_version_batch_delete_sorted(sequence<Point*> &P, shared_ptr<BaseNode> &old_root){
		auto new_root = node_copy(old_root);
		multi_version_roots.emplace_back(new_root);
		multi_version_batch_delete_sorted_node(new_root, P, 0, P.size(), 64);
		return new_root;
	}

	void Tree::batch_delete_sorted(sequence<Point*> &P){
		if (!P.size() || root == nullptr) return;
		else batch_delete_sorted_node(root, P, 0, P.size(), 64);
	}

	template <class HashTable>
	void Tree::batch_delete_sorted(sequence<Point*> &P, HashTable &ht){
		if (!P.size() || root == nullptr) return;
		else batch_delete_sorted_node(root, P, 0, P.size(), 64, ht);
	}

	void Tree::tree_compress(shared_ptr<BaseNode> &x){
		if (x->is_leaf()){
			return;
		}
		auto cur_inte = static_cast<InteNode*>(x.get());
		if (cur_inte->l_son != nullptr) tree_compress(cur_inte->l_son);
		if (cur_inte->r_son != nullptr) tree_compress(cur_inte->r_son);

		if (cur_inte->l_son == nullptr && cur_inte->r_son == nullptr) x.reset();
		else if (cur_inte->l_son == nullptr) x = move(cur_inte->r_son);
		else if (cur_inte->r_son == nullptr) x = move(cur_inte->l_son);
	}

	//	store all pointers of a root in a set
	void Tree::tree_node_count(shared_ptr<BaseNode> &x, set<uintptr_t> &s){
		if (!x) return;
		if (x->is_leaf()){
			auto cur_leaf = static_cast<LeafNode*>(x.get());
			auto node_addr = reinterpret_cast<uintptr_t>(cur_leaf); 
			if (s.find(node_addr) != s.end()) return;
			s.insert(node_addr);
			return;
		}
		auto cur_inte = static_cast<InteNode*>(x.get());
		auto node_addr = reinterpret_cast<uintptr_t>(cur_inte);
		if (s.find(node_addr) != s.end()) return;
		s.insert(node_addr);
		tree_node_count(cur_inte->l_son, s);
		tree_node_count(cur_inte->r_son, s);
	}

	size_t Tree::node_count(){
		set<uintptr_t > s;
		// cout << "seq set size = " << multi_version_roots.size() << endl;
		for (auto &root: multi_version_roots){
			tree_node_count(root, s);
		}
		return s.size();
	}


	int Tree::tree_hash(shared_ptr<BaseNode> &x){
		node_cnt++;
		if (x->is_leaf()){
			leaf_cnt++;
			int ret = 0;
			auto cur_leaf = static_cast<LeafNode*>(x.get());
			// cout << "leaf size = " << cur_leaf->records.size() << endl;
			record_cnt += cur_leaf->records.size();
			for (auto &pt: cur_leaf->records){
				int x_val = static_cast<int>(pt->x);
				int y_val = static_cast<int>(pt->y);
				ret ^= parlay::hash64(x_val ^ y_val);
			}
			return parlay::hash64(ret);
		}
		int L = 0, R = 0;
		auto cur_inte = static_cast<InteNode*>(x.get());
		if (cur_inte->l_son != nullptr) L = tree_hash(cur_inte->l_son);
		if (cur_inte->r_son != nullptr) R = tree_hash(cur_inte->r_son);
		return parlay::hash64(L ^ R);
	}

	void Tree::in_order_traverse(shared_ptr<BaseNode> &x){
		#ifdef USE_MBR
		print_mbr(x->mbr);
		#endif
		if (x->is_leaf()){
			auto cur_leaf = static_cast<LeafNode*>(x.get());
			if (cur_leaf->get_num_points() > leaf_size){
			// 	print_mbr(x->mbr);
			// 	cout << cur_leaf->get_num_points() << endl;
			}
			return;
		}
		auto cur_inte = static_cast<InteNode*>(x.get());
		if (cur_inte->l_son != nullptr) in_order_traverse(cur_inte->l_son);
		if (cur_inte->r_son != nullptr) in_order_traverse(cur_inte->r_son);
	}
	
	void Tree::in_order_traverse(shared_ptr<BaseNode> &x, unsigned x_prefix, unsigned y_prefix, size_t b, bool x_splitter){
		if (x->is_leaf()){
			auto cur_leaf = static_cast<LeafNode*>(x.get());
			cur_leaf->print_records();
			return;
		}
		auto cur_inte = static_cast<InteNode*>(x.get());
		if (x_splitter){
			if (cur_inte->l_son != nullptr) {
				in_order_traverse(cur_inte->l_son, x_prefix, y_prefix, b, !x_splitter);
			}
			if (cur_inte->r_son != nullptr){
				in_order_traverse(cur_inte->r_son, x_prefix + (1 << (b - 1)), y_prefix, b, !x_splitter);
			}
		}
		else{
			if (cur_inte->l_son != nullptr){
				in_order_traverse(cur_inte->l_son, x_prefix, y_prefix, b - 1, !x_splitter);
			}
			if (cur_inte->r_son != nullptr){
 				in_order_traverse(cur_inte->r_son, x_prefix, y_prefix + (1 << (b - 1)), b - 1, !x_splitter);
			}
		}
	}
	
// range count using mbr
#ifdef USE_MBR
	size_t Tree::range_count_node(shared_ptr<BaseNode> &x, Bounding_Box &query_mbr){
		// visited_inte++;
		int flag = mbr_mbr_relation(x->mbr, query_mbr);
		if (flag < 0) return 0;
		if (flag > 0) return x->get_num_points();
		if (x->is_leaf()){	// we have to scan the leaf to report the number of points;
			// visited_leaf++;
			size_t ret = 0;
			auto cur_leaf = static_cast<LeafNode*>(x.get());
			for (auto &p: cur_leaf->records){
				if (point_in_mbr(*p, query_mbr)){
					ret += 1;
				}
			}
			return ret;
		}
		else{
			auto cur_inte = static_cast<InteNode*>(x.get());
			size_t ret_L = 0, ret_R = 0;
			if (cur_inte->l_son != nullptr) ret_L = range_count_node(cur_inte->l_son, query_mbr);
			if (cur_inte->r_son != nullptr) ret_R = range_count_node(cur_inte->r_son, query_mbr);
			return ret_L + ret_R;
		}
		return -1; // unexpected error happens if the code runs to here.
	}

	size_t Tree::range_count(Bounding_Box &query_mbr){
		// visited_inte = 0, visited_leaf = 0;
		size_t ret = range_count_node(root, query_mbr);		
		// visited_inte -= visited_leaf;
		// cout << "visited nodes (interior/leaf): " << visited_inte << ", " << visited_leaf << endl;
		return ret;
	}

	size_t Tree::range_count_node(shared_ptr<BaseNode> &x, Bounding_Box &query_mbr, Bounding_Box &cur_mbr, bool x_splitter){
		int flag = mbr_mbr_relation(x->mbr, query_mbr);
		if (flag < 0) return 0;
		if (flag > 0) return x->get_num_points();
		if (x->is_leaf()){	// we have to scan the leaf to report the number of points;
			size_t ret = 0;
			auto cur_leaf = static_cast<LeafNode*>(x.get());
			for (auto &p: cur_leaf->records){
				if (point_in_mbr(*p, query_mbr)){
					ret += 1;
				}
			}
			return ret;
		}
		else{
			auto cur_inte = static_cast<InteNode*>(x.get());
			size_t ret_L = 0, ret_R = 0;
			if (cur_inte->l_son != nullptr) {
				auto L_box = cur_mbr;
				if (x_splitter) L_box.second.x = cur_inte->l_son->mbr.second.x;
				else L_box.second.y = cur_inte->l_son->mbr.second.y;
				ret_L = range_count_node(cur_inte->l_son, query_mbr, L_box, !x_splitter);
			}
			if (cur_inte->r_son != nullptr) {
				auto R_box = cur_mbr;
				if (x_splitter) R_box.first.x = cur_inte->r_son->mbr.first.x;
				else R_box.first.y = cur_inte->r_son->mbr.first.y;
				ret_R = range_count_node(cur_inte->r_son, query_mbr, R_box, !x_splitter);
			}
			return ret_L + ret_R;
		}
		return -1; // unexpected error happens if the code runs to here.
	}
#endif

	size_t Tree::range_count_node(shared_ptr<BaseNode> &x, Bounding_Box &query_mbr, Bounding_Box &cur_mbr, FT x_prefix, FT y_prefix, size_t b, bool x_splitter){
		int flag = mbr_mbr_relation(cur_mbr, query_mbr);
		if (flag < 0) return 0;
		if (flag > 0) return x->get_num_points();
		if (x->is_leaf()){	// we have to scan the leaf to report the number of points;
			size_t ret = 0;
			auto cur_leaf = static_cast<LeafNode*>(x.get());
			for (auto &p: cur_leaf->records){
				if (point_in_mbr(*p, query_mbr)){
					ret += 1;
				}
			}
			return ret;
		}
		else{
			auto cur_inte = static_cast<InteNode*>(x.get());
			FT splitter = 1.0 * (1u << (b - 1));
			size_t ret_L = 0, ret_R = 0;
			if (cur_inte->l_son != nullptr) {
				auto L_box = cur_mbr;
				if (x_splitter){
					L_box.second.x = min(x_prefix + splitter - FT_EPS, L_box.second.x);
					ret_L = range_count_node(cur_inte->l_son, query_mbr, L_box, x_prefix, y_prefix, b, !x_splitter);
				}
				else{
					L_box.second.y = min(y_prefix + splitter - FT_EPS, L_box.second.y);
					ret_L = range_count_node(cur_inte->l_son, query_mbr, L_box, x_prefix, y_prefix, b - 1, !x_splitter);
				}
			}
			if (cur_inte->r_son != nullptr) {
				auto R_box = cur_mbr;
				if (x_splitter){
					R_box.first.x = max(x_prefix + splitter, R_box.first.x);
					ret_R = range_count_node(cur_inte->r_son, query_mbr, R_box, x_prefix + splitter, y_prefix, b, !x_splitter);
				}
				else{
					R_box.first.y = max(y_prefix + splitter, R_box.first.y);
					ret_R = range_count_node(cur_inte->r_son, query_mbr, R_box, x_prefix, y_prefix + splitter, b - 1, !x_splitter);
				}
			}
			return ret_L + ret_R;
		}
		return -1; // unexpected error happens if the code runs to here.
	}

	size_t Tree::range_count(Bounding_Box &query_mbr, Bounding_Box &cur_mbr){
		// size_t ret = range_count_node(root, query_mbr, cur_mbr, true);		
		size_t ret = range_count_node(root, query_mbr, cur_mbr, 0.0, 0.0, 32, true);		
		return ret;
	}

	#ifdef USE_MBR
	template <class Out> 
	void Tree::range_report_node(shared_ptr<BaseNode> &x, Bounding_Box &query_mbr, size_t &cnt, Out &out){
		if (x->is_leaf()){
			auto cur_leaf = static_cast<LeafNode*>(x.get());
			for (auto &p: cur_leaf->records){
				if (point_in_mbr(*p, query_mbr)){
					out[cnt++] = p;
				}
			}
			return;
		}
		auto cur_inte = static_cast<InteNode*>(x.get());
		if (cur_inte->l_son != nullptr){
			if (mbr_mbr_relation(cur_inte->l_son->mbr, query_mbr) >= 0) range_report_node(cur_inte->l_son, query_mbr, cnt, out);
		}
		if (cur_inte->r_son != nullptr){
 			if (mbr_mbr_relation(cur_inte->r_son->mbr, query_mbr) >= 0) range_report_node(cur_inte->r_son, query_mbr, cnt, out);
		}
	}

	template <class Out>
	void Tree::range_report(Bounding_Box &query_mbr, size_t &cnt, Out &out){
		range_report_node(root, query_mbr, cnt, out);
	}
	#endif

	template <class Out> 
	void Tree::range_report_node(shared_ptr<BaseNode> &x, Bounding_Box &query_mbr, Bounding_Box &cur_mbr, FT x_prefix, FT y_prefix, size_t b, bool x_splitter, size_t &cnt, Out &out){
		// print_mbr(cur_mbr);
		// print_mbr(query_mbr);
		if (x->is_leaf()){
			auto cur_leaf = static_cast<LeafNode*>(x.get());
			for (auto &p: cur_leaf->records){
				if (point_in_mbr(*p, query_mbr)){
					out[cnt++] = p;
				}
			}
			return;
		}
		auto cur_inte = static_cast<InteNode*>(x.get());
		FT splitter = 1.0 * (1u << (b - 1));
		if (cur_inte->l_son != nullptr){
			auto L_box = cur_mbr;
			if (x_splitter){
				L_box.second.x = min(x_prefix + splitter - FT_EPS, L_box.second.x);
				if (mbr_mbr_relation(L_box, query_mbr) >= 0) range_report_node(cur_inte->l_son, query_mbr, L_box, x_prefix, y_prefix, b, !x_splitter, cnt, out);
			}
			else{
				L_box.second.y = min(y_prefix + splitter - FT_EPS, L_box.second.y);
				if (mbr_mbr_relation(L_box, query_mbr) >= 0) range_report_node(cur_inte->l_son, query_mbr, L_box, x_prefix, y_prefix, b - 1, !x_splitter, cnt, out);
			}
		}
		if (cur_inte->r_son != nullptr){
			auto R_box = cur_mbr;
			if (x_splitter){
				R_box.first.x = max(x_prefix + splitter, R_box.first.x);
 				if (mbr_mbr_relation(R_box, query_mbr) >= 0) range_report_node(cur_inte->r_son, query_mbr, R_box, x_prefix + splitter, y_prefix, b, !x_splitter, cnt, out);
			}
			else{
				R_box.first.y = max(y_prefix + splitter, R_box.first.y);
 				if (mbr_mbr_relation(R_box, query_mbr) >= 0) range_report_node(cur_inte->r_son, query_mbr, R_box, x_prefix, y_prefix + splitter, b - 1, !x_splitter, cnt, out);
			}
		}
	}

	template <class Out>
	void Tree::range_report(Bounding_Box &query_mbr, Bounding_Box &cur_mbr, size_t &cnt, Out &out){
		range_report_node(root, query_mbr, cur_mbr, 0.0, 0.0, 32, true, cnt, out);
	}

	#ifdef USE_MBR
	template <class T>
	void Tree::knn_report_node(shared_ptr<BaseNode> &x, size_t &k, Point* query_point, T &nn_res){
		if (x->is_leaf()){
			auto cur_leaf = static_cast<LeafNode*>(x.get());
			for (auto &p: cur_leaf->records){
				auto cur_sqrdis = point_point_sqrdis(p, query_point);
				if (nn_res.size() < k){
					nn_res.push({p, cur_sqrdis});
				}
				else if (cur_sqrdis < nn_res.top().second){
					nn_res.pop();
					nn_res.push({p, cur_sqrdis});
				}
			}
			return;
		}
		auto cur_inte = static_cast<InteNode*>(x.get());
		auto l_son_sqrdis = FT_INF_MAX, r_son_sqrdis = FT_INF_MAX;
		if (cur_inte->l_son != nullptr) l_son_sqrdis = point_mbr_sqrdis(query_point, cur_inte->l_son->mbr);
		if (cur_inte->r_son != nullptr) r_son_sqrdis = point_mbr_sqrdis(query_point, cur_inte->r_son->mbr);
		if (l_son_sqrdis <= r_son_sqrdis){ // first go left
			if (nn_res.size() < k || l_son_sqrdis < nn_res.top().second){
				knn_report_node(cur_inte->l_son, k, query_point, nn_res);
			}
			if (nn_res.size() < k || r_son_sqrdis < nn_res.top().second){
				knn_report_node(cur_inte->r_son, k, query_point, nn_res);
			}
		}
		else{	// first go right
			if (nn_res.size() < k || r_son_sqrdis < nn_res.top().second){
				knn_report_node(cur_inte->r_son, k, query_point, nn_res);
			}
			if (nn_res.size() < k || l_son_sqrdis < nn_res.top().second){
				knn_report_node(cur_inte->l_son, k, query_point, nn_res);
			}
		}
		return;
	}

	auto Tree::knn_report(size_t &k, Point* query_point){
		priority_queue<nn_pair, vector<nn_pair>, nn_pair_cmp> nn_res;
		knn_report_node(root, k, query_point, nn_res);
		return nn_res;
	}
	#endif

	auto Tree::knn_report(size_t &k, Point* query_point, Bounding_Box &cur_mbr){
		priority_queue<nn_pair, vector<nn_pair>, nn_pair_cmp> nn_res;
		knn_report_node(root, k, query_point, cur_mbr, 0.0, 0.0, 32, true, nn_res);
		return nn_res;
	}

	template <class T>
	void Tree::knn_report_node(shared_ptr<BaseNode> &x, size_t &k, Point* query_point, Bounding_Box &cur_mbr, FT x_prefix, FT y_prefix, size_t b, bool x_splitter, T &nn_res){
		if (x->is_leaf()){
			auto cur_leaf = static_cast<LeafNode*>(x.get());
			for (auto &p: cur_leaf->records){
				auto cur_sqrdis = point_point_sqrdis(p, query_point);
				if (nn_res.size() < k){
					nn_res.push({p, cur_sqrdis});
				}
				else if (cur_sqrdis < nn_res.top().second){
					nn_res.pop();
					nn_res.push({p, cur_sqrdis});
				}
			}
			return;
		}
		auto cur_inte = static_cast<InteNode*>(x.get());
		auto l_son_sqrdis = FT_INF_MAX, r_son_sqrdis = FT_INF_MAX;
		FT splitter = 1.0 * (1u << (b - 1));
		auto L_box = cur_mbr, R_box = cur_mbr;
		if (cur_inte->l_son != nullptr) {
			if (x_splitter){
				L_box.second.x = min(x_prefix + splitter - FT_EPS, L_box.second.x);
			}
			else{
				L_box.second.y = min(y_prefix + splitter - FT_EPS, L_box.second.y);
			}
			l_son_sqrdis = point_mbr_sqrdis(query_point, L_box);
		}
		if (cur_inte->r_son != nullptr) {
			if (x_splitter){
				R_box.first.x = max(x_prefix + splitter, R_box.first.x);
			}
			else{
				R_box.first.y = max(y_prefix + splitter, R_box.first.y);
			}
			r_son_sqrdis = point_mbr_sqrdis(query_point, R_box);
		}
		// auto go_left = [&](){
		// 	if (nn_res.size() < k || l_son_sqrdis < nn_res.top().second){
		// 		if (x_splitter) knn_report_node(cur_inte->l_son, k, query_point, L_box, x_prefix, y_prefix, b, !x_splitter, nn_res);
		// 		else knn_report_node(cur_inte->l_son, k, query_point, L_box, x_prefix, y_prefix, b - 1, !x_splitter, nn_res);
		// 	}
		// };
		// auto go_right = [&](){
		// 	if (nn_res.size() < k || r_son_sqrdis < nn_res.top().second){
		// 		if (x_splitter) knn_report_node(cur_inte->r_son, k, query_point, R_box, x_prefix + splitter, y_prefix, b, !x_splitter, nn_res);
		// 		else knn_report_node(cur_inte->r_son, k, query_point, R_box, x_prefix, y_prefix + splitter, b - 1, !x_splitter, nn_res);
		// 	}
		// };
		if (l_son_sqrdis <= r_son_sqrdis){ // first go left
			// go_left();
			// go_right();
			if (nn_res.size() < k || l_son_sqrdis < nn_res.top().second){
				if (x_splitter) knn_report_node(cur_inte->l_son, k, query_point, L_box, x_prefix, y_prefix, b, !x_splitter, nn_res);
				else knn_report_node(cur_inte->l_son, k, query_point, L_box, x_prefix, y_prefix, b - 1, !x_splitter, nn_res);
			}
			if (nn_res.size() < k || r_son_sqrdis < nn_res.top().second){
				if (x_splitter) knn_report_node(cur_inte->r_son, k, query_point, R_box, x_prefix + splitter, y_prefix, b, !x_splitter, nn_res);
				else knn_report_node(cur_inte->r_son, k, query_point, R_box, x_prefix, y_prefix + splitter, b - 1, !x_splitter, nn_res);
			}
		}
		else{	// first go right
			// go_right();
			// go_left();
			if (nn_res.size() < k || r_son_sqrdis < nn_res.top().second){
				if (x_splitter) knn_report_node(cur_inte->r_son, k, query_point, R_box, x_prefix + splitter, y_prefix, b, !x_splitter, nn_res);
				else knn_report_node(cur_inte->r_son, k, query_point, R_box, x_prefix, y_prefix + splitter, b - 1, !x_splitter, nn_res);
				
			}
			if (nn_res.size() < k || l_son_sqrdis < nn_res.top().second){
				if (x_splitter) knn_report_node(cur_inte->l_son, k, query_point, L_box, x_prefix, y_prefix, b, !x_splitter, nn_res);
				else knn_report_node(cur_inte->l_son, k, query_point, L_box, x_prefix, y_prefix, b - 1, !x_splitter, nn_res);
			}
		}
		return;
	}

}