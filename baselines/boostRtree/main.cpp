#include<bits/stdc++.h>
#include <boost/geometry.hpp>
#include <boost/geometry/index/rtree.hpp>
#include "helper/geobase.h"
#include "helper/parse_command_line.h"
#include "helper/time_loop.h"

#define TEST

using namespace std;
namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

// Define point and R-tree value type
typedef bg::model::point<double, 2, bg::cs::cartesian> Point;
typedef pair<Point, size_t> Value; // (point, ID)
typedef bgi::rtree<Value, bgi::quadratic<32>> RTree;

geobase::Bounding_Box largest_mbr;  //  largest mbr 
size_t maxSize = 100;   //  max points in a region

/* PT is a sequence of points */
template<typename PT>
auto convert_P(PT &P){
    vector<Value> ret(P.size());
    parlay::parallel_for(0, P.size(), [&](size_t i){
        ret[i] = make_pair(Point(P[i].x, P[i].y), P[i].id);
    });
    return ret;
}

/* T is a sequence of points */
template<typename T>
auto convert_Q(T &Q){
    vector<pair<Point, Point> > ret(Q.size());
    parlay::parallel_for(0, Q.size(), [&](size_t i){
        auto min_pt = Q[i].first;
        auto max_pt = Q[i].second;
        ret[i] = make_pair(Point(min_pt.x, min_pt.y), Point(max_pt.x, max_pt.y));
    });
    return ret;
}

template <typename T>
auto diff_by_id(T &a, T &b){
    auto id_cmp = [&](auto &lhs, auto &rhs) { return lhs.second < rhs.second; };
    
    auto sorted_a = parlay::sort(a, id_cmp);
    auto sorted_b = parlay::sort(b, id_cmp);

    vector<Value> added(b.size());
    vector<Value> removed(a.size());

    size_t cnt_add = 0, cnt_remove = 0;
    size_t i = 0, j = 0;

    while (i < a.size() && j < b.size()){
        if (sorted_a[i].second == sorted_b[j].second){ // same point, check coordinate
            if (bg::equals(sorted_a[i].first, sorted_b[j].first)){ //  same coordinate, do nothing
                i++, j++;
            }
            else{
                removed[cnt_remove++] = sorted_a[i++];
                added[cnt_add++] = sorted_b[j++];
            }
        }
        else if (sorted_a[i].second < sorted_b[j].second){ // first smaller, in A not in B
            removed[cnt_remove++] = sorted_a[i++];
        }
        else{ //  second smaller, in B not in A
            added[cnt_add++] = sorted_b[j++];
        }
    }
    while (i < sorted_a.size()){
        removed[cnt_remove++] = sorted_a[i++];
    }
    while (j < sorted_b.size()){
        added[cnt_add++] = sorted_b[j++];
    }
    removed.resize(cnt_remove);
    added.resize(cnt_add);

    return make_tuple(added, removed);
}

void print_values(const std::vector<Value>& values) {
    for (const auto& val : values) {
        double x = bg::get<0>(val.first);
        double y = bg::get<1>(val.first);
        size_t id = val.second;
        cout << fixed << setprecision(6) << "ID: " << id << " at (" << x << ", " << y << ")\n";
    }
}

template<typename PT, typename RQ>
void spatial_diff_test(PT &P, RQ &range_queries, parlay::sequence<size_t> &batch_sizes, size_t &insert_ratio){
    auto P_conv = convert_P(P);
    RTree rtree(P_conv.begin(), P_conv.end());  //  original version
    auto max_batch_size = batch_sizes[batch_sizes.size() - 1];

    /* get insert, delete points */
    auto P_test = geobase::shuffle_point(P, max_batch_size);
    auto [P_insert_set, P_delete_set] = geobase::split_insert_delete(P_test, insert_ratio, P.size());
    auto Q_conv = convert_Q(range_queries);
    
    for (auto &batch_size: batch_sizes){
        cout << "[INFO] Batch Size: " << batch_size << endl;
        auto insert_num = batch_size / 10 * insert_ratio;
        auto delete_num = batch_size / 10 * (10 - insert_ratio);
        
        auto P_insert = P_insert_set.substr(0, insert_num);
        auto P_delete = P_delete_set.substr(0, delete_num);

        /* apply delete, insert */
        auto P_newver = geobase::collect_newver_point(P, P_insert, P_delete);
        auto P_newver_conv = convert_P(P_newver);
        RTree rtree2(P_newver_conv.begin(), P_newver_conv.end());

        parlay::sequence<size_t> addCnt(range_queries.size());
        parlay::sequence<size_t> removeCnt(range_queries.size());
        vector<Value> res1, res2;

        auto avg_time = time_loop(
            3, 1.0,
            [&](){},
            [&](){
                for (size_t i = 0; i < Q_conv.size(); i++){
                    res1.clear();
                    res2.clear();
                    rtree.query(bgi::intersects(bg::model::box<Point>(Q_conv[i].first, Q_conv[i].second)), std::back_inserter(res1));
                    rtree2.query(bgi::intersects(bg::model::box<Point>(Q_conv[i].first, Q_conv[i].second)), std::back_inserter(res2));
                    auto [add, remove] = diff_by_id(res1, res2);
                    addCnt[i] = add.size();
                    removeCnt[i] = remove.size();
                }
            },
            [&]{}
        );

        cout << fixed << setprecision(6) << "[Boost-Rtree] spatial-diff time (avg): " << avg_time << endl;

        #ifdef TEST
        string file_name = "output/rtree_spatial_diff-" + to_string(batch_size); 
        ofstream spatialDiffOut(file_name);
        for (size_t i = 0; i < range_queries.size(); i++){
            spatialDiffOut<< addCnt[i] << " " << removeCnt[i] << endl;
        }
        #endif
    }
}

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

	string task = cmd.getOptionValue("-t");
	string input_file = cmd.getOptionValue("-i");
	int is_real = cmd.getOptionIntValue("-real", 0);

	/* read input file */
	ifstream fin(input_file);
	parlay::sequence<geobase::Point> P;
	geobase::read_pts(P, fin, is_real);	//	change to true if id is contained.

    auto P_conv = convert_P(P); // get point format for boost
    RTree rtree(P_conv.begin(), P_conv.end());
    cout << "bulk-load finished." << endl;

    if (task == "spatial-diff"){
		string query_file = cmd.getOptionValue("-r");
		auto [cnt, queries] = geobase::read_range_query(query_file, 8, maxSize);
		// size_t batch_size = 10000;
        // queries = queries.substr(0, 1);


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
        // parlay::sequence<size_t> ratios = {1, 2, 3, 4, 5, 6, 7, 8, 9};
		// queries = queries.substr(6148, 1);
        size_t ratio = 5;
        spatial_diff_test(P, queries, batch_sizes, ratio);
	}

}

int main(int argc, char** argv){
    run(argc, argv);
    cout << "test" << endl;
    return 0;
}