#include <iostream>
#include <vector>
#include <chrono>
#include <fstream>
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/index/rtree.hpp>
#include "src/geo/point.hpp"
#include <unordered_set>

using namespace std;
namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

typedef bg::model::point<double, 2, bg::cs::cartesian> BoostPoint;
typedef pair<BoostPoint, size_t> Value;
typedef bgi::rtree<Value, bgi::quadratic<32>> RTree;

int main() {
    ifstream fin_base("dataset/bhutan_evolution/base_upto_2020.txt");
    std::vector<geobase::Point> P_base_vec;
    size_t b_tid; double b_tx, b_ty;
    while (fin_base >> b_tid >> b_tx >> b_ty) {
        if (b_tx < 0 || b_ty < 0) continue;
        P_base_vec.push_back(geobase::Point(b_tid, b_tx*1000000, b_ty*1000000));
    }
    
    vector<Value> pts;
    for(size_t i=0; i<P_base_vec.size(); i++) {
        pts.push_back(make_pair(BoostPoint(P_base_vec[i].x, P_base_vec[i].y), P_base_vec[i].id));
    }
    RTree tree(pts.begin(), pts.end());
    
    // Simulate 1.5M random deletions
    unordered_set<size_t> removed_ids;
    for(size_t i=0; i<1500000; i++) {
        removed_ids.insert(P_base_vec[(i * 17) % P_base_vec.size()].id);
    }
    
    auto t1 = chrono::high_resolution_clock::now();
    for(int i=0; i<100; i++) {
        BoostPoint bg_q(P_base_vec[i].x, P_base_vec[i].y);
        size_t k = 100;
        size_t search_k = k * 2;
        vector<Value> result;
        while (true) {
            result.clear();
            for (auto it = tree.qbegin(bgi::nearest(bg_q, search_k)); it != tree.qend(); ++it) {
                if (removed_ids.find(it->second) == removed_ids.end()) {
                    result.push_back(*it);
                    if (result.size() == k) break;
                }
            }
            if (result.size() == k || search_k >= tree.size()) break;
            search_k = min(search_k * 2, (size_t)tree.size());
        }
    }
    auto t2 = chrono::high_resolution_clock::now();
    cout << "Doubling approach took " << chrono::duration_cast<chrono::microseconds>(t2-t1).count() / 1000.0 << " ms" << endl;
    return 0;
}
