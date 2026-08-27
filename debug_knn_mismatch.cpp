#include <iostream>
#include <fstream>
#include <vector>
#include "src/geo/point.hpp"
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/index/rtree.hpp>

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
    cout << "Base points: " << P_base_vec.size() << endl;
    
    ifstream fin_2021("dataset/bhutan_evolution/2021_merged.txt");
    while (fin_2021 >> b_tid >> b_tx >> b_ty) {
        if (b_tx < 0 || b_ty < 0) continue;
        P_base_vec.push_back(geobase::Point(b_tid, b_tx*1000000, b_ty*1000000));
    }
    cout << "Total points after 2021: " << P_base_vec.size() << endl;

    vector<Value> pts;
    for(size_t i=0; i<P_base_vec.size(); i++) {
        pts.push_back(make_pair(BoostPoint(P_base_vec[i].x, P_base_vec[i].y), P_base_vec[i].id));
    }
    RTree tree(pts.begin(), pts.end());
    
    // Reproduce the exact queries used in run_experiments.py
    // They are chosen randomly with seed 42 from base dataset!
    srand(42);
    parlay::sequence<geobase::Point> current_knn(100);
    for(int i=0; i<100; i++) current_knn[i] = P_base_vec[rand() % 650335]; // rand from base

    geobase::Point q = current_knn[8]; // KNN_10_8
    cout << "Query point 8: ID " << q.id << " (" << q.x << ", " << q.y << ")" << endl;
    
    BoostPoint bg_q(q.x, q.y);
    vector<Value> result;
    for (auto it = tree.qbegin(bgi::nearest(bg_q, 10u)); it != tree.qend(); ++it) {
        result.push_back(*it);
    }
    
    size_t hash_val = 0;
    cout << "Top 10 points (Boost):" << endl;
    for(const auto& v : result) {
        double dx = q.x - v.first.get<0>();
        double dy = q.y - v.first.get<1>();
        cout << "ID: " << v.second << " Dist: " << (dx*dx + dy*dy) << endl;
        hash_val += v.second;
    }
    cout << "Boost Hash: " << hash_val << endl;
    
    return 0;
}
