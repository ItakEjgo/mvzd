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
    
    ifstream fin_2021("dataset/bhutan_evolution/2021_merged.txt");
    while (fin_2021 >> b_tid >> b_tx >> b_ty) {
        if (b_tx < 0 || b_ty < 0) continue;
        P_base_vec.push_back(geobase::Point(b_tid, b_tx*1000000, b_ty*1000000));
    }

    vector<Value> pts;
    for(size_t i=0; i<P_base_vec.size(); i++) {
        pts.push_back(make_pair(BoostPoint(P_base_vec[i].x, P_base_vec[i].y), P_base_vec[i].id));
    }
    RTree tree(pts.begin(), pts.end());
    
    geobase::Point q(5400561020, 8.9976e+07, 2.75181e+07); // The actual query point
    
    BoostPoint bg_q(q.x, q.y);
    vector<Value> result;
    for (auto it = tree.qbegin(bgi::nearest(bg_q, 10u)); it != tree.qend(); ++it) {
        result.push_back(*it);
    }
    
    size_t hash_val = 0;
    for(const auto& v : result) {
        hash_val += v.second;
    }
    cout << "Boost Hash: " << hash_val << endl;
    return 0;
}
