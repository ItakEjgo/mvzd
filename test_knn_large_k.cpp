#include <iostream>
#include <vector>
#include <chrono>
#include <fstream>
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/index/rtree.hpp>
#include "src/geo/point.hpp"

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
    
    auto t1 = chrono::high_resolution_clock::now();
    int count = 0;
    BoostPoint bg_q(P_base_vec[0].x, P_base_vec[0].y);
    for (auto it = tree.qbegin(bgi::nearest(bg_q, 1500000u)); it != tree.qend(); ++it) {
        count++;
        if (count >= 1500000) break;
    }
    auto t2 = chrono::high_resolution_clock::now();
    cout << "Real data K=1.5M took " << chrono::duration_cast<chrono::microseconds>(t2-t1).count() / 1000.0 << " ms" << endl;
    return 0;
}
