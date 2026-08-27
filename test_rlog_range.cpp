#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <queue>
#include <unordered_set>
#include <chrono>

#include "parlay/internal/get_time.h"
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
    vector<Value> pts;
    for(size_t i=0; i<6500000; i++) {
        pts.push_back(make_pair(BoostPoint(i*10.0, i*10.0), i));
    }
    cout << "Bulk loading RTree..." << endl;
    RTree tree(pts.begin(), pts.end());
    
    cout << "Testing Range Query..." << endl;
    bg::model::box<BoostPoint> box(BoostPoint(100.0, 100.0), BoostPoint(200.0, 200.0));
    
    parlay::internal::timer t; t.start();
    for (int i=0; i<100; i++) {
        std::vector<Value> snap_res;
        tree.query(bgi::intersects(box), std::back_inserter(snap_res));
    }
    cout << "100 Range Queries took " << t.stop()*1000.0 << " ms" << endl;

    return 0;
}
