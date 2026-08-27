#include <iostream>
#include <vector>
#include <chrono>
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/index/rtree.hpp>

using namespace std;
namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

typedef bg::model::point<double, 2, bg::cs::cartesian> BoostPoint;
typedef pair<BoostPoint, size_t> Value;
typedef bgi::rtree<Value, bgi::quadratic<32>> RTree;

int main() {
    vector<Value> pts;
    for(size_t i=0; i<6500000; i++) {
        pts.push_back(make_pair(BoostPoint(i*1.0, i*1.0), i));
    }
    auto t1 = chrono::high_resolution_clock::now();
    cout << "Bulk loading..." << endl;
    RTree tree(pts.begin(), pts.end());
    auto t2 = chrono::high_resolution_clock::now();
    cout << "Bulk load took " << chrono::duration_cast<chrono::milliseconds>(t2-t1).count() << " ms" << endl;
    
    cout << "Extracting..." << endl;
    vector<Value> out(tree.begin(), tree.end());
    auto t3 = chrono::high_resolution_clock::now();
    cout << "Extraction took " << chrono::duration_cast<chrono::milliseconds>(t3-t2).count() << " ms" << endl;
    return 0;
}
