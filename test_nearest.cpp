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
    RTree tree(pts.begin(), pts.end());
    
    BoostPoint bg_q(0.0, 0.0);
    auto t1 = chrono::high_resolution_clock::now();
    cout << "Creating nearest iterator with 6500000..." << endl;
    auto it = tree.qbegin(bgi::nearest(bg_q, (unsigned int)tree.size()));
    auto t2 = chrono::high_resolution_clock::now();
    cout << "Iterator creation took " << chrono::duration_cast<chrono::milliseconds>(t2-t1).count() << " ms" << endl;
    
    cout << "Iterating 10 elements..." << endl;
    int count = 0;
    for(; it != tree.qend(); ++it) {
        count++;
        if (count >= 10) break;
    }
    auto t3 = chrono::high_resolution_clock::now();
    cout << "Iteration took " << chrono::duration_cast<chrono::milliseconds>(t3-t2).count() << " ms" << endl;
    return 0;
}
