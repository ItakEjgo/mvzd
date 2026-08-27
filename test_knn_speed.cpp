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
    for(size_t i=0; i<650000; i++) {
        pts.push_back(make_pair(BoostPoint(i*1.0, i*1.0), i));
    }
    RTree tree(pts.begin(), pts.end());
    BoostPoint bg_q(0.0, 0.0);
    
    auto t1 = chrono::high_resolution_clock::now();
    for(int i=0; i<100; i++) {
        int count = 0;
        for (auto it = tree.qbegin(bgi::nearest(bg_q, (unsigned)tree.size())); it != tree.qend(); ++it) {
            count++;
            if (count >= 100) break;
        }
    }
    auto t2 = chrono::high_resolution_clock::now();
    cout << "Tree size parameter took " << chrono::duration_cast<chrono::microseconds>(t2-t1).count() / 1000.0 << " ms" << endl;

    auto t3 = chrono::high_resolution_clock::now();
    for(int i=0; i<100; i++) {
        int count = 0;
        for (auto it = tree.qbegin(bgi::nearest(bg_q, 100u)); it != tree.qend(); ++it) {
            count++;
            if (count >= 100) break;
        }
    }
    auto t4 = chrono::high_resolution_clock::now();
    cout << "K parameter took " << chrono::duration_cast<chrono::microseconds>(t4-t3).count() / 1000.0 << " ms" << endl;
    return 0;
}
