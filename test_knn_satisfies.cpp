#include <iostream>
#include <vector>
#include <unordered_set>
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
    RTree tree;
    // Insert 10 points at x=1 to 10
    for(size_t i=1; i<=10; i++) {
        tree.insert(make_pair(BoostPoint(i, 0), i));
    }
    
    // We want 3 points nearest to origin, but ID must be >= 5 (so 1,2,3,4 are deleted)
    BoostPoint bg_q(0,0);
    vector<Value> result;
    
    auto is_valid = [](Value const& v) { return v.second >= 5; };
    
    tree.query(bgi::nearest(bg_q, 3) && bgi::satisfies(is_valid), std::back_inserter(result));
    
    cout << "Returned count: " << result.size() << endl;
    for(auto v : result) cout << "ID: " << v.second << " ";
    cout << endl;
    
    return 0;
}
