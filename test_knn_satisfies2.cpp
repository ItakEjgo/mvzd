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

class MyTree {
    RTree tree;
    unordered_set<size_t> removed_ids;
public:
    MyTree() {
        for(size_t i=1; i<=10; i++) tree.insert(make_pair(BoostPoint(i, 0), i));
        removed_ids.insert(1); removed_ids.insert(2); removed_ids.insert(3); removed_ids.insert(4);
    }
    
    void query() {
        BoostPoint bg_q(0,0);
        vector<Value> result;
        auto is_alive = [this](Value const& v) { return removed_ids.find(v.second) == removed_ids.end(); };
        tree.query(bgi::nearest(bg_q, 3) && bgi::satisfies(is_alive), std::back_inserter(result));
        for(auto v : result) cout << "ID: " << v.second << " ";
        cout << endl;
    }
};

int main() {
    MyTree t;
    t.query();
    return 0;
}
