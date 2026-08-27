#include <iostream>
#include <fstream>
#include <vector>
#include "parlay/internal/get_time.h"
#include "src/mvq.hpp"
#include "src/cpamz.hpp"
#include "src/global_config.hpp"
#include "src/cpambb.hpp"
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/index/rtree.hpp>

using namespace std;
namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;
typedef bg::model::point<double, 2, bg::cs::cartesian> BoostPoint;
typedef pair<BoostPoint, size_t> Value;
typedef bgi::rtree<Value, bgi::quadratic<32>, bgi::indexable<Value>, bgi::equal_to<Value>> RTree;

class RlogTree {
public:
    shared_ptr<RTree> snapshot;
    void build_base(const std::vector<Value>& base_data) {
        snapshot = make_shared<RTree>(base_data.begin(), base_data.end());
    }
};

int main() {
    ifstream fin_base("dataset/bhutan_evolution/base_upto_2020.txt");
    std::vector<geobase::Point> P_base_vec;
    size_t b_tid; double b_tx, b_ty;
    while (fin_base >> b_tid >> b_tx >> b_ty) {
        if (b_tx < 0 || b_ty < 0) continue;
        P_base_vec.push_back(geobase::Point(b_tid, b_tx*1000000, b_ty*1000000));
    }
    parlay::sequence<geobase::Point> P_base(P_base_vec.begin(), P_base_vec.end());
    auto P_base_set = geobase::get_sorted_points(P_base);
    
    vector<Value> P_base_conv(P_base.size());
    for(size_t i=0; i<P_base.size(); i++) P_base_conv[i] = make_pair(BoostPoint(P_base[i].x, P_base[i].y), P_base[i].id);

    parlay::internal::timer timer;

    timer.start();
    mvq::Tree mvzd_tree(mvq::Config::get().leaf_size);
    mvzd_tree.build(P_base_set);
    cout << "MVZD: " << timer.stop() << " s\n";

    timer.start();
    auto cpambb_master = CPAMBB::map_init(P_base_set, false);
    cout << "CPAMBB: " << timer.stop() << " s\n";

    timer.start();
    RlogTree rlog;
    rlog.build_base(P_base_conv);
    cout << "Rlog: " << timer.stop() << " s\n";

    return 0;
}
