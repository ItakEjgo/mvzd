#include <iostream>
#include <fstream>
#include "src/geo/point.hpp"
#include "src/cpambb.hpp"
using namespace std;
int main() {
    ifstream fin_base("dataset/bhutan_evolution/base_upto_2020.txt");
    std::vector<geobase::Point> pts;
    size_t b_tid; double b_tx, b_ty;
    while (fin_base >> b_tid >> b_tx >> b_ty) {
        if (b_tx < 0 || b_ty < 0) continue;
        pts.push_back(geobase::Point(b_tid, b_tx*1000000, b_ty*1000000));
    }
    parlay::sequence<geobase::Point> P(pts.begin(), pts.end());
    auto tree = CPAMBB::map_init(P, false);
    cout << "Nodes inserted: " << pts.size() << endl;
    cout << "Used nodes reported by GC: " << CPAMBB::zmap::GC::used_node() << endl;
    
    std::unordered_map<size_t, bool> mmp;
    auto f_noop = [&](const auto &et){ return 0; };
    double real_size = 1.0 * tree.size_in_bytes(f_noop, mmp);
    cout << "Real size_in_bytes: " << real_size / 1024.0 / 1024.0 << " MB" << endl;
    return 0;
}
