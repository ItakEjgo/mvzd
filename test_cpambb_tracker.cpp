#include <iostream>
#include <fstream>
#include <vector>
#include "src/geo/point.hpp"
#include "src/cpambb.hpp"

using namespace std;

std::unordered_map<size_t, bool> cpambb_mmp;
double cpambb_global_mem = 0;

double mem_cpambb(const CPAMBB::zmap& latest_branch) {
    auto f_noop = [&](const auto &et){ return 0; };
    cpambb_global_mem += latest_branch.size_in_bytes(f_noop, cpambb_mmp) / (1024.0 * 1024.0);
    return cpambb_global_mem;
}

int main() {
    ifstream fin_base("dataset/bhutan_evolution/base_upto_2020.txt");
    std::vector<geobase::Point> pts;
    size_t b_tid; double b_tx, b_ty;
    while (fin_base >> b_tid >> b_tx >> b_ty) {
        if (b_tx < 0 || b_ty < 0) continue;
        pts.push_back(geobase::Point(b_tid, b_tx*1000000, b_ty*1000000));
    }
    parlay::sequence<geobase::Point> P(pts.begin(), pts.end());
    
    auto tree_c = CPAMBB::map_init(P, false);
    cout << "Base mem: " << mem_cpambb(tree_c) << " MB" << endl;
    cout << "Second call mem: " << mem_cpambb(tree_c) << " MB" << endl;
    
    // Add 10 points
    auto P2 = P.substr(0, 10);
    for(int i=0; i<10; i++) P2[i].id += 9999999;
    tree_c = CPAMBB::map_insert(P2, tree_c, false);
    
    cout << "After small insert mem: " << mem_cpambb(tree_c) << " MB" << endl;
    
    return 0;
}
