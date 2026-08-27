#include <iostream>
#include "src/cpambb.hpp"
#include <parlay/alloc.h>
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
    
    size_t initial_bytes = parlay::num_used_bytes();
    cout << "Initial Parlay used bytes: " << initial_bytes / 1024.0 / 1024.0 << " MB" << endl;
    
    auto tree = CPAMBB::map_init(P, false);
    
    size_t after_bytes = parlay::num_used_bytes();
    cout << "After tree build Parlay used bytes: " << (after_bytes - initial_bytes) / 1024.0 / 1024.0 << " MB" << endl;
    
    cout << "Used nodes reported by GC: " << CPAMBB::zmap::GC::used_node() << endl;
    
    return 0;
}
