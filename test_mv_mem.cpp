#include <iostream>
#include <fstream>
#include <vector>
#include "src/geo/point.hpp"
#include "src/cpambb.hpp"
#include "src/mvq.hpp"
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
    
    ifstream fin_2021("dataset/bhutan_evolution/2021_merged.txt");
    std::vector<geobase::Point> pts21;
    while (fin_2021 >> b_tid >> b_tx >> b_ty) {
        if (b_tx < 0 || b_ty < 0) continue;
        pts21.push_back(geobase::Point(b_tid, b_tx*1000000, b_ty*1000000));
    }
    parlay::sequence<geobase::Point> P21(pts21.begin(), pts21.end());

    cout << "--- CPAMBB ---" << endl;
    auto tree_c = CPAMBB::map_init(P, false);
    std::unordered_map<size_t, bool> mmp;
    auto f_noop = [&](const auto &et){ return 0; };
    double base_mem_c = tree_c.size_in_bytes(f_noop, mmp);
    cout << "Base: " << base_mem_c / 1024.0 / 1024.0 << " MB" << endl;
    
    vector<CPAMBB::zmap> hist_c;
    hist_c.push_back(tree_c);
    
    // insert 2 batches from 2021
    size_t batch_size = P21.size() / 2;
    for(int b=0; b<2; b++) {
        auto batch = P21.substr(b * batch_size, batch_size);
        tree_c = CPAMBB::map_insert(batch, tree_c, false);
        hist_c.push_back(tree_c);
    }
    
    mmp.clear();
    double total_mem_c = 0;
    for(auto& t : hist_c) total_mem_c += t.size_in_bytes(f_noop, mmp);
    cout << "After 2 batches (1 year BP=50): " << total_mem_c / 1024.0 / 1024.0 << " MB" << endl;

    cout << "\n--- MVZD ---" << endl;
    mvq::Tree tree_m(mvq::Config::get().leaf_size);
    auto P_set = geobase::get_sorted_points(P);
    tree_m.build(P_set);
    auto master_m = tree_m.root;
    tree_m.multi_version_roots.push_back(master_m);
    
    auto stat_base = tree_m.get_tree_statistics();
    cout << "Base: " << stat_base.get_total() << " MB" << endl;
    
    parlay::sequence<geobase::Point> empty_pts;
    for(int b=0; b<2; b++) {
        auto batch = P21.substr(b * batch_size, batch_size);
        auto P_batch_set = geobase::get_sorted_points(batch);
        master_m = tree_m.multi_version_batch_insert_sorted(P_batch_set, master_m);
        tree_m.multi_version_roots.push_back(master_m);
    }
    
    auto stat_after = tree_m.get_tree_statistics();
    cout << "After 2 batches (1 year BP=50): " << stat_after.get_total() << " MB" << endl;

    return 0;
}
