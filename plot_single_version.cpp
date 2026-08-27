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
    
    ofstream fout("single_version_mem.csv");
    fout << "Year,Algo,Mem_MB\n";
    
    // CPAMBB
    auto tree_c = CPAMBB::map_init(P, false);
    std::unordered_map<size_t, bool> mmp;
    auto f_noop = [&](const auto &et){ return 0; };
    fout << "2020,CPAMBB," << tree_c.size_in_bytes(f_noop, mmp) / 1024.0 / 1024.0 << "\n";
    
    // MVZD
    mvq::Tree tree_m(mvq::Config::get().leaf_size);
    auto P_set = geobase::get_sorted_points(P);
    tree_m.build(P_set);
    auto master_m = tree_m.root;
    auto stat = tree_m.get_tree_statistics();
    fout << "2020,MVZD," << stat.get_total() << "\n";

    for (int year = 2021; year <= 2025; year++) {
        string merged_file = "dataset/bhutan_evolution/" + to_string(year) + "_merged.txt";
        ifstream fin_year(merged_file);
        if (!fin_year.is_open()) continue;
        
        std::vector<geobase::Point> temp_vec;
        while (fin_year >> b_tid >> b_tx >> b_ty) {
            if (b_tx < 0 || b_ty < 0) continue;
            temp_vec.push_back(geobase::Point(b_tid, b_tx*1000000, b_ty*1000000));
        }
        parlay::sequence<geobase::Point> P_year(temp_vec.begin(), temp_vec.end());
        
        // CPAMBB
        tree_c = CPAMBB::map_insert(P_year, tree_c, false);
        mmp.clear();
        fout << year << ",CPAMBB," << tree_c.size_in_bytes(f_noop, mmp) / 1024.0 / 1024.0 << "\n";
        
        // MVZD
        parlay::sequence<geobase::Point> empty_pts;
        master_m = tree_m.commit(master_m, P_year, empty_pts);
        // Note: we do NOT store the old master_m, so it gets GC'd. This is pure single-version tracking!
        stat = tree_m.get_tree_statistics();
        fout << year << ",MVZD," << stat.get_total() << "\n";
    }

    return 0;
}
