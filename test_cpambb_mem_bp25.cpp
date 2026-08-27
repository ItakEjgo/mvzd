#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include "src/geo/point.hpp"
#include "src/cpambb.hpp"
#include <parlay/sequence.h>

using namespace std;

// New Method: O(1) Amortized global memory tracker (including hidden leaf arrays)
std::unordered_map<size_t, bool> cpambb_mmp;
double cpambb_global_mem = 0;

double mem_cpambb_new(const CPAMBB::zmap& latest_branch) {
    auto f_noop = [&](const auto &et){ return 0; };
    cpambb_global_mem += latest_branch.size_in_bytes(f_noop, cpambb_mmp) / (1024.0 * 1024.0);
    return cpambb_global_mem;
}

// Old Method: Only counts regular interior nodes via PAM GC
double mem_cpambb_old() {
    return (sizeof(CPAMBB::zmap::GC::regular_node) * CPAMBB::zmap::GC::used_node()) / (1024.0 * 1024.0);
}

int main() {
    cout << "Loading dataset..." << endl;
    ifstream fin_base("dataset/bhutan_evolution/base_upto_2020.txt");
    std::vector<geobase::Point> pts;
    size_t b_tid; double b_tx, b_ty;
    while (fin_base >> b_tid >> b_tx >> b_ty) {
        if (b_tx < 0 || b_ty < 0) continue;
        pts.push_back(geobase::Point(b_tid, b_tx*1000000, b_ty*1000000));
    }
    parlay::sequence<geobase::Point> P_base(pts.begin(), pts.end());

    ifstream fin_2021("dataset/bhutan_evolution/2021_merged.txt");
    std::vector<geobase::Point> pts21;
    while (fin_2021 >> b_tid >> b_tx >> b_ty) {
        if (b_tx < 0 || b_ty < 0) continue;
        pts21.push_back(geobase::Point(b_tid, b_tx*1000000, b_ty*1000000));
    }
    parlay::sequence<geobase::Point> P21(pts21.begin(), pts21.end());

    vector<CPAMBB::zmap> global_history;

    cout << "Building Base Tree (2020)..." << endl;
    auto tree = CPAMBB::map_init(P_base, false);
    global_history.push_back(tree);

    cout << "================================================================" << endl;
    cout << "[Base Tree]  Old Method (Interior Only): " << mem_cpambb_old() << " MB" << endl;
    cout << "[Base Tree]  New Method (True Memory)  : " << mem_cpambb_new(tree) << " MB" << endl;
    cout << "================================================================" << endl;

    double bp = 25.0; // BP=25 means 4 batches
    size_t batch_size = max(1ul, (size_t)ceil(P21.size() * (bp / 100.0)));
    vector<parlay::sequence<geobase::Point>> batches;
    for (size_t i = 0; i < P21.size(); i += batch_size) {
        size_t end = min(P21.size(), i + batch_size);
        batches.push_back(P21.substr(i, end - i));
    }
    
    cout << "\nStarting 2021 Forward (BP=25%, " << batches.size() << " Batches)..." << endl;
    
    vector<CPAMBB::zmap> branches(batches.size());
    
    for(size_t b = 0; b < batches.size(); b++) {
        auto cv = CPAMBB::map_insert(batches[b], tree, false);
        branches[b] = cv;
        global_history.push_back(cv);
        
        cout << "[Batch " << b << " Fork] Old Method: " << mem_cpambb_old() << " MB | New Method: " << mem_cpambb_new(cv) << " MB" << endl;
    }
    
    auto master = tree;
    for(size_t b = 0; b < batches.size(); b++) {
        master = std::get<0>(CPAMBB::map_merge(tree, master, branches[b]));
        global_history.push_back(master);
        
        cout << "[Batch " << b << " Merge] Old Method: " << mem_cpambb_old() << " MB | New Method: " << mem_cpambb_new(master) << " MB" << endl;
    }

    cout << "================================================================" << endl;
    cout << "TEST COMPLETED." << endl;
    cout << "FINAL DIFFERENCE -> Old: " << mem_cpambb_old() << " MB vs New: " << mem_cpambb_new(master) << " MB" << endl;
    cout << "================================================================" << endl;

    return 0;
}
