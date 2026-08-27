#include <iostream>
#include <fstream>
#include <vector>
#include "src/geo/point.hpp"

using namespace std;

int main() {
    ifstream fin_base("dataset/bhutan_evolution/base_upto_2020.txt");
    std::vector<geobase::Point> P_base_vec;
    size_t b_tid; double b_tx, b_ty;
    while (fin_base >> b_tid >> b_tx >> b_ty) {
        if (b_tx < 0 || b_ty < 0) continue;
        P_base_vec.push_back(geobase::Point(b_tid, b_tx*1000000, b_ty*1000000));
    }
    
    srand(42);
    vector<int> K_list = {1, 10, 100};
    // BaseQuery (is_warmup=true)
    for(size_t K: K_list) for(int i=0; i<100; i++) rand();
    // BaseQuery (is_warmup=false)
    for(size_t K: K_list) for(int i=0; i<100; i++) rand();
    
    // 2021 MidQuery (is_warmup=true)
    for(size_t K: K_list) {
        if (K == 10) {
            for(int i=0; i<100; i++) {
                int r = rand() % P_base_vec.size();
                if (i == 8) cout << "MidQuery Warmup KNN_10_8 ID: " << P_base_vec[r].id << endl;
            }
        } else {
            for(int i=0; i<100; i++) rand();
        }
    }
    // 2021 MidQuery (is_warmup=false)
    for(size_t K: K_list) {
        if (K == 10) {
            for(int i=0; i<100; i++) {
                int r = rand() % P_base_vec.size();
                if (i == 8) {
                    cout << "MidQuery Actual KNN_10_8 ID: " << P_base_vec[r].id << " (" << P_base_vec[r].x << ", " << P_base_vec[r].y << ")" << endl;
                }
            }
        } else {
            for(int i=0; i<100; i++) rand();
        }
    }
    
    return 0;
}
