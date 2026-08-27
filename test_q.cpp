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
    cout << "Base points: " << P_base_vec.size() << endl;
    
    ifstream fq("dataset/bhutan_evolution/bhutan-med.qry");
    string dummy; getline(fq, dummy); getline(fq, dummy);
    int expected, c=0;
    double min_x, min_y, max_x, max_y;
    while(fq >> expected >> min_x >> min_y >> max_x >> max_y && c < 5) {
        int count = 0;
        for(auto& p : P_base_vec) {
            if (p.x >= min_x && p.x <= max_x && p.y >= min_y && p.y <= max_y) count++;
        }
        cout << "Query " << c << " expected: " << expected << " actual: " << count << endl;
        c++;
    }
    return 0;
}
