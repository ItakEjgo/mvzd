#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <cmath>
#include <algorithm>
#include "src/geo/point.hpp"

using namespace std;

int main() {
    ifstream fin_base("dataset/bhutan_evolution/base_upto_2020.txt");
    std::vector<geobase::Point> pts;
    size_t b_tid; double b_tx, b_ty;
    while (fin_base >> b_tid >> b_tx >> b_ty) {
        if (b_tx < 0 || b_ty < 0) continue;
        pts.push_back(geobase::Point(b_tid, b_tx*1000000, b_ty*1000000));
    }
    size_t N = pts.size();
    cout << "Loaded " << N << " points." << endl;
    
    double n14 = pow(N, 0.25);
    double n12 = pow(N, 0.5);
    cout << "N^1/4 = " << n14 << ", N^1/2 = " << n12 << endl;
    
    mt19937 rng(42);
    
    auto generate_file = [&](string filename, int min_c, int max_c) {
        ofstream fout(filename);
        fout << "42\n100 2\n"; // dummy seed and header
        
        uniform_int_distribution<int> count_dist(min_c, max_c);
        uniform_int_distribution<int> pt_dist(0, N-1);
        
        for(int i=0; i<100; i++) {
            int target_c = count_dist(rng);
            auto center = pts[pt_dist(rng)];
            
            double low_delta = 0.0;
            double high_delta = 10000000.0; // 10 degrees roughly
            double delta = 100.0;
            int actual_count = 0;
            
            // binary search for delta
            for(int iter=0; iter<50; iter++) {
                delta = (low_delta + high_delta) / 2.0;
                actual_count = 0;
                for(const auto& p : pts) {
                    if (p.x >= center.x - delta && p.x <= center.x + delta &&
                        p.y >= center.y - delta && p.y <= center.y + delta) {
                        actual_count++;
                    }
                }
                if (actual_count < target_c) {
                    low_delta = delta;
                } else if (actual_count > target_c) {
                    high_delta = delta;
                } else {
                    break;
                }
            }
            // Ensure integer coordinates for the bounding box
            long long minx = center.x - delta;
            long long miny = center.y - delta;
            long long maxx = center.x + delta;
            long long maxy = center.y + delta;
            fout << actual_count << " " << minx << " " << miny << " " << maxx << " " << maxy << "\n";
        }
        cout << "Generated " << filename << endl;
    };
    
    generate_file("dataset/bhutan_evolution/bhutan-small.qry", 5, max((int)n14, 5));
    generate_file("dataset/bhutan_evolution/bhutan-med.qry", (int)n14 + 1, (int)n12);
    generate_file("dataset/bhutan_evolution/bhutan-large.qry", (int)n12 + 1, (int)n12 * 4); // cap large at N^1/2 * 4 roughly 3200
    
    return 0;
}
