#include <iostream>
#include <fstream>
#include <vector>
#include "src/geo/point.hpp"
#include <queue>

using namespace std;

typedef pair<geobase::Point, double> nn_pair;

struct nn_pair_cmp {
    bool operator()(const nn_pair &lhs, const nn_pair &rhs) const {
        return lhs.second < rhs.second ||
               (lhs.second == rhs.second && lhs.first.id > rhs.first.id);
    }
};

int main() {
    ifstream fin_base("dataset/bhutan_evolution/base_upto_2020.txt");
    std::vector<geobase::Point> pts;
    size_t b_tid; double b_tx, b_ty;
    while (fin_base >> b_tid >> b_tx >> b_ty) {
        if (b_tx < 0 || b_ty < 0) continue;
        pts.push_back(geobase::Point(b_tid, b_tx*1000000, b_ty*1000000));
    }
    
    ifstream fin_2021("dataset/bhutan_evolution/2021_merged.txt");
    while (fin_2021 >> b_tid >> b_tx >> b_ty) {
        if (b_tx < 0 || b_ty < 0) continue;
        pts.push_back(geobase::Point(b_tid, b_tx*1000000, b_ty*1000000));
    }

    geobase::Point q(5400561020, 8.9976e+07, 2.75181e+07); // KNN_10_8 query point
    size_t K = 10;
    
    auto calc_dist = [](const geobase::Point& p1, const geobase::Point& p2) {
        double dx = p1.x - p2.x, dy = p1.y - p2.y;
        return dx*dx + dy*dy;
    };
    
    priority_queue<nn_pair, vector<nn_pair>, nn_pair_cmp> pq;
    for(auto& p : pts) {
        double d = calc_dist(p, q);
        nn_pair cur = {p, d};
        if(pq.size() < K) {
            pq.push(cur);
        } else {
            nn_pair_cmp cmp;
            if(cmp(cur, pq.top())) {
                pq.pop();
                pq.push(cur);
            }
        }
    }
    
    size_t hash_val = 0;
    while(!pq.empty()) {
        hash_val += pq.top().first.id;
        pq.pop();
    }
    
    cout << "Brute Force Hash: " << hash_val << endl;
    return 0;
}
