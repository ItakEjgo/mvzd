#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <queue>
#include <unordered_set>
#include <chrono>

#include "parlay/internal/get_time.h"
#include <boost/geometry.hpp>
#include "src/geo/point.hpp"
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/index/rtree.hpp>

using namespace std;
namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

typedef bg::model::point<double, 2, bg::cs::cartesian> BoostPoint;
typedef pair<BoostPoint, size_t> Value;
typedef bgi::rtree<Value, bgi::quadratic<32>> RTree;

struct RlogBranch {
    std::vector<Value> insert_log;
    std::vector<Value> remove_log;
};

class RlogTree {
public:
    struct MaxHeapCmp {
        bool operator()(const std::pair<double, Value>& a, const std::pair<double, Value>& b) const {
            return a.first < b.first;
        }
    };
    int threshold;
    int current_batches;
    shared_ptr<RTree> snapshot;
    std::vector<Value> insert_log;
    std::vector<Value> remove_log;
    
    RlogTree(int thresh) : threshold(thresh), current_batches(0) { snapshot = make_shared<RTree>(); }
    
    void build_base(const std::vector<Value>& base_data) {
        parlay::internal::timer t; t.start();
        snapshot = make_shared<RTree>(base_data.begin(), base_data.end());
        cout << "    [RlogTree] build_base took " << t.stop() * 1000.0 << " ms" << endl;
    }
    
    void compact() {
        if (insert_log.empty() && remove_log.empty()) return;
        parlay::internal::timer t_total; t_total.start();
        
        parlay::internal::timer t_step; t_step.start();
        std::vector<Value> all_pts(snapshot->begin(), snapshot->end());
        cout << "      [compact] Extract from snapshot took " << t_step.stop() * 1000.0 << " ms" << endl;
        
        t_step.start();
        all_pts.insert(all_pts.end(), insert_log.begin(), insert_log.end());
        cout << "      [compact] Append insert_log took " << t_step.stop() * 1000.0 << " ms" << endl;
        
        t_step.start();
        if (!remove_log.empty()) {
            std::unordered_set<size_t> removed_ids;
            for (const auto& val : remove_log) removed_ids.insert(val.second);
            std::vector<Value> filtered_pts;
            filtered_pts.reserve(all_pts.size());
            for (const auto& val : all_pts) {
                if (removed_ids.find(val.second) == removed_ids.end()) {
                    filtered_pts.push_back(val);
                }
            }
            cout << "      [compact] Filter remove_log took " << t_step.stop() * 1000.0 << " ms" << endl;
            t_step.start();
            snapshot = make_shared<RTree>(filtered_pts.begin(), filtered_pts.end());
            cout << "      [compact] Bulk load RTree (filtered) took " << t_step.stop() * 1000.0 << " ms" << endl;
        } else {
            cout << "      [compact] No removes to filter." << endl;
            snapshot = make_shared<RTree>(all_pts.begin(), all_pts.end());
            cout << "      [compact] Bulk load RTree took " << t_step.stop() * 1000.0 << " ms" << endl;
        }
        
        insert_log.clear();
        remove_log.clear();
        cout << "    [RlogTree] total compact took " << t_total.stop() * 1000.0 << " ms" << endl;
    }
    
    void merge(const RlogBranch& branch) {
        parlay::internal::timer t; t.start();
        if (!branch.insert_log.empty()) insert_log.insert(insert_log.end(), branch.insert_log.begin(), branch.insert_log.end());
        if (!branch.remove_log.empty()) remove_log.insert(remove_log.end(), branch.remove_log.begin(), branch.remove_log.end());
        cout << "    [RlogTree] merge logs took " << t.stop() * 1000.0 << " ms" << endl;
        
        current_batches++;
        if (current_batches >= threshold) {
            cout << "    [RlogTree] Threshold reached. Calling compact()..." << endl;
            compact();
            current_batches = 0;
        }
    }
    
    std::vector<Value> knn_report(const geobase::Point& q, size_t k) const {
        std::unordered_set<size_t> removed_ids;
        for (const auto& val : remove_log) removed_ids.insert(val.second);
        auto calc_sqr_dist = [](const geobase::Point& p1, const BoostPoint& p2) {
            double dx = p1.x - p2.get<0>(), dy = p1.y - p2.get<1>();
            return dx*dx + dy*dy;
        };
        std::priority_queue<std::pair<double, Value>, std::vector<std::pair<double, Value>>, MaxHeapCmp> max_heap;
        for (const auto& val : insert_log) {
            if (removed_ids.find(val.second) == removed_ids.end()) {
                double dist = calc_sqr_dist(q, val.first);
                if (max_heap.size() < k) max_heap.push({dist, val});
                else if (dist < max_heap.top().first) { max_heap.pop(); max_heap.push({dist, val}); }
            }
        }
        BoostPoint bg_q(q.x, q.y);
        for (auto it = snapshot->qbegin(bgi::nearest(bg_q, (unsigned int)snapshot->size())); it != snapshot->qend(); ++it) {
            if (removed_ids.find(it->second) == removed_ids.end()) {
                double dist = calc_sqr_dist(q, it->first);
                if (max_heap.size() < k) max_heap.push({dist, *it});
                else if (dist < max_heap.top().first) { max_heap.pop(); max_heap.push({dist, *it}); }
                else break;
            }
        }
        std::vector<Value> result;
        while (!max_heap.empty()) { result.push_back(max_heap.top().second); max_heap.pop(); }
        return result;
    }
};

int main() {
    cout << "=== Rlog_1yr Bottleneck Diagnostics ===" << endl;
    
    parlay::internal::timer t_main; t_main.start();
    
    cout << "1. Loading base dataset..." << endl;
    ifstream fin_base("dataset/bhutan_evolution/base_upto_2020.txt");
    vector<Value> P_base_conv;
    size_t b_tid; double b_tx, b_ty;
    while (fin_base >> b_tid >> b_tx >> b_ty) {
        if (b_tx < 0 || b_ty < 0) continue;
        P_base_conv.push_back(make_pair(BoostPoint(b_tx*1000000, b_ty*1000000), b_tid));
    }
    cout << "   Base points: " << P_base_conv.size() << endl;

    cout << "2. Initializing Rlog_1yr (threshold=2)..." << endl;
    RlogTree rlog(2);
    rlog.build_base(P_base_conv);

    cout << "3. Loading 2021 batch..." << endl;
    string merged_file = "dataset/bhutan_evolution/2021_merged.txt";
    ifstream fin_year(merged_file);
    vector<Value> P_year_conv;
    while (fin_year >> b_tid >> b_tx >> b_ty) {
        if (b_tx < 0 || b_ty < 0) continue;
        P_year_conv.push_back(make_pair(BoostPoint(b_tx*1000000, b_ty*1000000), b_tid));
    }
    cout << "   2021 total points: " << P_year_conv.size() << endl;
    
    size_t batch_size = P_year_conv.size() / 2;
    vector<Value> batch0(P_year_conv.begin(), P_year_conv.begin() + batch_size);
    vector<Value> batch1(P_year_conv.begin() + batch_size, P_year_conv.end());

    cout << "4. Simulating Forward Fork..." << endl;
    RlogBranch branch0, branch1;
    branch0.insert_log = batch0;
    branch1.insert_log = batch1;

    cout << "5. Simulating Forward Merge Batch 0..." << endl;
    parlay::internal::timer tm; tm.start();
    rlog.merge(branch0);
    cout << "   Merge Batch 0 took " << tm.stop()*1000.0 << " ms" << endl;

    cout << "6. Simulating Forward Merge Batch 1..." << endl;
    tm.start();
    rlog.merge(branch1);
    cout << "   Merge Batch 1 took " << tm.stop()*1000.0 << " ms" << endl;

    cout << "7. Testing KNN Queries (k=100)..." << endl;
    geobase::Point qp(0, P_base_conv[0].first.get<0>(), P_base_conv[0].first.get<1>());
    tm.start();
    for (int i=0; i<100; i++) {
        rlog.knn_report(qp, 100);
    }
    cout << "   100 KNN queries took " << tm.stop()*1000.0 << " ms" << endl;

    cout << "=== Diagnostics Complete in " << t_main.stop() << " sec ===" << endl;
    return 0;
}
