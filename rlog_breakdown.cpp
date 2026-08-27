#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <queue>
#include <unordered_set>
#include <chrono>
#include <atomic>

#include "parlay/internal/get_time.h"
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/index/rtree.hpp>
#include "src/geo/point.hpp"

using namespace std;
namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

typedef bg::model::point<double, 2, bg::cs::cartesian> BoostPoint;
typedef pair<BoostPoint, size_t> Value;

inline std::atomic<size_t> boost_live_mem(0);

template <typename T>
class TrackingAllocator {
public:
    typedef T value_type;
    TrackingAllocator() = default;
    template <typename U> TrackingAllocator(const TrackingAllocator<U>&) {}
    
    T* allocate(std::size_t n) {
        boost_live_mem.fetch_add(n * sizeof(T), std::memory_order_relaxed);
        return static_cast<T*>(::operator new(n * sizeof(T)));
    }
    void deallocate(T* p, std::size_t n) {
        boost_live_mem.fetch_sub(n * sizeof(T), std::memory_order_relaxed);
        ::operator delete(p);
    }
};

template <typename T, typename U>
bool operator==(const TrackingAllocator<T>&, const TrackingAllocator<U>&) { return true; }
template <typename T, typename U>
bool operator!=(const TrackingAllocator<T>&, const TrackingAllocator<U>&) { return false; }

typedef bgi::rtree<Value, bgi::quadratic<32>, bgi::indexable<Value>, bgi::equal_to<Value>, TrackingAllocator<Value>> RTree;

struct RlogBranch {
    std::vector<Value> insert_log;
    std::vector<Value> remove_log;
};

ofstream debug_log;

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
    
    RlogTree() : threshold(9999999), current_batches(0) { snapshot = make_shared<RTree>(); }
    RlogTree(int thresh) : threshold(thresh), current_batches(0) { snapshot = make_shared<RTree>(); }
    
    void build_base(const std::vector<Value>& base_data) {
        parlay::internal::timer t; t.start();
        snapshot = make_shared<RTree>(base_data.begin(), base_data.end());
        debug_log << "      [RlogTree::build_base] took " << t.stop()*1000.0 << " ms\n";
    }
    void compact() {
        parlay::internal::timer t; t.start();
        if (insert_log.empty() && remove_log.empty()) return;
        
        std::vector<Value> all_pts(snapshot->begin(), snapshot->end());
        all_pts.insert(all_pts.end(), insert_log.begin(), insert_log.end());
        
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
            snapshot = make_shared<RTree>(filtered_pts.begin(), filtered_pts.end());
        } else {
            snapshot = make_shared<RTree>(all_pts.begin(), all_pts.end());
        }
        
        insert_log.clear();
        remove_log.clear();
        debug_log << "      [RlogTree::compact] took " << t.stop()*1000.0 << " ms\n";
    }
    void merge(const RlogBranch& branch) {
        parlay::internal::timer t; t.start();
        if (!branch.insert_log.empty()) insert_log.insert(insert_log.end(), branch.insert_log.begin(), branch.insert_log.end());
        if (!branch.remove_log.empty()) remove_log.insert(remove_log.end(), branch.remove_log.begin(), branch.remove_log.end());
        debug_log << "      [RlogTree::merge logs] took " << t.stop()*1000.0 << " ms\n";
        current_batches++;
        if (current_batches >= threshold) {
            debug_log << "      [RlogTree::merge] threshold triggered, compacting...\n";
            compact();
            current_batches = 0;
        }
    }
};

int main() {
    debug_log.open("rlog_breakdown.log");
    debug_log << "=== RlogTree Breakdown Test (Rlog_1yr, bp=50) ===\n";
    
    parlay::internal::timer t_main; t_main.start();
    
    debug_log << "[Step 1] Loading base dataset...\n";
    ifstream fin_base("dataset/bhutan_evolution/base_upto_2020.txt");
    vector<Value> P_base_conv;
    size_t b_tid; double b_tx, b_ty;
    while (fin_base >> b_tid >> b_tx >> b_ty) {
        if (b_tx < 0 || b_ty < 0) continue;
        P_base_conv.push_back(make_pair(BoostPoint(b_tx*1000000, b_ty*1000000), b_tid));
    }
    debug_log << "Loaded " << P_base_conv.size() << " points.\n";

    debug_log << "[Step 2] Initializing Rlog_1yr (threshold=2)...\n";
    RlogTree rlog(2);
    rlog.build_base(P_base_conv);

    debug_log << "[Step 3] Loading 2021 batch...\n";
    string merged_file = "dataset/bhutan_evolution/2021_merged.txt";
    ifstream fin_year(merged_file);
    vector<Value> P_year_conv;
    while (fin_year >> b_tid >> b_tx >> b_ty) {
        if (b_tx < 0 || b_ty < 0) continue;
        P_year_conv.push_back(make_pair(BoostPoint(b_tx*1000000, b_ty*1000000), b_tid));
    }
    
    size_t batch_size = P_year_conv.size() / 2; // bp=50
    vector<Value> batch0(P_year_conv.begin(), P_year_conv.begin() + batch_size);
    vector<Value> batch1(P_year_conv.begin() + batch_size, P_year_conv.end());

    vector<RlogBranch> history;

    debug_log << "[Step 4] Simulating Forward Fork Batch 0...\n";
    parlay::internal::timer t_fork; t_fork.start();
    RlogBranch branch0;
    branch0.insert_log = batch0;
    history.push_back(branch0);
    debug_log << "Fork Batch 0 took " << t_fork.stop()*1000.0 << " ms\n";

    debug_log << "[Step 5] Simulating Forward Fork Batch 1...\n";
    t_fork.start();
    RlogBranch branch1;
    branch1.insert_log = batch1;
    history.push_back(branch1);
    debug_log << "Fork Batch 1 took " << t_fork.stop()*1000.0 << " ms\n";

    debug_log << "[Step 6] Simulating Forward Merge Batch 0...\n";
    parlay::internal::timer t_merge; t_merge.start();
    rlog.merge(history[0]);
    debug_log << "Merge Batch 0 overall took " << t_merge.stop()*1000.0 << " ms\n";

    debug_log << "[Step 7] Simulating Forward Merge Batch 1...\n";
    t_merge.start();
    rlog.merge(history[1]);
    debug_log << "Merge Batch 1 overall took " << t_merge.stop()*1000.0 << " ms\n";

    debug_log << "=== Breakdown Complete in " << t_main.stop() << " sec ===\n";
    debug_log.close();
    return 0;
}
