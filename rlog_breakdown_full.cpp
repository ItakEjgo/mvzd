#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <queue>
#include <unordered_set>
#include <atomic>

#include "parlay/internal/get_time.h"
#include "helper/time_loop.h"
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/index/rtree.hpp>
#include "src/geo/point.hpp"
#include "src/mvq.hpp"

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
template <typename T, typename U> bool operator==(const TrackingAllocator<T>&, const TrackingAllocator<U>&) { return true; }
template <typename T, typename U> bool operator!=(const TrackingAllocator<T>&, const TrackingAllocator<U>&) { return false; }

typedef bgi::rtree<Value, bgi::quadratic<32>, bgi::indexable<Value>, bgi::equal_to<Value>, TrackingAllocator<Value>> RTree;

struct RlogBranch {
    std::vector<Value> insert_log;
    std::vector<Value> remove_log;
};

ofstream blog;
#define LOG(msg) blog << msg << std::endl; cout << msg << std::endl;

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
        snapshot = make_shared<RTree>(base_data.begin(), base_data.end());
    }
    void compact() {
        if (insert_log.empty() && remove_log.empty()) return;
        parlay::internal::timer t; t.start();
        std::vector<Value> all_pts(snapshot->begin(), snapshot->end());
        all_pts.insert(all_pts.end(), insert_log.begin(), insert_log.end());
        
        if (!remove_log.empty()) {
            std::unordered_set<size_t> removed_ids;
            for (const auto& val : remove_log) removed_ids.insert(val.second);
            std::vector<Value> filtered_pts;
            filtered_pts.reserve(all_pts.size());
            for (const auto& val : all_pts) {
                if (removed_ids.find(val.second) == removed_ids.end()) filtered_pts.push_back(val);
            }
            snapshot = make_shared<RTree>(filtered_pts.begin(), filtered_pts.end());
        } else {
            snapshot = make_shared<RTree>(all_pts.begin(), all_pts.end());
        }
        
        insert_log.clear();
        remove_log.clear();
        LOG("        [Compact] took " << t.stop() * 1000.0 << " ms");
    }
    void merge(const RlogBranch& branch) {
        if (!branch.insert_log.empty()) insert_log.insert(insert_log.end(), branch.insert_log.begin(), branch.insert_log.end());
        if (!branch.remove_log.empty()) remove_log.insert(remove_log.end(), branch.remove_log.begin(), branch.remove_log.end());
        current_batches++;
        if (current_batches >= threshold) {
            compact();
            current_batches = 0;
        }
    }
    std::vector<Value> range_report(const geobase::Bounding_Box& q) const {
        std::vector<Value> snap_res;
        bg::model::box<BoostPoint> box(BoostPoint(q.first.x, q.first.y), BoostPoint(q.second.x, q.second.y));
        snapshot->query(bgi::intersects(box), std::back_inserter(snap_res));
        std::unordered_set<size_t> removed_ids;
        for (const auto& val : remove_log) removed_ids.insert(val.second);
        std::vector<Value> result;
        for (const auto& val : snap_res) {
            if (removed_ids.find(val.second) == removed_ids.end()) result.push_back(val);
        }
        for (const auto& val : insert_log) {
            if (removed_ids.find(val.second) == removed_ids.end()) {
                if (bg::within(val.first, box)) result.push_back(val);
            }
        }
        return result;
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
        for (auto it = snapshot->qbegin(bgi::nearest(bg_q, (unsigned)snapshot->size())); it != snapshot->qend(); ++it) {
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

template<typename F_Init, typename F_Run, typename F_End>
double run_batch_time_loop(F_Init initf, F_Run runf, F_End endf) {
    return time_loop(3, 1.0, initf, runf, endf) * 1000.0;
}

int main() {
    blog.open("rlog_full_breakdown.log");
    LOG("=== Rlog_1yr Detailed Breakdown ===");
    parlay::internal::timer global_t; global_t.start();
    
    ifstream fin_base("dataset/bhutan_evolution/base_upto_2020.txt");
    std::vector<geobase::Point> P_base_vec;
    size_t b_tid; double b_tx, b_ty;
    while (fin_base >> b_tid >> b_tx >> b_ty) {
        if (b_tx < 0 || b_ty < 0) continue;
        P_base_vec.push_back(geobase::Point(b_tid, b_tx*1000000, b_ty*1000000));
    }
    vector<Value> P_base_conv(P_base_vec.size());
    for(size_t i=0; i<P_base_vec.size(); i++) P_base_conv[i] = make_pair(BoostPoint(P_base_vec[i].x, P_base_vec[i].y), P_base_vec[i].id);
    LOG("[Init] Loaded Base Data: " << P_base_conv.size() << " pts");

    auto load_q = [](string f) {
        parlay::sequence<geobase::Bounding_Box> queries;
        size_t dummy_maxSize;
        auto res = geobase::read_range_query(f, 8, dummy_maxSize);
        auto q = std::get<1>(res);
        queries = q.substr(0, std::min((size_t)100, q.size()));
        return queries;
    };
    auto queries_large = load_q("dataset/bhutan_evolution/bhutan-large.qry");
    
    parlay::sequence<geobase::Point> current_knn(100);
    srand(42);
    for(int i=0; i<100; i++) current_knn[i] = P_base_vec[rand() % P_base_vec.size()];

    RlogTree master(2);
    master.build_base(P_base_conv);
    vector<RlogBranch> global_history;
    LOG("[Init] Built Base Tree");

    auto run_queries = [&](string phase) {
        parlay::internal::timer qt; qt.start();
        LOG("  [Query] Starting " << phase << " Range Queries...");
        auto runf_range = [&]() {
            for(size_t i=0; i<queries_large.size(); i++) {
                auto res = master.range_report(queries_large[i]);
            }
        };
        double t_range = run_batch_time_loop([](){}, runf_range, [](){});
        LOG("  [Query] Range Queries took " << t_range << " ms");

        LOG("  [Query] Starting " << phase << " KNN Queries (k=100)...");
        auto runf_knn = [&]() {
            for(size_t i=0; i<100; i++) {
                auto res = master.knn_report(current_knn[i], 100);
            }
        };
        double t_knn = run_batch_time_loop([](){}, runf_knn, [](){});
        LOG("  [Query] KNN Queries took " << t_knn << " ms");
        LOG("  [Query] Total Phase took " << qt.stop()*1000.0 << " ms");
    };

    run_queries("BaseQuery");

    // Load 2021 Data
    string merged_file = "dataset/bhutan_evolution/2021_merged.txt";
    ifstream fin_year(merged_file);
    vector<Value> P_year_conv;
    while (fin_year >> b_tid >> b_tx >> b_ty) {
        if (b_tx < 0 || b_ty < 0) continue;
        P_year_conv.push_back(make_pair(BoostPoint(b_tx*1000000, b_ty*1000000), b_tid));
    }
    size_t batch_size = P_year_conv.size() / 2;
    vector<vector<Value>> batches;
    batches.push_back(vector<Value>(P_year_conv.begin(), P_year_conv.begin() + batch_size));
    batches.push_back(vector<Value>(P_year_conv.begin() + batch_size, P_year_conv.end()));
    LOG("\n--- YEAR 2021 (FORWARD) ---");
    LOG("Loaded 2021 data, total points: " << P_year_conv.size());

    parlay::internal::timer t_year; t_year.start();
    for (size_t b = 0; b < 2; b++) {
        parlay::internal::timer t; t.start();
        RlogBranch br;
        br.insert_log = batches[b];
        global_history.push_back(br);
        LOG("  [Fork] Batch " << b << " took " << t.stop()*1000.0 << " ms");
    }

    parlay::internal::timer tm; tm.start();
    for (size_t b = 0; b < 2; b++) {
        parlay::internal::timer t; t.start();
        master.merge(global_history[global_history.size() - 2 + b]);
        LOG("  [Merge] Batch " << b << " step took " << t.stop()*1000.0 << " ms");
    }
    LOG("  [Forward] Total Merge took " << tm.stop()*1000.0 << " ms");
    
    run_queries("MidQuery 2021");

    LOG("=== Total Execution Time: " << global_t.stop() << " sec ===");
    return 0;
}
