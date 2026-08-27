#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <queue>
#include <sys/stat.h>
#include <unordered_set>
#include <atomic>

#include "parlay/internal/get_time.h"
#include <cpam/parse_command_line.h>
#include "helper/time_loop.h"

#include "src/mvq.hpp"
#include "src/cpamz.hpp"
#include "src/global_config.hpp"
#include "src/cpambb.hpp"

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/index/rtree.hpp>

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

typedef bgi::rtree<Value, bgi::quadratic<32>, bgi::indexable<Value>, bgi::equal_to<Value>> RTree;


size_t mvzd_nodes_touched = 0;

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
    mutable bool cache_valid = false;
    mutable std::unordered_set<size_t> removed_ids;
    
    void update_cache() const {
        if (!cache_valid) {
            removed_ids.clear();
            removed_ids.reserve(remove_log.size()); // 优化1：防 rehash 扩容
            for (const auto& val : remove_log) removed_ids.insert(val.second);
            cache_valid = true;
        }
    }
    
    RlogTree() : threshold(9999999), current_batches(0) { snapshot = make_shared<RTree>(); }
    RlogTree(int thresh) : threshold(thresh), current_batches(0) { snapshot = make_shared<RTree>(); }
    
    void build_base(const std::vector<Value>& base_data) {
        snapshot = make_shared<RTree>(base_data.begin(), base_data.end());
    }
    void compact() {
        if (insert_log.empty() && remove_log.empty()) return;
        
        std::vector<Value> next_pts;
        next_pts.reserve(snapshot->size() + insert_log.size());
        
        if (!remove_log.empty()) {
            update_cache();
            auto is_alive = [this](Value const& v) { return removed_ids.find(v.second) == removed_ids.end(); };
            
            // 直接在从树里取点的迭代器层级，使用 is_alive 过滤 (std::copy_if)
            std::copy_if(snapshot->begin(), snapshot->end(), std::back_inserter(next_pts), is_alive);
            
            // 同样对 insert_log 进行过滤提取
            std::copy_if(insert_log.begin(), insert_log.end(), std::back_inserter(next_pts), is_alive);
        } else {
            next_pts.assign(snapshot->begin(), snapshot->end());
            next_pts.insert(next_pts.end(), insert_log.begin(), insert_log.end());
        }
        
        snapshot = make_shared<RTree>(next_pts.begin(), next_pts.end());
        
        // 优化2：真正释放内存，防止假性 OOM 和内存指标虚高
        std::vector<Value>().swap(insert_log); 
        std::vector<Value>().swap(remove_log);
        removed_ids = std::unordered_set<size_t>(); // 释放哈希表底层的 bucket 内存
        cache_valid = true; // 空表即为有效状态
    }
    void commit_inserts(const std::vector<Value>& new_pts) {
        insert_log.insert(insert_log.end(), new_pts.begin(), new_pts.end());
        
    }
    void commit_removes(const std::vector<Value>& del_pts) {
        remove_log.insert(remove_log.end(), del_pts.begin(), del_pts.end());
        cache_valid = false;
    }
    void merge(const RlogBranch& branch) {
        // 优化3：只插入数据时，根本不需要废弃关于“死亡名单”的哈希表缓存！
        if (!branch.insert_log.empty()) { insert_log.insert(insert_log.end(), branch.insert_log.begin(), branch.insert_log.end()); }
        if (!branch.remove_log.empty()) { remove_log.insert(remove_log.end(), branch.remove_log.begin(), branch.remove_log.end()); cache_valid = false; }
        current_batches++;
        if (current_batches >= threshold) {
            compact();
            current_batches = 0;
        }
    }
    std::vector<Value> range_report(const geobase::Bounding_Box& q) const {
        update_cache();
        std::vector<Value> result;
        bg::model::box<BoostPoint> box(BoostPoint(q.first.x, q.first.y), BoostPoint(q.second.x, q.second.y));
        auto is_alive = [this](Value const& v) { return removed_ids.find(v.second) == removed_ids.end(); };
        
        // 优雅：直接在树的遍历底层完成过滤，连中间临时数组 snap_res 都省了
        snapshot->query(bgi::intersects(box) && bgi::satisfies(is_alive), std::back_inserter(result));
        
        for (const auto& val : insert_log) {
            if (removed_ids.find(val.second) == removed_ids.end()) {
                // 优化4：修正语义不一致。RTree 用的是 intersects(包含边界), 这里之前用 within(不包含边界) 是有 Bug 的。
                if (bg::intersects(val.first, box)) result.push_back(val);
            }
        }
        return result;
    }
    std::vector<Value> knn_report(const geobase::Point& q, size_t k) const {
        update_cache();
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
        auto is_alive = [this](Value const& v) { return removed_ids.find(v.second) == removed_ids.end(); };
        
        for (auto it = snapshot->qbegin(bgi::nearest(bg_q, (unsigned)k) && bgi::satisfies(is_alive)); it != snapshot->qend(); ++it) {
            double dist = calc_sqr_dist(q, it->first);
            if (max_heap.size() < k) {
                max_heap.push({dist, *it});
            } else if (dist < max_heap.top().first) { 
                max_heap.pop(); max_heap.push({dist, *it}); 
            } else {
                break; // 优化5：提前终止！如果树里找出的点已经比 max_heap 里的点更远，后续的树节点只会更远，直接打断 RTree 遍历！
            }
        }
        
        std::vector<Value> result;
        while (!max_heap.empty()) { result.push_back(max_heap.top().second); max_heap.pop(); }
        return result;
    }
    size_t size() const { return snapshot->size() + insert_log.size() + remove_log.size(); }
};

struct BranchRes { double fork_ms; double commit_ms; double merge_ms; double mem_mb; };








std::unordered_map<size_t, bool> cpambb_mmp;
double cpambb_global_mem = 0;

double mem_cpambb(const CPAMBB::zmap& latest_branch) {
    auto f_noop = [&](const auto &et){ return 0; };
    cpambb_global_mem += latest_branch.size_in_bytes(f_noop, cpambb_mmp) / (1024.0 * 1024.0);
    return cpambb_global_mem;
}

#include <unistd.h>
size_t get_rss_bytes() {
    std::ifstream statm("/proc/self/statm");
    size_t size, resident;
    if (statm >> size >> resident) {
        return resident * sysconf(_SC_PAGESIZE);
    }
    return 0;
}
size_t initial_rss = 0;

double mem_mvzd() { return mvq::global_live_mem.load(std::memory_order_relaxed) / (1024.0 * 1024.0); }

double mem_boost() { 
    size_t current_rss = get_rss_bytes();
    return (current_rss > initial_rss) ? (current_rss - initial_rss) / (1024.0 * 1024.0) : 0.0;
}
double mem_rlog(const RlogTree& master, const vector<RlogBranch>& history) {
    return mem_boost();
}

void log_branch(ofstream& fout, const string& algo, int batch_idx, BranchRes res) {
    fout << std::flush;
    fout << algo << " | " << batch_idx << " | " << res.fork_ms << " | " << res.commit_ms << " | " << res.merge_ms << " | " << res.mem_mb << "\n";
}

struct YearData { int year; vector<parlay::sequence<geobase::Point>> batches; };

template<typename F_Init, typename F_Run, typename F_End>
double run_batch_time_loop(F_Init initf, F_Run runf, F_End endf) {
    return time_loop(3, 1.0, initf, runf, endf) * 1000.0;
}

int main(int argc, char** argv) {
    srand(42);
    cpam::commandLine cmd(argc, argv, "");
    double bp = cmd.getOptionDoubleValue("-bp", 50.0);
    string run_algo = cmd.getOptionValue("-algo", "all");
    
    size_t batches_per_year = max(1ul, (size_t)ceil(100.0 / bp));
    
    std::cout << "[Init] Loading datasets..." << std::endl;
    ifstream fin_base("dataset/bhutan_evolution/base_upto_2020.txt");
    std::vector<geobase::Point> P_base_vec;
    size_t b_tid; double b_tx, b_ty;
    while (fin_base >> b_tid >> b_tx >> b_ty) {
        if (b_tx < 0 || b_ty < 0) continue;
        P_base_vec.push_back(geobase::Point(b_tid, b_tx*1000000, b_ty*1000000));
    }
    parlay::sequence<geobase::Point> P_base(P_base_vec.begin(), P_base_vec.end());

    mvq::Config::get().largest_mbr = geobase::Bounding_Box(geobase::Point(75000000.0, 26000000.0), geobase::Point(95000000.0, 35000000.0));
    auto P_base_set = geobase::get_sorted_points(P_base);

    auto load_q = [](string f) {
        parlay::sequence<geobase::Bounding_Box> queries;
        auto [cnt, q] = geobase::read_range_query(f, 8, mvq::Config::get().maxSize);
        queries = q.substr(0, std::min((size_t)100, q.size()));
        return queries;
    };
    auto queries_small = load_q("dataset/bhutan_evolution/bhutan-small.qry");
    auto queries_med = load_q("dataset/bhutan_evolution/bhutan-med.qry");
    auto queries_large = load_q("dataset/bhutan_evolution/bhutan-large.qry");
    vector<size_t> K_list = {1, 10, 100};

    parlay::internal::timer timer;
    
    vector<shared_ptr<mvq::BaseNode>> mvzd_global_history;
    vector<CPAMBB::zmap> cpambb_global_history;
    vector<shared_ptr<RTree>> boost_global_history;
    vector<vector<RlogBranch>> rlog_global_history(6);
    vector<vector<shared_ptr<RTree>>> rlog_master_history(6);
    
    mvq::Tree* mvzd_tree = nullptr;
    shared_ptr<mvq::BaseNode> mvzd_master;
    CPAMBB::zmap cpambb_master;
    vector<Value> P_base_conv;
    RTree* boost_master = nullptr;
    RlogTree *rlog_master[6];
    
    std::cout << "[Warmup] Warming up allocators..." << std::endl;
    auto sub_pts = P_base_set.substr(0, min((size_t)10000, P_base_set.size()));
    if (run_algo == "MVZD" || run_algo == "all") {
        mvq::Tree dummy(mvq::Config::get().leaf_size);
        dummy.build(sub_pts);
    }
    if (run_algo == "CPAMBB" || run_algo == "all") {
        auto dummy = CPAMBB::map_init(sub_pts, false);
    }
    vector<YearData> all_years;

    for (int year = 2021; year <= 2025; year++) {
        string merged_file = "dataset/bhutan_evolution/" + to_string(year) + "_merged.txt";
        ifstream fin_year(merged_file);
        if (!fin_year.is_open()) continue;
        
        std::vector<geobase::Point> temp_vec;
        size_t tid; double tx, ty;
        while (fin_year >> tid >> tx >> ty) {
            if (tx < 0 || ty < 0) continue;
            temp_vec.push_back(geobase::Point(tid, tx*1000000, ty*1000000));
        }
        parlay::sequence<geobase::Point> P_year(temp_vec.begin(), temp_vec.end());
        
        size_t batch_size = max(1ul, (size_t)ceil(P_year.size() * (bp / 100.0)));
        vector<parlay::sequence<geobase::Point>> batches;
        for (size_t i = 0; i < P_year.size(); i += batch_size) {
            size_t end = min(P_year.size(), i + batch_size);
            batches.push_back(P_year.substr(i, end - i));
        }
        all_years.push_back({year, batches});
    }

    mvq::global_live_mem = 0; 
    std::cout << "[Warmup] Completed." << std::endl;
    initial_rss = get_rss_bytes();

    if (run_algo == "MVZD" || run_algo == "all") {
        mvzd_tree = new mvq::Tree(mvq::Config::get().leaf_size);
        mvzd_tree->build(P_base_set);
        mvzd_master = mvzd_tree->root;
        mvzd_global_history.push_back(mvzd_master); // Initial state
    }
    if (run_algo == "CPAMBB" || run_algo == "all") {
        cpambb_master = CPAMBB::map_init(P_base_set, false);
        cpambb_global_history.push_back(cpambb_master); // Initial state
    }
    if (run_algo == "BoostR" || run_algo.find("Rlog") != string::npos || run_algo == "all") {
        if (P_base_conv.empty()) {
            P_base_conv.resize(P_base.size());
            for(size_t i=0; i<P_base.size(); i++) P_base_conv[i] = make_pair(BoostPoint(P_base[i].x, P_base[i].y), P_base[i].id);
        }
    }
    if (run_algo == "BoostR" || run_algo == "all") {
        boost_master = new RTree(P_base_conv.begin(), P_base_conv.end());
        boost_global_history.push_back(make_shared<RTree>(*boost_master)); // Initial state
    }
    if (run_algo.find("Rlog") != string::npos || run_algo == "all") {
        string names[] = {"Rlog_1yr", "Rlog_2yr", "Rlog_3yr", "Rlog_4yr", "Rlog_5yr", "Rlog_NoSnap"};
        int thresholds[6] = {
            1 * (int)batches_per_year, 2 * (int)batches_per_year, 3 * (int)batches_per_year, 
            4 * (int)batches_per_year, 5 * (int)batches_per_year, 9999999
        };
        for(int i=0; i<6; i++) {
            if (run_algo == names[i] || run_algo == "all") {
                rlog_master[i] = new RlogTree(thresholds[i]);
                rlog_master[i]->build_base(P_base_conv);
                rlog_master_history[i].push_back(rlog_master[i]->snapshot); // Track history
            }
        }
    }
    
    parlay::sequence<geobase::Point> empty_pts;
    mkdir("verification_results", 0777);

    auto run_queries_impl = [&](string prefix, bool is_warmup = false) {
        string fname = "verification_results/" + prefix + "_" + run_algo + ".txt";
        ofstream fout;
        if (!is_warmup) {
            fout.open(fname);
            fout << "Algo | QType | QID | Count | Hash | Time_ms | Mem_MB | Nodes\n";
        }
        
        auto run_range_set = [&](const parlay::sequence<geobase::Bounding_Box>& qs, string q_label) {
            size_t total_res_size = 0;
            if (run_algo == "MVZD" || run_algo == "all") {
                parlay::sequence<size_t> cnts(qs.size(), 0);
                parlay::sequence<parlay::sequence<geobase::Point>> outs(qs.size());
                for(size_t i=0; i<qs.size(); i++) outs[i] = parlay::sequence<geobase::Point>(1000000);
                auto initf = [&]() { ::mvzd_nodes_touched = 0; };
                auto runf = [&]() {
                    for(size_t i=0; i<qs.size(); i++) {
                        cnts[i] = 0;
                        geobase::Bounding_Box q_copy = qs[i];
                        mvzd_tree->range_report(mvzd_master, q_copy, mvq::Config::get().largest_mbr, cnts[i], outs[i]);
                    }
                };
                double total_ms = run_batch_time_loop(initf, runf, [](){});
                double avg_ms = total_ms / qs.size();
                size_t avg_nodes = ::mvzd_nodes_touched / qs.size();
                for(size_t i=0; i<qs.size(); i++) {
                    total_res_size += cnts[i];
                    size_t hash_val = 0;
                    for(size_t j=0; j<cnts[i]; j++) hash_val += outs[i][j].id;
                    if (!is_warmup) fout << "MVZD | " << q_label << " | " << i << " | " << cnts[i] << " | " << hash_val << " | " << avg_ms << " | " << mem_mvzd() << " | " << avg_nodes << "\n";
                }
            }
            if (run_algo == "CPAMBB" || run_algo == "all") {
                parlay::sequence<int64_t> cnts(qs.size(), 0);
                parlay::sequence<parlay::sequence<geobase::Point>> outs(qs.size());
                for(size_t i=0; i<qs.size(); i++) outs[i] = parlay::sequence<geobase::Point>(1000000);
                auto runf = [&]() {
                    for(size_t i=0; i<qs.size(); i++) {
                        geobase::Bounding_Box q_copy = qs[i];
                        cnts[i] = CPAMBB::range_report(cpambb_master, q_copy, outs[i], false);
                    }
                };
                double total_ms = run_batch_time_loop([](){}, runf, [](){});
                double avg_ms = total_ms / qs.size();
                for(size_t i=0; i<qs.size(); i++) {
                    size_t hash_val = 0;
                    for(int64_t j=0; j<cnts[i]; j++) hash_val += outs[i][j].id;
                    if (!is_warmup) fout << "CPAMBB | " << q_label << " | " << i << " | " << cnts[i] << " | " << hash_val << " | " << avg_ms << " | " << mem_cpambb(cpambb_master) << " | 0\n";
                }
            }
            if (run_algo == "BoostR" || run_algo == "all") {
                parlay::sequence<size_t> cnts(qs.size(), 0);
                parlay::sequence<size_t> hashes(qs.size(), 0);
                auto runf = [&]() {
                    for(size_t i=0; i<qs.size(); i++) {
                        std::vector<Value> result;
                        bg::model::box<BoostPoint> box(BoostPoint(qs[i].first.x, qs[i].first.y), BoostPoint(qs[i].second.x, qs[i].second.y));
                        boost_master->query(bgi::intersects(box), std::back_inserter(result));
                        cnts[i] = result.size();
                        size_t h = 0;
                        for(const auto& v : result) h += v.second;
                        hashes[i] = h;
                    }
                };
                double total_ms = run_batch_time_loop([](){}, runf, [](){});
                double avg_ms = total_ms / qs.size();
                for(size_t i=0; i<qs.size(); i++) {
                    if (!is_warmup) fout << "BoostR | " << q_label << " | " << i << " | " << cnts[i] << " | " << hashes[i] << " | " << avg_ms << " | " << mem_boost() << " | 0\n";
                }
            }
            if (run_algo.find("Rlog") != string::npos || run_algo == "all") {
                string names[] = {"Rlog_1yr", "Rlog_2yr", "Rlog_3yr", "Rlog_4yr", "Rlog_5yr", "Rlog_NoSnap"};
                for(int j=0; j<6; j++) {
                    if (run_algo == names[j] || run_algo == "all") {
                        parlay::sequence<size_t> cnts(qs.size(), 0);
                        parlay::sequence<size_t> hashes(qs.size(), 0);
                        auto runf = [&]() {
                            for(size_t i=0; i<qs.size(); i++) {
                                auto res = rlog_master[j]->range_report(qs[i]);
                                cnts[i] = res.size();
                                size_t h = 0;
                                for(const auto& v : res) h += v.second;
                                hashes[i] = h;
                            }
                        };
                        double total_ms = run_batch_time_loop([](){}, runf, [](){});
                        double avg_ms = total_ms / qs.size();
                        for(size_t i=0; i<qs.size(); i++) {
                            if (!is_warmup) fout << names[j] << " | " << q_label << " | " << i << " | " << cnts[i] << " | " << hashes[i] << " | " << avg_ms << " | " << mem_rlog(*rlog_master[j], rlog_global_history[j]) << " | " << rlog_master[j]->size() << "\n";
                        }
                    }
                }
            }
            if (!is_warmup) fout << std::flush;
            if (!is_warmup && (run_algo == "MVZD" || run_algo == "all")) {
                std::cout << "    [" << q_label << "] Avg Result Size: " << (double)total_res_size / qs.size() << " pts/query\n";
            }
        };

        run_range_set(queries_small, "Range_Small");
        run_range_set(queries_med, "Range_Med");
        run_range_set(queries_large, "Range_Large");

        for(size_t K : K_list) {
            parlay::sequence<geobase::Point> current_knn(100);
            for(int i=0; i<100; i++) current_knn[i] = P_base_vec[rand() % P_base_vec.size()];
            string qk = "KNN_" + to_string(K);
            
            if (run_algo == "MVZD" || run_algo == "all") {
                parlay::sequence<size_t> cnts(100, 0);
                parlay::sequence<size_t> hashes(100, 0);
                auto initf = [&]() { ::mvzd_nodes_touched = 0; };
                auto runf = [&]() {
                    for(size_t i=0; i<100; i++) {
                        auto res = mvzd_tree->knn_report(K, current_knn[i], mvq::Config::get().largest_mbr);
                        cnts[i] = res.size();
                        size_t h = 0;
                        while(!res.empty()) { h += res.top().first.id; res.pop(); }
                        hashes[i] = h;
                    }
                };
                double total_ms = run_batch_time_loop(initf, runf, [](){});
                double avg_ms = total_ms / 100.0;
                size_t avg_nodes = ::mvzd_nodes_touched / 100;
                for(size_t i=0; i<100; i++) {
                    if (!is_warmup) fout << "MVZD | " << qk << " | " << i << " | " << cnts[i] << " | " << hashes[i] << " | " << avg_ms << " | " << mem_mvzd() << " | " << avg_nodes << "\n";
                }
            }
            if (run_algo == "CPAMBB" || run_algo == "all") {
                parlay::sequence<size_t> cnts(100, 0);
                parlay::sequence<size_t> hashes(100, 0);
                auto runf = [&]() {
                    for(size_t i=0; i<100; i++) {
                        auto res = CPAMBB::knn(cpambb_master, current_knn[i], K);
                        cnts[i] = res.size();
                        size_t h = 0;
                        while(!res.empty()) { h += res.top().first.id; res.pop(); }
                        hashes[i] = h;
                    }
                };
                double total_ms = run_batch_time_loop([](){}, runf, [](){});
                double avg_ms = total_ms / 100.0;
                for(size_t i=0; i<100; i++) {
                    if (!is_warmup) fout << "CPAMBB | " << qk << " | " << i << " | " << cnts[i] << " | " << hashes[i] << " | " << avg_ms << " | " << mem_cpambb(cpambb_master) << " | 0\n";
                }
            }
            if (run_algo == "BoostR" || run_algo == "all") {
                parlay::sequence<size_t> cnts(100, 0);
                parlay::sequence<size_t> hashes(100, 0);
                auto runf = [&]() {
                    for(size_t i=0; i<100; i++) {
                        std::vector<Value> result;
                        BoostPoint pt(current_knn[i].x, current_knn[i].y);
                        boost_master->query(bgi::nearest(pt, K), std::back_inserter(result));
                        cnts[i] = result.size();
                        size_t h = 0;
                        for(const auto& v : result) h += v.second;
                        hashes[i] = h;
                    }
                };
                double total_ms = run_batch_time_loop([](){}, runf, [](){});
                double avg_ms = total_ms / 100.0;
                for(size_t i=0; i<100; i++) {
                    if (!is_warmup) fout << "BoostR | " << qk << " | " << i << " | " << cnts[i] << " | " << hashes[i] << " | " << avg_ms << " | " << mem_boost() << " | 0\n";
                }
            }
            if (run_algo.find("Rlog") != string::npos || run_algo == "all") {
                string names[] = {"Rlog_1yr", "Rlog_2yr", "Rlog_3yr", "Rlog_4yr", "Rlog_5yr", "Rlog_NoSnap"};
                for(int j=0; j<6; j++) {
                    if (run_algo == names[j] || run_algo == "all") {
                        parlay::sequence<size_t> cnts(100, 0);
                        parlay::sequence<size_t> hashes(100, 0);
                        auto runf = [&]() {
                            for(size_t i=0; i<100; i++) {
                                auto res = rlog_master[j]->knn_report(current_knn[i], K);
                                cnts[i] = res.size();
                                size_t h = 0;
                                for(const auto& v : res) h += v.second;
                                hashes[i] = h;
                            }
                        };
                        double total_ms = run_batch_time_loop([](){}, runf, [](){});
                        double avg_ms = total_ms / 100.0;
                        for(size_t i=0; i<100; i++) {
                            if (!is_warmup) fout << names[j] << " | " << qk << " | " << i << " | " << cnts[i] << " | " << hashes[i] << " | " << avg_ms << " | " << mem_rlog(*rlog_master[j], rlog_global_history[j]) << " | " << rlog_master[j]->size() << "\n";
                        }
                    }
                }
            }
            if (!is_warmup) fout << std::flush;
        }
    };

    run_queries_impl("BaseQuery", true);
    run_queries_impl("BaseQuery", false);
    cout << "    [BaseQuery] completed\n";


    for (size_t y_idx = 0; y_idx < all_years.size(); y_idx++) {
        auto& yd = all_years[y_idx];
        int year = yd.year;
        auto& batches = yd.batches;
        cout << "--- YEAR " << year << " (FORWARD) ---" << endl;
        
        ofstream fwd_out("verification_results/" + to_string(year) + "_Forward_Batch_" + run_algo + ".txt");
        fwd_out << "Algo | Batch | Fork_ms | Commit_ms | Merge_ms | Mem_MB\n";

        vector<shared_ptr<mvq::BaseNode>> mvzd_branches(batches.size());
        vector<CPAMBB::zmap> cpambb_branches(batches.size());

        for (size_t b = 0; b < batches.size(); b++) {
            auto& batch = batches[b];
            
            if (run_algo == "MVZD" || run_algo == "all") {
                timer.start(); auto v_mvzd = mvzd_master; double tf = timer.stop()*1000;
                timer.start(); auto cv_mvzd = mvzd_tree->commit(v_mvzd, batch, empty_pts); double tc = timer.stop()*1000;
                mvzd_branches[b] = cv_mvzd;
                mvzd_global_history.push_back(cv_mvzd); 
                log_branch(fwd_out, "MVZD", b, {tf, tc, 0, mem_mvzd()});
            }
            if (run_algo == "CPAMBB" || run_algo == "all") {
                timer.start(); auto v_cpambb = cpambb_master; double tf = timer.stop()*1000;
                timer.start(); auto cv_cpambb = CPAMBB::map_insert(batch, v_cpambb); double tc = timer.stop()*1000;
                cpambb_branches[b] = cv_cpambb;
                cpambb_global_history.push_back(cv_cpambb);
                log_branch(fwd_out, "CPAMBB", b, {tf, tc, 0, mem_cpambb(cv_cpambb)});
            }
            if (run_algo == "BoostR" || run_algo == "all") {
                timer.start(); RTree v_boost = *boost_master; double tf = timer.stop()*1000;
                vector<Value> P_branch_conv(batch.size());
                for(size_t k=0; k<batch.size(); k++) P_branch_conv[k] = make_pair(BoostPoint(batch[k].x, batch[k].y), batch[k].id);
                timer.start(); v_boost.insert(P_branch_conv.begin(), P_branch_conv.end()); double tc = timer.stop()*1000;
                boost_global_history.push_back(make_shared<RTree>(v_boost));
                log_branch(fwd_out, "BoostR", b, {tf, tc, 0, mem_boost()});
            }
            if (run_algo.find("Rlog") != string::npos || run_algo == "all") {
                vector<Value> P_branch_conv(batch.size());
                for(size_t k=0; k<batch.size(); k++) P_branch_conv[k] = make_pair(BoostPoint(batch[k].x, batch[k].y), batch[k].id);
                string names[] = {"Rlog_1yr", "Rlog_2yr", "Rlog_3yr", "Rlog_4yr", "Rlog_5yr", "Rlog_NoSnap"};
                for(int j=0; j<6; j++) {
                    if (run_algo == names[j] || run_algo == "all") {
                        timer.start(); RlogBranch v_rlog; double l_tf = timer.stop()*1000;
                        timer.start(); v_rlog.insert_log = P_branch_conv; double l_tc = timer.stop()*1000;
                        rlog_global_history[j].push_back(v_rlog);
                        log_branch(fwd_out, names[j], b, {l_tf, l_tc, 0, mem_rlog(*rlog_master[j], rlog_global_history[j])});
                    }
                }
            }
        }

        if (run_algo == "MVZD" || run_algo == "all") {
            timer.start();
            auto original_mvzd_master = mvzd_master;
            parlay::sequence<geobase::Point> total_add, total_remove;
            for (size_t b = 0; b < batches.size(); b++) {
                auto [add, remove] = mvzd_tree->two_version_diff(original_mvzd_master, mvzd_branches[b]);
                auto [insert_pts, delete_pts, update_pts] = geobase::filter_diff_results(add, remove);
                total_add.append(insert_pts);
                total_remove.append(delete_pts);
            }
            mvzd_master = mvzd_tree->commit(original_mvzd_master, total_add, total_remove);
            mvzd_global_history.push_back(mvzd_master);
            double tm = timer.stop()*1000;
            fwd_out << "MVZD | SUMMARY | 0 | 0 | " << tm << " | " << mem_mvzd() << "\n";
            cout << "    [Forward] Merge MVZD took " << tm << " ms\n";
        }
        if (run_algo == "CPAMBB" || run_algo == "all") {
            timer.start();
            auto original_cpambb_master = cpambb_master;
            parlay::sequence<geobase::Point> total_add, total_remove;
            for (size_t b = 0; b < batches.size(); b++) {
                auto [add, remove] = CPAMBB::map_diff(original_cpambb_master, cpambb_branches[b]);
                auto [insert_pts, delete_pts, update_pts] = geobase::filter_diff_results(add, remove);
                total_add.append(insert_pts);
                total_remove.append(delete_pts);
            }
            cpambb_master = CPAMBB::map_commit(original_cpambb_master, total_add, total_remove);
            cpambb_global_history.push_back(cpambb_master);
            double tm = timer.stop()*1000;
            fwd_out << "CPAMBB | SUMMARY | 0 | 0 | " << tm << " | " << mem_cpambb(cpambb_master) << "\n";
            cout << "    [Forward] Merge CPAMBB took " << tm << " ms\n";
        }
        if (run_algo == "BoostR" || run_algo == "all") {
            timer.start();
            for (size_t b = 0; b < batches.size(); b++) {
                vector<Value> P_branch_conv(batches[b].size());
                for(size_t k=0; k<batches[b].size(); k++) P_branch_conv[k] = make_pair(BoostPoint(batches[b][k].x, batches[b][k].y), batches[b][k].id);
                boost_master->insert(P_branch_conv.begin(), P_branch_conv.end());
            }
            boost_global_history.push_back(make_shared<RTree>(*boost_master));
            double tm = timer.stop()*1000;
            fwd_out << "BoostR | SUMMARY | 0 | 0 | " << tm << " | " << mem_boost() << "\n";
            cout << "    [Forward] Merge BoostR took " << tm << " ms\n";
        }
        if (run_algo.find("Rlog") != string::npos || run_algo == "all") {
            string names[] = {"Rlog_1yr", "Rlog_2yr", "Rlog_3yr", "Rlog_4yr", "Rlog_5yr", "Rlog_NoSnap"};
            for(int j=0; j<6; j++) {
                if (run_algo == names[j] || run_algo == "all") {
                    timer.start();
                    for (size_t b = 0; b < batches.size(); b++) {
                        rlog_master[j]->merge(rlog_global_history[j][rlog_global_history[j].size() - batches.size() + b]);
                    }
                    rlog_master_history[j].push_back(rlog_master[j]->snapshot);
                    double tm = timer.stop()*1000;
                    fwd_out << names[j] << " | SUMMARY | 0 | 0 | " << tm << " | " << mem_rlog(*rlog_master[j], rlog_global_history[j]) << "\n";
                    cout << "    [Forward] Merge " << names[j] << " took " << tm << " ms\n";
                }
            }
        }
        
        run_queries_impl(to_string(year) + "_MidQuery", true);
        run_queries_impl(to_string(year) + "_MidQuery", false);
        cout << "    [MidQuery] completed\n";
    }

    for (int y_idx = (int)all_years.size() - 1; y_idx >= 0; y_idx--) {
        auto& yd = all_years[y_idx];
        int year = yd.year;
        auto& batches = yd.batches;
        cout << "--- YEAR " << year << " (BACKWARD) ---" << endl;
        
        ofstream bwd_out("verification_results/" + to_string(year) + "_Backward_Batch_" + run_algo + ".txt");
        bwd_out << "Algo | Batch | Fork_ms | Commit_ms | Merge_ms | Mem_MB\n";

        vector<shared_ptr<mvq::BaseNode>> mvzd_branches(batches.size());
        vector<CPAMBB::zmap> cpambb_branches(batches.size());

        for (int b = (int)batches.size() - 1; b >= 0; b--) {
            auto& batch = batches[b];
            
            if (run_algo == "MVZD" || run_algo == "all") {
                timer.start(); auto v_mvzd = mvzd_master; double tf = timer.stop()*1000;
                timer.start(); auto cv_mvzd = mvzd_tree->commit(v_mvzd, empty_pts, batch); double tc = timer.stop()*1000;
                mvzd_branches[b] = cv_mvzd;
                mvzd_global_history.push_back(cv_mvzd);
                log_branch(bwd_out, "MVZD", b, {tf, tc, 0, mem_mvzd()});
            }
            if (run_algo == "CPAMBB" || run_algo == "all") {
                timer.start(); auto v_cpambb = cpambb_master; double tf = timer.stop()*1000;
                timer.start(); auto cv_cpambb = CPAMBB::map_delete(batch, v_cpambb); double tc = timer.stop()*1000;
                cpambb_branches[b] = cv_cpambb;
                cpambb_global_history.push_back(cv_cpambb);
                log_branch(bwd_out, "CPAMBB", b, {tf, tc, 0, mem_cpambb(cv_cpambb)});
            }
            if (run_algo == "BoostR" || run_algo == "all") {
                timer.start(); RTree v_boost = *boost_master; double tf = timer.stop()*1000;
                vector<Value> P_branch_conv(batch.size());
                for(size_t k=0; k<batch.size(); k++) P_branch_conv[k] = make_pair(BoostPoint(batch[k].x, batch[k].y), batch[k].id);
                timer.start(); v_boost.remove(P_branch_conv.begin(), P_branch_conv.end()); double tc = timer.stop()*1000;
                boost_global_history.push_back(make_shared<RTree>(v_boost));
                log_branch(bwd_out, "BoostR", b, {tf, tc, 0, mem_boost()});
            }
            if (run_algo.find("Rlog") != string::npos || run_algo == "all") {
                vector<Value> P_branch_conv(batch.size());
                for(size_t k=0; k<batch.size(); k++) P_branch_conv[k] = make_pair(BoostPoint(batch[k].x, batch[k].y), batch[k].id);
                string names[] = {"Rlog_1yr", "Rlog_2yr", "Rlog_3yr", "Rlog_4yr", "Rlog_5yr", "Rlog_NoSnap"};
                for(int j=0; j<6; j++) {
                    if (run_algo == names[j] || run_algo == "all") {
                        timer.start(); RlogBranch v_rlog; double l_tf = timer.stop()*1000;
                        timer.start(); v_rlog.remove_log = P_branch_conv; double l_tc = timer.stop()*1000;
                        rlog_global_history[j].push_back(v_rlog);
                        log_branch(bwd_out, names[j], b, {l_tf, l_tc, 0, mem_rlog(*rlog_master[j], rlog_global_history[j])});
                    }
                }
            }
        }

        if (run_algo == "MVZD" || run_algo == "all") {
            timer.start();
            auto original_mvzd_master = mvzd_master;
            parlay::sequence<geobase::Point> total_add, total_remove;
            for (int b = (int)batches.size() - 1; b >= 0; b--) {
                auto [add, remove] = mvzd_tree->two_version_diff(original_mvzd_master, mvzd_branches[b]);
                auto [insert_pts, delete_pts, update_pts] = geobase::filter_diff_results(add, remove);
                total_add.append(insert_pts);
                total_remove.append(delete_pts);
            }
            mvzd_master = mvzd_tree->commit(original_mvzd_master, total_add, total_remove);
            mvzd_global_history.push_back(mvzd_master);
            double tm = timer.stop()*1000;
            bwd_out << "MVZD | SUMMARY | 0 | 0 | " << tm << " | " << mem_mvzd() << "\n";
            cout << "    [Backward] Merge MVZD took " << tm << " ms\n";
        }
        if (run_algo == "CPAMBB" || run_algo == "all") {
            timer.start();
            auto original_cpambb_master = cpambb_master;
            parlay::sequence<geobase::Point> total_add, total_remove;
            for (int b = (int)batches.size() - 1; b >= 0; b--) {
                auto [add, remove] = CPAMBB::map_diff(original_cpambb_master, cpambb_branches[b]);
                auto [insert_pts, delete_pts, update_pts] = geobase::filter_diff_results(add, remove);
                total_add.append(insert_pts);
                total_remove.append(delete_pts);
            }
            cpambb_master = CPAMBB::map_commit(original_cpambb_master, total_add, total_remove);
            cpambb_global_history.push_back(cpambb_master);
            double tm = timer.stop()*1000;
            bwd_out << "CPAMBB | SUMMARY | 0 | 0 | " << tm << " | " << mem_cpambb(cpambb_master) << "\n";
            cout << "    [Backward] Merge CPAMBB took " << tm << " ms\n";
        }
        if (run_algo == "BoostR" || run_algo == "all") {
            timer.start();
            for (int b = (int)batches.size() - 1; b >= 0; b--) {
                vector<Value> P_branch_conv(batches[b].size());
                for(size_t k=0; k<batches[b].size(); k++) P_branch_conv[k] = make_pair(BoostPoint(batches[b][k].x, batches[b][k].y), batches[b][k].id);
                boost_master->remove(P_branch_conv.begin(), P_branch_conv.end());
            }
            boost_global_history.push_back(make_shared<RTree>(*boost_master));
            double tm = timer.stop()*1000;
            bwd_out << "BoostR | SUMMARY | 0 | 0 | " << tm << " | " << mem_boost() << "\n";
            cout << "    [Backward] Merge BoostR took " << tm << " ms\n";
        }
        if (run_algo.find("Rlog") != string::npos || run_algo == "all") {
            string names[] = {"Rlog_1yr", "Rlog_2yr", "Rlog_3yr", "Rlog_4yr", "Rlog_5yr", "Rlog_NoSnap"};
            for(int j=0; j<6; j++) {
                if (run_algo == names[j] || run_algo == "all") {
                    timer.start();
                    for (int b = (int)batches.size() - 1; b >= 0; b--) {
                        rlog_master[j]->merge(rlog_global_history[j][rlog_global_history[j].size() - batches.size() + ((batches.size()-1)-b)]);
                    }
                    rlog_master_history[j].push_back(rlog_master[j]->snapshot);
                    double tm = timer.stop()*1000;
                    bwd_out << names[j] << " | SUMMARY | 0 | 0 | " << tm << " | " << mem_rlog(*rlog_master[j], rlog_global_history[j]) << "\n";
                    cout << "    [Backward] Merge " << names[j] << " took " << tm << " ms\n";
                }
            }
        }
        
        run_queries_impl(to_string(year) + "_EndQuery", true);
        run_queries_impl(to_string(year) + "_EndQuery", false);
        cout << "    [EndQuery] completed\n";
    }

    return 0;
}
