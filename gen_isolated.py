import os

cpp_code = r"""
#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <sys/stat.h>
#include <unistd.h>
#include <unordered_set>
#include <boost/geometry.hpp>
#include <boost/geometry/index/rtree.hpp>
#include <numeric>
#include <cmath>
#include "geo/point.hpp"
#include "geo/operations.hpp"
#include "geo/io.hpp"
#include "mvq.hpp"
#include "cpambb.hpp"
#include "parlay/internal/get_time.h"
#include <cpam/parse_command_line.h>

using namespace std;
namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;
typedef bg::model::point<double, 2, bg::cs::cartesian> Point;
typedef pair<Point, size_t> Value;
typedef bgi::rtree<Value, bgi::quadratic<32>> RTree;

class RlogTree {
public:
    shared_ptr<RTree> snapshot;
    vector<Value> insert_log;
    vector<Value> remove_log;
    int version_depth;
    int threshold;

    RlogTree(int t = 4) : version_depth(0), threshold(t) { snapshot = make_shared<RTree>(); }
    void build_base(const vector<Value>& pts) {
        snapshot = make_shared<RTree>(pts.begin(), pts.end());
        insert_log.clear(); remove_log.clear(); version_depth = 0;
    }
    RlogTree(const RlogTree& other) {
        snapshot = other.snapshot; insert_log = other.insert_log; remove_log = other.remove_log;
        version_depth = other.version_depth; threshold = other.threshold;
    }
    RlogTree& operator=(const RlogTree& other) {
        if (this != &other) {
            snapshot = other.snapshot; insert_log = other.insert_log; remove_log = other.remove_log;
            version_depth = other.version_depth; threshold = other.threshold;
        }
        return *this;
    }
    void commit_inserts(const vector<Value>& pts) { insert_log.insert(insert_log.end(), pts.begin(), pts.end()); }
    void commit_removes(const vector<Value>& pts) { remove_log.insert(remove_log.end(), pts.begin(), pts.end()); }
    void compact() {
        auto new_snap = make_shared<RTree>(*snapshot);
        for(const auto& val : remove_log) {
            new_snap->remove(val);
        }
        new_snap->insert(insert_log.begin(), insert_log.end());
        snapshot = new_snap; insert_log.clear(); remove_log.clear(); version_depth = 0;
    }
    void merge(const RlogTree& branch, size_t original_log_size, size_t original_remove_size = 0) {
        size_t new_inserts = branch.insert_log.size() - original_log_size;
        if (new_inserts > 0) {
            insert_log.insert(insert_log.end(), branch.insert_log.begin() + original_log_size, branch.insert_log.end());
        }
        size_t new_removes = branch.remove_log.size() - original_remove_size;
        if (new_removes > 0) {
            remove_log.insert(remove_log.end(), branch.remove_log.begin() + original_remove_size, branch.remove_log.end());
        }
        version_depth++;
        if (version_depth >= threshold) compact();
    }
    struct MaxHeapCmp {
        bool operator()(const std::pair<double, Value>& a, const std::pair<double, Value>& b) const {
            return a.first < b.first;
        }
    };

    vector<Value> range_report(const geobase::Bounding_Box& q) const {
        std::unordered_set<size_t> removed_ids;
        for (const auto& val : remove_log) removed_ids.insert(val.second);

        vector<Value> snap_res;
        bg::model::box<Point> box(Point(q.first.x, q.first.y), Point(q.second.x, q.second.y));
        snapshot->query(bgi::intersects(box), std::back_inserter(snap_res));
        
        vector<Value> result;
        result.reserve(snap_res.size() + insert_log.size());
        
        for (const auto& val : snap_res) {
            if (removed_ids.find(val.second) == removed_ids.end()) {
                result.push_back(val);
            }
        }
        for (const auto& val : insert_log) {
            if (removed_ids.find(val.second) == removed_ids.end()) {
                if (bg::within(val.first, box)) result.push_back(val);
            }
        }
        return result;
    }

    vector<Value> knn_report(const geobase::Point& q, size_t k) const {
        std::unordered_set<size_t> removed_ids;
        for (const auto& val : remove_log) removed_ids.insert(val.second);

        auto calc_sqr_dist = [](const geobase::Point& p1, const Point& p2) {
            double dx = p1.x - p2.get<0>();
            double dy = p1.y - p2.get<1>();
            return dx * dx + dy * dy;
        };

        std::priority_queue<std::pair<double, Value>, std::vector<std::pair<double, Value>>, MaxHeapCmp> max_heap;

        for (const auto& val : insert_log) {
            if (removed_ids.find(val.second) == removed_ids.end()) {
                double dist = calc_sqr_dist(q, val.first);
                max_heap.push({dist, val});
                if (max_heap.size() > k) max_heap.pop();
            }
        }

        Point bg_q(q.x, q.y);
        for (auto it = snapshot->qbegin(bgi::nearest(bg_q, snapshot->size())); it != snapshot->qend(); ++it) {
            if (removed_ids.find(it->second) != removed_ids.end()) continue;
            
            double dist = calc_sqr_dist(q, it->first);
            if (max_heap.size() == k && dist >= max_heap.top().first) break;
            
            max_heap.push({dist, *it});
            if (max_heap.size() > k) max_heap.pop();
        }

        vector<Value> result;
        result.reserve(max_heap.size());
        while (!max_heap.empty()) {
            result.push_back(max_heap.top().second);
            max_heap.pop();
        }
        std::reverse(result.begin(), result.end());
        return result;
    }
    size_t size() const { return snapshot->size() + insert_log.size() + remove_log.size(); }
};

struct QRes { size_t count; size_t hash_val; double time_ms; double mem_mb; };
struct BranchRes { double fork_ms; double commit_ms; double merge_ms; double mem_mb; };

// Logical Memory Tracking Functions
double mem_mvzd(mvq::Tree& tree) {
    auto stat = tree.get_tree_statistics();
    return (stat.mem_inte_nodes + stat.mem_leaf_nodes) / (1024.0 * 1024.0);
}
double mem_cpambb() {
    size_t inte_used = CPAMBB::zmap::GC::used_node();
    size_t internal_nodes_space = sizeof(CPAMBB::zmap::GC::regular_node) * inte_used;
    return (internal_nodes_space) / (1024.0 * 1024.0);
}
double mem_boost(const RTree& master, size_t acc_branch_elements) {
    return ((master.size() + acc_branch_elements) * sizeof(Value)) / (1024.0 * 1024.0);
}
double mem_rlog(const RlogTree& master, size_t acc_branch_elements) {
    return ((master.snapshot->size() + master.insert_log.size() + master.remove_log.size() + acc_branch_elements) * sizeof(Value)) / (1024.0 * 1024.0);
}

// Queries
QRes time_query_mvzd(mvq::Tree& tree, auto m_cur, geobase::Bounding_Box q) {
    parlay::internal::timer t; t.start();
    size_t cnt=0; parlay::sequence<geobase::Point> out(1000000);
    tree.range_report(m_cur, q, mvq::Config::get().largest_mbr, cnt, out);
    double ms = t.stop() * 1000.0;
    size_t hash_val = 0;
    for(size_t i=0; i<cnt; i++) hash_val += out[i].id;
    return {cnt, hash_val, ms, mem_mvzd(tree)};
}

QRes time_query_cpambb(auto m_cur, geobase::Bounding_Box q) {
    parlay::internal::timer t; t.start();
    parlay::sequence<geobase::Point> out(1000000);
    int64_t cnt = CPAMBB::range_report(m_cur, q, out, false);
    double ms = t.stop() * 1000.0;
    size_t hash_val = 0;
    for(int64_t i=0; i<cnt; i++) hash_val += out[i].id;
    return {(size_t)cnt, hash_val, ms, mem_cpambb()}; 
}

QRes time_query_boost(RTree& tree, geobase::Bounding_Box q, size_t acc_elements) {
    parlay::internal::timer t; t.start();
    std::vector<Value> result;
    bg::model::box<Point> box(Point(q.first.x, q.first.y), Point(q.second.x, q.second.y));
    tree.query(bgi::intersects(box), std::back_inserter(result));
    double ms = t.stop() * 1000.0;
    size_t hash_val = 0;
    for(const auto& v : result) hash_val += v.second;
    return {result.size(), hash_val, ms, mem_boost(tree, acc_elements)};
}

QRes time_query_rlog(RlogTree& tree, geobase::Bounding_Box q, size_t acc_elements) {
    parlay::internal::timer t; t.start();
    auto result = tree.range_report(q);
    double ms = t.stop() * 1000.0;
    size_t hash_val = 0;
    for(const auto& v : result) hash_val += v.second;
    return {result.size(), hash_val, ms, mem_rlog(tree, acc_elements)};
}

// KNN Queries
QRes time_knn_mvzd(mvq::Tree& tree, auto m_cur, geobase::Point q, size_t k) {
    parlay::internal::timer t; t.start();
    tree.root = m_cur;
    auto res = tree.knn_report(k, q, mvq::Config::get().largest_mbr);
    double ms = t.stop() * 1000.0;
    size_t cnt = res.size();
    size_t hash_val = 0;
    while(!res.empty()) { hash_val += res.top().first.id; res.pop(); }
    return {cnt, hash_val, ms, mem_mvzd(tree)};
}

QRes time_knn_cpambb(auto m_cur, geobase::Point q, size_t k) {
    parlay::internal::timer t; t.start();
    auto res = CPAMBB::knn(m_cur, q, k);
    double ms = t.stop() * 1000.0;
    size_t cnt = res.size();
    size_t hash_val = 0;
    while(!res.empty()) { hash_val += res.top().first.id; res.pop(); }
    return {cnt, hash_val, ms, mem_cpambb()};
}

QRes time_knn_boost(RTree& tree, geobase::Point q, size_t k, size_t acc_elements) {
    parlay::internal::timer t; t.start();
    std::vector<Value> result;
    Point pt(q.x, q.y);
    tree.query(bgi::nearest(pt, k), std::back_inserter(result));
    double ms = t.stop() * 1000.0;
    size_t hash_val = 0;
    for(const auto& v : result) hash_val += v.second;
    return {result.size(), hash_val, ms, mem_boost(tree, acc_elements)};
}

QRes time_knn_rlog(RlogTree& tree, geobase::Point q, size_t k, size_t acc_elements) {
    parlay::internal::timer t; t.start();
    auto result = tree.knn_report(q, k);
    double ms = t.stop() * 1000.0;
    size_t hash_val = 0;
    for(const auto& v : result) hash_val += v.second;
    return {result.size(), hash_val, ms, mem_rlog(tree, acc_elements)};
}

void log_branch(ofstream& fout, const string& algo, int batch_idx, BranchRes res) {
    fout << std::flush;
    fout << algo << " | " << batch_idx << " | " << res.fork_ms << " | " << res.commit_ms << " | " << res.merge_ms << " | " << res.mem_mb << "\n";
}

struct YearData { int year; vector<parlay::sequence<geobase::Point>> batches; };

int main(int argc, char** argv) {
    cpam::commandLine cmd(argc, argv, "");
    double bp = cmd.getOptionDoubleValue("-bp", 10.0);
    string run_algo = cmd.getOptionValue("-algo", "all");
    
    size_t batches_per_year = max(1ul, (size_t)ceil(100.0 / bp));
    
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
        queries = q.substr(0, std::min((size_t)10, q.size()));
        return queries;
    };
    auto queries_small = load_q("dataset/bhutan_evolution/bhutan-small.qry");
    auto queries_med = load_q("dataset/bhutan_evolution/bhutan-med.qry");
    auto queries_large = load_q("dataset/bhutan_evolution/bhutan-large.qry");
    auto queries_knn = queries_med;
    vector<size_t> K_list = {1, 10, 100};

    parlay::internal::timer timer;
    
    size_t boost_acc_elements = 0;
    vector<size_t> rlog_acc_logs(6, 0);
    
    mvq::Tree* mvzd_tree = nullptr;
    shared_ptr<mvq::BaseNode> mvzd_master;
    
    CPAMBB::zmap cpambb_master;
    
    vector<Value> P_base_conv;
    RTree* boost_master = nullptr;
    
    RlogTree *rlog_master[6];
    
    if (run_algo == "MVZD" || run_algo == "all") {
        mvzd_tree = new mvq::Tree(mvq::Config::get().leaf_size);
        mvzd_tree->build(P_base_set);
        mvzd_master = mvzd_tree->root;
        mvzd_tree->multi_version_roots.push_back(mvzd_master);
    }
    if (run_algo == "CPAMBB" || run_algo == "all") {
        cpambb_master = CPAMBB::map_init(P_base_set, false);
    }
    if (run_algo == "BoostR" || run_algo.find("Rlog") != string::npos || run_algo == "all") {
        P_base_conv.resize(P_base.size());
        for(size_t i=0; i<P_base.size(); i++) P_base_conv[i] = make_pair(Point(P_base[i].x, P_base[i].y), P_base[i].id);
    }
    if (run_algo == "BoostR" || run_algo == "all") {
        boost_master = new RTree(P_base_conv.begin(), P_base_conv.end());
    }
    if (run_algo.find("Rlog") != string::npos || run_algo == "all") {
        int thresholds[6] = {
            1 * (int)batches_per_year, 2 * (int)batches_per_year, 3 * (int)batches_per_year, 
            4 * (int)batches_per_year, 5 * (int)batches_per_year, 9999999
        };
        for(int i=0; i<6; i++) {
            rlog_master[i] = new RlogTree(thresholds[i]);
            rlog_master[i]->build_base(P_base_conv);
        }
    }
    
    parlay::sequence<geobase::Point> empty_pts;
    mkdir("verification_results", 0777);

    auto run_queries_impl = [&](string prefix) {
        string fname = "verification_results/" + prefix + "_" + run_algo + ".txt";
        ofstream fout(fname);
        fout << "Algo | QType | QID | Count | Hash | Time_ms | Mem_MB\n";
        
        auto run_range_set = [&](const parlay::sequence<geobase::Bounding_Box>& qs, string q_label) {
            for(size_t i=0; i<qs.size(); i++) {
                auto q = qs[i];
                if (run_algo == "MVZD" || run_algo == "all") {
                    auto r = time_query_mvzd(*mvzd_tree, mvzd_master, q);
                    fout << "MVZD | " << q_label << " | " << i << " | " << r.count << " | " << r.hash_val << " | " << r.time_ms << " | " << r.mem_mb << "\n";
                }
                if (run_algo == "CPAMBB" || run_algo == "all") {
                    auto r = time_query_cpambb(cpambb_master, q);
                    fout << "CPAMBB | " << q_label << " | " << i << " | " << r.count << " | " << r.hash_val << " | " << r.time_ms << " | " << r.mem_mb << "\n";
                }
                if (run_algo == "BoostR" || run_algo == "all") {
                    auto r = time_query_boost(*boost_master, q, boost_acc_elements);
                    fout << "BoostR | " << q_label << " | " << i << " | " << r.count << " | " << r.hash_val << " | " << r.time_ms << " | " << r.mem_mb << "\n";
                }
                if (run_algo.find("Rlog") != string::npos || run_algo == "all") {
                    string names[] = {"Rlog_1yr", "Rlog_2yr", "Rlog_3yr", "Rlog_4yr", "Rlog_5yr", "Rlog_NoSnap"};
                    for(int j=0; j<6; j++) {
                        if (run_algo == names[j] || run_algo == "all") {
                            auto r = time_query_rlog(*rlog_master[j], q, rlog_acc_logs[j]);
                            fout << names[j] << " | " << q_label << " | " << i << " | " << r.count << " | " << r.hash_val << " | " << r.time_ms << " | " << r.mem_mb << "\n";
                        }
                    }
                }
                fout << std::flush;
            }
        };

        run_range_set(queries_small, "Range_Small");
        run_range_set(queries_med, "Range_Med");
        run_range_set(queries_large, "Range_Large");

        for(size_t K : K_list) {
            for(size_t i=0; i<queries_knn.size(); i++) {
                geobase::Point qp = queries_knn[i].first;
                string qk = "KNN_" + to_string(K);
                if (run_algo == "MVZD" || run_algo == "all") {
                    auto r = time_knn_mvzd(*mvzd_tree, mvzd_master, qp, K);
                    fout << "MVZD | " << qk << " | " << i << " | " << r.count << " | " << r.hash_val << " | " << r.time_ms << " | " << r.mem_mb << "\n";
                }
                if (run_algo == "CPAMBB" || run_algo == "all") {
                    auto r = time_knn_cpambb(cpambb_master, qp, K);
                    fout << "CPAMBB | " << qk << " | " << i << " | " << r.count << " | " << r.hash_val << " | " << r.time_ms << " | " << r.mem_mb << "\n";
                }
                if (run_algo == "BoostR" || run_algo == "all") {
                    auto r = time_knn_boost(*boost_master, qp, K, boost_acc_elements);
                    fout << "BoostR | " << qk << " | " << i << " | " << r.count << " | " << r.hash_val << " | " << r.time_ms << " | " << r.mem_mb << "\n";
                }
                if (run_algo.find("Rlog") != string::npos || run_algo == "all") {
                    string names[] = {"Rlog_1yr", "Rlog_2yr", "Rlog_3yr", "Rlog_4yr", "Rlog_5yr", "Rlog_NoSnap"};
                    for(int j=0; j<6; j++) {
                        if (run_algo == names[j] || run_algo == "all") {
                            auto r = time_knn_rlog(*rlog_master[j], qp, K, rlog_acc_logs[j]);
                            fout << names[j] << " | " << qk << " | " << i << " | " << r.count << " | " << r.hash_val << " | " << r.time_ms << " | " << r.mem_mb << "\n";
                        }
                    }
                }
                fout << std::flush;
            }
        }
    };

    run_queries_impl("BaseQuery");
    cout << "    [BaseQuery] completed\n";

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

    // ---------------------------------------------------------
    // Phase 1: Forward Loop (All Years)
    // STAR TOPOLOGY (Concurrent Forks, Merge at end of year)
    // ---------------------------------------------------------
    for (size_t y_idx = 0; y_idx < all_years.size(); y_idx++) {
        auto& yd = all_years[y_idx];
        int year = yd.year;
        auto& batches = yd.batches;
        cout << "--- YEAR " << year << " (FORWARD) ---" << endl;
        
        ofstream fwd_out("verification_results/" + to_string(year) + "_Forward_Batch_" + run_algo + ".txt");
        fwd_out << "Algo | Batch | Fork_ms | Commit_ms | Merge_ms | Mem_MB\n";

        vector<shared_ptr<mvq::BaseNode>> mvzd_branches(batches.size());
        vector<CPAMBB::zmap> cpambb_branches(batches.size());
        vector<RTree> boost_branches(batches.size());
        vector<vector<RlogTree>> rlog_branches(6, vector<RlogTree>(batches.size()));

        // Step 1: Independent Fork & Commit
        for (size_t b = 0; b < batches.size(); b++) {
            auto& batch = batches[b];
            
            if (run_algo == "MVZD" || run_algo == "all") {
                timer.start(); auto v_mvzd = mvzd_master; double tf = timer.stop()*1000;
                timer.start(); auto cv_mvzd = mvzd_tree->commit(v_mvzd, batch, empty_pts); double tc = timer.stop()*1000;
                mvzd_branches[b] = cv_mvzd;
                mvzd_tree->multi_version_roots.push_back(cv_mvzd); 
                log_branch(fwd_out, "MVZD", b, {tf, tc, 0, mem_mvzd(*mvzd_tree)});
            }
            if (run_algo == "CPAMBB" || run_algo == "all") {
                timer.start(); auto v_cpambb = cpambb_master; double tf = timer.stop()*1000;
                timer.start(); auto cv_cpambb = CPAMBB::map_insert(batch, v_cpambb); double tc = timer.stop()*1000;
                cpambb_branches[b] = cv_cpambb;
                log_branch(fwd_out, "CPAMBB", b, {tf, tc, 0, mem_cpambb()});
            }
            if (run_algo == "BoostR" || run_algo.find("Rlog") != string::npos || run_algo == "all") {
                vector<Value> P_branch_conv(batch.size());
                for(size_t k=0; k<batch.size(); k++) P_branch_conv[k] = make_pair(Point(batch[k].x, batch[k].y), batch[k].id);
                if (run_algo == "BoostR" || run_algo == "all") {
                    timer.start(); RTree v_boost = *boost_master; double tf = timer.stop()*1000;
                    timer.start(); v_boost.insert(P_branch_conv.begin(), P_branch_conv.end()); double tc = timer.stop()*1000;
                    boost_acc_elements += v_boost.size(); 
                    boost_branches[b] = v_boost;
                    log_branch(fwd_out, "BoostR", b, {tf, tc, 0, mem_boost(*boost_master, boost_acc_elements)});
                }
                if (run_algo.find("Rlog") != string::npos || run_algo == "all") {
                    string names[] = {"Rlog_1yr", "Rlog_2yr", "Rlog_3yr", "Rlog_4yr", "Rlog_5yr", "Rlog_NoSnap"};
                    for(int j=0; j<6; j++) {
                        if (run_algo == names[j] || run_algo == "all") {
                            timer.start(); RlogTree v_rlog = *rlog_master[j]; double l_tf = timer.stop()*1000;
                            timer.start(); v_rlog.commit_inserts(P_branch_conv); double l_tc = timer.stop()*1000;
                            rlog_acc_logs[j] += v_rlog.insert_log.size() + v_rlog.remove_log.size();
                            rlog_branches[j][b] = v_rlog;
                            log_branch(fwd_out, names[j], b, {l_tf, l_tc, 0, mem_rlog(*rlog_master[j], rlog_acc_logs[j])});
                        }
                    }
                }
            }
            cout << "    [Forward] Fork+Commit Batch [" << b + 1 << "/" << batches.size() << "] completed\n";
        }

        // Step 2: Merge Back to Master
        if (run_algo == "MVZD" || run_algo == "all") {
            timer.start();
            for (size_t b = 0; b < batches.size(); b++) {
                mvzd_master = std::get<0>(mvzd_tree->merge(mvzd_master, mvzd_master, mvzd_branches[b]));
                mvzd_tree->multi_version_roots.push_back(mvzd_master);
            }
            cout << "    [Forward] Merge MVZD took " << timer.stop()*1000 << " ms\n";
        }
        if (run_algo == "CPAMBB" || run_algo == "all") {
            timer.start();
            for (size_t b = 0; b < batches.size(); b++) {
                cpambb_master = std::get<0>(CPAMBB::map_merge(cpambb_master, cpambb_master, cpambb_branches[b]));
            }
            cout << "    [Forward] Merge CPAMBB took " << timer.stop()*1000 << " ms\n";
        }
        if (run_algo == "BoostR" || run_algo == "all") {
            timer.start();
            for (size_t b = 0; b < batches.size(); b++) {
                vector<Value> P_branch_conv(batches[b].size());
                for(size_t k=0; k<batches[b].size(); k++) P_branch_conv[k] = make_pair(Point(batches[b][k].x, batches[b][k].y), batches[b][k].id);
                boost_master->insert(P_branch_conv.begin(), P_branch_conv.end());
            }
            cout << "    [Forward] Merge BoostR took " << timer.stop()*1000 << " ms\n";
        }
        if (run_algo.find("Rlog") != string::npos || run_algo == "all") {
            string names[] = {"Rlog_1yr", "Rlog_2yr", "Rlog_3yr", "Rlog_4yr", "Rlog_5yr", "Rlog_NoSnap"};
            for(int j=0; j<6; j++) {
                if (run_algo == names[j] || run_algo == "all") {
                    timer.start();
                    for (size_t b = 0; b < batches.size(); b++) {
                        rlog_master[j]->merge(rlog_branches[j][b], rlog_master[j]->insert_log.size(), rlog_master[j]->remove_log.size());
                    }
                    cout << "    [Forward] Merge " << names[j] << " took " << timer.stop()*1000 << " ms\n";
                }
            }
        }
        
        run_queries_impl(to_string(year) + "_MidQuery");
        cout << "    [MidQuery] completed\n";
    }

    // ---------------------------------------------------------
    // Phase 2: Backward Loop (All Years in Reverse)
    // STAR TOPOLOGY
    // ---------------------------------------------------------
    for (int y_idx = (int)all_years.size() - 1; y_idx >= 0; y_idx--) {
        auto& yd = all_years[y_idx];
        int year = yd.year;
        auto& batches = yd.batches;
        cout << "--- YEAR " << year << " (BACKWARD) ---" << endl;
        
        ofstream bwd_out("verification_results/" + to_string(year) + "_Backward_Batch_" + run_algo + ".txt");
        bwd_out << "Algo | Batch | Fork_ms | Commit_ms | Merge_ms | Mem_MB\n";

        vector<shared_ptr<mvq::BaseNode>> mvzd_branches(batches.size());
        vector<CPAMBB::zmap> cpambb_branches(batches.size());
        vector<RTree> boost_branches(batches.size());
        vector<vector<RlogTree>> rlog_branches(6, vector<RlogTree>(batches.size()));

        // Step 1: Independent Fork & Delete
        for (int b = (int)batches.size() - 1; b >= 0; b--) {
            auto& batch = batches[b];
            
            if (run_algo == "MVZD" || run_algo == "all") {
                timer.start(); auto v_mvzd = mvzd_master; double tf = timer.stop()*1000;
                timer.start(); auto cv_mvzd = mvzd_tree->commit(v_mvzd, empty_pts, batch); double tc = timer.stop()*1000;
                mvzd_branches[b] = cv_mvzd;
                mvzd_tree->multi_version_roots.push_back(cv_mvzd);
                log_branch(bwd_out, "MVZD", b, {tf, tc, 0, mem_mvzd(*mvzd_tree)});
            }
            if (run_algo == "CPAMBB" || run_algo == "all") {
                timer.start(); auto v_cpambb = cpambb_master; double tf = timer.stop()*1000;
                timer.start(); auto cv_cpambb = CPAMBB::map_delete(batch, v_cpambb); double tc = timer.stop()*1000;
                cpambb_branches[b] = cv_cpambb;
                log_branch(bwd_out, "CPAMBB", b, {tf, tc, 0, mem_cpambb()});
            }
            if (run_algo == "BoostR" || run_algo.find("Rlog") != string::npos || run_algo == "all") {
                vector<Value> P_branch_conv(batch.size());
                for(size_t k=0; k<batch.size(); k++) P_branch_conv[k] = make_pair(Point(batch[k].x, batch[k].y), batch[k].id);
                if (run_algo == "BoostR" || run_algo == "all") {
                    timer.start(); RTree v_boost = *boost_master; double tf = timer.stop()*1000;
                    timer.start(); v_boost.remove(P_branch_conv.begin(), P_branch_conv.end()); double tc = timer.stop()*1000;
                    boost_acc_elements += v_boost.size();
                    boost_branches[b] = v_boost;
                    log_branch(bwd_out, "BoostR", b, {tf, tc, 0, mem_boost(*boost_master, boost_acc_elements)});
                }
                if (run_algo.find("Rlog") != string::npos || run_algo == "all") {
                    string names[] = {"Rlog_1yr", "Rlog_2yr", "Rlog_3yr", "Rlog_4yr", "Rlog_5yr", "Rlog_NoSnap"};
                    for(int j=0; j<6; j++) {
                        if (run_algo == names[j] || run_algo == "all") {
                            timer.start(); RlogTree v_rlog = *rlog_master[j]; double l_tf = timer.stop()*1000;
                            timer.start(); v_rlog.commit_removes(P_branch_conv); double l_tc = timer.stop()*1000;
                            rlog_acc_logs[j] += v_rlog.insert_log.size() + v_rlog.remove_log.size();
                            rlog_branches[j][b] = v_rlog;
                            log_branch(bwd_out, names[j], b, {l_tf, l_tc, 0, mem_rlog(*rlog_master[j], rlog_acc_logs[j])});
                        }
                    }
                }
            }
            int completed = batches.size() - b;
            cout << "    [Backward] Fork+Commit Batch [" << completed << "/" << batches.size() << "] completed\n";
        }

        // Merge Back Phase
        if (run_algo == "MVZD" || run_algo == "all") {
            timer.start();
            for (int b = (int)batches.size() - 1; b >= 0; b--) {
                mvzd_master = std::get<0>(mvzd_tree->merge(mvzd_master, mvzd_master, mvzd_branches[b]));
                mvzd_tree->multi_version_roots.push_back(mvzd_master);
            }
            cout << "    [Backward] Merge MVZD took " << timer.stop()*1000 << " ms\n";
        }
        if (run_algo == "CPAMBB" || run_algo == "all") {
            timer.start();
            for (int b = (int)batches.size() - 1; b >= 0; b--) {
                cpambb_master = std::get<0>(CPAMBB::map_merge(cpambb_master, cpambb_master, cpambb_branches[b]));
            }
            cout << "    [Backward] Merge CPAMBB took " << timer.stop()*1000 << " ms\n";
        }
        if (run_algo == "BoostR" || run_algo == "all") {
            timer.start();
            for (int b = (int)batches.size() - 1; b >= 0; b--) {
                vector<Value> P_branch_conv(batches[b].size());
                for(size_t k=0; k<batches[b].size(); k++) P_branch_conv[k] = make_pair(Point(batches[b][k].x, batches[b][k].y), batches[b][k].id);
                boost_master->remove(P_branch_conv.begin(), P_branch_conv.end());
            }
            cout << "    [Backward] Merge BoostR took " << timer.stop()*1000 << " ms\n";
        }
        if (run_algo.find("Rlog") != string::npos || run_algo == "all") {
            string names[] = {"Rlog_1yr", "Rlog_2yr", "Rlog_3yr", "Rlog_4yr", "Rlog_5yr", "Rlog_NoSnap"};
            for(int j=0; j<6; j++) {
                if (run_algo == names[j] || run_algo == "all") {
                    timer.start();
                    for (int b = (int)batches.size() - 1; b >= 0; b--) {
                        rlog_master[j]->merge(rlog_branches[j][b], rlog_master[j]->insert_log.size(), rlog_master[j]->remove_log.size());
                    }
                    cout << "    [Backward] Merge " << names[j] << " took " << timer.stop()*1000 << " ms\n";
                }
            }
        }
        
        run_queries_impl(to_string(year) + "_EndQuery");
        cout << "    [EndQuery] completed\n";
    }

    return 0;
}
"""

with open("verify_bench.cpp", "w") as f:
    f.write(cpp_code)
