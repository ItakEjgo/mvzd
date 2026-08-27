import re
with open("verify_bench.cpp", "r") as f:
    text = f.read()

# 1. Inject TrackingAllocator
allocator_code = """
#include <atomic>
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
"""
text = re.sub(r"typedef bgi::rtree<Value, bgi::quadratic<32>> RTree;", allocator_code, text)

# 2. Refactor RlogTree to RlogBranch
rlog_branch_code = """struct RlogBranch {
    std::vector<Value> insert_log;
    std::vector<Value> remove_log;
};

class RlogTree {"""
text = text.replace("class RlogTree {", rlog_branch_code)

# 3. Fix RlogTree merge
merge_old = """    void merge(const RlogTree& branch, size_t original_ins_size, size_t original_rem_size) {
        size_t new_inserts = branch.insert_log.size() - original_ins_size;
        if (new_inserts > 0) insert_log.insert(insert_log.end(), branch.insert_log.begin() + original_ins_size, branch.insert_log.end());
        size_t new_removes = branch.remove_log.size() - original_rem_size;
        if (new_removes > 0) remove_log.insert(remove_log.end(), branch.remove_log.begin() + original_rem_size, branch.remove_log.end());
        if (insert_log.size() + remove_log.size() >= (size_t)threshold) compact();
    }"""
merge_new = """    void merge(const RlogBranch& branch) {
        if (!branch.insert_log.empty()) insert_log.insert(insert_log.end(), branch.insert_log.begin(), branch.insert_log.end());
        if (!branch.remove_log.empty()) remove_log.insert(remove_log.end(), branch.remove_log.begin(), branch.remove_log.end());
        if (insert_log.size() + remove_log.size() >= (size_t)threshold) compact();
    }"""
text = text.replace(merge_old, merge_new)

# 4. Fix mem_boost
text = text.replace("double mem_boost(const RTree& master, size_t acc_branch_elements) {\n    return ((master.size() + acc_branch_elements) * sizeof(Value)) / (1024.0 * 1024.0);\n}", "double mem_boost() { return boost_live_mem.load(std::memory_order_relaxed) / (1024.0 * 1024.0); }")

# 5. Fix mem_rlog
mem_rlog_old = "double mem_rlog(const RlogTree& master, size_t acc_branch_elements) { return ((master.size() + acc_branch_elements) * sizeof(Value)) / (1024.0 * 1024.0); }"
mem_rlog_new = """double mem_rlog(const RlogTree& master, const vector<RlogBranch>& history) {
    size_t log_elements = master.insert_log.size() + master.remove_log.size();
    for(const auto& br : history) log_elements += br.insert_log.size() + br.remove_log.size();
    return mem_boost() + (log_elements * sizeof(Value)) / (1024.0 * 1024.0);
}"""
text = text.replace(mem_rlog_old, mem_rlog_new)

# 6. Change mem_boost() call signatures
text = re.sub(r"mem_boost\(\*boost_master, boost_acc_elements\)", "mem_boost()", text)
text = re.sub(r"size_t boost_acc_elements = 0;\n", "", text)

# 7. Add rlog_master_history & rlog_global_history, and patch init loop
init_old = """    vector<size_t> rlog_acc_logs(6, 0);
    
    mvq::Tree* mvzd_tree = nullptr;
    shared_ptr<mvq::BaseNode> mvzd_master;
    CPAMBB::zmap cpambb_master;
    vector<CPAMBB::zmap> cpambb_global_history;
    
    vector<Value> P_base_conv;
    RTree* boost_master = nullptr;
    RlogTree *rlog_master[6];"""
init_new = """    vector<shared_ptr<mvq::BaseNode>> mvzd_global_history;
    vector<shared_ptr<RTree>> boost_global_history;
    vector<vector<RlogBranch>> rlog_global_history(6);
    vector<vector<shared_ptr<RTree>>> rlog_master_history(6);
    
    mvq::Tree* mvzd_tree = nullptr;
    shared_ptr<mvq::BaseNode> mvzd_master;
    CPAMBB::zmap cpambb_master;
    vector<CPAMBB::zmap> cpambb_global_history;
    
    vector<Value> P_base_conv;
    RTree* boost_master = nullptr;
    RlogTree *rlog_master[6];"""
text = text.replace(init_old, init_new)

# Fix BoostR global history init
text = text.replace("boost_master = new RTree(P_base_conv.begin(), P_base_conv.end());", "boost_master = new RTree(P_base_conv.begin(), P_base_conv.end());\n        boost_global_history.push_back(make_shared<RTree>(*boost_master));")

# Fix Rlog init
rlog_init_old = """        for(int i=0; i<6; i++) {
            rlog_master[i] = new RlogTree(thresholds[i]);
            rlog_master[i]->build_base(P_base_conv);
        }"""
rlog_init_new = """        string names[] = {"Rlog_1yr", "Rlog_2yr", "Rlog_3yr", "Rlog_4yr", "Rlog_5yr", "Rlog_NoSnap"};
        for(int i=0; i<6; i++) {
            if (run_algo == names[i] || run_algo == "all") {
                rlog_master[i] = new RlogTree(thresholds[i]);
                rlog_master[i]->build_base(P_base_conv);
                rlog_master_history[i].push_back(rlog_master[i]->snapshot);
            }
        }"""
text = text.replace(rlog_init_old, rlog_init_new)

# Fix mem_rlog calls inside run_queries_impl
text = re.sub(r"mem_rlog\(\*rlog_master\[j\], rlog_acc_logs\[j\]\)", r"mem_rlog(*rlog_master[j], rlog_global_history[j])", text)

# Fix Forward branch creation
fwd_rlog_old = """                            timer.start(); RlogTree v_rlog = *rlog_master[j]; double l_tf = timer.stop()*1000;
                            timer.start(); v_rlog.commit_inserts(P_branch_conv); double l_tc = timer.stop()*1000;
                            rlog_acc_logs[j] += v_rlog.insert_log.size() + v_rlog.remove_log.size();
                            rlog_branches[j][b] = v_rlog;
                            log_branch(fwd_out, names[j], b, {l_tf, l_tc, 0, mem_rlog(*rlog_master[j], rlog_acc_logs[j])});"""
fwd_rlog_new = """                            timer.start(); RlogBranch v_rlog; double l_tf = timer.stop()*1000;
                            timer.start(); v_rlog.insert_log = P_branch_conv; double l_tc = timer.stop()*1000;
                            rlog_global_history[j].push_back(v_rlog);
                            log_branch(fwd_out, names[j], b, {l_tf, l_tc, 0, mem_rlog(*rlog_master[j], rlog_global_history[j])});"""
text = text.replace(fwd_rlog_old, fwd_rlog_new)

# Fix Forward BoostR Branch creation
fwd_boost_old = """                    timer.start(); RTree v_boost = *boost_master; double tf = timer.stop()*1000;
                    timer.start(); v_boost.insert(P_branch_conv.begin(), P_branch_conv.end()); double tc = timer.stop()*1000;
                    boost_acc_elements += v_boost.size(); 
                    boost_branches[b] = v_boost;
                    log_branch(fwd_out, "BoostR", b, {tf, tc, 0, mem_boost()});"""
fwd_boost_new = """                    timer.start(); RTree v_boost = *boost_master; double tf = timer.stop()*1000;
                    timer.start(); v_boost.insert(P_branch_conv.begin(), P_branch_conv.end()); double tc = timer.stop()*1000;
                    boost_global_history.push_back(make_shared<RTree>(v_boost));
                    log_branch(fwd_out, "BoostR", b, {tf, tc, 0, mem_boost()});"""
text = text.replace(fwd_boost_old, fwd_boost_new)


# Fix Forward Merge
fwd_rlog_merge_old = """                    for (size_t b = 0; b < batches.size(); b++) {
                        rlog_master[j]->merge(rlog_branches[j][b], rlog_master[j]->insert_log.size(), rlog_master[j]->remove_log.size());
                    }
                    double tm = timer.stop()*1000;
                    fwd_out << names[j] << " | SUMMARY | 0 | 0 | " << tm << " | " << mem_rlog(*rlog_master[j], rlog_acc_logs[j]) << "\\n";"""
fwd_rlog_merge_new = """                    for (size_t b = 0; b < batches.size(); b++) {
                        rlog_master[j]->merge(rlog_global_history[j][rlog_global_history[j].size() - batches.size() + b]);
                    }
                    rlog_master_history[j].push_back(rlog_master[j]->snapshot);
                    double tm = timer.stop()*1000;
                    fwd_out << names[j] << " | SUMMARY | 0 | 0 | " << tm << " | " << mem_rlog(*rlog_master[j], rlog_global_history[j]) << "\\n";"""
text = text.replace(fwd_rlog_merge_old, fwd_rlog_merge_new)

fwd_boost_merge_old = """            for (size_t b = 0; b < batches.size(); b++) {
                vector<Value> P_branch_conv(batches[b].size());
                for(size_t k=0; k<batches[b].size(); k++) P_branch_conv[k] = make_pair(BoostPoint(batches[b][k].x, batches[b][k].y), batches[b][k].id);
                boost_master->insert(P_branch_conv.begin(), P_branch_conv.end());
            }
            double tm = timer.stop()*1000;"""
fwd_boost_merge_new = """            for (size_t b = 0; b < batches.size(); b++) {
                vector<Value> P_branch_conv(batches[b].size());
                for(size_t k=0; k<batches[b].size(); k++) P_branch_conv[k] = make_pair(BoostPoint(batches[b][k].x, batches[b][k].y), batches[b][k].id);
                boost_master->insert(P_branch_conv.begin(), P_branch_conv.end());
            }
            boost_global_history.push_back(make_shared<RTree>(*boost_master));
            double tm = timer.stop()*1000;"""
text = text.replace(fwd_boost_merge_old, fwd_boost_merge_new)


# Fix Backward branch creation
bwd_rlog_old = """                            timer.start(); RlogTree v_rlog = *rlog_master[j]; double l_tf = timer.stop()*1000;
                            timer.start(); v_rlog.commit_removes(P_branch_conv); double l_tc = timer.stop()*1000;
                            rlog_acc_logs[j] += v_rlog.insert_log.size() + v_rlog.remove_log.size();
                            rlog_branches[j][b] = v_rlog;
                            log_branch(bwd_out, names[j], b, {l_tf, l_tc, 0, mem_rlog(*rlog_master[j], rlog_acc_logs[j])});"""
bwd_rlog_new = """                            timer.start(); RlogBranch v_rlog; double l_tf = timer.stop()*1000;
                            timer.start(); v_rlog.remove_log = P_branch_conv; double l_tc = timer.stop()*1000;
                            rlog_global_history[j].push_back(v_rlog);
                            log_branch(bwd_out, names[j], b, {l_tf, l_tc, 0, mem_rlog(*rlog_master[j], rlog_global_history[j])});"""
text = text.replace(bwd_rlog_old, bwd_rlog_new)

# Fix Backward BoostR Branch creation
bwd_boost_old = """                    timer.start(); RTree v_boost = *boost_master; double tf = timer.stop()*1000;
                    timer.start(); v_boost.remove(P_branch_conv.begin(), P_branch_conv.end()); double tc = timer.stop()*1000;
                    boost_acc_elements += v_boost.size();
                    boost_branches[b] = v_boost;
                    log_branch(bwd_out, "BoostR", b, {tf, tc, 0, mem_boost()});"""
bwd_boost_new = """                    timer.start(); RTree v_boost = *boost_master; double tf = timer.stop()*1000;
                    timer.start(); v_boost.remove(P_branch_conv.begin(), P_branch_conv.end()); double tc = timer.stop()*1000;
                    boost_global_history.push_back(make_shared<RTree>(v_boost));
                    log_branch(bwd_out, "BoostR", b, {tf, tc, 0, mem_boost()});"""
text = text.replace(bwd_boost_old, bwd_boost_new)

# Fix Backward Merge
bwd_rlog_merge_old = """                    for (int b = (int)batches.size() - 1; b >= 0; b--) {
                        rlog_master[j]->merge(rlog_branches[j][b], rlog_master[j]->insert_log.size(), rlog_master[j]->remove_log.size());
                    }
                    double tm = timer.stop()*1000;
                    bwd_out << names[j] << " | SUMMARY | 0 | 0 | " << tm << " | " << mem_rlog(*rlog_master[j], rlog_acc_logs[j]) << "\\n";"""
bwd_rlog_merge_new = """                    for (int b = (int)batches.size() - 1; b >= 0; b--) {
                        rlog_master[j]->merge(rlog_global_history[j][rlog_global_history[j].size() - batches.size() + ((batches.size()-1)-b)]);
                    }
                    rlog_master_history[j].push_back(rlog_master[j]->snapshot);
                    double tm = timer.stop()*1000;
                    bwd_out << names[j] << " | SUMMARY | 0 | 0 | " << tm << " | " << mem_rlog(*rlog_master[j], rlog_global_history[j]) << "\\n";"""
text = text.replace(bwd_rlog_merge_old, bwd_rlog_merge_new)

bwd_boost_merge_old = """            for (int b = (int)batches.size() - 1; b >= 0; b--) {
                vector<Value> P_branch_conv(batches[b].size());
                for(size_t k=0; k<batches[b].size(); k++) P_branch_conv[k] = make_pair(BoostPoint(batches[b][k].x, batches[b][k].y), batches[b][k].id);
                boost_master->remove(P_branch_conv.begin(), P_branch_conv.end());
            }
            double tm = timer.stop()*1000;"""
bwd_boost_merge_new = """            for (int b = (int)batches.size() - 1; b >= 0; b--) {
                vector<Value> P_branch_conv(batches[b].size());
                for(size_t k=0; k<batches[b].size(); k++) P_branch_conv[k] = make_pair(BoostPoint(batches[b][k].x, batches[b][k].y), batches[b][k].id);
                boost_master->remove(P_branch_conv.begin(), P_branch_conv.end());
            }
            boost_global_history.push_back(make_shared<RTree>(*boost_master));
            double tm = timer.stop()*1000;"""
text = text.replace(bwd_boost_merge_old, bwd_boost_merge_new)

# Remove unused rlog_branches and boost_branches declarations
text = re.sub(r"\s*vector<RTree> boost_branches\(batches\.size\(\)\);\s*", "\n        ", text)
text = re.sub(r"\s*vector<vector<RlogTree>> rlog_branches\(6, vector<RlogTree>\(batches\.size\(\)\)\);\s*", "\n        ", text)

with open("verify_bench.cpp", "w") as f:
    f.write(text)

