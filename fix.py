import re

with open("verify_bench.cpp", "r") as f:
    text = f.read()

# 1. Fix compact
compact_old = """    void compact() {
        if (insert_log.empty() && remove_log.empty()) return;
        auto new_snap = make_shared<RTree>(*snapshot);
        for (const auto& val : insert_log) new_snap->insert(val);
        for (const auto& val : remove_log) {
            std::vector<Value> res;
            new_snap->query(bgi::intersects(val.first), std::back_inserter(res));
            for(auto& r : res) {
                if (r.second == val.second) {
                    new_snap->remove(r);
                    break;
                }
            }
        }
        snapshot = new_snap;
        insert_log.clear();
        remove_log.clear();
    }"""
compact_new = """    void compact() {
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
    }"""
text = text.replace(compact_old, compact_new)

# 2. Fix threshold
rlog_tree_def_old = """    int threshold;
    shared_ptr<RTree> snapshot;
    std::vector<Value> insert_log;
    std::vector<Value> remove_log;
    
    RlogTree() : threshold(9999999) { snapshot = make_shared<RTree>(); }
    RlogTree(int thresh) : threshold(thresh) { snapshot = make_shared<RTree>(); }"""

rlog_tree_def_new = """    int threshold;
    int current_batches;
    shared_ptr<RTree> snapshot;
    std::vector<Value> insert_log;
    std::vector<Value> remove_log;
    
    RlogTree() : threshold(9999999), current_batches(0) { snapshot = make_shared<RTree>(); }
    RlogTree(int thresh) : threshold(thresh), current_batches(0) { snapshot = make_shared<RTree>(); }"""
text = text.replace(rlog_tree_def_old, rlog_tree_def_new)

rlog_merge_old = """    void merge(const RlogBranch& branch) {
        if (!branch.insert_log.empty()) insert_log.insert(insert_log.end(), branch.insert_log.begin(), branch.insert_log.end());
        if (!branch.remove_log.empty()) remove_log.insert(remove_log.end(), branch.remove_log.begin(), branch.remove_log.end());
        if (insert_log.size() + remove_log.size() >= (size_t)threshold) compact();
    }"""
rlog_merge_new = """    void merge(const RlogBranch& branch) {
        if (!branch.insert_log.empty()) insert_log.insert(insert_log.end(), branch.insert_log.begin(), branch.insert_log.end());
        if (!branch.remove_log.empty()) remove_log.insert(remove_log.end(), branch.remove_log.begin(), branch.remove_log.end());
        current_batches++;
        if (current_batches >= threshold) {
            compact();
            current_batches = 0;
        }
    }"""
text = text.replace(rlog_merge_old, rlog_merge_new)
text = text.replace("if (insert_log.size() + remove_log.size() >= (size_t)threshold) compact();", "")

with open("verify_bench.cpp", "w") as f:
    f.write(text)

