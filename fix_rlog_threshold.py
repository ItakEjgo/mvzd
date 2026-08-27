import re

with open("verify_bench.cpp", "r") as f:
    text = f.read()

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

