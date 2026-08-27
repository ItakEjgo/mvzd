import re

with open("verify_bench.cpp", "r") as f:
    text = f.read()

start_idx = text.find("    void compact() {")
end_idx = text.find("    void commit_inserts(const std::vector<Value>& new_pts) {")

if start_idx != -1 and end_idx != -1:
    old_compact = text[start_idx:end_idx]
    new_compact = """    void compact() {
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
        insert_log.clear();
        remove_log.clear();
        cache_valid = false;
    }
"""
    text = text[:start_idx] + new_compact + text[end_idx:]
    with open("verify_bench.cpp", "w") as f:
        f.write(text)
    print("Replaced compact with std::copy_if successfully.")
else:
    print("Failed to find compact bounds.")
