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
        // 预分配足够的空间，避免 push_back 时的多次动态扩容
        next_pts.reserve(snapshot->size() + insert_log.size());
        
        if (!remove_log.empty()) {
            // 完美复用我们之前写好的全局哈希表，干掉局部的冗余建表
            update_cache();
            
            // 直接遍历树结构，跳过中间数组，直接把存活点塞进 next_pts
            for (const auto& val : *snapshot) {
                if (removed_ids.find(val.second) == removed_ids.end()) {
                    next_pts.push_back(val);
                }
            }
            // 同样处理 insert_log
            for (const auto& val : insert_log) {
                if (removed_ids.find(val.second) == removed_ids.end()) {
                    next_pts.push_back(val);
                }
            }
        } else {
            // 如果没有删除操作，直接用最快的方式拼接
            next_pts.assign(snapshot->begin(), snapshot->end());
            next_pts.insert(next_pts.end(), insert_log.begin(), insert_log.end());
        }
        
        // 批量构建出新树 (Boost RTree bulk-loading)
        snapshot = make_shared<RTree>(next_pts.begin(), next_pts.end());
        
        insert_log.clear();
        remove_log.clear();
        cache_valid = false;
    }
"""
    text = text[:start_idx] + new_compact + text[end_idx:]
    with open("verify_bench.cpp", "w") as f:
        f.write(text)
    print("Replaced compact successfully.")
else:
    print("Failed to find compact bounds.")
