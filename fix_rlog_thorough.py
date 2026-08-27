import re

with open("verify_bench.cpp", "r") as f:
    text = f.read()

# 1. update_cache rehashing
text = text.replace(
"""    void update_cache() const {
        if (!cache_valid) {
            removed_ids.clear();
            for (const auto& val : remove_log) removed_ids.insert(val.second);
            cache_valid = true;
        }
    }""",
"""    void update_cache() const {
        if (!cache_valid) {
            removed_ids.clear();
            removed_ids.reserve(remove_log.size()); // 优化1：防 rehash 扩容
            for (const auto& val : remove_log) removed_ids.insert(val.second);
            cache_valid = true;
        }
    }""")

# 2. compact memory leak
text = text.replace(
"""        snapshot = make_shared<RTree>(next_pts.begin(), next_pts.end());
        insert_log.clear();
        remove_log.clear();
        cache_valid = false;
    }""",
"""        snapshot = make_shared<RTree>(next_pts.begin(), next_pts.end());
        
        // 优化2：真正释放内存，防止假性 OOM 和内存指标虚高
        std::vector<Value>().swap(insert_log); 
        std::vector<Value>().swap(remove_log);
        removed_ids = std::unordered_set<size_t>(); // 释放哈希表底层的 bucket 内存
        cache_valid = true; // 空表即为有效状态
    }""")

# 3. merge cache invalidation bug
text = text.replace(
"""    void merge(const RlogBranch& branch) {
        if (!branch.insert_log.empty()) { insert_log.insert(insert_log.end(), branch.insert_log.begin(), branch.insert_log.end()); cache_valid = false; }
        if (!branch.remove_log.empty()) { remove_log.insert(remove_log.end(), branch.remove_log.begin(), branch.remove_log.end()); cache_valid = false; }""",
"""    void merge(const RlogBranch& branch) {
        // 优化3：只插入数据时，根本不需要废弃关于“死亡名单”的哈希表缓存！
        if (!branch.insert_log.empty()) { insert_log.insert(insert_log.end(), branch.insert_log.begin(), branch.insert_log.end()); }
        if (!branch.remove_log.empty()) { remove_log.insert(remove_log.end(), branch.remove_log.begin(), branch.remove_log.end()); cache_valid = false; }""")

# 4. range_report semantic inconsistency
text = text.replace(
"""        for (const auto& val : insert_log) {
            if (removed_ids.find(val.second) == removed_ids.end()) {
                if (bg::within(val.first, box)) result.push_back(val);
            }
        }""",
"""        for (const auto& val : insert_log) {
            if (removed_ids.find(val.second) == removed_ids.end()) {
                // 优化4：修正语义不一致。RTree 用的是 intersects(包含边界), 这里之前用 within(不包含边界) 是有 Bug 的。
                if (bg::intersects(val.first, box)) result.push_back(val);
            }
        }""")

# 5. knn early break
text = text.replace(
"""        for (auto it = snapshot->qbegin(bgi::nearest(bg_q, (unsigned)k) && bgi::satisfies(is_alive)); it != snapshot->qend(); ++it) {
            double dist = calc_sqr_dist(q, it->first);
            if (max_heap.size() < k) max_heap.push({dist, *it});
            else if (dist < max_heap.top().first) { max_heap.pop(); max_heap.push({dist, *it}); }
        }""",
"""        for (auto it = snapshot->qbegin(bgi::nearest(bg_q, (unsigned)k) && bgi::satisfies(is_alive)); it != snapshot->qend(); ++it) {
            double dist = calc_sqr_dist(q, it->first);
            if (max_heap.size() < k) {
                max_heap.push({dist, *it});
            } else if (dist < max_heap.top().first) { 
                max_heap.pop(); max_heap.push({dist, *it}); 
            } else {
                break; // 优化5：提前终止！如果树里找出的点已经比 max_heap 里的点更远，后续的树节点只会更远，直接打断 RTree 遍历！
            }
        }""")

with open("verify_bench.cpp", "w") as f:
    f.write(text)
print("Applied all thorough optimizations.")
