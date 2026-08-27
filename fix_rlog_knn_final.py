import re

with open("verify_bench.cpp", "r") as f:
    text = f.read()

text = text.replace(
"""    void commit_removes(const std::vector<Value>& del_pts) {
        remove_log.insert(remove_log.end(), del_pts.begin(), del_pts.end());
        
    }""",
"""    void commit_removes(const std::vector<Value>& del_pts) {
        remove_log.insert(remove_log.end(), del_pts.begin(), del_pts.end());
        cache_valid = false;
    }""")

old_knn = """    std::vector<Value> knn_report(const geobase::Point& q, size_t k) const {
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
        for (auto it = snapshot->qbegin(bgi::nearest(bg_q, (unsigned)(k + remove_log.size()))); it != snapshot->qend(); ++it) {
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
    }"""

new_knn = """    std::vector<Value> knn_report(const geobase::Point& q, size_t k) const {
        update_cache();
        auto calc_sqr_dist = [](const geobase::Point& p1, const BoostPoint& p2) {
            double dx = p1.x - p2.get<0>(), dy = p1.y - p2.get<1>();
            return dx*dx + dy*dy;
        };
        std::vector<Value> result;
        size_t search_k = k * 2;
        if (search_k == 0) search_k = 1;
        while (true) {
            std::priority_queue<std::pair<double, Value>, std::vector<std::pair<double, Value>>, MaxHeapCmp> max_heap;
            for (const auto& val : insert_log) {
                if (removed_ids.find(val.second) == removed_ids.end()) {
                    double dist = calc_sqr_dist(q, val.first);
                    if (max_heap.size() < k) max_heap.push({dist, val});
                    else if (dist < max_heap.top().first) { max_heap.pop(); max_heap.push({dist, val}); }
                }
            }
            
            size_t tree_count = 0;
            double last_tree_dist = -1;
            BoostPoint bg_q(q.x, q.y);
            for (auto it = snapshot->qbegin(bgi::nearest(bg_q, (unsigned)search_k)); it != snapshot->qend(); ++it) {
                tree_count++;
                double dist = calc_sqr_dist(q, it->first);
                last_tree_dist = dist;
                if (removed_ids.find(it->second) == removed_ids.end()) {
                    if (max_heap.size() < k) max_heap.push({dist, *it});
                    else if (dist < max_heap.top().first) { max_heap.pop(); max_heap.push({dist, *it}); }
                }
            }
            
            if (tree_count < search_k) {
                while (!max_heap.empty()) { result.push_back(max_heap.top().second); max_heap.pop(); }
                break;
            }
            if (max_heap.size() == k && last_tree_dist >= max_heap.top().first) {
                while (!max_heap.empty()) { result.push_back(max_heap.top().second); max_heap.pop(); }
                break;
            }
            if (snapshot->size() == 0 || search_k >= snapshot->size()) {
                while (!max_heap.empty()) { result.push_back(max_heap.top().second); max_heap.pop(); }
                break;
            }
            search_k = min(search_k * 2, (size_t)snapshot->size());
        }
        return result;
    }"""

text = text.replace(old_knn, new_knn)

with open("verify_bench.cpp", "w") as f:
    f.write(text)

