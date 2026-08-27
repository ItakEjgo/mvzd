import re

with open("verify_bench.cpp", "r") as f:
    text = f.read()

# We need to replace the entire knn_report function
# Find the start of knn_report
start_idx = text.find("    std::vector<Value> knn_report(const geobase::Point& q, size_t k) const {")
# Find the end of knn_report
# It ends right before `    size_t size() const {`
end_idx = text.find("    size_t size() const {")

if start_idx != -1 and end_idx != -1:
    old_knn = text[start_idx:end_idx]
    new_knn = """    std::vector<Value> knn_report(const geobase::Point& q, size_t k) const {
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
            if (max_heap.size() < k) max_heap.push({dist, *it});
            else if (dist < max_heap.top().first) { max_heap.pop(); max_heap.push({dist, *it}); }
        }
        
        std::vector<Value> result;
        while (!max_heap.empty()) { result.push_back(max_heap.top().second); max_heap.pop(); }
        return result;
    }
"""
    text = text[:start_idx] + new_knn + text[end_idx:]
    with open("verify_bench.cpp", "w") as f:
        f.write(text)
    print("Replaced knn_report successfully.")
else:
    print("Failed to find knn_report bounds.")
