import re

with open("verify_bench.cpp", "r") as f:
    text = f.read()

start_idx = text.find("    std::vector<Value> range_report(const geobase::Bounding_Box& q) const {")
end_idx = text.find("    std::vector<Value> knn_report(const geobase::Point& q, size_t k) const {")

if start_idx != -1 and end_idx != -1:
    new_range = """    std::vector<Value> range_report(const geobase::Bounding_Box& q) const {
        update_cache();
        std::vector<Value> result;
        bg::model::box<BoostPoint> box(BoostPoint(q.first.x, q.first.y), BoostPoint(q.second.x, q.second.y));
        auto is_alive = [this](Value const& v) { return removed_ids.find(v.second) == removed_ids.end(); };
        
        // 优雅：直接在树的遍历底层完成过滤，连中间临时数组 snap_res 都省了
        snapshot->query(bgi::intersects(box) && bgi::satisfies(is_alive), std::back_inserter(result));
        
        for (const auto& val : insert_log) {
            if (removed_ids.find(val.second) == removed_ids.end()) {
                if (bg::within(val.first, box)) result.push_back(val);
            }
        }
        return result;
    }
"""
    text = text[:start_idx] + new_range + text[end_idx:]
    with open("verify_bench.cpp", "w") as f:
        f.write(text)
    print("Replaced range_report successfully.")
else:
    print("Failed to find range_report bounds.")
