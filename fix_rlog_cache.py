import re

with open("verify_bench.cpp", "r") as f:
    text = f.read()

rlog_vars_old = """    int threshold;
    int current_batches;
    shared_ptr<RTree> snapshot;
    std::vector<Value> insert_log;
    std::vector<Value> remove_log;"""

rlog_vars_new = """    int threshold;
    int current_batches;
    shared_ptr<RTree> snapshot;
    std::vector<Value> insert_log;
    std::vector<Value> remove_log;
    mutable bool cache_valid = false;
    mutable std::unordered_set<size_t> removed_ids;
    
    void update_cache() const {
        if (!cache_valid) {
            removed_ids.clear();
            for (const auto& val : remove_log) removed_ids.insert(val.second);
            cache_valid = true;
        }
    }"""
text = text.replace(rlog_vars_old, rlog_vars_new)

text = text.replace("insert_log.clear();\n        remove_log.clear();", "insert_log.clear();\n        remove_log.clear();\n        cache_valid = false;")
text = text.replace("if (!branch.insert_log.empty()) insert_log.insert(insert_log.end(), branch.insert_log.begin(), branch.insert_log.end());", "if (!branch.insert_log.empty()) { insert_log.insert(insert_log.end(), branch.insert_log.begin(), branch.insert_log.end()); cache_valid = false; }")
text = text.replace("if (!branch.remove_log.empty()) remove_log.insert(remove_log.end(), branch.remove_log.begin(), branch.remove_log.end());", "if (!branch.remove_log.empty()) { remove_log.insert(remove_log.end(), branch.remove_log.begin(), branch.remove_log.end()); cache_valid = false; }")

report_old_range = """    std::vector<Value> range_report(const geobase::Bounding_Box& q) const {
        std::vector<Value> snap_res;
        bg::model::box<BoostPoint> box(BoostPoint(q.first.x, q.first.y), BoostPoint(q.second.x, q.second.y));
        snapshot->query(bgi::intersects(box), std::back_inserter(snap_res));
        std::unordered_set<size_t> removed_ids;
        for (const auto& val : remove_log) removed_ids.insert(val.second);"""

report_new_range = """    std::vector<Value> range_report(const geobase::Bounding_Box& q) const {
        update_cache();
        std::vector<Value> snap_res;
        bg::model::box<BoostPoint> box(BoostPoint(q.first.x, q.first.y), BoostPoint(q.second.x, q.second.y));
        snapshot->query(bgi::intersects(box), std::back_inserter(snap_res));"""

text = text.replace(report_old_range, report_new_range)

report_old_knn = """    std::vector<Value> knn_report(const geobase::Point& q, size_t k) const {
        std::unordered_set<size_t> removed_ids;
        for (const auto& val : remove_log) removed_ids.insert(val.second);"""

report_new_knn = """    std::vector<Value> knn_report(const geobase::Point& q, size_t k) const {
        update_cache();"""

text = text.replace(report_old_knn, report_new_knn)

with open("verify_bench.cpp", "w") as f:
    f.write(text)

