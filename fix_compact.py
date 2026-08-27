import re
with open("verify_bench.cpp", "r") as f:
    text = f.read()

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
with open("verify_bench.cpp", "w") as f:
    f.write(text)
