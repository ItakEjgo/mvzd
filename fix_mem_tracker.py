import re

with open("verify_bench.cpp", "r") as f:
    text = f.read()

new_mem_funcs = """
double mem_mvzd() {
    if (mvzd_global_history.empty()) return 0;
    mvq::tree_stat stat;
    mvq::NodeMap node_map;
    for (auto rt: mvzd_global_history){
        if (rt) mvzd_tree->get_tree_statistics(rt, stat, node_map);
    }
    return (stat.mem_inte_nodes + stat.mem_leaf_nodes) / (1024.0 * 1024.0);
}

double mem_cpambb() {
    if (cpambb_global_history.empty()) return 0;
    std::unordered_map<size_t, bool> mmp;
    auto f_noop = [&](const auto &et){ return 0; };
    double total_mem = 0;
    for(auto t : cpambb_global_history) {
        total_mem += 1.0 * t.size_in_bytes(f_noop, mmp);
    }
    return total_mem / 1024.0 / 1024.0;
}
"""

text = re.sub(r'double mem_mvzd\(\) \{.*?\}\n', '', text, flags=re.DOTALL)
text = re.sub(r'double mem_cpambb\(\) \{.*?\}\n', '', text, flags=re.DOTALL)

# Insert the new functions before double mem_boost
text = text.replace("double mem_boost() {", new_mem_funcs + "\ndouble mem_boost() {")

with open("verify_bench.cpp", "w") as f:
    f.write(text)
