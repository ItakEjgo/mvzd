import re

with open("verify_bench.cpp", "r") as f:
    text = f.read()

text = text.replace(
    "double mem_mvzd() { return mvq::global_live_mem.load(std::memory_order_relaxed) / (1024.0 * 1024.0); }",
    """double mem_mvzd() {
    if (!mvzd_tree || !mvzd_master) return mvq::global_live_mem.load(std::memory_order_relaxed) / (1024.0 * 1024.0);
    mvq::tree_stat stat;
    mvq::NodeMap node_map;
    mvzd_tree->get_tree_statistics(mvzd_master, stat, node_map);
    return (stat.mem_inte_nodes + stat.mem_leaf_nodes) / (1024.0 * 1024.0);
}"""
)

with open("verify_bench.cpp", "w") as f:
    f.write(text)

