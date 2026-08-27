import re

with open("verify_bench.cpp", "r") as f:
    text = f.read()

text = re.sub(r'double mem_mvzd\(.*?\) \{.*?\n\}', '', text, flags=re.DOTALL)
text = re.sub(r'double mem_cpambb\(.*?\) \{.*?\n\}', '', text, flags=re.DOTALL)
text = text.replace("mem_mvzd(mvzd_global_history, mvzd_tree)", "mem_mvzd()")
text = text.replace("mem_cpambb(cpambb_global_history)", "mem_cpambb()")

new_mem = """
double mem_mvzd() { return mvq::global_live_mem.load(std::memory_order_relaxed) / (1024.0 * 1024.0); }
double mem_cpambb() { return (sizeof(CPAMBB::zmap::GC::regular_node) * CPAMBB::zmap::GC::used_node()) / (1024.0 * 1024.0); }
"""

text = text.replace("double mem_boost() {", new_mem + "\ndouble mem_boost() {")

with open("verify_bench.cpp", "w") as f:
    f.write(text)

