import re

with open("verify_bench.cpp", "r") as f:
    text = f.read()

# Replace double mem_cpambb()
text = re.sub(r'double mem_cpambb\(\) \{.*?\}\n', '', text)

new_cpambb_mem = """
std::unordered_map<size_t, bool> cpambb_mmp;
double cpambb_global_mem = 0;

double mem_cpambb(const CPAMBB::zmap& latest_branch) {
    auto f_noop = [&](const auto &et){ return 0; };
    cpambb_global_mem += latest_branch.size_in_bytes(f_noop, cpambb_mmp) / (1024.0 * 1024.0);
    return cpambb_global_mem;
}
"""

text = text.replace("double mem_mvzd()", new_cpambb_mem + "\ndouble mem_mvzd()")

# Replace calls to mem_cpambb()
# Note: we need to pass the appropriate branch or master
# In range_set for CPAMBB:
text = text.replace("mem_cpambb()", "mem_cpambb(cpambb_master)")

# In CPAMBB batch loops:
# log_branch(fwd_out, "CPAMBB", b, {tf, tc, 0, mem_cpambb()});
# -> we must pass cv_cpambb to log_branch because it's the branch we just created!
text = text.replace("log_branch(fwd_out, \"CPAMBB\", b, {tf, tc, 0, mem_cpambb(cpambb_master)});", "log_branch(fwd_out, \"CPAMBB\", b, {tf, tc, 0, mem_cpambb(cv_cpambb)});")
text = text.replace("log_branch(bwd_out, \"CPAMBB\", b, {tf, tc, 0, mem_cpambb(cpambb_master)});", "log_branch(bwd_out, \"CPAMBB\", b, {tf, tc, 0, mem_cpambb(cv_cpambb)});")

with open("verify_bench.cpp", "w") as f:
    f.write(text)

