import re

with open("verify_bench.cpp", "r") as f:
    text = f.read()

old_mem = """double mem_rlog(const RlogTree& master, const vector<RlogBranch>& history) {
    size_t log_elements = master.insert_log.size() + master.remove_log.size();
    for(const auto& br : history) log_elements += br.insert_log.size() + br.remove_log.size();
    return mem_boost() + (log_elements * sizeof(Value)) / (1024.0 * 1024.0);
}"""

new_mem = """double mem_rlog(const RlogTree& master, const vector<RlogBranch>& history) {
    size_t bytes = 0;
    bytes += sizeof(master.insert_log) + master.insert_log.capacity() * sizeof(Value);
    bytes += sizeof(master.remove_log) + master.remove_log.capacity() * sizeof(Value);
    for(const auto& br : history) {
        bytes += sizeof(RlogBranch);
        bytes += br.insert_log.capacity() * sizeof(Value);
        bytes += br.remove_log.capacity() * sizeof(Value);
    }
    return mem_boost() + (double)bytes / (1024.0 * 1024.0);
}"""

if old_mem in text:
    text = text.replace(old_mem, new_mem)
    with open("verify_bench.cpp", "w") as f:
        f.write(text)
    print("Patched RlogTree memory calculation successfully.")
else:
    print("Could not find the target code to patch.")

