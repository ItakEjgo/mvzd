import re

with open("verify_bench.cpp", "r") as f:
    text = f.read()

text = text.replace("double mem_mvzd() {", "double mem_mvzd(vector<shared_ptr<mvq::BaseNode>>& mvzd_global_history, mvq::Tree* mvzd_tree) {")
text = text.replace("double mem_cpambb() {", "double mem_cpambb(vector<CPAMBB::zmap>& cpambb_global_history) {")

text = text.replace("mem_mvzd()", "mem_mvzd(mvzd_global_history, mvzd_tree)")
text = text.replace("mem_cpambb()", "mem_cpambb(cpambb_global_history)")

with open("verify_bench.cpp", "w") as f:
    f.write(text)

