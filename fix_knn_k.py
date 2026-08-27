import re
with open("verify_bench.cpp", "r") as f:
    text = f.read()

text = text.replace("snapshot->qbegin(bgi::nearest(bg_q, snapshot->size()))", "snapshot->qbegin(bgi::nearest(bg_q, (unsigned)(k + remove_log.size())))")

with open("verify_bench.cpp", "w") as f:
    f.write(text)

