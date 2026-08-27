import re

with open("verify_bench.cpp", "r") as f:
    content = f.read()

bad_block = """        for(int i=0; i<6; i++) {
            rlog_master[i] = new RlogTree(thresholds[i]);
            rlog_master[i]->build_base(P_base_conv);
            rlog_master_history[i].push_back(rlog_master[i]->snapshot);
        }"""

good_block = """        string names[] = {"Rlog_1yr", "Rlog_2yr", "Rlog_3yr", "Rlog_4yr", "Rlog_5yr", "Rlog_NoSnap"};
        for(int i=0; i<6; i++) {
            if (run_algo == names[i] || run_algo == "all") {
                rlog_master[i] = new RlogTree(thresholds[i]);
                rlog_master[i]->build_base(P_base_conv);
                rlog_master_history[i].push_back(rlog_master[i]->snapshot);
            }
        }"""

if bad_block in content:
    content = content.replace(bad_block, good_block)
    with open("verify_bench.cpp", "w") as f:
        f.write(content)
    print("Patched Rlog init.")
else:
    print("Block not found!")

