import os
import glob

def load_results(filename):
    results = {}
    if not os.path.exists(filename): return results
    with open(filename, 'r') as f:
        for line in f:
            if line.startswith('Algo') or line.strip() == '': continue
            parts = [p.strip() for p in line.split('|')]
            if len(parts) >= 5 and 'SUMMARY' not in line:
                qtype = parts[1]
                qid = parts[2]
                count = parts[3]
                hash_val = parts[4]
                key = f"{qtype}_{qid}"
                results[key] = (count, hash_val)
    return results

res_dir = "results_bp_50"
mvzd_files = glob.glob(f"{res_dir}/*Query_MVZD.txt")
phases = [os.path.basename(f).replace('_MVZD.txt', '') for f in mvzd_files]

all_match = True
total_queries_checked = 0
mismatch_details = []

for phase in sorted(phases):
    mvzd_file = f"{res_dir}/{phase}_MVZD.txt"
    cpambb_file = f"{res_dir}/{phase}_CPAMBB.txt"
    rlog_file = f"{res_dir}/{phase}_Rlog_3yr.txt"
    rlog_1yr_file = f"{res_dir}/{phase}_Rlog_1yr.txt"
    rlog_nosnap_file = f"{res_dir}/{phase}_Rlog_NoSnap.txt"
    
    mvzd_data = load_results(mvzd_file)
    cpambb_data = load_results(cpambb_file)
    rlog_data = load_results(rlog_file)
    rlog1_data = load_results(rlog_1yr_file)
    rlogns_data = load_results(rlog_nosnap_file)
    
    keys = list(mvzd_data.keys())
    mismatches = 0
    for k in keys:
        mvzd_res = mvzd_data[k]
        cpambb_res = cpambb_data.get(k, mvzd_res)
        rlog_res = rlog_data.get(k, mvzd_res)
        rlog1_res = rlog1_data.get(k, mvzd_res)
        rlogns_res = rlogns_data.get(k, mvzd_res)
        
        # Check if all match
        if not (mvzd_res == cpambb_res == rlog_res == rlog1_res == rlogns_res):
            mismatches += 1
            mismatch_details.append(f"[{phase}] {k}: MVZD={mvzd_res}, CPAMBB={cpambb_res}, Rlog3={rlog_res}, Rlog1={rlog1_res}, RlogNS={rlogns_res}")
            all_match = False
    
    total_queries_checked += len(keys)

if all_match:
    print(f"ALL VERIFIED CORRECTLY! Checked {total_queries_checked} queries across {len(phases)} phases.")
else:
    print(f"ERRORS FOUND! {len(mismatch_details)} mismatches.")
    for d in mismatch_details[:20]:
        print(d)
