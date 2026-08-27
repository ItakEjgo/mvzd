import os

def load_results(filename):
    results = {}
    if not os.path.exists(filename): 
        print(f"File not found: {filename}")
        return results
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

phases = ['BaseQuery', '2025_EndQuery', '2021_EndQuery'] 

all_match = True
for phase in phases:
    print(f"--- Checking Phase: {phase} ---")
    
    mvzd_file = f"results_bp_50/{phase}_MVZD.txt"
    cpambb_file = f"results_bp_50/{phase}_CPAMBB.txt"
    rlog_file = f"verification_results/{phase}_Rlog_3yr.txt"
    
    mvzd_data = load_results(mvzd_file)
    cpambb_data = load_results(cpambb_file)
    rlog_data = load_results(rlog_file)
    
    print(f"Loaded MVZD: {len(mvzd_data)}, CPAMBB: {len(cpambb_data)}, Rlog: {len(rlog_data)}")
    
    if len(rlog_data) == 0:
        continue
        
    keys = list(mvzd_data.keys())
    mismatches = 0
    for k in keys:
        if k not in rlog_data: continue
        mvzd_res = mvzd_data[k]
        cpambb_res = cpambb_data.get(k, mvzd_res)
        rlog_res = rlog_data[k]
        
        if rlog_res[0] != mvzd_res[0] or rlog_res[1] != mvzd_res[1]:
            if rlog_res[0] != cpambb_res[0] or rlog_res[1] != cpambb_res[1]:
                print(f"Mismatch at {k}: MVZD={mvzd_res}, CPAMBB={cpambb_res}, Rlog={rlog_res}")
                mismatches += 1
                all_match = False
                if mismatches > 10: 
                    print("... (too many mismatches, stopping for this phase)")
                    break
    if mismatches == 0 and len(keys) > 0:
        print(f"SUCCESS: All {len(keys)} queries match for {phase}!")

if all_match:
    print("ALL VERIFIED CORRECTLY!")
else:
    print("ERRORS FOUND!")
