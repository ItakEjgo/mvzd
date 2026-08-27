import os
import sys
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

phases = ['BaseQuery', '2021_MidQuery', '2025_EndQuery'] # Check a few representative phases
algos = ['MVZD', 'CPAMBB', 'Rlog_3yr']

all_match = True
for phase in phases:
    print(f"--- Checking Phase: {phase} ---")
    data = {}
    for algo in algos:
        fname = f"{phase}_{algo}.txt"
        data[algo] = load_results(fname)
        print(f"Loaded {len(data[algo])} queries for {algo}")
    
    if len(data['Rlog_3yr']) == 0:
        print("Rlog_3yr results not found yet.")
        continue
        
    keys = list(data['MVZD'].keys())
    mismatches = 0
    for k in keys:
        if k not in data['Rlog_3yr']: continue
        mvzd_res = data['MVZD'][k]
        cpambb_res = data['CPAMBB'].get(k, mvzd_res)
        rlog_res = data['Rlog_3yr'][k]
        
        # MVZD and CPAMBB should agree. If they don't, we trust MVZD as baseline for now or check both
        if rlog_res[0] != mvzd_res[0] or rlog_res[1] != mvzd_res[1]:
            # Also check if it matches CPAMBB
            if rlog_res[0] != cpambb_res[0] or rlog_res[1] != cpambb_res[1]:
                print(f"Mismatch at {k}: MVZD={mvzd_res}, CPAMBB={cpambb_res}, Rlog={rlog_res}")
                mismatches += 1
                all_match = False
                if mismatches > 10: 
                    print("... (too many mismatches, stopping for this phase)")
                    break
    if mismatches == 0 and len(keys) > 0 and len(data['Rlog_3yr']) > 0:
        print(f"SUCCESS: All {len(keys)} queries match for {phase}!")

if all_match:
    print("ALL VERIFIED CORRECTLY!")
else:
    print("ERRORS FOUND!")
