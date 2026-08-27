import pandas as pd
import glob

for year in [2021, 2025]:
    print(f"--- YEAR {year} MEMORY (MB) ---")
    files = glob.glob(f"results_bp_50/{year}_Forward_Batch_*.txt")
    for f in sorted(files):
        algo = f.split('_')[-1].split('.')[0]
        # skip numbers
        if '1yr' in f or '2yr' in f or '4yr' in f or '5yr' in f: continue
        try:
            df = pd.read_csv(f, sep="|", skipinitialspace=True)
            df.columns = [c.strip() for c in df.columns]
            df['Batch'] = df['Batch'].astype(str).str.strip()
            mem = df[df['Batch'] == 'SUMMARY']['Mem_MB'].values[0]
            print(f"{algo:15}: {mem:.2f} MB")
        except Exception as e:
            pass
