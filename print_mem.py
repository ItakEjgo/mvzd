import pandas as pd
for year in [2021, 2025]:
    print(f"--- YEAR {year} MEMORY (MB) ---")
    for algo in ["MVZD", "CPAMBB", "Rlog_3yr"]:
        try:
            df_fwd = pd.read_csv(f"results_bp_50/{year}_Forward_Batch_{algo}.txt", sep="|", skipinitialspace=True)
            df_fwd.columns = [c.strip() for c in df_fwd.columns]
            df_fwd['Batch'] = df_fwd['Batch'].astype(str).str.strip()
            mem = df_fwd[df_fwd['Batch'] == 'SUMMARY']['Mem_MB'].values[0]
            print(f"{algo}: {mem:.2f} MB")
        except Exception as e:
            print(f"{algo}: Error {e}")
