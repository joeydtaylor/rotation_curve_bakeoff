import pandas as pd

# Load full summaries with needed fields
egr  = pd.read_csv("egr_out/tables/summary.csv")
lcdm = pd.read_csv("lcdm_out/tables/summary.csv")

# Keep OK fits only (optional but recommended)
egr  = egr[egr["fit_status"] == "OK"]
lcdm = lcdm[lcdm["fit_status"] == "OK"]

# One row per ID: keep the best (lowest) BIC within each model
egr_best  = egr.sort_values(["ID","BIC"]).groupby("ID", as_index=False).first()
lcdm_best = lcdm.sort_values(["ID","BIC"]).groupby("ID", as_index=False).first()

# Select the columns we compare on
egr_best  = egr_best[["ID","BIC","AICc","s_frac","rho_AR1"]].rename(
    columns={"BIC":"BIC_EGR","AICc":"AICc_EGR","s_frac":"s_frac_EGR","rho_AR1":"rho_EGR"}
)
lcdm_best = lcdm_best[["ID","BIC","AICc","s_frac","rho_AR1"]].rename(
    columns={"BIC":"BIC_LCDM","AICc":"AICc_LCDM","s_frac":"s_frac_LCDM","rho_AR1":"rho_LCDM"}
)

# Inner join and deltas
m = pd.merge(egr_best, lcdm_best, on="ID", how="inner")
m["dBIC"]    = m["BIC_LCDM"]  - m["BIC_EGR"]
m["dAICc"]   = m["AICc_LCDM"] - m["AICc_EGR"]
m["d_sfrac"] = m["s_frac_LCDM"] - m["s_frac_EGR"]
m["d_rho"]   = m["rho_LCDM"] - m["rho_EGR"]

# Wins
wins_bic_egr   = int((m["dBIC"]  >  0).sum())
wins_bic_lcdm  = int((m["dBIC"]  <  0).sum())
wins_aicc_egr  = int((m["dAICc"] >  0).sum())
wins_aicc_lcdm = int((m["dAICc"] <  0).sum())

print(f"n galaxies: {len(m)}")
print(f"EGR wins (BIC):  {wins_bic_egr}/{len(m)}")
print(f"LCDM wins (BIC): {wins_bic_lcdm}/{len(m)}")
print(f"EGR wins (AICc): {wins_aicc_egr}/{len(m)}")
print(f"LCDM wins (AICc):{wins_aicc_lcdm}/{len(m)}")
print("\nQuartiles:")
print(m[["dBIC","dAICc","d_sfrac","d_rho"]].describe(percentiles=[.25,.5,.75]))

# Evidence bins (Kass–Raftery thresholds on ΔBIC)
bins = pd.cut(m["dBIC"], bins=[-1e9,-10,-6,-2,2,6,10,1e9],
              labels=["LCDM≫EGR","LCDM≫EGR(~6–10)","LCDM>EGR(~2–6)","~tie(±2)",
                      "EGR>LCDM(~2–6)","EGR≫LCDM(~6–10)","EGR≫LCDM"])
print("\nΔBIC evidence bins:\n", bins.value_counts().to_string())
m.sort_values("dBIC").to_csv("bakeoff_dedup.csv", index=False)
