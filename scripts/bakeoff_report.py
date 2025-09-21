import pandas as pd, numpy as np

df = pd.read_csv("bakeoff.csv")  # must contain float columns: dBIC, dAICc
n = len(df)

# basic hygiene
for c in ("dBIC","dAICc"):
    df[c] = pd.to_numeric(df[c], errors="coerce")

edges  = [-np.inf,-10,-6,-2,2,6,10,np.inf]
labels = ["LCDM≫EGR (≥10)","LCDM>EGR (6–10)","LCDM>EGR (2–6)","~tie (±2)",
          "EGR>LCDM (2–6)","EGR>LCDM (6–10)","EGR≫LCDM (≥10)"]

cats = pd.cut(df["dBIC"], bins=edges, labels=labels, include_lowest=True)
tbl  = cats.value_counts(dropna=True).reindex(labels).fillna(0).astype(int)
missing = n - int(tbl.sum())

# wins
wins_bic  = int((df["dBIC"]  > 0).sum())
wins_aicc = int((df["dAICc"] > 0).sum())

# quartiles
qB = df["dBIC"].quantile([.25,.5,.75]).to_dict()
qA = df["dAICc"].quantile([.25,.5,.75]).to_dict()

# print table with percentages and missing
print("dBIC")
for k in labels:
    v = int(tbl[k]); print(f"{k:20s} {v:4d}  ({v/n:.1%})")
if missing:
    n_nan = int(df["dBIC"].isna().sum())
    n_inf = int(np.isinf(df["dBIC"]).sum())
    print(f"\nUnbinned (NaN/±inf): {missing}  [NaN={n_nan}, ±inf={n_inf}]")

print(f"\nEGR wins (BIC):   {wins_bic}/{n}  ({wins_bic/n:.1%})")
print(f"EGR wins (AICc):  {wins_aicc}/{n}  ({wins_aicc/n:.1%})")
print(f"\nQuartiles dBIC :  Q1={qB[0.25]:.2f}, Med={qB[0.5]:.2f}, Q3={qB[0.75]:.2f}")
print(f"Quartiles dAICc:  Q1={qA[0.25]:.2f}, Med={qA[0.5]:.2f}, Q3={qA[0.75]:.2f}")

# after running the block above
out = pd.DataFrame({
    "Bin": labels,
    "Count": [int(tbl[k]) for k in labels],
    "Percent": [f"{(tbl[k]/n):.1%}" for k in labels],
})
out.to_csv("bakeoff_evidence_table.csv", index=False)
print("bakeoff_evidence_table.csv")
