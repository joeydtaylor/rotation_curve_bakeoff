import pandas as pd, numpy as np
from itertools import combinations_with_replacement

df = pd.read_csv("bakeoff.csv")

def hodges_lehmann(x):
    x = np.asarray(x, float)
    x = x[np.isfinite(x) & (x!=0)]
    # Walsh averages via vectorized add.outer (avoid Python loop)
    W = (x[:,None] + x[None,:]) / 2.0
    return float(np.median(W))

def hl_boot_ci(x, B=1500, alpha=0.05, seed=0):
    x = pd.to_numeric(x, errors="coerce").dropna().values
    x = x[x!=0]
    rng = np.random.default_rng(seed)
    hl0 = hodges_lehmann(x)
    hls = np.empty(B, float)
    n = x.size
    for b in range(B):
        xb = rng.choice(x, size=n, replace=True)
        hls[b] = hodges_lehmann(xb)
    lo, hi = np.percentile(hls, [100*alpha/2, 100*(1-alpha/2)])
    return hl0, float(lo), float(hi), int(n)

for col in ("dBIC","dAICc"):
    hl, lo, hi, n = hl_boot_ci(df[col], B=1500, alpha=0.05, seed=42)
    print(f"{col}: HL={hl:.3f}, 95% HL-CI≈[{lo:.3f}, {hi:.3f}], n={n}")
