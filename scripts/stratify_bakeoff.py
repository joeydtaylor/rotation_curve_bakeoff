#!/usr/bin/env python3
import argparse, json, numpy as np, pandas as pd
from pathlib import Path

def load_bakeoff(p):
    df = pd.read_csv(p)
    # join key
    if "galaxy" in df.columns:
        df["__KEY__"] = df["galaxy"].astype(str)
    elif "ID" in df.columns:
        df["__KEY__"] = df["ID"].astype(str)
    else:
        raise SystemExit("bakeoff.csv must have column 'galaxy' or 'ID'")
    # numeric deltas
    for c in ("dBIC","dAICc"):
        df[c] = pd.to_numeric(df[c], errors="coerce")
    return df

def load_wilsons(csv_p, json_p):
    df = pd.read_csv(csv_p)
    df["__KEY__"] = df["ID"].astype(str)
    # numeric where present
    for c in ("n","n_outer","r2_outer","frac_used","R_max","finite_outer"):
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")

    # augment from JSON (R_max etc. live there)
    meta = json.loads(Path(json_p).read_text())

    def pull(col):
        out = {}
        for k, v in meta.items():
            if isinstance(v, dict) and col in v:
                try:
                    out[k] = float(v[col])
                except Exception:
                    out[k] = np.nan
            else:
                out[k] = np.nan
        return out

    for col in ("R_max","finite_outer","n_outer","r2_outer","frac_used"):
        if col not in df.columns or df[col].isna().all():
            lk = pull(col)
            df[col] = df["__KEY__"].map(lambda k: lk.get(k, np.nan))

    return df

def qtercile_bins(series: pd.Series):
    """Return a category Series with up to tercile bins; None if degenerate."""
    s = pd.to_numeric(series, errors="coerce")
    s = s[np.isfinite(s)]
    if s.size < 6:
        return None, None
    try:
        bins = pd.qcut(series, q=3, labels=["Q1","Q2","Q3"], duplicates="drop")
        # need at least 2 distinct bins to be meaningful
        if hasattr(bins, "cat") and len(bins.cat.categories) >= 2:
            return bins, list(bins.cat.categories)
        return None, None
    except Exception:
        return None, None

def summarize(sub: pd.DataFrame, key="dBIC"):
    s = pd.to_numeric(sub[key], errors="coerce")
    return dict(
        n   = int(s.notna().sum()),
        win = int((s > 0).sum()),
        pct = float((s > 0).sum() / max(1, len(s))),
        q1  = float(s.quantile(0.25)),
        med = float(s.quantile(0.50)),
        q3  = float(s.quantile(0.75)),
    )

def main():
    ap = argparse.ArgumentParser(description="Stratify ΔBIC/ΔAICc by Wilsons covariates")
    ap.add_argument("--bakeoff", default="bakeoff.csv")
    ap.add_argument("--wilsons-csv",  required=True)
    ap.add_argument("--wilsons-json", required=True)
    ap.add_argument("--out", default="out/tables/report/strata.csv")
    args = ap.parse_args()

    bo = load_bakeoff(args.bakeoff)
    w  = load_wilsons(args.wilsons_csv, args.wilsons_json)

    # candidate covariates; keep those with data
    candidates = ["n_outer","r2_outer","R_max","finite_outer"]
    covs = [c for c in candidates if c in w.columns and np.isfinite(w[c]).sum() >= 6]
    if not covs:
        raise SystemExit("No usable covariates found (need any of n_outer, r2_outer, R_max, finite_outer)")

    # merge once
    m = bo.merge(w[["__KEY__"] + covs], on="__KEY__", how="left")

    rows = []
    for cov in covs:
        bins, labels = qtercile_bins(m[cov])
        if bins is None:
            # skip degenerate covariate (e.g., all values equal)
            continue
        m[cov+"_bin"] = bins
        for b in labels:
            sub = m[m[cov+"_bin"] == b]
            if len(sub) == 0:
                continue
            sb = summarize(sub, "dBIC")
            sa = summarize(sub, "dAICc")
            rows.append([cov, str(b), sb["n"], sb["win"], sb["pct"], sb["q1"], sb["med"], sb["q3"],
                                   sa["n"], sa["win"], sa["pct"], sa["q1"], sa["med"], sa["q3"]])

    out = pd.DataFrame(rows, columns=[
        "covariate","bin",
        "n_BIC","wins_BIC","wins_pct_BIC","q1_BIC","med_BIC","q3_BIC",
        "n_AICc","wins_AICc","wins_pct_AICc","q1_AICc","med_AICc","q3_AICc"
    ])
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out, index=False)
    print(args.out)

if __name__ == "__main__":
    main()
