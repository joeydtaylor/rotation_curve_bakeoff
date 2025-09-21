#!/usr/bin/env python3
import argparse, json, numpy as np, pandas as pd
from pathlib import Path
from textwrap import dedent

def load_bakeoff(p):
    df = pd.read_csv(p)
    for c in ("dBIC","dAICc"):
        df[c] = pd.to_numeric(df[c], errors="coerce")
    return df

def load_wilsons(csv_p, json_p):
    df = pd.read_csv(csv_p)
    meta = json.loads(Path(json_p).read_text())
    pass_ids = [k for k,v in meta.items() if v.get("status")=="PASS"]
    df = df[df["ID"].isin(pass_ids)].copy()
    # force numeric
    num = ["n","n_outer","r2_outer","frac_used","A","A_err","B","B_err",
           "alpha_log","alpha_log_err","cC2","cC2_err","alpha_Ric","alpha_Ric_err"]
    for c in num:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")
    return df

def evidence_table(df):
    n = len(df)
    edges  = [-np.inf,-10,-6,-2,2,6,10,np.inf]
    labels = ["LCDM≫EGR (≥10)","LCDM>EGR (6–10)","LCDM>EGR (2–6)","~tie (±2)",
              "EGR>LCDM (2–6)","EGR>LCDM (6–10)","EGR≫LCDM (≥10)"]
    cats = pd.cut(df["dBIC"], bins=edges, labels=labels, include_lowest=True)
    tbl = cats.value_counts(dropna=True).reindex(labels).fillna(0).astype(int)
    missing = int(n - tbl.sum())
    return tbl, missing

def robust(x):
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    if x.size == 0:
        return dict(n=0, median=np.nan, mad=np.nan, p16=np.nan, p84=np.nan)
    med = float(np.median(x))
    mad = float(1.4826 * np.median(np.abs(x - med)))
    p16, p84 = np.percentile(x, [16,84])
    return dict(n=int(x.size), median=med, mad=mad, p16=float(p16), p84=float(p84))

def signed_rank(x):
    # Wilcoxon signed-rank on ΔBIC and ΔAICc (approx, no ties handling beyond drop zeros)
    try:
        from scipy.stats import wilcoxon
        m = np.isfinite(x) & (x != 0)
        if m.sum() < 10: return np.nan
        return float(wilcoxon(x[m]).pvalue)
    except Exception:
        return np.nan

def main():
    ap = argparse.ArgumentParser(description="Fuse bakeoff with Wilsons; emit report artifacts")
    ap.add_argument("--bakeoff", default="bakeoff.csv")
    ap.add_argument("--wilsons-csv",  default="out/tables/wilsons_fixedmu.csv")
    ap.add_argument("--wilsons-json", default="out/tables/wilsons_fixedmu.json")
    ap.add_argument("--outdir", default="out/tables/report")
    args = ap.parse_args()

    outdir = Path(args.outdir); outdir.mkdir(parents=True, exist_ok=True)

    # Load
    bo = load_bakeoff(args.bakeoff)
    w  = load_wilsons(args.wilsons_csv, args.wilsons_json)

    # Evidence
    n = len(bo)
    tbl, missing = evidence_table(bo)
    wins_bic  = int((bo["dBIC"]  > 0).sum())
    wins_aicc = int((bo["dAICc"] > 0).sum())
    qB = bo["dBIC"].quantile([.25,.5,.75]).to_dict()
    qA = bo["dAICc"].quantile([.25,.5,.75]).to_dict()
    p_wilcoxon_bic  = signed_rank(bo["dBIC"].to_numpy())
    p_wilcoxon_aicc = signed_rank(bo["dAICc"].to_numpy())

    # Wilson aggregates (robust)
    R_cC2  = robust(w["cC2"].values)
    R_aRic = robust(w["alpha_Ric"].values)
    R_B    = robust(w["B"].values)

    # Save tables
    ev = pd.DataFrame({"Bin": tbl.index, "Count": tbl.values,
                       "Percent": [f"{(v/n):.1%}" for v in tbl.values]})
    ev.to_csv(outdir/"bakeoff_evidence_table.csv", index=False)
    agg = pd.DataFrame([
        ["cC2",       R_cC2["n"],  R_cC2["median"],  R_cC2["mad"],  R_cC2["p16"],  R_cC2["p84"]],
        ["alpha_Ric", R_aRic["n"], R_aRic["median"], R_aRic["mad"], R_aRic["p16"], R_aRic["p84"]],
        ["B",         R_B["n"],    R_B["median"],    R_B["mad"],    R_B["p16"],    R_B["p84"]],
    ], columns=["param","n","median","MAD","p16","p84"])
    agg.to_csv(outdir/"wilsons_summary_compact.csv", index=False)

    # Markdown block
    md = dedent(f"""
    ### Model comparison (bakeoff)
    - **EGR wins (BIC):** {wins_bic}/{n} ({wins_bic/n:.1%}); **EGR wins (AICc):** {wins_aicc}/{n} ({wins_aicc/n:.1%}).
    - Evidence bins (Jeffreys, ΔBIC):  
      {", ".join([f"{row.Bin}: {row.Count} ({row.Percent})" for _,row in ev.iterrows()])}.
    - ΔBIC quartiles: Q1={qB[0.25]:.2f}, Med={qB[0.5]:.2f}, Q3={qB[0.75]:.2f}; ΔAICc quartiles: Q1={qA[0.25]:.2f}, Med={qA[0.5]:.2f}, Q3={qA[0.75]:.2f}.
    - Signed-rank (Wilcoxon) vs 0: p_BIC={p_wilcoxon_bic:.2e}, p_AICc={p_wilcoxon_aicc:.2e}.

    ### IR Wilsons (fixed-μ extraction; medians ± MAD; units M_Pl^2)
    - c_C^2 = {R_cC2["median"]:.3g} ± {R_cC2["mad"]:.3g}  [p16={R_cC2["p16"]:.3g}, p84={R_cC2["p84"]:.3g}] (n={R_cC2["n"]})
    - α_Ric = {R_aRic["median"]:.3g} ± {R_aRic["mad"]:.3g}  [p16={R_aRic["p16"]:.3g}, p84={R_aRic["p84"]:.3g}] (n={R_aRic["n"]})
    - Scalar combo bound: use |B| ≤ {np.percentile(np.abs(w["B"].dropna()),84):.3g}; median(B)={R_B["median"]:.3g}, MAD={R_B["mad"]:.3g} (n={R_B["n"]})
    """).strip()
    (outdir/"report.md").write_text(md)

    # LaTeX block
    tex = dedent(f"""
    % --- Bakeoff summary ---
    \\textbf{{EGR vs. $\\Lambda$CDM (BIC)}}: {wins_bic}/{n} ({wins_bic/n:.1%}).\\quad
    \\textbf{{EGR vs. $\\Lambda$CDM (AICc)}}: {wins_aicc}/{n} ({wins_aicc/n:.1%}).\\\\
    $\\Delta\\mathrm{{BIC}}$ quartiles: Q1={qB[0.25]:.2f}, med={qB[0.5]:.2f}, Q3={qB[0.75]:.2f}.\\quad
    $\\Delta\\mathrm{{AICc}}$ quartiles: Q1={qA[0.25]:.2f}, med={qA[0.5]:.2f}, Q3={qA[0.75]:.2f}.\\\\
    Wilcoxon vs 0: $p_{{\\rm BIC}}={p_wilcoxon_bic:.2e}$, $p_{{\\rm AICc}}={p_wilcoxon_aicc:.2e}$.\\\\[4pt]
    % --- Wilsons ---
    IR Wilsons (medians $\\pm$ MAD, units $M_\\mathrm{{Pl}}^2$):\\\\
    $c_{{C^2}} = {R_cC2["median"]:.3g} \\pm {R_cC2["mad"]:.3g}$\\,,\\qquad
    $\\alpha_{{\\rm Ric}} = {R_aRic["median"]:.3g} \\pm {R_aRic["mad"]:.3g}$.\\\\
    Scalar-sector bound: $|B| \\lesssim {np.percentile(np.abs(w["B"].dropna()),84):.3g}$ (84th percentile of $|B|$).
    """).strip()
    (outdir/"report.tex").write_text(tex)

    print(f"wrote: {outdir}/bakeoff_evidence_table.csv")
    print(f"wrote: {outdir}/wilsons_summary_compact.csv")
    print(f"wrote: {outdir}/report.md")
    print(f"wrote: {outdir}/report.tex")

if __name__ == "__main__":
    main()
