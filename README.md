# Rotation-Curve Bakeoff — EGR screen-channel vs. ΛCDM

*Reproducible model-selection + IR EFT Wilson extraction*

This repo contains everything needed to reproduce the results we cite in the paper:

* head-to-head fits of **EGR** (causal screen-channel) vs **ΛCDM** halo on galaxy rotation curves,
* information-criterion comparisons (BIC/AICc) with evidence bins and effect sizes,
* **IR EFT Wilsons** $\{c_{C^2},\alpha_{\rm Ric}\}$ extracted from the EGR small-$k$ limit,
* robustness checks stratified by data-quality covariates.

The repo ships with CSV artifacts so anyone can regenerate the published tables exactly. Scripts also support full re-runs if you replace inputs.

---

## Contents

```
rotation_curve_bakeoff/
├─ bakeoff.csv                 # per-galaxy ΔBIC, ΔAICc (same rows used for the paper)
├─ wilsons_fixedmu.csv         # per-galaxy small-k regression + mapped Wilsons (EGR)
├─ wilsons_fixedmu.json        # metadata for stratification (R_max, r2_outer, etc.)
├─ scripts/
│  ├─ bakeoff_report.py        # evidence bins, win rates, quartiles → CSV
│  ├─ bakeoff_wilsons_report.py# fuse bakeoff + Wilsons → summary CSV + MD/TeX
│  ├─ bootstrap_ci_for_hl.py   # Hodges–Lehmann effect sizes + bootstrap CIs
│  └─ stratify_bakeoff.py      # tercile stratification by data-quality covariates
└─ out/tables/report/          # generated outputs (created on first run)
```

---

## Environment

* Python ≥ 3.10
* `pandas`, `numpy`, `scipy` (for Wilcoxon), optionally `matplotlib` if you plot extras.

Install:

```bash
pip install pandas numpy scipy
```

---

## Definitions (used throughout)

* **ΔBIC** := `BIC_EGR − BIC_LCDM`. Positive favors EGR.
  Same sign convention for **ΔAICc**.
* **Jeffreys bins** for ΔBIC:
  `(-∞,−10], (−10,−6], (−6,−2], (−2,2], (2,6], (6,10], (10,∞)`
  Reported as: `LCDM≫EGR`, `LCDM>`, `LCDM>`, `~tie`, `EGR>`, `EGR>`, `EGR≫`.
* **Small-k Wilson extraction (EGR only):**
  For each galaxy we compute $\mu_{\rm eff} = V_{\rm fit}^2/V_{\rm bar}^2$, map $R\to k=\pi/R$, select an outer-R window, and regress $\delta\mu=\mu_{\rm eff}-1$ onto $\{k^2, k^4, k^4\log(k^2/\mu^2)\}$ with a **global $\mu=0.3001615492546175$** and **ridge $=10^{-4}$**. The coefficients map to $\{c_{C^2},\alpha_{\rm Ric}\}$ (units $M_{\rm Pl}^2$). The scalar combo $B$ is reported as a **bound**, not a point estimate.

---

## Quick start (reproduce the paper tables)

From repo root:

```bash
# 1) Evidence bins, win rates, quartiles
python scripts/bakeoff_report.py
# -> out/tables/report/bakeoff_evidence_table.csv

# 2) Fuse bakeoff + Wilsons; emit compact summary + Markdown/LaTeX blocks
python scripts/bakeoff_wilsons_report.py \
  --bakeoff bakeoff.csv \
  --wilsons-csv  wilsons_fixedmu.csv \
  --wilsons-json wilsons_fixedmu.json \
  --outdir out/tables/report
# -> out/tables/report/{bakeoff_evidence_table.csv,wilsons_summary_compact.csv,report.md,report.tex}

# 3) Effect sizes (Hodges–Lehmann) with bootstrap CIs
python scripts/bootstrap_ci_for_hl.py
# -> prints HL and 95% HL-CIs for ΔBIC and ΔAICc (n excludes NaNs)

# 4) Robustness: stratify wins and Δ-quartiles by data-quality terciles
python scripts/stratify_bakeoff.py \
  --bakeoff bakeoff.csv \
  --wilsons-csv  wilsons_fixedmu.csv \
  --wilsons-json wilsons_fixedmu.json \
  --out out/tables/report/strata.csv
# -> out/tables/report/strata.csv
```

**What to cite:**

* `out/tables/report/bakeoff_evidence_table.csv` — counts and percentages per Jeffreys bin; win rates; quartiles.
* `out/tables/report/wilsons_summary_compact.csv` — medians ± MAD for $c_{C^2}$ and $\alpha_{\rm Ric}$ (units $M_{\rm Pl}^2$); include 84th percentile of $|B|$ as the scalar bound.
* `out/tables/report/strata.csv` — BIC/AICc win-rates and Δ-quartiles across terciles of `n_outer`, `r2_outer`, `R_max`, `finite_outer`.
* `out/tables/report/report.md` / `report.tex` — drop-in text blocks for Results/Methods.

---

## Regenerating inputs (optional)

The repo includes `bakeoff.csv` and `wilsons_fixedmu.{csv,json}` to ensure exact reproducibility.
To regenerate them from scratch:

1. **EGR fits + NPZ export** in your rotation-curve fitter repo (not this one):
   run the EGR fitter, ensure it writes `*_rc_arrays.npz` per galaxy.

2. **Extract small-k Wilsons** there using the same fixed-μ settings:

```bash
python scripts/extract_wilsons.py \
  --tables out/tables \
  --mpl 1.0 \
  --mu 0.3001615492546175 \
  --ridge 1e-4 \
  --outcsv  wilsons_fixedmu.csv \
  --outjson wilsons_fixedmu.json
```

3. Copy `wilsons_fixedmu.csv` and `wilsons_fixedmu.json` into this repo root and rerun the steps in **Quick start**.

4. To create or update `bakeoff.csv`, compute per-galaxy `dBIC` and `dAICc` as defined above with consistent sign and write a CSV with columns at least:

```
galaxy, dBIC, dAICc
```

---

## Column schemas

### `bakeoff.csv`

* `galaxy` or `ID` (string)
* `dBIC` (float) = `BIC_EGR − BIC_LCDM`
* `dAICc` (float) = `AICc_EGR − AICc_LCDM`

NaNs are allowed; they’re excluded from counts and effect sizes.

### `wilsons_fixedmu.csv`

Per galaxy (EGR only):
`ID, n, n_outer, r2_outer, frac_used, A, A_err, B, B_err, alpha_log, alpha_log_err, cC2, cC2_err, alpha_Ric, alpha_Ric_err`
Medians ± MAD for the mapped Wilsons are produced by `bakeoff_wilsons_report.py`.

### `wilsons_fixedmu.json`

Metadata keyed by `ID`. Expected fields include `R_max`, `finite_outer`, `n_outer`, `r2_outer`, `mu_scale` (used for stratification and QA).

---

## Troubleshooting

* **Bins don’t sum to N:** you have NaNs or ±inf in `dBIC`. The report prints the unbinned count; fix upstream or accept exclusion.
* **Join errors in stratification:** ensure `wilsons_fixedmu.json` is the mate of your CSV; the script falls back to JSON for missing covariates.
* **Different numbers than the paper:** verify sign convention for Δ’s, the fixed global `--mu 0.3001615492546175`, and `--ridge 1e-4`.
