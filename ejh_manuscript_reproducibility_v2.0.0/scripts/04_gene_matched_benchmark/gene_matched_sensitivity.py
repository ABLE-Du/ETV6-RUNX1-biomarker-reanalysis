# -*- coding: utf-8 -*-
"""
(1) Gene-matched TARGET sensitivity: recompute every pairwise correspondence on the
    exact 113 genes measurable in TARGET, so the benchmark is no longer gene-set-mismatched.
    This does NOT replace the 219-gene primary result.
(2) Li95 published clinical phenotype (age / WBC) - Fisher tests on the official source table.
(3) COG552 PAM C1/C2 split re-tabulated from the published source data.
"""
import pandas as pd, numpy as np, json, sys, io, os
from pathlib import Path
from scipy.stats import spearmanr, fisher_exact
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
RNG = np.random.default_rng(20260919)

b = pd.read_csv(Path(os.environ.get("EJH_TARGET_BENCHMARK", "results/B6_TARGET_benchmark.tsv")), sep="\t",
                keep_default_na=False)
for c in ["oksa_coefficient_slow_vs_fast", "li_C1_minus_C2",
          "target_primary_log2FC", "target_protocol_log2FC"]:
    b[c] = pd.to_numeric(b[c].replace("", np.nan), errors="coerce")

def boot_rho(x, y, nboot=10000, seed=3):
    x = np.asarray(x, float); y = np.asarray(y, float); n = len(x)
    rng = np.random.default_rng(seed); out = np.empty(nboot)
    for i in range(nboot):
        idx = rng.integers(0, n, n)
        if len(set(x[idx])) < 2 or len(set(y[idx])) < 2:
            out[i] = np.nan; continue
        out[i] = spearmanr(x[idx], y[idx]).statistic
    return np.nanpercentile(out, 2.5), np.nanpercentile(out, 97.5)

def perm_rho(x, y, nperm=20000, seed=4):
    x = np.asarray(x, float); y = np.asarray(y, float)
    obs = spearmanr(x, y).statistic
    rng = np.random.default_rng(seed); null = np.empty(nperm)
    for i in range(nperm):
        null[i] = spearmanr(x, rng.permutation(y)).statistic
    p = (np.sum(np.abs(null) >= abs(obs)) + 1) / (nperm + 1)
    return obs, null.mean(), null.std(ddof=1), p

def concord(x, y):
    x = np.asarray(x, float); y = np.asarray(y, float)
    m = np.sign(x) == np.sign(y)
    return 100.0 * m.mean()

rows = []
# ---------- primary (frozen, for reference only) ----------
p219 = b.dropna(subset=["oksa_coefficient_slow_vs_fast", "li_C1_minus_C2"])
rows.append(dict(comparison="Oksa vs Li194 (PRIMARY, frozen)", gene_universe="219 (Oksa x Li intersection)",
                 n=len(p219), rho=spearmanr(p219.oksa_coefficient_slow_vs_fast, p219.li_C1_minus_C2).statistic,
                 concordance_pct=concord(p219.oksa_coefficient_slow_vs_fast, p219.li_C1_minus_C2),
                 rho_CI_lo=np.nan, rho_CI_hi=np.nan, P=np.nan, note="unchanged primary result"))

# ---------- gene-matched 113 ----------
tg = b.dropna(subset=["target_primary_log2FC"]).copy()
print("genes measurable in TARGET:", len(tg))
sets = {
    "Oksa vs Li194": ("oksa_coefficient_slow_vs_fast", "li_C1_minus_C2"),
    "Oksa vs TARGET": ("oksa_coefficient_slow_vs_fast", "target_primary_log2FC"),
    "Li194 vs TARGET": ("li_C1_minus_C2", "target_primary_log2FC"),
    "Oksa vs TARGET (protocol-adjusted)": ("oksa_coefficient_slow_vs_fast", "target_protocol_log2FC"),
}
for label, (cx, cy) in sets.items():
    d = tg.dropna(subset=[cx, cy])
    x = d[cx].values; y = d[cy].values
    rho, nm, ns, p = perm_rho(x, y)
    lo, hi = boot_rho(x, y)
    rows.append(dict(comparison=label, gene_universe="113 (gene-matched to TARGET)", n=len(d),
                     rho=rho, concordance_pct=concord(x, y), rho_CI_lo=lo, rho_CI_hi=hi,
                     P=p, note="null: gene-identity permutation, 20,000"))
G = pd.DataFrame(rows)
G.to_csv("TARGET_113_GENE_MATCHED_SENSITIVITY.tsv", sep="\t", index=False)
print()
print(G.to_string(index=False))

# ---------------------------------------------------------------- Li95 clinical
print()
print("=" * 70)
print("Li95 validation-cohort clinical phenotype (re-tabulated from official Source Data)")
print("=" * 70)
tab = [("age >= 5 y", 33, 10, 28, 24, "C1 enriched (adverse = older)"),
       ("WBC >= 50 x10^9/L", 7, 11, 53, 23, "C2 enriched (favourable = higher WBC)")]
cl = []
for name, a, bb, c, d, note in tab:
    orv, p = fisher_exact([[a, bb], [c, d]])
    import statsmodels.api as sm
    lo, hi = sm.stats.Table2x2(np.array([[a, bb], [c, d]])).oddsratio_confint()
    cl.append(dict(cohort="Li95 validation (n=95)", feature=name, C1_positive=a, C1_negative=bb,
                   C2_positive=c, C2_negative=d, OR_C1_vs_C2=orv, CI_lo=lo, CI_hi=hi, P=p,
                   direction=note))
CL = pd.DataFrame(cl)
print(CL.to_string(index=False))

# ---------------------------------------------------------------- COG552
print()
print("=" * 70)
print("COG552 PAM classification (re-tabulated from published source data)")
print("=" * 70)
c5 = pd.read_csv("li95_out/COG552_pam_calls.tsv", sep="\t", keep_default_na=False)
c5["call"] = np.where(pd.to_numeric(c5.PAM_C1) > pd.to_numeric(c5.PAM_C2), "C1", "C2")
vc = c5["call"].value_counts().to_dict()
print("n =", len(c5), " calls:", vc)
print("NOTE: source data contain UMAP coordinates + PAM posterior probabilities only.")
print("      No expression matrix and no clinical variables are available ->")
print("      COG552 is used as PUBLISHED_INDEPENDENT_CONTEXT only.")

json.dump({
    "gene_matched_113": G.to_dict("records"),
    "li95_clinical": CL.to_dict("records"),
    "cog552": {"n": int(len(c5)), "calls": {k: int(v) for k, v in vc.items()},
               "data_available": "UMAP + PAM posteriors only; no expression, no clinical"},
}, open("_gene_matched_status.json", "w"), indent=1, default=str)
print("\nwritten TARGET_113_GENE_MATCHED_SENSITIVITY.tsv")
