# -*- coding: utf-8 -*-
"""Paired bootstrap of the DIFFERENCE in rho on the same 113 genes:
   (Oksa-Li194)  minus  (Oksa-TARGET)  and  minus (Li194-TARGET)."""
import pandas as pd, numpy as np, json, sys, io, os
from pathlib import Path
from scipy.stats import spearmanr
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

b = pd.read_csv(Path(os.environ.get("EJH_TARGET_BENCHMARK", "results/B6_TARGET_benchmark.tsv")), sep="\t", keep_default_na=False)
for c in ["oksa_coefficient_slow_vs_fast", "li_C1_minus_C2", "target_primary_log2FC"]:
    b[c] = pd.to_numeric(b[c].replace("", np.nan), errors="coerce")
d = b.dropna(subset=["target_primary_log2FC", "oksa_coefficient_slow_vs_fast", "li_C1_minus_C2"])
X = d[["oksa_coefficient_slow_vs_fast", "li_C1_minus_C2", "target_primary_log2FC"]].values
n = len(X)
print("n genes =", n)

def rhos(M):
    return (spearmanr(M[:, 0], M[:, 1]).statistic,
            spearmanr(M[:, 0], M[:, 2]).statistic,
            spearmanr(M[:, 1], M[:, 2]).statistic)

obs = rhos(X)
print("observed  Oksa-Li = %.4f   Oksa-TARGET = %.4f   Li-TARGET = %.4f" % obs)

rng = np.random.default_rng(20260919)
NB = 20000
d1 = np.empty(NB); d2 = np.empty(NB)
for i in range(NB):
    idx = rng.integers(0, n, n)
    M = X[idx]
    if min(len(set(M[:, j])) for j in range(3)) < 2:
        d1[i] = np.nan; d2[i] = np.nan; continue
    a, bb, cc = rhos(M)
    d1[i] = a - bb
    d2[i] = a - cc
out = []
for lab, arr, o in [("rho(Oksa-Li194) - rho(Oksa-TARGET)", d1, obs[0] - obs[1]),
                    ("rho(Oksa-Li194) - rho(Li194-TARGET)", d2, obs[0] - obs[2])]:
    lo, hi = np.nanpercentile(arr, [2.5, 97.5])
    out.append(dict(comparison=lab, delta_rho=o, CI_lo=lo, CI_hi=hi,
                    frac_null_ge=float(np.nanmean(np.abs(arr) >= abs(o)))))
    print("%-42s delta = %.4f  95%% CI [%.4f, %.4f]" % (lab, o, lo, hi))
pd.DataFrame(out).to_csv("RHO_DIFFERENCE_BOOTSTRAP.tsv", sep="\t", index=False)
print("written RHO_DIFFERENCE_BOOTSTRAP.tsv")
