# -*- coding: utf-8 -*-
"""
Module B - local maturation / immunophenotypic triangulation.

Pre-specified (from Li et al., NOT chosen here):
  adverse / C1 state  ->  HSC-like : higher CD34, higher GATA2, older age
  favourable / C2 state -> pre/pro-B-like : higher CD38, higher WBC

Local primary endpoint (frozen): day-19 MRD >= 0.1%  vs  < 0.1%.
No post-hoc cut-offs. No threshold search. No new score unless both markers are
measured in a substantial fraction (they are not - see output).
"""
import pandas as pd, numpy as np, sys, io, json, re, os
from pathlib import Path
from scipy.stats import mannwhitneyu, fisher_exact, rankdata
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

RNG = np.random.default_rng(20260919)
OUT = "."

# ---------------------------------------------------------------- load
coh = pd.read_csv(Path(os.environ["EJH_DEIDENTIFIED_COHORT"]),
                  sep="\t", keep_default_na=False, dtype=str)
flow = pd.read_csv("LOCAL_FLOW_2020_RAW.tsv", sep="\t", keep_default_na=False, dtype=str)

def f(x):
    x = ("" if x is None else str(x)).strip()
    if x in ("", "/", "-", "—", "无", "NA", "nan"):
        return np.nan
    try:
        return float(x)
    except Exception:
        return np.nan

# ---------------------------------------------------------------- 1. field audit
print("=" * 72)
print("[1] LOCAL_FLOW_FIELD_AUDIT")
print("=" * 72)
flow_ids = set(flow["study_id"])
audit_rows = []
for _, r in coh.iterrows():
    sid = r["study_id"]
    if sid in flow_ids:
        fr = flow[flow["study_id"] == sid].iloc[0]
        c34r, c38r = fr["CD34_raw"].strip(), fr["CD38_raw"].strip()
        c34_avail = c34r not in ("", "/", "-", "无")
        c38_avail = c38r not in ("", "/", "-", "无")
        audit_rows.append({
            "study_id": sid, "treatment_era": r["treatment_era"],
            "CD34_raw_available": c34_avail,
            "CD34_data_type": "percent_positive_blasts" if c34_avail else "not_recorded",
            "CD34_curated_value": f(c34r),
            "CD34_source": "2020 sheet col46 'CD34'" if c34_avail else "",
            "CD38_raw_available": c38_avail,
            "CD38_data_type": "percent_positive_blasts" if c38_avail else "not_recorded",
            "CD38_curated_value": f(c38r),
            "CD38_source": "2020 sheet col47 'CD38'" if c38_avail else "",
            "manual_review_required": False,
            "usable_for_analysis": c34_avail or c38_avail,
        })
    else:
        # eras 2005/2009, 2015, 2025 carry only a lineage string; no antigen columns exist
        audit_rows.append({
            "study_id": sid, "treatment_era": r["treatment_era"],
            "CD34_raw_available": False, "CD34_data_type": "not_recorded",
            "CD34_curated_value": np.nan, "CD34_source": "",
            "CD38_raw_available": False, "CD38_data_type": "not_recorded",
            "CD38_curated_value": np.nan, "CD38_source": "",
            "manual_review_required": False, "usable_for_analysis": False,
        })
aud = pd.DataFrame(audit_rows)
aud.to_csv("LOCAL_FLOW_FIELD_AUDIT.tsv", sep="\t", index=False)
print("rows:", len(aud))
print("CD34 available:", int(aud.CD34_raw_available.sum()),
      " CD38 available:", int(aud.CD38_raw_available.sum()))
print("by era:")
print(aud.groupby("treatment_era")[["CD34_raw_available", "CD38_raw_available"]].sum())

# ---------------------------------------------------------------- helpers
def bh(p):
    p = np.asarray(p, float); n = len(p); o = np.argsort(p)
    q = np.empty(n); prev = 1.0
    for i in range(n - 1, -1, -1):
        idx = o[i]
        prev = min(prev, p[idx] * n / (i + 1))
        q[idx] = prev
    return q

def cliffs_delta(a, b):
    """P(a>b) - P(a<b).  a = poor-response group, b = good-response group."""
    a = np.asarray(a, float); b = np.asarray(b, float)
    n1, n2 = len(a), len(b)
    if n1 == 0 or n2 == 0:
        return np.nan
    comb = np.concatenate([a, b])
    r = rankdata(comb)
    R1 = r[:n1].sum()
    U = R1 - n1 * (n1 + 1) / 2.0
    return 2.0 * U / (n1 * n2) - 1.0

def hl_shift(a, b, nboot=20000, alpha=0.05, seed=1):
    """Hodges-Lehmann location shift (median of all pairwise a_i - b_j) with
    bootstrap percentile CI."""
    a = np.asarray(a, float); b = np.asarray(b, float)
    d = (a[:, None] - b[None, :]).ravel()
    hl = np.median(d)
    rng = np.random.default_rng(seed)
    bs = np.empty(nboot)
    for i in range(nboot):
        aa = rng.choice(a, len(a), replace=True)
        bb = rng.choice(b, len(b), replace=True)
        bs[i] = np.median((aa[:, None] - bb[None, :]).ravel())
    return hl, np.percentile(bs, 100 * alpha / 2), np.percentile(bs, 100 * (1 - alpha / 2))

def cliffs_ci(a, b, nboot=20000, alpha=0.05, seed=2):
    rng = np.random.default_rng(seed)
    a = np.asarray(a, float); b = np.asarray(b, float)
    bs = np.empty(nboot)
    for i in range(nboot):
        bs[i] = cliffs_delta(rng.choice(a, len(a), replace=True),
                             rng.choice(b, len(b), replace=True))
    return np.percentile(bs, 100 * alpha / 2), np.percentile(bs, 100 * (1 - alpha / 2))

# ---------------------------------------------------------------- 2. CD34 / CD38
print()
print("=" * 72)
print("[2] CD34 / CD38  vs  day-19 MRD >= 0.1%")
print("=" * 72)
m = aud.merge(coh[["study_id", "D19_MRD_ge_0p1"]], on="study_id")
m["poor"] = m["D19_MRD_ge_0p1"].str.upper() == "TRUE"
sub = m[m.usable_for_analysis].copy()
print("flow-evaluable subset n =", len(sub), " eras:", dict(sub.treatment_era.value_counts()))
print("D19 MRD evaluable within subset:", int((sub.D19_MRD_ge_0p1 != "").sum()))

res = []
for marker, expected in [("CD34", "higher in poor response"), ("CD38", "lower in poor response")]:
    colv = "CD34_curated_value" if marker == "CD34" else "CD38_curated_value"
    d = sub[sub[colv].notna() & (sub.D19_MRD_ge_0p1 != "")]
    poor = d.loc[d.poor, colv].astype(float).values
    good = d.loc[~d.poor, colv].astype(float).values
    n_eval = len(d)
    n_poor, n_good = len(poor), len(good)
    missing = int(len(aud) - d[colv].notna().sum())
    if n_poor >= 1 and n_good >= 1:
        U, p = mannwhitneyu(poor, good, alternative="two-sided")
        hl, lo, hi = hl_shift(poor, good)
        dl = cliffs_delta(poor, good)
        dlo, dhi = cliffs_ci(poor, good)
    else:
        U, p, hl, lo, hi, dl, dlo, dhi = [np.nan] * 8
    res.append(dict(marker=marker, n_evaluable=n_eval, n_poor=n_poor, n_good=n_good,
                    n_missing_in_cohort=missing,
                    median_poor=np.median(poor) if n_poor else np.nan,
                    median_good=np.median(good) if n_good else np.nan,
                    HL_shift=hl, HL_CI_lo=lo, HL_CI_hi=hi,
                    cliffs_delta=dl, cliffs_CI_lo=dlo, cliffs_CI_hi=dhi,
                    U=U, P=p, expected_direction=expected,
                    observed_direction=("higher in poor" if np.median(poor) > np.median(good) else "lower in poor") if n_poor and n_good else "NA",
                    era_strata=dict(d.treatment_era.value_counts()).__str__()))
R = pd.DataFrame(res)
R["BH_q_2markers"] = bh(R["P"].values)
R.to_csv("LOCAL_CD34_CD38_RESULTS.tsv", sep="\t", index=False)
print(R[["marker", "n_evaluable", "n_poor", "n_good", "median_poor", "median_good",
         "HL_shift", "HL_CI_lo", "HL_CI_hi", "cliffs_delta", "cliffs_CI_lo", "cliffs_CI_hi",
         "P", "BH_q_2markers", "observed_direction"]].to_string(index=False))

# composite gate
n_both = int((aud.CD34_raw_available & aud.CD38_raw_available).sum())
print("\n[n_both markers measured] =", n_both, "of", len(aud),
      "-> composite gate:", "BUILD" if n_both >= 0.7 * len(aud) else "NOT BUILD (insufficient joint coverage)")

# ---------------------------------------------------------------- 3. age / WBC
print()
print("=" * 72)
print("[3] Age >= 5 y   and   WBC >= 50 x 10^9/L   vs  day-19 MRD >= 0.1%")
print("=" * 72)
c = coh.copy()
c["poor"] = c["D19_MRD_ge_0p1"].str.upper() == "TRUE"
c["age"] = c["age_at_diagnosis_years_num"].map(f)
c["wbc"] = c["WBC_num"].map(f)
c["mrd_eval"] = c["D19_MRD_ge_0p1"] != ""

def cmh(strata, exp, out):
    """Cochran-Mantel-Haenszel OR with Robins-Breslow-Greenland CI."""
    num = den = 0.0
    s_num = s_den = 0.0
    for s, g in strata:
        a = int(((g[exp]) & (g[out])).sum()); b = int(((g[exp]) & (~g[out])).sum())
        cc = int(((~g[exp]) & (g[out])).sum()); d = int(((~g[exp]) & (~g[out])).sum())
        n = a + b + cc + d
        if n == 0:
            continue
        num += a * d / n; den += b * cc / n
        if n > 1:
            s_num += (a + d) * a * d / n ** 2
            s_den += (a + d) ** 2 / n ** 2 * 0  # placeholder, replaced below
    if den == 0:
        return np.nan, np.nan, np.nan
    orv = num / den
    # RBG variance
    v = 0.0
    P = Q = R_ = S = 0.0
    for s, g in strata:
        a = int(((g[exp]) & (g[out])).sum()); b = int(((g[exp]) & (~g[out])).sum())
        cc = int(((~g[exp]) & (g[out])).sum()); d = int(((~g[exp]) & (~g[out])).sum())
        n = a + b + cc + d
        if n <= 1:
            continue
        P += (a + d) / n * (a * d / n)
        Q += (b + cc) / n * (b * cc / n)
        R_ += (a + d) / n * (b * cc / n)
        S += (b + cc) / n * (a * d / n)
    if P == 0 or Q == 0:
        return orv, np.nan, np.nan
    var = P / (2 * P ** 2) + (R_ + S) / (2 * P * Q) + Q / (2 * Q ** 2)
    lo = np.exp(np.log(orv) - 1.96 * np.sqrt(var))
    hi = np.exp(np.log(orv) + 1.96 * np.sqrt(var))
    return orv, lo, hi

aw = []
for name, exp_series, expected in [
        ("age_ge_5y", c["age"] >= 5.0, "more frequent in poor response"),
        ("WBC_ge_50", c["wbc"] >= 50.0, "less frequent in poor response")]:
    d = c[c["mrd_eval"] & exp_series.notna()].copy()
    d["exp"] = exp_series[d.index].values
    e = d["exp"].values.astype(bool); o = d["poor"].values.astype(bool)
    a = int((e & o).sum()); b = int((e & ~o).sum())
    cc = int((~e & o).sum()); dd = int((~e & ~o).sum())
    orv, p = fisher_exact([[a, b], [cc, dd]])
    import statsmodels.api as sm
    zero_cell = min(a, b, cc, dd) == 0
    if zero_cell:
        # conditional MLE is not identifiable at the boundary; do not print a CI
        lo, hi = np.nan, np.nan
        orv = 0.0
    else:
        tbl = sm.stats.Table2x2(np.array([[a, b], [cc, dd]]))
        lo, hi = tbl.oddsratio_confint()
    # era-aware CMH (strata with no exposure or no outcome contribute zero and are dropped by n>0)
    strata = [(s, g) for s, g in d.groupby("treatment_era")]
    cmh_or, cmh_lo, cmh_hi = cmh(strata, "exp", "poor")
    aw.append(dict(feature=name, n_evaluable=len(d), a=a, b=b, c=cc, d=dd,
                   OR=orv, CI_lo=lo, CI_hi=hi, P=p, zero_cell=zero_cell,
                   CMH_OR_era_adjusted=cmh_or, CMH_CI_lo=cmh_lo, CMH_CI_hi=cmh_hi,
                   expected_direction=expected,
                   observed_direction="more frequent in poor" if (a / (a + b)) > (cc / (cc + dd)) else "less frequent in poor",
                   pct_poor_exposed=a / (a + b) if a + b else np.nan,
                   pct_good_exposed=cc / (cc + dd) if cc + dd else np.nan))
A = pd.DataFrame(aw)
A["BH_q_2features"] = bh(A["P"].values)
A.to_csv("LOCAL_AGE_WBC_TRIANGULATION.tsv", sep="\t", index=False)
print(A[["feature", "n_evaluable", "a", "b", "c", "d", "OR", "CI_lo", "CI_hi", "P",
         "BH_q_2features", "CMH_OR_era_adjusted", "CMH_CI_lo", "CMH_CI_hi",
         "pct_poor_exposed", "pct_good_exposed", "observed_direction"]].to_string(index=False))

# ---------------------------------------------------------------- 4. era sensitivity
print()
print("=" * 72)
print("[4] Era sensitivity")
print("=" * 72)
sens = []
for _, r in R.iterrows():
    mk = r["marker"]
    sens.append(dict(analysis=f"{mk} vs D19 MRD", n=r["n_evaluable"],
                     era_strata=r["era_strata"],
                     era_adjustment="not applicable - single era (all 2020)",
                     note="MRD reporting threshold constant within the flow-evaluable subset; "
                          "no era stratification possible and none needed. Generalisability "
                          "to other eras untested."))
for _, r in A.iterrows():
    sens.append(dict(analysis=f"{r['feature']} vs D19 MRD", n=r["n_evaluable"],
                     era_strata="4 eras (2005/2009 n=32, 2015 n=41, 2020 n=32, 2025 n=10)",
                     era_adjustment="Cochran-Mantel-Haenszel",
                     note=f"unadjusted OR {r['OR']:.2f} ({r['CI_lo']:.2f}-{r['CI_hi']:.2f}); "
                          f"era-adjusted OR {r['CMH_OR_era_adjusted']:.2f} ({r['CMH_CI_lo']:.2f}-{r['CMH_CI_hi']:.2f})"))
S = pd.DataFrame(sens)
S.to_csv("LOCAL_MATURATION_ERA_SENSITIVITY.tsv", sep="\t", index=False)
print(S.to_string(index=False))

# ---------------------------------------------------------------- 5. gate
print()
print("=" * 72)
print("[5] LOCAL_MATURATION_TRIANGULATION gate")
print("=" * 72)
c34 = R[R.marker == "CD34"].iloc[0]
c38 = R[R.marker == "CD38"].iloc[0]
pred_cd34 = c34["median_poor"] > c34["median_good"]
pred_cd38 = c38["median_poor"] < c38["median_good"]
print("CD34 direction as predicted (higher in poor):", bool(pred_cd34),
      " medians %.1f vs %.1f" % (c34["median_poor"], c34["median_good"]),
      " P=%.4f q=%.4f" % (c34["P"], c34["BH_q_2markers"]))
print("CD38 direction as predicted (lower in poor):", bool(pred_cd38),
      " medians %.1f vs %.1f" % (c38["median_poor"], c38["median_good"]),
      " P=%.4f q=%.4f" % (c38["P"], c38["BH_q_2markers"]))
n_cd38 = int(c38["n_evaluable"])
gate = "NOT_SUPPORTED"
if (pred_cd34 and c34["BH_q_2markers"] < 0.05) or (pred_cd38 and c38["BH_q_2markers"] < 0.05 and n_cd38 >= 20):
    gate = "STRONG"
elif pred_cd34 or pred_cd38:
    gate = "SUPPORTIVE"
print("\nLOCAL_MATURATION_TRIANGULATION =", gate)
print("LOCAL_CD34_ANALYSABLE =", bool(c34["n_evaluable"] >= 10))
print("LOCAL_CD38_ANALYSABLE =", bool(n_cd38 >= 10), "(n=%d, missingness %.0f%%)" % (n_cd38, 100 * (1 - n_cd38 / 32)))

json.dump(dict(LOCAL_CD34_ANALYSABLE=bool(c34["n_evaluable"] >= 10),
               LOCAL_CD38_ANALYSABLE=bool(n_cd38 >= 10),
               LOCAL_MATURATION_TRIANGULATION=gate,
               n_flow_evaluable=len(sub), n_cd34=int(c34["n_evaluable"]), n_cd38=n_cd38),
          open("_local_maturation_status.json", "w"), indent=1)
print("\nwritten: LOCAL_FLOW_FIELD_AUDIT.tsv, LOCAL_CD34_CD38_RESULTS.tsv,")
print("         LOCAL_AGE_WBC_TRIANGULATION.tsv, LOCAL_MATURATION_ERA_SENSITIVITY.tsv")
