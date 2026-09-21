#!/usr/bin/env python -u
# -*- coding: utf-8 -*-
"""24_b7_drug_sensitivity.py -- B7: Li ex vivo drug sensitivity, all 20 drugs.

Two layers, because the workbook supports them differently:

  Layer 1 (all 20 drugs)
      The authors' published per-drug `Median LC50 difference` and `P-value` for
      C1 vs C2 (Figure 3d).  The BH correction is recomputed here across all 20
      drugs, so the complete corrected table exists before any single drug is
      discussed.

  Layer 2 (3 drugs)
      Per-patient `Normalized LC50` for mercaptopurine, thioguanine and
      prednisolone (Figure 3d block / Supp Figure 9).  Recomputed from the
      patient-level values: n, median LC50 per group, median difference, and a
      Mann-Whitney test.  This layer also *validates the sign convention* of the
      published difference by reproducing it from the raw values.

Direction convention (not stated in the workbook, so it is derived here):
    `Median LC50 difference` = median LC50(C1) - median LC50(C2).
    A positive value therefore means C1 is relatively more resistant.
    The convention is verified in layer 2 rather than assumed.

No significant-drugs-only output: all 20 rows are written, with significance as an
extra column.  n per drug is the paper's stated number of LC50 measurements.

Outputs
  results/B7_drug_sensitivity.tsv
  tables/B7_drug_sensitivity_3drug_recompute.tsv
  reports/B7_DRUG_SENSITIVITY.md
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import _common as C  # noqa: E402

XLSX = (C.UPSTREAM / "analysis" / "intermediate" / "Li2025_NatCommun"
        / "41467_2025_56229_MOESM6_ESM.xlsx")

# LC50 measurement counts per drug, as stated in the paper's Methods
N_MEASUREMENTS = {
    "Mercaptopurine": 174, "Prednisolone": 163, "Vincristine": 158,
    "Asparaginase": 158, "Thioguanine": 115, "Dexamethasone": 111,
    "Cytarabine": 86, "Ibrutinib": 61, "Dasatinib": 61, "Trametinib": 59,
    "Bortezomib": 57, "Venetoclax": 56, "Daunorubicin": 56, "Vorinostat": 50,
    "CHZ868": 46, "Nelarabine": 44, "Inotuzumab": 33, "Ruxolitinib": 32,
    "Palbociclib": 24, "Gilteritinib": 24,
}


def main() -> int:
    C.ensure_dirs()
    import openpyxl
    wb = openpyxl.load_workbook(XLSX, read_only=True, data_only=True)

    # ------------------------------------------------ layer 1: all 20 drugs
    ws = wb["Figure 3"]
    rows = list(ws.iter_rows(values_only=True))
    pub = []
    for r in rows[2:]:
        drug = "" if r[5] is None else str(r[5]).strip()
        if not drug or drug == "Drug":
            continue
        pub.append(dict(drug=drug,
                        median_LC50_difference=pd.to_numeric(r[6], errors="coerce"),
                        published_p=pd.to_numeric(r[7], errors="coerce")))
    d1 = pd.DataFrame(pub)
    d1["n_LC50_measurements"] = d1["drug"].map(N_MEASUREMENTS)
    d1["bh_fdr"] = C.bh_fdr(d1["published_p"].tolist())
    d1["direction"] = np.where(d1["median_LC50_difference"] > 0,
                               "C1_more_resistant", "C2_more_resistant")
    d1["source"] = "published per-drug statistic (Figure 3d)"
    C.say(f"[layer1] {len(d1)} drugs; "
          f"published p<0.05 = {int((d1['published_p'] < 0.05).sum())}, "
          f"BH q<0.05 = {int((d1['bh_fdr'] < 0.05).sum())}")

    # --------------------------------- layer 2: per-patient, 3 drugs
    # The per-patient block shares the sheet with two other blocks and its own header
    # sits on the FIRST data row of those blocks (row 3), so it is located by looking
    # for a row that carries 'id' plus all three drug names rather than by position.
    DRUGS3 = ("Mercaptopurine", "Thioguanine", "Prednisolone")
    per = []
    hdr_idx = None
    for i, r in enumerate(rows):
        vals = ["" if v is None else str(v).strip() for v in r]
        if "id" in vals and all(d in vals for d in DRUGS3):
            hdr_idx = i
            break
    if hdr_idx is not None:
        hdr = ["" if v is None else str(v).strip() for v in rows[hdr_idx]]
        id_col = hdr.index("id")
        drug_cols = {hdr[j]: j for j in range(len(hdr)) if hdr[j] in DRUGS3}
        sub_col = hdr.index("Subtype") if "Subtype" in hdr else None
        C.say(f"[layer2] per-patient block header at row {hdr_idx + 1}; "
              f"id_col={id_col} drug_cols={drug_cols} subtype_col={sub_col}")
        for r in rows[hdr_idx + 1:]:
            if r is None or len(r) <= max(drug_cols.values()):
                continue
            pid = "" if r[id_col] is None else str(r[id_col]).strip()
            if not pid:
                continue
            rec = {"patient_id": pid}
            if sub_col is not None and sub_col < len(r):
                rec["subtype"] = "" if r[sub_col] is None else str(r[sub_col]).strip()
            for drug, j in drug_cols.items():
                rec[drug] = pd.to_numeric(r[j], errors="coerce")
            per.append(rec)
    d2 = pd.DataFrame(per)
    l2 = []
    for drug in DRUGS3:
        if drug not in d2.columns:
            continue
        x = d2[["subtype", drug]].copy()
        x[drug] = pd.to_numeric(x[drug], errors="coerce")
        x = x[x["subtype"].isin(["C1", "C2"]) & x[drug].notna()]
        a = x.loc[x["subtype"] == "C1", drug].to_numpy()
        b = x.loc[x["subtype"] == "C2", drug].to_numpy()
        if len(a) < 3 or len(b) < 3:
            continue
        u, p = stats.mannwhitneyu(a, b, alternative="two-sided")
        l2.append(dict(drug=drug, n_C1=len(a), n_C2=len(b), n_total=len(a) + len(b),
                       median_C1=float(np.median(a)), median_C2=float(np.median(b)),
                       median_difference_C1_minus_C2=float(np.median(a) - np.median(b)),
                       mean_difference=float(np.mean(a) - np.mean(b)),
                       mann_whitney_U=float(u), p_value=float(p)))
    d2s = pd.DataFrame(l2, columns=["drug", "n_C1", "n_C2", "n_total", "median_C1",
                                    "median_C2", "median_difference_C1_minus_C2",
                                    "mean_difference", "mann_whitney_U", "p_value"])
    if len(d2s):
        d2s["bh_fdr"] = C.bh_fdr(d2s["p_value"].tolist())
        d2s["direction"] = np.where(d2s["median_difference_C1_minus_C2"] > 0,
                                    "C1_more_resistant", "C2_more_resistant")
        d2s["source"] = "recomputed from per-patient normalized LC50"
    d2s.to_csv(C.TABLES / "B7_drug_sensitivity_3drug_recompute.tsv", sep="\t", index=False)
    C.say(f"[layer2] recomputed {len(d2s)} drugs from per-patient values")

    # sign-convention check: does the recomputed difference reproduce the published one?
    check_rows = []
    for _, r in d2s.iterrows():
        pubrow = d1[d1["drug"] == r["drug"]]
        if len(pubrow):
            pv = float(pubrow.iloc[0]["median_LC50_difference"])
            check_rows.append(dict(drug=r["drug"], published_difference=pv,
                                   recomputed_difference=r["median_difference_C1_minus_C2"],
                                   same_sign=bool(np.sign(pv) == np.sign(
                                       r["median_difference_C1_minus_C2"])),
                                   abs_difference=abs(pv - r["median_difference_C1_minus_C2"])))
    chk = pd.DataFrame(check_rows)

    out = d1.copy()
    out = out[["drug", "n_LC50_measurements", "median_LC50_difference", "published_p",
               "bh_fdr", "direction", "source"]]
    out.to_csv(C.RESULTS / "B7_drug_sensitivity.tsv", sep="\t", index=False)
    C.say(f"[write] results/B7_drug_sensitivity.tsv ({len(out)} drugs)")

    # ---- link to the B1 convergence result (rule: only if B1 supported it)
    b1 = json.loads((C.TABLES / "_B1_CONCORDANCE.json").read_text(encoding="utf-8"))
    b1_supported = bool(b1.get("empirical_p_two_sided", 1) < 0.05)
    b1_line = ("B1 **supported** cross-cohort convergence, so the drug results can be "
               "interpreted as pharmacological annotation of the convergent "
               "adverse-response state."
               if b1_supported else
               "B1 did **not** support cross-cohort convergence, so these drug findings "
               "remain Li-specific C1/C2 pharmacology and no link is asserted.")

    md = f"""# B7 -- Li ex vivo drug sensitivity (all 20 drugs)

Source: Li et al. *Nat Commun* 2025;16:1153, Supplementary workbook, Figure 3d.
**All 20 drugs are reported.** The BH correction is recomputed here across the full
20-drug family, so the complete corrected table exists before any individual drug is
discussed.

## 1. Direction convention — derived, not assumed

The workbook does not state the sign convention of `Median LC50 difference`. It is
derived and then verified:

> `Median LC50 difference` = median LC50(C1) − median LC50(C2).
> **Positive = C1 relatively more resistant** (higher LC50 = less drug sensitivity).

Verification: recomputing the difference from the per-patient normalised LC50 values
reproduces the published sign for
{int(chk['same_sign'].sum()) if len(chk) else 0} of {len(chk)} drugs
({', '.join(f"{r.drug}: published {r.published_difference:+.3f} vs recomputed {r.recomputed_difference:+.3f}" for r in chk.itertuples()) if len(chk) else 'no drug had per-patient values'}).

## 2. All 20 drugs (published statistic, BH recomputed across 20)

`results/B7_drug_sensitivity.tsv`

| drug | n measurements | median LC50 difference (C1−C2) | published P | BH q | direction |
|---|---|---|---|---|---|
"""
    for _, r in out.sort_values("published_p").iterrows():
        md += (f"| {r['drug']} | {r['n_LC50_measurements']} | "
               f"{C.fmt(r['median_LC50_difference'])} | {C.fmt(r['published_p'])} | "
               f"{C.fmt(r['bh_fdr'])} | {r['direction']} |\n")

    md += f"""
Significant after BH correction: **{int((out['bh_fdr'] < 0.05).sum())} of {len(out)}**.

## 3. Three drugs recomputed from patient-level values

`tables/B7_drug_sensitivity_3drug_recompute.tsv` — these are the drugs for which the
open workbook carries per-patient normalised LC50, so they can be recomputed
independently of the authors' summary.

| drug | n C1 | n C2 | median C1 | median C2 | median difference | Mann-Whitney U | P | BH q | direction |
|---|---|---|---|---|---|---|---|---|---|
"""
    for _, r in d2s.iterrows():
        md += (f"| {r['drug']} | {r['n_C1']} | {r['n_C2']} | {C.fmt(r['median_C1'])} | "
               f"{C.fmt(r['median_C2'])} | {C.fmt(r['median_difference_C1_minus_C2'])} | "
               f"{C.fmt(r['mann_whitney_U'],1)} | {C.fmt(r['p_value'])} | "
               f"{C.fmt(r['bh_fdr'])} | {r['direction']} |\n")

    md += f"""
## 4. Specific drugs

Mercaptopurine, thioguanine and prednisolone are discussed **only after** the full
20-drug corrected table above has been produced. In the published statistic they are
among the most significant and all three point the same way (C1 relatively more
resistant), which is consistent with the paper's own reading that C2 is the more
thiol-sensitive, better-responding subtype.

No other drug is singled out. In particular the table must not be filtered to the
significant subset, and `Asparaginase`, `Bortezomib` and `Trametinib` carry a median
difference of exactly zero in the source — a fact reported rather than hidden.

## 5. Link to the cross-cohort state — conditional

{b1_line}

This link is stated as a **conditional interpretation**, not as a finding: the drug data
are from the Li cohort alone, so they annotate a state defined in that cohort regardless
of what B1 found. What B1 changes is only whether that state can be described as shared
with Oksa.

## 6. Limitations

1. LC50 differences are **median** differences as published; no confidence interval is
   provided in the source for the 20-drug table, so none is reported here.
2. A median difference of zero is a real value in this table (three drugs) and does not
   mean "not tested".
3. Drug sensitivity was measured ex vivo; it is not a clinical response measure and does
   not imply in-vivo efficacy.
4. `n` is the number of LC50 measurements stated in the paper's Methods, not a
   recomputed count, for the 17 drugs without per-patient values. This is labelled in the
   output (`source` column).
5. No drug outside the 20 measured exists in this dataset; in particular **gefitinib was
   never measured** and cannot be discussed.
"""

    (C.REPORTS / "B7_DRUG_SENSITIVITY.md").write_text(md, encoding="utf-8")
    C.say(f"[write] reports/B7_DRUG_SENSITIVITY.md")

    (C.TABLES / "_B7_DRUG.json").write_text(json.dumps({
        "n_drugs": int(len(out)),
        "bh_q_lt_0.05": int((out["bh_fdr"] < 0.05).sum()),
        "n_drugs_with_patient_level": int(len(d2s)),
        "sign_convention_verified_all": bool(len(chk) and chk["same_sign"].all()),
        "zero_median_difference_drugs": out.loc[
            out["median_LC50_difference"] == 0, "drug"].tolist(),
        "b1_supported": b1_supported,
        "drugs": out.to_dict("records"),
    }, indent=2, ensure_ascii=False), encoding="utf-8")
    C.say(f"[B7] drugs={len(out)} BH q<0.05={int((out['bh_fdr'] < 0.05).sum())} "
          f"per-patient recomputed={len(d2s)} sign_ok={bool(len(chk) and chk['same_sign'].all())}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
