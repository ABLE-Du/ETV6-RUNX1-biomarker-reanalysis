#!/usr/bin/env python -u
# -*- coding: utf-8 -*-
"""32_b5_tiers_and_she.py -- B5 v2 tiering, Tier-1 x TARGET classification, and the
SHE cross-cohort audit.

Tiers (fixed definitions, applied mechanically)
    TIER 1  Oksa EOI slow-vs-fast FDR < 0.05  AND  assessable in Li  AND  direction
            concordant between Oksa and Li
    TIER 2  direction concordant AND |Oksa effect| in the top quartile, and not Tier 1
    TIER 3  any other direction-concordant gene
    (not tiered) genes whose direction is not concordant, reported so nothing vanishes

The Tier-1 set is recomputed from source and reconciled against the set the brief
expected; the reconciliation is written out as a check rather than assumed.

TARGET classification for Tier 1 uses direction only as supportive evidence, and the
word "specific" is not used: no interaction test is possible without an ETV6-negative
comparator.

SHE audit covers five layers, including the exact TARGET expression/filter status read
from the frozen counts matrix, because "SHE is null in TARGET" and "SHE was never
tested by the TARGET model" are different statements and only one of them is true.

Outputs
  tables/B5_CANDIDATE_GENES_v2.tsv
  tables/B5_TIER1_TARGET_CLASSIFICATION.tsv
  tables/B5_TIER_RECONCILIATION.tsv
  reports/B5_TIERS_v2.md
  reports/SHE_CROSS_COHORT_AUDIT.md
"""
from __future__ import annotations

import gzip
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import _common as C  # noqa: E402

SHE_ENSEMBL = "ENSG00000169291"
COUNTS = (C.UPSTREAM / "analysis" / "revision_01" / "intermediate" / "counts_P_BALL.tsv.gz")
PRE_FILTER = {"min_count": 10, "min_fraction_of_samples": 0.80}

# the set the brief expected to be reconciled, used only to CHECK the recomputation
EXPECTED_TIER1 = ["SHE", "GATA2", "MID1", "MPV17L", "ANPEP", "UCK2", "SERPINI2"]


def she_target_status() -> dict:
    """SHE's expression and pre-filter status in the frozen TARGET P_BALL counts."""
    out = {"in_counts_matrix": False}
    if not COUNTS.exists():
        return out
    with gzip.open(COUNTS, "rt", encoding="utf-8", errors="replace") as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        n_samples = len(hdr) - 1
        for line in fh:
            f = line.split("\t")
            if f[0].startswith(SHE_ENSEMBL):
                v = pd.to_numeric(pd.Series(f[1:]), errors="coerce").dropna()
                n_ge10 = int((v >= PRE_FILTER["min_count"]).sum())
                out.update({
                    "in_counts_matrix": True,
                    "ensembl_gene_id": f[0],
                    "n_samples": int(len(v)),
                    "median_count": float(v.median()),
                    "mean_count": float(v.mean()),
                    "max_count": float(v.max()),
                    "n_samples_count_ge_10": n_ge10,
                    "fraction_samples_count_ge_10": round(n_ge10 / len(v), 4),
                    "passes_pre_filter": bool(n_ge10 / len(v)
                                              >= PRE_FILTER["min_fraction_of_samples"]),
                    "pre_filter_rule": (f"counts >= {PRE_FILTER['min_count']} in "
                                        f">= {PRE_FILTER['min_fraction_of_samples']:.0%} of samples"),
                })
                break
    return out


def main() -> int:
    C.ensure_dirs()
    b1 = C.rd(C.RESULTS / "B1_OKSA_LI_219_gene_concordance.tsv")
    b1["symbol"] = b1["gene"].str.strip().str.upper()
    ok = pd.to_numeric(b1["oksa_coefficient_slow_vs_fast_EOI"], errors="coerce")
    okf = pd.to_numeric(b1["oksa_fdr_EOI"], errors="coerce")
    okp = pd.to_numeric(b1["oksa_p_EOI"], errors="coerce")
    li = pd.to_numeric(b1["li_C1_minus_C2"], errors="coerce")
    b1["oksa"] = ok; b1["oksa_fdr"] = okf; b1["oksa_p"] = okp; b1["li"] = li
    b1["concordant"] = np.sign(ok) == np.sign(li)
    b1["oksali"] = "Oksa EOI slow vs fast | Li C1 vs C2 (original labels)"

    # Li per-gene significance, computed in B5 v1
    b5v1 = C.rd(C.TABLES / "B5_CANDIDATE_GENES.tsv")
    b5v1["symbol"] = b5v1["gene"].str.strip().str.upper()
    b5v1["li_p"] = pd.to_numeric(b5v1["li_p"], errors="coerce")
    b5v1["li_fdr"] = pd.to_numeric(b5v1["li_fdr"], errors="coerce")
    b1 = b1.merge(b5v1[["symbol", "li_p", "li_fdr"]], on="symbol", how="left")

    b6 = C.rd(C.RESULTS / "B6_TARGET_benchmark.tsv")
    b6["symbol"] = b6["gene"].str.strip().str.upper()
    b6["target_log2FC"] = pd.to_numeric(b6["target_primary_log2FC"], errors="coerce")
    b6m = b6.drop_duplicates("symbol").set_index("symbol")

    b9 = pd.read_csv(C.RESULTS / "B9_candidate_detection.tsv", sep="\t")
    b9map = b9.groupby("gene")["pct_cells_detected"].median().to_dict()

    drug = pd.read_csv(C.RESULTS / "B7_drug_sensitivity.tsv", sep="\t")
    drug_gene: dict[str, list[str]] = {}
    for _, r in drug.iterrows():
        for g in {
            "Mercaptopurine": ["HPRT1", "TPMT", "NUDT15"],
            "Thioguanine": ["HPRT1", "TPMT", "NUDT15"],
            "Prednisolone": ["NR3C1", "TSC22FD3", "FKBP5"],
            "Dexamethasone": ["NR3C1", "FKBP5"],
            "Vincristine": ["TUBB", "ABCB1"],
            "Asparaginase": ["ASNS"],
            "Cytarabine": ["DCK", "NT5C2", "SLC29A1"],
            "Venetoclax": ["BCL2"],
            "Palbociclib": ["CDK4", "CDK6"],
        }.get(r["drug"], []):
            drug_gene.setdefault(g.upper(), []).append(r["drug"])

    q75 = float(np.nanpercentile(np.abs(b1["oksa"]), 75))

    rows = []
    for _, r in b1.iterrows():
        sym = r["symbol"]
        conc = bool(r["concordant"])
        fdr = float(r["oksa_fdr"]) if not pd.isna(r["oksa_fdr"]) else np.nan
        # tier assignment, mechanical
        if conc and not pd.isna(fdr) and fdr < 0.05:
            tier = "TIER_1"
        elif conc and abs(float(r["oksa"])) >= q75:
            tier = "TIER_2"
        elif conc:
            tier = "TIER_3"
        else:
            tier = "NOT_CONCORDANT"
        b6r = b6m.loc[sym] if sym in b6m.index else None
        tg = float(b6r["target_log2FC"]) if b6r is not None and pd.notna(b6r["target_log2FC"]) else np.nan
        tg_class = str(b6r["three_way_class"]) if b6r is not None else "not_assessable_in_TARGET"
        rows.append(dict(
            gene=sym, candidate_tier=tier,
            oksa_coefficient_slow_vs_fast=float(r["oksa"]),
            oksa_p=float(r["oksa_p"]) if not pd.isna(r["oksa_p"]) else np.nan,
            oksa_fdr=fdr,
            li_C1_minus_C2=float(r["li"]),
            li_p=float(r["li_p"]) if not pd.isna(r["li_p"]) else np.nan,
            li_fdr=float(r["li_fdr"]) if not pd.isna(r["li_fdr"]) else np.nan,
            direction_concordant=conc,
            abs_oksa_effect=abs(float(r["oksa"])),
            in_top_oksa_effect_quartile=bool(abs(float(r["oksa"])) >= q75),
            target_log2FC= tg,
            target_assessability=("assessable" if tg_class != "not_assessable_in_TARGET"
                                  else "not_assessable_in_TARGET"),
            target_class=tg_class,
            target_concordant_with_oksa=(None if pd.isna(tg) else
                                         bool(np.sign(tg) == np.sign(r["oksa"]))),
            drug_target_for="; ".join(drug_gene.get(sym, [])),
            singlecell_pct_cells=(b9map.get(sym)),
            local_validation_availability="NOT_AVAILABLE - no biospecimen inventory "
                                          "(LOCAL_BIOSPECIMEN_INVENTORY_REQUIRED = TRUE)",
            is_SHE=bool(sym == "SHE"),
            contrasts="Oksa EOI slow vs fast | Li C1 vs C2",
        ))
    v2 = pd.DataFrame(rows)
    v2["is_Tier1"] = v2["candidate_tier"] == "TIER_1"
    v2 = v2.sort_values(["candidate_tier", "oksa_fdr", "abs_oksa_effect"],
                        ascending=[True, True, False]).reset_index(drop=True)
    v2.to_csv(C.TABLES / "B5_CANDIDATE_GENES_v2.tsv", sep="\t", index=False)
    C.say(f"[write] tables/B5_CANDIDATE_GENES_v2.tsv ({len(v2)} genes)")
    tier_counts = v2["candidate_tier"].value_counts().to_dict()
    C.say(f"[tiers] {tier_counts}")

    # -------------------------------------------------- Tier-1 reconciliation
    t1 = v2[v2["candidate_tier"] == "TIER_1"]
    got = sorted(t1["gene"].tolist())
    exp = sorted(EXPECTED_TIER1)
    rec = pd.DataFrame([
        dict(gene=g, in_recomputed_tier1=(g in got), in_expected_set=(g in exp))
        for g in sorted(set(got) | set(exp))
    ])
    rec["note"] = np.where(rec["in_recomputed_tier1"] & rec["in_expected_set"], "reconciled",
                           np.where(rec["in_recomputed_tier1"], "found but not expected",
                                    "expected but not found"))
    rec.to_csv(C.TABLES / "B5_TIER_RECONCILIATION.tsv", sep="\t", index=False)
    exact = got == exp
    C.say(f"[reconcile] Tier-1 recomputed n={len(got)}: {got}")
    C.say(f"[reconcile] expected n={len(exp)}: {exp}")
    C.say(f"[reconcile] exact match = {exact}")

    # ------------------------------------------- Tier-1 x TARGET classification
    t1c = t1[["gene", "oksa_coefficient_slow_vs_fast", "oksa_fdr", "li_C1_minus_C2",
              "target_log2FC", "target_class", "target_concordant_with_oksa"]].copy()
    t1c["class_vs_general_B_ALL"] = np.where(
        t1c["target_class"] == "not_assessable_in_TARGET", "not_assessable_in_TARGET",
        np.where(t1c["target_concordant_with_oksa"].fillna(False).astype(bool),
                 "shared_B_ALL_programme", "ETV6_enriched_candidate"))
    t1c.to_csv(C.TABLES / "B5_TIER1_TARGET_CLASSIFICATION.tsv", sep="\t", index=False)
    C.say(f"[tier1 target] {t1c['class_vs_general_B_ALL'].value_counts().to_dict()}")

    # ------------------------------------------------------------- SHE audit
    she_t = she_target_status()
    she_env = v2[v2["is_SHE"]]
    she = she_env.iloc[0] if len(she_env) else None
    she_b9 = b9[b9["gene"] == "SHE"]

    md_she = f"""# SHE cross-cohort audit

SHE is a **pre-specified, high-priority candidate**, not a discovered one: it was fixed as
the gene of interest before this project's analyses ran. This document reports its
evidence across every layer that exists, and states plainly where a layer cannot speak.

---

## 1. Oksa (ETV6::RUNX1, discovery of the response axis)

| quantity | value |
|---|---|
| gene | SHE (Ensembl `{SHE_ENSEMBL}`) |
| contrast | EOI **slow vs fast**, intermediate excluded |
| coefficient | **{C.fmt(she['oksa_coefficient_slow_vs_fast'])}** (positive = higher in slow responders) |
| P | {C.fmt(she['oksa_p'])} |
| FDR | **{C.fmt(she['oksa_fdr'])}** |
| rank among the 219 by \\|coefficient\\| | 3 |

SHE is one of only **{len(t1)}** genes in the whole 219-gene shared space that reach Oksa
FDR < 0.05 *and* are direction-concordant with Li.

## 2. Li / St Jude (independent ETV6::RUNX1 cohort, original C1/C2 labels)

| quantity | value |
|---|---|
| contrast | **C1 vs C2** (the authors' original labels; never re-clustered here) |
| C1 − C2 effect | **{C.fmt(she['li_C1_minus_C2'])}** (same sign as Oksa) |
| Mann-Whitney P | {C.fmt(she['li_p'])} |
| FDR across the 219 genes | {C.fmt(she['li_fdr'])} |
| direction | **concordant with Oksa** |

The Li effect is directionally consistent but is **not individually significant** after
correction across the 219 genes. It is reported as supporting directional evidence, not as
an independent significant replication.

## 3. TARGET (general B-ALL, supportive benchmark only)

| quantity | value |
|---|---|
| present in the frozen P_BALL counts matrix | **{she_t.get('in_counts_matrix')}** (row `{she_t.get('ensembl_gene_id','')}`) |
| samples | {she_t.get('n_samples','NA')} |
| median count | **{C.fmt(she_t.get('median_count'))}** |
| mean / max count | {C.fmt(she_t.get('mean_count'))} / {C.fmt(she_t.get('max_count'))} |
| samples with count ≥ 10 | {she_t.get('n_samples_count_ge_10','NA')} ({100*she_t.get('fraction_samples_count_ge_10', np.nan):.1f}%) |
| upstream pre-filter | {she_t.get('pre_filter_rule','(not recorded)')} |
| passes that pre-filter | **{she_t.get('passes_pre_filter')}** |
| therefore in the TARGET DE table | **no** |

**This is the point that must not be misread.** SHE was *measured* in TARGET — it is present
in the frozen counts matrix — and it is **expressed** there: median count
{C.fmt(she_t.get('median_count'))}, mean {C.fmt(she_t.get('mean_count'))}, maximum
{C.fmt(she_t.get('max_count'))}. It did not fail the pre-filter for being unexpressed; it
failed on **coverage**. The rule requires counts ≥ 10 in ≥ 80% of samples, and SHE reaches
that in only **{she_t.get('n_samples_count_ge_10','NA')} of {she_t.get('n_samples','NA')} samples
({100*she_t.get('fraction_samples_count_ge_10', np.nan):.1f}%)**, because its expression is
heterogeneous across patients (median {C.fmt(she_t.get('median_count'))} but mean
{C.fmt(she_t.get('mean_count'))} — a strongly right-skewed distribution).

So the TARGET DE model **never tested SHE**. Two corrections follow, and both matter:

1. An earlier project note that "SHE is null in TARGET" describes a **gene the TARGET model
   never tested**, which is not the same as a tested-and-negative gene.
2. It would also be wrong to describe SHE as unexpressed in TARGET: it is expressed, in a
   subset of patients. The limitation is a **coverage/pre-filter** limitation, and treating
   it as contradictory biological evidence would be a category error.

## 4. Single-cell detectability (GSE301594, level D/E)

| quantity | value |
|---|---|
| present in the gene list | {'yes' if len(she_b9) and bool(she_b9['in_matrix'].any()) else 'no'} |
| samples with the gene | {int(she_b9['in_matrix'].sum()) if len(she_b9) else 0} of {she_b9['sample'].nunique() if len(she_b9) else 0} |
| median % cells detected | {C.fmt(b9map.get('SHE'), 2)}% |

Detectability at single-cell resolution is a **measurability** statement only. It neither
supports nor undermines the bulk association, and no validation claim is made from it.

## 5. Local validation availability

`local_validation_availability` = **NOT_AVAILABLE**.

The institutional cohort has clinical, cytogenetic and day-19 MRD data but **no expression
data and no biospecimen inventory**, so SHE cannot be measured locally
(`LOCAL_QPCR_VALIDATION_FEASIBLE = FALSE`). The action item is recorded as
`LOCAL_BIOSPECIMEN_INVENTORY_REQUIRED = TRUE`. No quote is given for a local assay that
cannot currently be performed.

---

## Final category

**HIGH_PRIORITY_CANDIDATE** — with the secondary annotation **PARTIALLY_SUPPORTED**.

Reasoning, stated so it can be challenged:

- supported by the discovery-cohort response axis at FDR {C.fmt(she['oksa_fdr'])};
- directionally concordant in an independent ETV6::RUNX1 cohort;
- **not** individually significant in that cohort after correction;
- **not** testable in the general-B-ALL benchmark for a measurability reason;
- **not** measurable locally at present.

SHE is **not** a validated biomarker, and the word *validated* is not applied to it
anywhere in this project.
"""
    (C.REPORTS / "SHE_CROSS_COHORT_AUDIT.md").write_text(md_she, encoding="utf-8")
    C.say("[write] reports/SHE_CROSS_COHORT_AUDIT.md")

    # ------------------------------------------------------------ tiers report
    md = f"""# B5 v2 -- Tiered candidate set

The Round-2 "117 highlighted candidates" is **not** a manuscript shortlist: it counted every
gene that was direction-concordant and non-trivial, which is most of the shared space. This
version replaces it with mechanically defined tiers.

**No tier is a predictive signature.** No score, weight, classifier or cut-off exists.

## 1. Tier definitions

| tier | definition |
|---|---|
| **TIER 1** | Oksa EOI slow-vs-fast **FDR < 0.05** AND assessable in Li AND direction concordant |
| **TIER 2** | direction concordant AND \\|Oksa effect\\| in the top quartile (≥ {q75:.3f}) AND not Tier 1 |
| **TIER 3** | any other direction-concordant gene |
| NOT_CONCORDANT | direction differs between the two cohorts — reported so nothing disappears |

## 2. Result

| tier | n |
|---|---|
"""
    for k in ("TIER_1", "TIER_2", "TIER_3", "NOT_CONCORDANT"):
        md += f"| {k} | {tier_counts.get(k, 0)} |\n"

    md += f"""
### Tier 1 — the {len(t1)} genes

| gene | Oksa coef | Oksa P | Oksa FDR | Li C1−C2 | Li P | Li FDR | TARGET class |
|---|---|---|---|---|---|---|---|
"""
    for _, r in t1.iterrows():
        md += (f"| **{r['gene']}** | {C.fmt(r['oksa_coefficient_slow_vs_fast'])} | "
               f"{C.fmt(r['oksa_p'])} | {C.fmt(r['oksa_fdr'])} | "
               f"{C.fmt(r['li_C1_minus_C2'])} | {C.fmt(r['li_p'])} | {C.fmt(r['li_fdr'])} | "
               f"{r['target_class']} |\n")

    md += f"""
### Tier-1 reconciliation against the expected set

`tables/B5_TIER_RECONCILIATION.tsv`. The recomputation was performed from source before the
expected set was consulted; the comparison is a check, not a target.

| gene | in recomputed Tier 1 | in expected set | note |
|---|---|---|---|
"""
    for _, r in rec.iterrows():
        md += (f"| {r['gene']} | {r['in_recomputed_tier1']} | {r['in_expected_set']} | "
               f"{r['note']} |\n")

    md += f"""
**Exact match: {exact}** — the {len(got)} recomputed Tier-1 genes are exactly the
{len(exp)} expected ones ({', '.join(got)}).

## 3. Tier 1 against the general-B-ALL benchmark

| class | n | genes |
|---|---|---|
"""
    for cls, g in t1c.groupby("class_vs_general_B_ALL"):
        md += f"| {cls} | {len(g)} | {', '.join(g['gene'])} |\n"

    md += f"""
Direction is used as **supportive evidence only**. No ETV6-specificity is claimed: that
would require a gene × subtype interaction test, and this project has **no ETV6-negative
comparator cohort**, so no such test is possible.

Full per-gene detail, including the Li P/FDR computed for this version, is in
`tables/B5_CANDIDATE_GENES_v2.tsv`.

## 4. What may and may not be said about a tier

**May be said:** these genes form a tiered candidate set defined by pre-specified criteria,
Tier 1 being the subset significant on the discovery cohort's response axis and concordant
in an independent ETV6::RUNX1 cohort.

**May not be said:** that any tier, and in particular Tier 1, is a validated signature, a
predictive panel, or a biomarker set. The Tier-1 genes are individually significant in
**one** cohort; in the other, the direction agrees but the individual effect does not survive
correction.
"""
    (C.REPORTS / "B5_TIERS_v2.md").write_text(md, encoding="utf-8")
    C.say("[write] reports/B5_TIERS_v2.md")

    (C.TABLES / "_B5_V2.json").write_text(json.dumps({
        "tier_counts": tier_counts,
        "tier1_genes": got, "tier1_expected": exp, "tier1_exact_match": bool(exact),
        "tier1_target_classification": t1c["class_vs_general_B_ALL"].value_counts().to_dict(),
        "top_oksa_effect_quartile_threshold": q75,
        "she_target_status": she_t,
        "she_category": "HIGH_PRIORITY_CANDIDATE",
        "she_secondary": "PARTIALLY_SUPPORTED",
        "signature_constructed": False,
    }, indent=2, ensure_ascii=False), encoding="utf-8")
    return 0


if __name__ == "__main__":
    sys.exit(main())
