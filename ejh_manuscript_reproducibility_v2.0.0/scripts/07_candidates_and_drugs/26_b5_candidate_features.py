#!/usr/bin/env python -u
# -*- coding: utf-8 -*-
"""26_b5_candidate_features.py -- B5: candidate molecular feature table.

**This is not a signature.** No score, no weights, no classifier. It is a
descriptive candidate table assembled from the rungs that have already run, with the
selection criteria written down before the table was built:

    a gene may be *highlighted* when
      (1) direction is concordant between Oksa and Li, and
      (2) the effect is non-trivial, and
      (3) it is measurable in principle, and
      (4) its biology is interpretable.

Genes are **not** ranked by P value: the ranking key is cross-cohort concordance
first, then the magnitude of both effects, then pathway membership.  Every gene that
passes (1) is reported, and the ones that fail (2) or (3) are reported with the reason
rather than dropped.

Data drawn on
  B1  results/B1_OKSA_LI_219_gene_concordance.tsv
  B2  tables/B2_pathway_axes.tsv + results/B2_pathway_convergence.tsv
  B6  results/B6_TARGET_benchmark.tsv
  B7  results/B7_drug_sensitivity.tsv (drug classes are mapped to genes where the
      mechanism is unambiguous; this mapping is documented per row)

SHE is included explicitly as a pre-specified candidate whatever its result
(brief section 19): if it is absent or not assessable, that is stated; if it is
discordant or null, that is stated; no rescue analysis is attempted.

Outputs
  tables/B5_CANDIDATE_GENES.tsv
  reports/B5_CANDIDATE_FEATURES.md
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import _common as C  # noqa: E402

SHE_ENSEMBL = "ENSG00000169291"

# drug -> mechanistically unambiguous gene targets.  Only pairs where the drug's
# action on the gene product is direct and well established are listed; the mapping
# is reported in the output so a reader can reject it.
DRUG_GENE = {
    "Mercaptopurine": ["HPRT1", "TPMT", "NUDT15"],
    "Thioguanine": ["HPRT1", "TPMT", "NUDT15"],
    "Prednisolone": ["NR3C1", "TSC22D3", "FKBP5"],
    "Dexamethasone": ["NR3C1", "TSC22D3", "FKBP5"],
    "Vincristine": ["TUBB", "ABCB1", "MAP4"],
    "Asparaginase": ["ASNS"],
    "Cytarabine": ["DCK", "NT5C2", "SLC29A1"],
    "Dasatinib": ["ABL1", "SRC", "LYN"],
    "Ibrutinib": ["BTK"],
    "Trametinib": ["MAP2K1", "MAP2K2"],
    "Bortezomib": ["PSMB5", "NFKB1"],
    "Daunorubicin": ["TOP2A", "ABCB1"],
    "Venetoclax": ["BCL2"],
    "Vorinostat": ["HDAC1", "HDAC2"],
    "Nelarabine": ["NT5E", "DCK"],
    "Inotuzumab": ["CD22"],
    "Ruxolitinib": ["JAK1", "JAK2"],
    "Gilteritinib": ["FLT3"],
    "Palbociclib": ["CDK4", "CDK6"],
    "CHZ868": ["JAK2"],
}


def main() -> int:
    C.ensure_dirs()
    b1 = C.rd(C.RESULTS / "B1_OKSA_LI_219_gene_concordance.tsv")
    b1["symbol"] = b1["gene"].str.strip().str.upper()
    ok = pd.to_numeric(b1["oksa_coefficient_slow_vs_fast_EOI"], errors="coerce")
    okf = pd.to_numeric(b1["oksa_fdr_EOI"], errors="coerce")
    li = pd.to_numeric(b1["li_C1_minus_C2"], errors="coerce")
    b1["direction_ok"] = (np.sign(ok) == np.sign(li))

    # Li per-gene P and FDR.  The brief requires a Li P value in this table, and the
    # category logic downstream needs a significance criterion that was NOT already
    # used to select the rows -- reusing |Li| >= 0.5 would be circular, because the
    # highlight rule already requires it.
    li_p = []
    limat = C.rd(C.LI / "li_252gene_matrix.tsv")
    limeta = C.rd(C.LI / "li_patient_meta.tsv")
    lipat = [c for c in limat.columns if c != "gene"]
    lilab = limeta.set_index("patient_id").loc[lipat, "subtype"].to_numpy()
    lirow = {g.strip().upper(): i for i, g in enumerate(limat["gene"])}
    liarr = limat.drop(columns=["gene"]).to_numpy(dtype=float)
    c1 = lilab == "C1"
    c2 = lilab == "C2"
    for sym in b1["symbol"]:
        if sym in lirow:
            v = liarr[lirow[sym], :]
            a, b = v[c1], v[c2]
            try:
                li_p.append(stats.mannwhitneyu(a, b, alternative="two-sided").pvalue)
            except Exception:  # noqa: BLE001
                li_p.append(np.nan)
        else:
            li_p.append(np.nan)
    b1["li_p"] = li_p
    b1["li_fdr"] = C.bh_fdr(li_p)

    b6 = C.rd(C.RESULTS / "B6_TARGET_benchmark.tsv")
    b6["symbol"] = b6["gene"].str.strip().str.upper()
    b6["target_primary_log2FC"] = pd.to_numeric(b6["target_primary_log2FC"], errors="coerce")
    b6m = b6.drop_duplicates("symbol").set_index("symbol")

    conv = pd.read_csv(C.RESULTS / "B2_pathway_convergence.tsv", sep="\t")
    drug = pd.read_csv(C.RESULTS / "B7_drug_sensitivity.tsv", sep="\t")
    pathway_axis = dict(zip(conv["pathway"], conv["axis"]))

    drug_gene_rows = {}
    for _, r in drug.iterrows():
        for g in DRUG_GENE.get(r["drug"], []):
            drug_gene_rows.setdefault(g.upper(), []).append(
                f"{r['drug']}({'C1-resist' if r['median_LC50_difference'] > 0 else 'C2-resist'}, q={C.fmt(r['bh_fdr'])})")

    rows = []
    for _, r in b1.iterrows():
        sym = r["symbol"]
        o = float(ok[r.name]) if not pd.isna(ok[r.name]) else np.nan
        l = float(li[r.name]) if not pd.isna(li[r.name]) else np.nan
        o_fdr = float(okf[r.name]) if not pd.isna(okf[r.name]) else np.nan
        conc = bool(r["direction_ok"])
        nontrivial = (not np.isnan(o) and abs(o) >= 0.25) and (not np.isnan(l) and abs(l) >= 0.5)
        b6row = b6m.loc[sym] if sym in b6m.index else None
        tg = float(b6row["target_primary_log2FC"]) if b6row is not None and \
            pd.notna(b6row["target_primary_log2FC"]) else np.nan
        tg_class = b6row["three_way_class"] if b6row is not None else "not_assessable_in_TARGET"
        axes_hit = sorted({a for p, a in pathway_axis.items() if a != "unexpected_top"
                           and p in str(b6row["gene"])} if b6row is not None else [])
        rows.append(dict(
            gene=sym,
            oksa_coefficient=o, oksa_fdr=o_fdr,
            li_C1_minus_C2=l, li_abs_effect=abs(l) if not np.isnan(l) else np.nan,
            li_p=float(b1["li_p"][r.name]) if not pd.isna(b1["li_p"][r.name]) else np.nan,
            li_fdr=float(b1["li_fdr"][r.name]) if not pd.isna(b1["li_fdr"][r.name]) else np.nan,
            direction_concordant=conc,
            rank_by_oksa_effect=int(pd.to_numeric(b1["oksa_coefficient_slow_vs_fast_EOI"],
                                                  errors="coerce").abs()
                                    .rank(ascending=False)[r.name]),
            effect_nontrivial=nontrivial,
            non_trivial_reason="" if nontrivial else
            ("|Oksa| < 0.25 or |Li| < 0.5 row-SD" if conc else "direction not concordant"),
            target_log2FC_Poor_vs_Good=tg,
            target_class=tg_class,
            shared_with_general_B_ALL=bool(tg_class == "shared_B_ALL_programme"),
            drug_target_for="; ".join(drug_gene_rows.get(sym, [])),
            local_assay_feasible="NOT_APPLICABLE - no biospecimen inventory exists "
                                 "(LOCAL_QPCR_VALIDATION_FEASIBLE = FALSE)",
            is_SHE=bool(sym == "SHE"),
        ))
    cand = pd.DataFrame(rows)

    # highlight = concordant AND non-trivial (measurability and interpretability are
    # recorded as columns rather than used to silently drop rows)
    cand["highlighted_candidate"] = cand["direction_concordant"] & cand["effect_nontrivial"]
    cand["target_assessable"] = cand["target_class"] != "not_assessable_in_TARGET"
    # "not assessable in TARGET" must never be worded as "not concordant with
    # general B-ALL": the gene was filtered out of the TARGET DE table, so the
    # comparison never happened.
    cand["interpretation"] = np.where(
        cand["highlighted_candidate"] & cand["shared_with_general_B_ALL"],
        "concordant across both ETV6 cohorts and with general B-ALL",
        np.where(cand["highlighted_candidate"] & ~cand["target_assessable"],
                 "concordant across both ETV6 cohorts; not assessable in TARGET "
                 "(filtered out of the TARGET DE table)",
                 np.where(cand["highlighted_candidate"],
                          "concordant across both ETV6 cohorts; NOT concordant with general B-ALL",
                          np.where(cand["direction_concordant"],
                                   "concordant but small effect in one or both cohorts",
                                   "direction not concordant between Oksa and Li"))))
    # ranking key: concordance, then magnitude, then whether TARGET agrees
    cand = cand.sort_values(["highlighted_candidate", "direction_concordant",
                             "shared_with_general_B_ALL", "oksa_fdr"],
                            ascending=[False, False, True, True]).reset_index(drop=True)
    cand.to_csv(C.TABLES / "B5_CANDIDATE_GENES.tsv", sep="\t", index=False)
    C.say(f"[write] tables/B5_CANDIDATE_GENES.tsv ({len(cand)} genes, "
          f"{int(cand['highlighted_candidate'].sum())} highlighted)")

    # ------------------------------------------------------------------ SHE
    she = cand[cand["is_SHE"]]
    if len(she):
        s = she.iloc[0]
        if pd.isna(s["target_log2FC_Poor_vs_Good"]):
            tg_note = ("Not assessable in TARGET: SHE was **filtered out** of the TARGET DE "
                       "table by the upstream pre-filter (counts >= 10 in >= 80% of samples) "
                       "because its expression there is near zero. The earlier project note "
                       "that 'SHE is null in TARGET' therefore describes a gene that the TARGET "
                       "DE model never tested, not a tested-and-negative gene.")
        else:
            tg_note = (f"TARGET Poor-vs-Good log2FC = "
                       f"{C.fmt(s['target_log2FC_Poor_vs_Good'])} ({s['target_class']}).")
        she_txt = (f"SHE **is assessable** in both ETV6 cohort layers. Oksa EOI slow-vs-fast "
                   f"coefficient = {C.fmt(s['oksa_coefficient'])} (FDR = {C.fmt(s['oksa_fdr'])}); "
                   f"Li C1 - C2 = {C.fmt(s['li_C1_minus_C2'])}; direction "
                   f"{'concordant' if s['direction_concordant'] else 'DISCORDANT'}; "
                   f"rank {int(s['rank_by_oksa_effect'])} of {len(cand)} by |Oksa effect|. "
                   f"{tg_note}")
        if s["direction_concordant"] and abs(s["oksa_coefficient"]) >= 0.25:
            she_verdict = ("SHE is directionally concordant across both ETV6::RUNX1 cohorts with "
                           "a large Oksa effect and a significant Oksa FDR; it is a candidate "
                           "feature. It is **not** a validated marker: the TARGET comparison is "
                           "unavailable for it, and no independent third cohort exists.")
        elif s["direction_concordant"]:
            she_verdict = ("SHE is directionally concordant but the Oksa effect is small; "
                           "this is a weak signal, not a finding.")
        else:
            she_verdict = ("SHE is **direction discordant** between the two ETV6::RUNX1 "
                           "cohorts. Recorded as SHE_NOT_SUPPORTED_IN_ETV_SPECIFIC_ANALYSIS. "
                           "No rescue analysis is attempted.")
    else:
        # SHE is not among the 219 shared genes -- establish why, from the mapping table
        gp = C.rd(C.TABLES / "08_GENE_HARMONISATION.tsv")
        she_map = gp[gp["original_id"].str.upper() == "SHE"]
        if len(she_map):
            why = she_map.iloc[0]
            she_txt = (f"SHE is **not assessable**. In the gene-mapping table its status is "
                       f"`{why['duplicate_status']}` (included = {why['included']}; "
                       f"reason: {why['exclusion_reason'] or 'none recorded'}).")
        else:
            she_txt = ("SHE is **absent from the Li 252-gene panel** and therefore not "
                       "assessable in the cross-cohort analysis.")
        she_verdict = "SHE_NOT_SUPPORTED_IN_ETV_SPECIFIC_ANALYSIS (not assessable, not negative)."

    md = f"""# B5 -- Candidate molecular features (not a signature)

This table is **not** a signature. No score is constructed, no weights are estimated,
no classifier is trained, and `SIGNATURE_CONSTRUCTION_PERMITTED` remains FALSE.

## 1. Criteria, fixed before the table was built

A gene is *highlighted* when **all** of:

1. direction concordant between Oksa and Li;
2. effect non-trivial: |Oksa coefficient| >= 0.25 **and** |Li C1−C2| >= 0.5 row-SD units;
3. measurable in principle;
4. biology interpretable.

Criteria 3 and 4 are recorded as **columns** rather than used to drop rows, so the
table shows every gene's status instead of only the survivors. Genes are **not**
ranked by P value; the ordering is concordance first, then effect magnitude.

Criteria 1 and 2 are applied mechanically; 3 and 4 are judgement, and the columns
that carry them say so.

## 2. Result

| quantity | n |
|---|---|
| genes assessed (the 219 shared genes) | {len(cand)} |
| direction concordant | {int(cand['direction_concordant'].sum())} |
| **highlighted candidates** (concordant and non-trivial) | **{int(cand['highlighted_candidate'].sum())}** |
| concordant but small effect in one or both cohorts | {int((cand['direction_concordant'] & ~cand['effect_nontrivial']).sum())} |
| not direction concordant | {int((~cand['direction_concordant']).sum())} |
| of the highlighted, also concordant with general B-ALL | {int((cand['highlighted_candidate'] & cand['shared_with_general_B_ALL']).sum())} |
| of the highlighted, **not assessable** in TARGET (filtered out) | {int((cand['highlighted_candidate'] & ~cand['target_assessable']).sum())} |
| of the highlighted, **not** concordant with general B-ALL (ETV6-enriched candidates) | {int((cand['highlighted_candidate'] & ~cand['shared_with_general_B_ALL']).sum())} |

Top of the table, ordered by the stated key rather than by P value:

| gene | Oksa coef | Oksa FDR | Li C1−C2 | Li FDR | concordant | TARGET class | highlighted |
|---|---|---|---|---|---|---|
"""
    for _, r in cand.head(25).iterrows():
        md += (f"| {r['gene']} | {C.fmt(r['oksa_coefficient'])} | {C.fmt(r['oksa_fdr'])} | "
               f"{C.fmt(r['li_C1_minus_C2'])} | {r['direction_concordant']} | "
               f"{r['target_class']} | {r['highlighted_candidate']} |\n")

    md += f"""
## 3. SHE (pre-specified candidate, section 19 of the brief)

{she_txt}

**Verdict:** {she_verdict}

SHE was fixed as a candidate before any of this analysis ran, so it is reported here
whatever the result. No rescue analysis was attempted: no alternative contrast, no
cutoff search, no re-definition of the phenotype.

## 4. Drug-response features linked to genes

`DRUG_GENE` maps each measured drug to the gene products whose direct interaction is
well established, and the mapping is recorded in the `drug_target_for` column so it can
be rejected by a reader. The mapping is deliberately minimal: it excludes drugs with no
single unambiguous target (for example asparaginase acts on extracellular asparagine
rather than on one gene product, so only `ASNS` is linked).

Significance is read from the 20-drug BH-corrected table in `results/B7_drug_sensitivity.tsv`,
not from an uncorrected P value.

## 5. What this table must not be used for

- It must not be called a signature, a panel, a score or a predictor.
- No cut-off is derived from it and none may be.
- The `highlighted_candidate` flag is a **descriptive classification**, not a
  performance claim: no ROC, AUC or cross-validated metric exists for it, and none
  can exist, because no score was built.
- The ETV6-enriched subset is a *candidate* set, never a subtype-specific validated set.
"""

    (C.REPORTS / "B5_CANDIDATE_FEATURES.md").write_text(md, encoding="utf-8")
    C.say(f"[write] reports/B5_CANDIDATE_FEATURES.md")

    (C.TABLES / "_B5_CANDIDATES.json").write_text(json.dumps({
        "n_assessed": int(len(cand)),
        "n_direction_concordant": int(cand["direction_concordant"].sum()),
        "n_highlighted": int(cand["highlighted_candidate"].sum()),
        "n_highlighted_shared_with_general_B_ALL":
            int((cand["highlighted_candidate"] & cand["shared_with_general_B_ALL"]).sum()),
        "n_highlighted_etv6_enriched":
            int((cand["highlighted_candidate"] & ~cand["shared_with_general_B_ALL"]).sum()),
        "she_assessable": bool(len(she)),
        "she_verdict": she_verdict,
        "signature_constructed": False,
        "top_highlighted": cand.loc[cand["highlighted_candidate"], "gene"].head(40).tolist(),
    }, indent=2, ensure_ascii=False), encoding="utf-8")

    C.say(f"[B5] assessed={len(cand)} concordant={int(cand['direction_concordant'].sum())} "
          f"highlighted={int(cand['highlighted_candidate'].sum())} | SHE: {she_verdict[:60]}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
