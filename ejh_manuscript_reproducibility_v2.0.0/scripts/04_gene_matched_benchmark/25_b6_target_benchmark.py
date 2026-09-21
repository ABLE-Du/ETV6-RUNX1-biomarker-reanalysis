#!/usr/bin/env python -u
# -*- coding: utf-8 -*-
"""25_b6_target_benchmark.py -- B6: is the Oksa/Li signal ETV6-enriched or shared
general-B-ALL biology?

The question is decided by *relative* concordance, not by a significance test in
isolation.  Three effect vectors, all with the same sign convention
(positive = higher in the poor-response group):

    Oksa     EOI slow vs fast coefficient
    Li       C1 - C2
    TARGET   DESeq2 log2FoldChange, contrast Poor vs Good, P_BALL n = 191
             (frozen upstream output; contrast verified in the upstream script as
              results(dds, contrast = c("MRD_group", "Poor", "Good")))

Design
  1  concordance of Oksa with Li and with TARGET over the same 219 genes;
  2  a paired comparison (McNemar) of per-gene directional agreement, because the
     two concordance rates are computed on the same genes and are therefore
     dependent;
  3  per-gene classification into shared versus ETV6-enriched candidate.

Naming rule, enforced by the output: a gene concordant with TARGET is labelled
`shared_B_ALL_programme`; a gene concordant with Li but not TARGET is an
`ETV6_enriched_candidate`.  The word **ETV6-specific is forbidden** in this project's
output unless a formal subtype-by-gene interaction test supports it, which the
available data cannot provide.

The upstream protocol-adjusted TARGET design is reported as a sensitivity, because
protocol is the dominant prognostic factor in that cohort and the direction of a gene
can depend on whether it is adjusted for.

Outputs
  results/B6_TARGET_benchmark.tsv
  tables/B6_concordance_comparison.tsv
  reports/B6_TARGET_BENCHMARK.md
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

DE_DIR = C.UPSTREAM / "analysis" / "revision_01" / "results" / "deseq2"


def load_de(name: str) -> pd.DataFrame:
    d = pd.read_csv(DE_DIR / name, sep="\t")
    d["symbol"] = d["gene_symbol"].astype(str).str.strip().str.upper()
    return d


def main() -> int:
    C.ensure_dirs()
    primary = load_de("DEG_MRD_primary_ALL_GENES.tsv.gz")
    protocol = load_de("DEG_MRD_sens_protocol_ALL_GENES.tsv.gz")
    C.say(f"[target] primary DE genes={len(primary)}  protocol-adjusted genes={len(protocol)}")

    b1 = C.rd(C.RESULTS / "B1_OKSA_LI_219_gene_concordance.tsv")
    b1["symbol"] = b1["gene"].str.strip().str.upper()
    b1["oksa"] = pd.to_numeric(b1["oksa_coefficient_slow_vs_fast_EOI"], errors="coerce")
    b1["li"] = pd.to_numeric(b1["li_C1_minus_C2"], errors="coerce")

    def attach(de: pd.DataFrame, tag: str) -> pd.DataFrame:
        m = de.drop_duplicates("symbol").set_index("symbol")
        out = b1.copy()
        out[f"{tag}_log2FC"] = out["symbol"].map(m["log2FoldChange"])
        out[f"{tag}_padj"] = out["symbol"].map(m["padj"])
        return out

    df = attach(primary, "target_primary")
    df["target_protocol_log2FC"] = df["symbol"].map(
        protocol.drop_duplicates("symbol").set_index("symbol")["log2FoldChange"])

    n_no_target = int(df["target_primary_log2FC"].isna().sum())

    def concord(mask_a: pd.Series, mask_b: pd.Series) -> dict:
        m = mask_a.notna() & mask_b.notna()
        if m.sum() < 5:
            return dict(n=int(m.sum()), rho=np.nan, pct=np.nan)
        rho = stats.spearmanr(mask_a[m], mask_b[m]).statistic
        agree = (np.sign(mask_a[m]) == np.sign(mask_b[m]))
        return dict(n=int(m.sum()), rho=float(rho), pct=float(100 * agree.mean()),
                    agree=agree, idx=m)

    ok_li = concord(df["oksa"], df["li"])
    ok_tg = concord(df["oksa"], df["target_primary_log2FC"])
    li_tg = concord(df["li"], df["target_primary_log2FC"])
    ok_tgp = concord(df["oksa"], df["target_protocol_log2FC"])

    # paired McNemar on the genes where all three are available
    m3 = (df["oksa"].notna() & df["li"].notna()
          & df["target_primary_log2FC"].notna()).to_numpy()
    a_li = (np.sign(df["oksa"][m3]) == np.sign(df["li"][m3])).to_numpy()
    a_tg = (np.sign(df["oksa"][m3]) == np.sign(df["target_primary_log2FC"][m3])).to_numpy()
    n_both = int((a_li & a_tg).sum())
    n_li_only = int((a_li & ~a_tg).sum())
    n_tg_only = int((~a_li & a_tg).sum())
    n_neither = int((~a_li & ~a_tg).sum())
    try:
        mcn = stats.binomtest(n_li_only, n_li_only + n_tg_only, 0.5, alternative="two-sided")
        mcn_p = float(mcn.pvalue)
    except Exception:  # noqa: BLE001
        mcn_p = float("nan")

    # per-gene classification on the three-way complete cases.
    # Genes absent from the TARGET filtered table must be labelled not-assessable --
    # a missing comparison is not a discordant one, and letting np.where fall through
    # would silently have classified 106 genes as "discordant_with_both".
    cls = np.full(len(df), "not_assessable_in_TARGET", dtype=object)
    cls[m3] = np.where(
        a_li & a_tg, "shared_B_ALL_programme",
        np.where(a_li & ~a_tg, "ETV6_enriched_candidate",
                 np.where(~a_li & a_tg, "TARGET_only_concordant",
                          "discordant_with_both")))
    df["three_way_class"] = cls

    out = df[["symbol", "oksa", "li", "li_abs_effect", "target_primary_log2FC",
              "target_primary_padj", "target_protocol_log2FC", "three_way_class"]].copy()
    out = out.rename(columns={"symbol": "gene", "oksa": "oksa_coefficient_slow_vs_fast",
                              "li": "li_C1_minus_C2"})
    out["oksa_target_concordant"] = np.sign(out["oksa_coefficient_slow_vs_fast"]) == \
        np.sign(out["target_primary_log2FC"])
    out["oksa_li_concordant"] = np.sign(out["oksa_coefficient_slow_vs_fast"]) == np.sign(out["li_C1_minus_C2"])
    out.to_csv(C.RESULTS / "B6_TARGET_benchmark.tsv", sep="\t", index=False)
    C.say(f"[write] results/B6_TARGET_benchmark.tsv ({len(out)} rows)")

    comp = pd.DataFrame([
        dict(pair="Oksa_vs_Li", contrast_1="EOI slow vs fast", contrast_2="C1 vs C2",
             n=ok_li["n"], spearman_rho=ok_li["rho"], direction_concordant_pct=ok_li["pct"],
             cohort_specificity="both ETV6::RUNX1"),
        dict(pair="Oksa_vs_TARGET_primary", contrast_1="EOI slow vs fast",
             contrast_2="Poor vs Good (P_BALL, ~MRD_group)", n=ok_tg["n"],
             spearman_rho=ok_tg["rho"], direction_concordant_pct=ok_tg["pct"],
             cohort_specificity="TARGET is general B-ALL"),
        dict(pair="Oksa_vs_TARGET_protocol_adjusted", contrast_1="EOI slow vs fast",
             contrast_2="Poor vs Good (P_BALL, protocol + MRD_group)", n=ok_tgp["n"],
             spearman_rho=ok_tgp["rho"], direction_concordant_pct=ok_tgp["pct"],
             cohort_specificity="TARGET general B-ALL, protocol-adjusted (sensitivity)"),
        dict(pair="Li_vs_TARGET_primary", contrast_1="C1 vs C2",
             contrast_2="Poor vs Good (P_BALL)", n=li_tg["n"],
             spearman_rho=li_tg["rho"], direction_concordant_pct=li_tg["pct"],
             cohort_specificity="cross-dataset"),
    ])
    comp.to_csv(C.TABLES / "B6_concordance_comparison.tsv", sep="\t", index=False)
    C.say(f"[write] tables/B6_concordance_comparison.tsv")

    na_target = int((df["three_way_class"] == "not_assessable_in_TARGET").sum())
    etv6_enriched = int((df["three_way_class"] == "ETV6_enriched_candidate").sum())
    shared = int((df["three_way_class"] == "shared_B_ALL_programme").sum())
    tg_only = int((df["three_way_class"] == "TARGET_only_concordant").sum())
    disc = int((df["three_way_class"] == "discordant_with_both").sum())

    stronger_than_target = bool(ok_li["rho"] > ok_tg["rho"])
    verdict = ("ETV6-enriched candidate programme" if (stronger_than_target and mcn_p < 0.05)
               else "shared B-ALL poor-response programme")

    md = f"""# B6 -- TARGET supportive benchmark: ETV6-enriched or shared B-ALL biology?

TARGET is **not** an ETV6::RUNX1 cohort here. Its only permitted role is as a
general-B-ALL benchmark, and its phenotype is a different endpoint (day-29 MRD at
0.01%, versus Oksa day-29 at 0.1% and Li early-induction MRD). Nothing below treats
TARGET as validation of an ETV6 finding.

**Effect vectors compared** (all with the same sign convention, positive = higher in the
poor-response group):

| source | effect | n genes |
|---|---|---|
| Oksa | EOI slow vs fast coefficient | {len(b1)} (the 219 shared genes) |
| Li | C1 − C2 | {len(b1)} |
| TARGET | DESeq2 log2FoldChange, contrast `c("MRD_group","Poor","Good")` | {len(primary):,} tested, {n_no_target} of the 219 absent |

The TARGET contrast was verified in the upstream script
(`analysis/scripts/22_DE_MRD_BALL.R` line 124) rather than assumed, and the frozen
output was reused unchanged.

---

## 1. Relative concordance — the number that decides the question

| comparison | n genes | Spearman rho | direction concordant |
|---|---|---|---|
"""
    for _, r in comp.iterrows():
        md += (f"| {r['pair']} | {r['n']} | {C.fmt(r['spearman_rho'])} | "
               f"{C.fmt(r['direction_concordant_pct'],1)}% |\n")

    md += f"""
Oksa–Li concordance (rho = {C.fmt(ok_li['rho'])}) is
**{'stronger than' if stronger_than_target else 'not stronger than'}** Oksa–TARGET
concordance on the same genes (rho = {C.fmt(ok_tg['rho'])}). Under the
protocol-adjusted TARGET design the comparison is rho = {C.fmt(ok_tgp['rho'])}.

## 2. Paired test (McNemar) on the three-way complete cases

Because both concordance rates are computed on the same genes, they are dependent and
cannot be compared with two independent tests. The paired breakdown is:

| pattern | n |
|---|---|
| concordant with both Li and TARGET | {n_both} |
| concordant with Li only | {n_li_only} |
| concordant with TARGET only | {n_tg_only} |
| concordant with neither | {n_neither} |

Exact McNemar (binomial) on the discordant pairs ({n_li_only} vs {n_tg_only}):
**P = {C.fmt(mcn_p)}**.

## 3. Per-gene classification

| class | n | how it is read |
|---|---|---|
| `shared_B_ALL_programme` | {shared} | agrees with Oksa and with general B-ALL — cannot be attributed to ETV6::RUNX1 |
| `ETV6_enriched_candidate` | {etv6_enriched} | agrees with Oksa and Li but not with general B-ALL — the only class a subtype-enriched claim could be built on |
| `TARGET_only_concordant` | {tg_only} | agrees with Oksa and TARGET but not Li |
| `discordant_with_both` | {disc} | agrees with neither |
| `not_assessable_in_TARGET` | {na_target} | absent from the TARGET filtered table — reported as not assessable, never as discordant |

## 4. Verdict

**{verdict}.**

{"The Oksa–Li concordance is significantly stronger than Oksa–TARGET concordance on the same genes (McNemar P = " + C.fmt(mcn_p) + "), so the shared signal is not fully explained by general B-ALL biology." if (stronger_than_target and mcn_p < 0.05) else "Oksa–Li concordance is not demonstrably stronger than Oksa–TARGET concordance on the same genes (McNemar P = " + C.fmt(mcn_p) + "). The concordant genes must therefore be described as a **shared B-ALL poor-response programme**, and no ETV6-enriched claim is made."}

Permitted labels in the manuscript: *shared B-ALL poor-response programme* and
*ETV6-enriched candidate programme*.

Forbidden label: **ETV6-specific**. No formal gene-by-subtype interaction test is
possible here — there is no ETV6-negative comparator cohort in this project — so the
word is not used anywhere.

## 5. Limitations

1. TARGET's endpoint is day-29 MRD at **0.01%**; Oksa's is day-29 at **0.1%**; Li's is
   early-induction MRD. The three phenotypes are related but not identical, and no
   endpoint was altered to improve agreement.
2. TARGET is a general B-ALL cohort and its composition differs from an
   ETV6::RUNX1-enriched series; a gene can differ in direction for that reason alone.
3. The 219 genes are a targeted-panel/transcriptome intersection, not a random sample.
4. Protocol adjustment changes some TARGET directions, which is why both designs are
   shown; the primary reading uses the unadjusted design, matching the unadjusted
   Oksa and Li contrasts.
5. {n_no_target} of the 219 genes are absent from the TARGET filtered table and are
   reported as not assessable there, not as discordant.
"""

    (C.REPORTS / "B6_TARGET_BENCHMARK.md").write_text(md, encoding="utf-8")
    C.say(f"[write] reports/B6_TARGET_BENCHMARK.md")

    (C.TABLES / "_B6_BENCHMARK.json").write_text(json.dumps({
        "oksali_rho": ok_li["rho"], "oksali_pct": ok_li["pct"], "oksali_n": ok_li["n"],
        "oksatarget_rho": ok_tg["rho"], "oksatarget_pct": ok_tg["pct"], "oksatarget_n": ok_tg["n"],
        "oksatarget_protocol_rho": ok_tgp["rho"],
        "mcnemar_p": mcn_p, "n_li_only": n_li_only, "n_tg_only": n_tg_only,
        "n_both": n_both, "n_neither": n_neither,
        "shared_B_ALL_programme_n": shared,
        "ETV6_enriched_candidate_n": etv6_enriched,
        "target_only_concordant_n": tg_only, "discordant_with_both_n": disc, "not_assessable_in_target_n": na_target,
        "n_genes_absent_from_target": n_no_target,
        "verdict": verdict,
        "etv6_specific_claim_made": False,
    }, indent=2, ensure_ascii=False), encoding="utf-8")

    C.say(f"[B6] Oksa-Li rho={ok_li['rho']:.3f} ({ok_li['pct']:.1f}%) vs "
          f"Oksa-TARGET rho={ok_tg['rho']:.3f} ({ok_tg['pct']:.1f}%) | "
          f"McNemar P={mcn_p:.4f} | shared={shared} ETV6_enriched={etv6_enriched} | {verdict}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
