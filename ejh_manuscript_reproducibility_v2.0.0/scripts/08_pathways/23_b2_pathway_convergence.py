#!/usr/bin/env python -u
# -*- coding: utf-8 -*-
"""23_b2_pathway_convergence.py -- B2 step 2: pathway-level *convergence*.

Wording discipline: the Li cohort has a 252-gene panel, so it cannot support a
genome-wide GSEA.  Nothing here is called "independent pathway validation".  What is
tested is whether pathway membership is associated with cross-cohort agreement,
and whether the Oksa-enriched pathways are over-represented among the Li genes that
agree.  That is *pathway-level convergence*.

Three questions, each with an explicit universe (the universe is stated because a
hypergeometric test is meaningless without one):

  A  are the Li panel genes over-represented among Oksa-enriched pathway members?
     universe = the 15,928 Oksa-ranked genes; sample = the 219 shared genes
  B  within the 219 shared genes, is pathway membership associated with
     direction concordance?
     universe = 219; sample = the concordant genes
  C  within the 219 shared genes, is pathway membership associated with a strong
     Li C1-vs-C2 effect?
     universe = 219; sample = top-quartile |Li effect|

Also reports every pre-specified biological axis and every unexpected top pathway,
so the axis table is not a significance-filtered view.

Outputs
  results/B2_pathway_convergence.tsv
  tables/B2_pathway_axes.tsv
  reports/B2_PATHWAY_CONVERGENCE.md
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

# pre-specified axes -> pathway-name patterns (MSigDB name substrings)
AXES = {
    "B_cell_developmental_state": ["B_CELL", "LYMPHOCYTE", "VDJ", "IMMUNOGLOBULIN",
                                  "GERMINAL_CENTER", "B_CELL_RECEPTOR", "PRE_B"],
    "cell_cycle_E2F": ["CELL_CYCLE", "E2F", "DNA_REPLICATION", "MITOTIC", "G2M",
                       "S_PHASE", "MITOSIS", "CHECKPOINT"],
    "MYC": ["MYC"],
    "TNFA_NFKB": ["TNFA", "NFKB", "TNF_", "_TNF"],
    "JAK_STAT": ["JAK", "STAT"],
    "hypoxia": ["HYPOXIA", "HIF"],
    "apoptosis": ["APOPTOSIS", "BCL2", "CASPASE", "DEATH", "BAX"],
    "DNA_repair": ["DNA_REPAIR", "MISMATCH_REPAIR", "HOMOLOGOUS_RECOMBINATION",
                   "NUCLEOTIDE_EXCISION", "BASE_EXCISION", "DOUBLE_STRAND_BREAK",
                   "FANCONI", "TRANSLESION"],
    "RAS_MAPK": ["RAS", "MAPK", "ERK", "RAF", "MEK"],
    "BCR_signalling": ["BCR_SIGNALING", "BCR_PATHWAY", "SIGNALING_BY_BCR"],
    "metabolic_programmes": ["OXIDATIVE_PHOSPHORYLATION", "GLYCOLYSIS", "METABOLISM",
                             "LIPID", "AMINO_ACID", "FATTY_ACID", "TRANSLATION",
                             "RIBOSOME", "CHOLESTEROL", "TCA_CYCLE"],
    "unexpected_top": [],   # filled with the top pathways not matching any axis
}


def axis_of(pathway: str) -> str:
    p = pathway.upper()
    for ax, pats in AXES.items():
        if ax == "unexpected_top":
            continue
        if any(pat in p for pat in pats):
            return ax
    return "unexpected_top"


def main() -> int:
    C.ensure_dirs()
    gsea = pd.read_csv(C.RESULTS / "B2_oksa_fgsea_all_pathways.tsv.gz", sep="\t")
    b1 = C.rd(C.RESULTS / "B1_OKSA_LI_219_gene_concordance.tsv")
    genes = [g.strip().upper() for g in b1["gene"]]
    concordant = [g.strip().upper() for g in
                  b1.loc[b1["direction_concordant"].isin(["positive_positive",
                                                          "negative_negative"]), "gene"]]
    li_abs = pd.to_numeric(b1["li_abs_effect"], errors="coerce").to_numpy()
    q75 = float(np.nanpercentile(li_abs, 75))
    strong_li = [g for g, v in zip(genes, li_abs) if not np.isnan(v) and v >= q75]

    oksa_universe = int(pd.read_csv(C.ROOT / "intermediate" / "B2_oksa_fgsea_qc.tsv",
                                    sep="\t").set_index("item").loc["unique_symbols", "value"])

    sets = {}
    for _, r in gsea.iterrows():
        le = [x.strip().upper() for x in str(r["leadingEdge"]).split(";") if x.strip()]
        sets[r["pathway"]] = le

    gset = set(genes)
    cons_set = set(concordant)
    strong_set = set(strong_li)

    rows = []
    for _, r in gsea.iterrows():
        le = sets[r["pathway"]]
        le_shared = [g for g in le if g in gset]
        le_oksa = [g for g in le if g]
        if not le_shared and not le_oksa:
            continue
        # A: universe = all Oksa-ranked genes, sample = the 219 shared genes
        a_a = len(set(le_oksa))
        A = stats.hypergeom.sf(len(set(le_shared)) - 1, oksa_universe, a_a,
                               len(gset)) if a_a else np.nan
        # B: universe = 219, sample = concordant
        b_hit = len(set(le_shared) & cons_set)
        B = stats.hypergeom.sf(b_hit - 1, len(gset), len(le_shared), len(cons_set)) \
            if le_shared else np.nan
        # C: universe = 219, sample = strong Li effect
        c_hit = len(set(le_shared) & strong_set)
        Cp = stats.hypergeom.sf(c_hit - 1, len(gset), len(le_shared), len(strong_set)) \
            if le_shared else np.nan
        rows.append(dict(
            pathway=r["pathway"], collection=r["collection"], oksa_NES=r["NES"],
            oksa_padj=r["padj"], oksa_direction=r["direction"],
            pathway_size=len(set(le_oksa)), n_in_shared_219=len(le_shared),
            axis=axis_of(r["pathway"]),
            A_shared_genes_in_pathway=len(set(le_shared)),
            A_panel_enrichment_p=A,
            B_concordant_in_pathway=b_hit, B_concordant_p=B,
            C_strongLi_in_pathway=c_hit, C_strongLi_p=Cp,
        ))
    res = pd.DataFrame(rows)

    # BH within each question's family (all pathways with a defined test)
    for col, pcol, qcol in [("A", "A_panel_enrichment_p", "A_panel_q"),
                            ("B", "B_concordant_p", "B_concordant_q"),
                            ("C", "C_strongLi_p", "C_strongLi_q")]:
        m = res[pcol].notna()
        res[qcol] = np.nan
        if m.any():
            res.loc[m, qcol] = C.bh_fdr(res.loc[m, pcol].tolist())
    res = res.sort_values(["B_concordant_q", "oksa_padj"]).reset_index(drop=True)
    res.to_csv(C.RESULTS / "B2_pathway_convergence.tsv", sep="\t", index=False)
    C.say(f"[write] results/B2_pathway_convergence.tsv ({len(res)} pathways)")

    # ------------------------------------------------------------ axis table
    axis_rows = []
    for ax in AXES:
        g = res[res["axis"] == ax]
        if not len(g):
            axis_rows.append(dict(axis=ax, n_pathways=0, n_padj_lt_0_05=0,
                                  best_pathway="", best_NES=np.nan, best_padj=np.nan,
                                  best_direction="", best_B_q=np.nan))
            continue
        b = g.sort_values("oksa_padj").iloc[0]
        axis_rows.append(dict(
            axis=ax, n_pathways=len(g), n_padj_lt_0_05=int((g["oksa_padj"] < 0.05).sum()),
            best_pathway=b["pathway"], best_NES=b["oksa_NES"], best_padj=b["oksa_padj"],
            best_direction=b["oksa_direction"],
            best_B_q=b.get("B_concordant_q", np.nan)))
    axdf = pd.DataFrame(axis_rows).sort_values(["n_padj_lt_0_05", "best_padj"],
                                              ascending=[False, True])
    axdf.to_csv(C.TABLES / "B2_pathway_axes.tsv", sep="\t", index=False)
    C.say(f"[write] tables/B2_pathway_axes.tsv ({len(axdf)} pre-specified axes)")

    n_sig = int((gsea["padj"] < 0.05).sum())
    n_shared_pathways = int(res["n_in_shared_219"].gt(0).sum())
    top_shared = res[res["n_in_shared_219"] > 0].sort_values("B_concordant_q").head(15)
    b_sig = int((res["B_concordant_q"] < 0.05).sum())
    direct = res[res["oksa_NES"] > 0]
    direct_pct = 100 * (direct["B_concordant_p"] < 0.5).mean() if len(direct) else np.nan

    md = f"""# B2 -- Pathway-level convergence (Oksa x Li)

Ranking: the **EOI slow-vs-fast** coefficient for all **{oksa_universe:,}** genes in Oksa
Supplementary Table 22 — the full genome-scale ranking, **not** a DEG list. Pre-ranked
fgsea (Hallmark, Reactome, GO:BP, KEGG). All {len(gsea):,} tested pathways are saved in
`results/B2_oksa_fgsea_all_pathways.tsv.gz` with NES, p, padj and leadingEdge;
**{n_sig}** reach padj < 0.05, and significance is an extra column rather than a filter.

## What is and is not claimed here

The Li cohort has a **252-gene panel**, so a genome-wide GSEA cannot be run on it.
Accordingly:

- the Oksa pathway result is a genuine genome-scale enrichment result;
- the Li contribution is an **over-representation test**, i.e. are Oksa-enriched pathways
  disproportionately represented among the Li genes that agree;
- this is called **pathway-level convergence**, and the phrase
  *"Li independently validates the pathways"* is not used anywhere, because a 252-gene
  panel cannot support the claim.

## 1. The Oksa pathway landscape

| collection | pathways tested | padj < 0.05 |
|---|---|---|
"""
    for coll, g in gsea.groupby("collection"):
        md += f"| {coll} | {len(g)} | {int((g['padj'] < 0.05).sum())} |\n"

    md += f"""
Top 12 pathways by padj. Every one of them is enriched in **slow** responders
(NES > 0):

| pathway | collection | NES | p | padj | axis |
|---|---|---|---|---|---|
"""
    for _, r in gsea.sort_values("padj").head(12).iterrows():
        md += (f"| {r['pathway']} | {r['collection']} | {C.fmt(r['NES'])} | "
               f"{C.fmt(r['pval'])} | {C.fmt(r['padj'])} | {axis_of(r['pathway'])} |\n")

    md += f"""
The dominant theme is a **TNFA/NFKB inflammatory programme together with a broad
translational / ribosomal programme** (eukaryotic translation initiation and elongation,
EIF2AK4-GCN2 amino-acid deficiency response, rRNA processing, SRP-dependent targeting,
ribosome-associated quality control, nonsense-mediated decay, KEGG ribosome) plus
erythroid and myeloid homeostasis gene sets. Whether the haematopoietic-lineage signal
reflects a less mature or more multilineage transcriptome in slow responders is a
hypothesis, not a finding; nothing here establishes a mechanism.

## 2. Pre-specified biological axes

`tables/B2_pathway_axes.tsv` — all twelve axes are reported, including those with no
significant pathway, so the axis table is not a significance-filtered view.

| axis | pathways | padj<0.05 | best pathway | NES | padj | direction |
|---|---|---|---|---|---|---|
"""
    for _, r in axdf.iterrows():
        md += (f"| {r['axis']} | {r['n_pathways']} | {r['n_padj_lt_0_05']} | "
               f"{r['best_pathway']} | {C.fmt(r['best_NES'])} | {C.fmt(r['best_padj'])} | "
               f"{r['best_direction']} |\n")

    md += f"""
## 3. Convergence tests

Universes are stated because a hypergeometric test is meaningless without one.

| question | universe | sample | n pathways tested | q < 0.05 |
|---|---|---|---|---|
| A: are the Li panel genes concentrated in Oksa pathway members? | {oksa_universe:,} Oksa-ranked genes | the {len(gset)} shared genes | {int(res['A_panel_enrichment_p'].notna().sum())} | {int((res['A_panel_q'] < 0.05).sum())} |
| B: within the {len(gset)} shared genes, does pathway membership predict direction concordance? | {len(gset)} shared genes | {len(cons_set)} concordant genes | {int(res['B_concordant_p'].notna().sum())} | {b_sig} |
| C: within the {len(gset)} shared genes, does pathway membership predict a strong Li effect? | {len(gset)} shared genes | {len(strong_set)} top-quartile genes | {int(res['C_strongLi_p'].notna().sum())} | {int((res['C_strongLi_q'] < 0.05).sum())} |

Pathways with at least one shared gene: {n_shared_pathways}.

### Strongest convergence signals (question B)

| pathway | collection | Oksa NES | Oksa padj | shared genes | concordant | q |
|---|---|---|---|---|---|---|
"""
    for _, r in top_shared.iterrows():
        md += (f"| {r['pathway']} | {r['collection']} | {C.fmt(r['oksa_NES'])} | "
               f"{C.fmt(r['oksa_padj'])} | {r['n_in_shared_219']} | {r['B_concordant_in_pathway']} | "
               f"{C.fmt(r['B_concordant_q'])} |\n")

    md += f"""
## 4. Interpretation

{"**Supported.** " + f"{b_sig} pathways show over-representation of concordant genes beyond chance within the shared gene space (question B), and the Oksa-enriched pathways are enriched for shared genes (question A, {int((res['A_panel_q'] < 0.05).sum())} pathways). This supports a statement of *pathway-level convergence* between the two ETV6::RUNX1 cohorts." if (b_sig or int((res['A_panel_q'] < 0.05).sum())) else "**Not supported.** No pathway shows significant over-representation of concordant genes after correction. Reported as a null result; the analysis was not redesigned."}

The following are explicitly **not** claimed:

- that Li independently validates any pathway (a 252-gene panel cannot);
- that any pathway is ETV6::RUNX1-specific (no ETV6-negative comparator was analysed
  here; rung B6 addresses the general-B-ALL alternative);
- causality, or that a pathway drives poor response.

## 5. Limitations

1. The Oksa ranking is a published effect vector from the authors' model, not our fit.
2. Duplicate gene symbols in the Oksa table (33 symbols) were collapsed by the mean
   coefficient; the count is recorded in `intermediate/B2_oksa_fgsea_qc.tsv`.
3. GO:BP contributes most pathways and therefore most of the multiple-testing burden;
   collection-level results are reported separately above rather than only pooled.
4. Pathway membership is correlated within a collection, so BH across 6,009 pathways is
   anti-conservative in the usual way for gene-set testing.
"""

    (C.REPORTS / "B2_PATHWAY_CONVERGENCE.md").write_text(md, encoding="utf-8")
    C.say(f"[write] reports/B2_PATHWAY_CONVERGENCE.md")

    (C.TABLES / "_B2_CONVERGENCE.json").write_text(json.dumps({
        "n_pathways_tested": int(len(gsea)), "n_padj_lt_0_05": n_sig,
        "per_collection": {c: int(len(g)) for c, g in gsea.groupby("collection")},
        "n_pathways_with_shared_gene": n_shared_pathways,
        "B_q_lt_0.05": b_sig, "A_q_lt_0.05": int((res["A_panel_q"] < 0.05).sum()),
        "C_q_lt_0.05": int((res["C_strongLi_q"] < 0.05).sum()),
        "oksa_universe": oksa_universe, "n_shared_genes": len(gset),
        "n_concordant": len(cons_set), "n_strong_li": len(strong_set),
        "axes": axdf.to_dict("records"),
    }, indent=2, ensure_ascii=False), encoding="utf-8")

    C.say(f"[B2] pathways={len(gsea)} sig={n_sig} | with shared gene={n_shared_pathways} | "
          f"A_q<0.05={int((res['A_panel_q'] < 0.05).sum())} B_q<0.05={b_sig} "
          f"C_q<0.05={int((res['C_strongLi_q'] < 0.05).sum())}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
