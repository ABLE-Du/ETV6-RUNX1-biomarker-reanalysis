#!/usr/bin/env python -u
# -*- coding: utf-8 -*-
"""08_gene_harmonisation.py -> tables/08_GENE_HARMONISATION.tsv
                            -> reports/08_GENE_HARMONISATION_PLAN.md

Round-1 deliverable 8: the gene space.

Two cohorts have to be joined: Oksa (Ensembl id + symbol, 19,588 genes, GRCh38 /
Ensembl 103) and Li (symbol only, 252 genes).  The join is symbol-mediated, and
symbols are ambiguous, so the honest artefact is a per-gene mapping table that
records every decision -- including the ones that end in "not mapped".

Only the overlap can carry cross-cohort claims, so this script computes the size of
that overlap rather than assuming it is 252.
"""
from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OKSA = ROOT / "intermediate" / "oksa"
LI = ROOT / "intermediate" / "li"
TABLES = ROOT / "tables"
REPORTS = ROOT / "reports"

COLUMNS = ["dataset", "original_id", "ensembl", "symbol", "mapping_source",
           "duplicate_status", "included", "exclusion_reason"]

# the Oksa EOI contrast, read by its composed block name
OKSA_EFFECT = "EOI slow VS fast | Coefficient slow VS fast"
OKSA_P = "EOI slow VS fast | P-value"
OKSA_FDR = "EOI slow VS fast | Adjusted p-value"

# Explicit, curated old -> current symbol renames.
#
# Only entries that were verified against BOTH sides are listed: the old symbol is
# absent from Oksa Supplementary Table 22 AND the new symbol is present.  Nothing
# here is inferred by fuzzy matching -- a fuzzy probe over the 38 unmapped symbols
# produced obviously wrong hits (e.g. "C3orf56" -> "C3", "HRC" -> "CTHRC1"), which
# is exactly why automatic similarity matching is forbidden in this project.
SYMBOL_ALIASES = {
    "DFNA5": "GSDME",
    "PVRL2": "NECTIN2",
    "INADL": "PATJ",
    "LHFP": "LHFPL6",
    "KIAA1211": "CRACD",
}


def say(*a):
    print(*a, flush=True)


def main() -> int:
    TABLES.mkdir(parents=True, exist_ok=True)
    oksa = pd.read_csv(OKSA / "SuppTable22_DE_analyses.tsv", sep="\t",
                       dtype=str, keep_default_na=False)
    li = pd.read_csv(LI / "li_252gene_matrix.tsv", sep="\t", dtype=str,
                     keep_default_na=False)
    li_genes = [g.strip() for g in li["gene"] if g.strip()]

    o_ens = oksa["Gene Ensembl ID"].str.strip()
    o_sym = oksa["Gene symbol"].str.strip()
    o_sym_upper = o_sym.str.upper()
    li_upper = [g.upper() for g in li_genes]

    # duplicate symbols inside Oksa
    sym_counts = o_sym_upper.value_counts()
    dup_syms = set(sym_counts[sym_counts > 1].index)

    rows = []
    # ---- Oksa side, all 19,588 genes ---------------------------------------
    for ens, sym, symu in zip(o_ens, o_sym, o_sym_upper):
        dup = "duplicate_symbol_in_oksa" if symu in dup_syms else "unique"
        if not sym:
            rows.append(["OKSA", ens, ens, "", "ensembl_from_source",
                         "no_symbol", "False", "oksa_row_has_no_gene_symbol"])
        else:
            rows.append(["OKSA", ens, ens, sym, "ensembl_from_source", dup, "True", ""])

    # ---- the 252 Li genes mapped onto Oksa --------------------------------
    oksa_by_sym: dict[str, list[str]] = {}
    for ens, symu in zip(o_ens, o_sym_upper):
        if symu:
            oksa_by_sym.setdefault(symu, []).append(ens)

    mapped = not_mapped = ambiguous = aliased = 0
    pair_rows = []
    for g, gu in zip(li_genes, li_upper):
        via = "symbol_exact"
        lookup = gu
        if gu not in oksa_by_sym and gu in SYMBOL_ALIASES:
            new = SYMBOL_ALIASES[gu].upper()
            if new in oksa_by_sym:
                lookup = new
                via = f"symbol_alias:{gu}->{new}"
                aliased += 1
        hits = oksa_by_sym.get(lookup, [])
        if len(hits) == 1:
            mapped += 1
            status, ens, incl, why = via, hits[0], "True", ""
        elif len(hits) > 1:
            ambiguous += 1
            status, ens, incl, why = ("multiple_ensembl_for_symbol", ";".join(hits),
                                      "False", "ambiguous_symbol_requires_gene_id")
        else:
            not_mapped += 1
            status, ens, incl, why = ("not_in_oksa", "", "False",
                                      "absent_from_oksa_table_22_filtered_gene_universe")
        pair_rows.append([g, ens, status, incl])
        rows.append(["LI", g, ens, g, "symbol_from_supp_figure_1", status, incl, why])

    n_li = len(li_genes)
    n_oksasym_unique = int(o_sym_upper.replace("", pd.NA).dropna().nunique())
    n_dup_rows = int(o_sym_upper.duplicated().sum())
    out = TABLES / "08_GENE_HARMONISATION.tsv"
    with out.open("w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(COLUMNS)
        w.writerows(rows)
    say(f"[write] {out.relative_to(ROOT)}  ({len(rows)} rows)")

    # ---- the join table that rung B1 will actually use ---------------------
    # Both symbols are carried: the Li side must be looked up in the Li matrix by the
    # ORIGINAL symbol, while the Oksa side is keyed by the Oksa symbol.  For the five
    # alias-recovered genes these differ (e.g. Li 'DFNA5' vs Oksa 'GSDME'), and keeping
    # only one of them silently dropped those genes from the concordance analysis.
    keep = {g: e for g, e, st, inc in pair_rows if inc == "True"}
    join = oksa[o_ens.isin(set(keep.values()))][
        ["Gene Ensembl ID", "Gene symbol", OKSA_EFFECT, OKSA_P, OKSA_FDR]].copy()
    join = join.rename(columns={OKSA_EFFECT: "oksa_coefficient_slow_vs_fast_EOI",
                                OKSA_P: "oksa_p_EOI",
                                OKSA_FDR: "oksa_fdr_EOI",
                                "Gene symbol": "oksa_symbol"})
    # map back to the Li original symbol for every included gene
    li_orig = {e: g for g, e, st, inc in pair_rows if inc == "True"}
    join["li_original_symbol"] = join["Gene Ensembl ID"].map(li_orig)
    join = join[["Gene Ensembl ID", "oksa_symbol", "li_original_symbol",
                 "oksa_coefficient_slow_vs_fast_EOI", "oksa_p_EOI", "oksa_fdr_EOI"]]
    join.to_csv(TABLES / "08_OKSA_LI_gene_overlap.tsv", sep="\t", index=False)
    say(f"[write] tables/08_OKSA_LI_gene_overlap.tsv ({len(join)} genes, "
        f"{int((join['oksa_symbol'] != join['li_original_symbol']).sum())} with differing symbols)")

    ae = pd.to_numeric(oksa["Average expression"], errors="coerce")
    ae_min, ae_max = float(ae.min()), float(ae.max())
    n_le0 = int((ae <= 0).sum())

    md = f"""# 08 -- Gene harmonisation plan

Machine-readable: `tables/08_GENE_HARMONISATION.tsv` (every gene, every decision) and
`tables/08_OKSA_LI_gene_overlap.tsv` (the {len(join)} genes that can carry a cross-cohort claim).

---

## The problem

| cohort | identifier given | genes | build / annotation |
|---|---|---|---|
| OKSA | Ensembl gene id **and** symbol | {len(oksa):,} | GRCh38, Ensembl 103 (declared by the submitter) |
| LI | **symbol only** | {n_li} | not declared in the supplementary matrix |
| TARGET | Ensembl id + symbol | per `gene_annotation_aligned.tsv` | GRCh38 asserted from MAF contigs only — never machine-verified for the expression layer |

The Oksa↔Li join is therefore symbol-mediated, and symbols are not unique.

## Result of the join (computed, not assumed)

| outcome | genes | consequence |
|---|---|------|
| unique symbol, mapped to exactly one Oksa Ensembl id | **{mapped}** | usable for rung B1 |
| of which recovered only via an explicit old->current symbol rename | {aliased} | usable, but flagged with the alias used |
| symbol absent from Oksa Supplementary Table 22 | {not_mapped} | dropped with an explicit reason, never silently |
| symbol maps to several Oksa Ensembl ids | {ambiguous} | held out pending a gene-id resolution |

**The overlap is {mapped}, not 252.** Any cross-cohort statement must therefore be phrased
against {mapped} genes, and the {not_mapped + ambiguous} non-mapped genes must be listed as
not-assessable rather than as negative findings.

### Why {not_mapped} genes are absent — and why that is not a negative result

Oksa Supplementary Table 22 is **not a whole-transcriptome table**. It contains {len(oksa):,} rows
but only **{n_oksasym_unique:,} unique gene symbols** ({n_dup_rows:,} rows are duplicate symbols), and it
omits genes that a whole-genome DE would normally carry (for example `FGF6`, `DSCAM`, `IL5`,
`EYA1` are all absent while `ACTB` and `GAPDH` are present). Its `Average expression` column
runs from {ae_min:.2f} to {ae_max:.2f} with {n_le0:,} rows at or below zero, i.e. it is a centred log
scale over a **pre-filtered gene universe**.

The consequence is a wording rule with teeth:

> A Li gene that cannot be found in Oksa Supplementary Table 22 is **not assessable** in that
> cohort. It must never be reported as "not differentially expressed", "no effect" or
> "not replicated".

### Symbol-version drift, and why alias matching is explicit

{aliased} of the {mapped} mapped genes were only recoverable after applying an **explicit,
individually verified** old->current rename:

| old symbol | current symbol |
|---|---|
"""
    for old, new in SYMBOL_ALIASES.items():
        md += f"| `{old}` | `{new}` |\n"

    md += f"""
Automatic similarity matching is **forbidden** here. A fuzzy probe over the unmapped symbols
produced obviously wrong suggestions (`C3orf56` -> `C3`, `HRC` -> `CTHRC1`, `IL5` -> `IL5RA`),
so any future alias addition must be verified the same way these five were: the old symbol
absent from Oksa **and** the new symbol present.

Oksa-side duplicate symbols: {len(dup_syms)} symbols are used by more than one Ensembl id
and are flagged `duplicate_symbol_in_oksa` in the mapping table.

## Rules (frozen in `config/09_ANALYSIS_LOCK.yaml`)

1. **Ensembl gene id is the primary key** wherever it exists. Symbols are a display and
   join convenience only.
2. **A duplicate symbol is never silently dropped.** Each occurrence is emitted with
   `duplicate_status` and `included`, so the decision is auditable; ambiguous symbols are
   held out, not collapsed to a "representative" gene.
3. **Cross-platform comparison uses direction, rank and z-score only.** Absolute expression
   is never merged: the Li matrix is row-standardised (so it carries no absolute scale),
   Oksa's table carries coefficients from its own model, and TARGET uses a different pipeline.
4. **Unmapped genes are recorded, never imputed**, and are reported as *not assessable*.
5. Counts and TPM are never mixed: integer counts go to count models, TPM only to
   visualisation and gene-set scoring.

## The Oksa contrast that must be read

Supplementary Table 22 is a **five-block wide sheet**: `Mid-induction slow VS fast`,
`EOI slow VS fast`, `EOI slow VS fast + intermediate`, `EOI slow + intermediate VS fast`,
`EOC slow VS fast`, each with its own coefficient / P / adjusted P. The extraction composes
the two header rows so each block is explicitly named; the project's pre-specified primary
contrast is **`EOI slow VS fast`**, and the overlap table uses exactly that block.

Note the intermediate group appears on both sides in the two middle blocks. Those two blocks
are sensitivity framings of the same ordinal data and must never be treated as independent
replications of the primary contrast.
"""

    (REPORTS / "08_GENE_HARMONISATION_PLAN.md").write_text(md, encoding="utf-8")
    say(f"[write] reports/08_GENE_HARMONISATION_PLAN.md")
    say(f"[join] mapped={mapped} not_mapped={not_mapped} ambiguous={ambiguous} "
        f"oksa_dup_symbols={len(dup_syms)}")

    (TABLES / "_GENE_HARMONISATION.json").write_text(json.dumps({
        "oksa_genes": len(oksa), "li_genes": n_li,
        "mapped_unique": mapped, "mapped_via_alias": aliased, "oksa_unique_symbols": n_oksasym_unique, "oksa_duplicate_rows": n_dup_rows, "symbol_aliases": SYMBOL_ALIASES, "not_mapped": not_mapped, "ambiguous": ambiguous,
        "oksa_duplicate_symbols": len(dup_syms),
        "join_table_genes": len(join),
        "primary_oksa_contrast": "EOI slow VS fast",
    }, indent=2, ensure_ascii=False), encoding="utf-8")
    return 0


if __name__ == "__main__":
    sys.exit(main())
