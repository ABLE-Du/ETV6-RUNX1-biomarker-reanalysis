# Reproducibility deposit — cross-cohort convergence of adverse early-treatment-response states in ETV6::RUNX1-positive childhood ALL

Manuscript: *Cross-cohort convergence of transcriptional states associated with adverse early treatment
response in ETV6::RUNX1-positive childhood acute lymphoblastic leukaemia* (submitted to
*European Journal of Haematology*).

This audited public deposit contains the **analysis scripts, environment specifications and non-identifiable derived
tables** that reproduce every number reported in the paper. It does **not** contain any source data that
must be obtained from the original publishers, and it contains **no individual-level institutional data**.

Repository release: **v2.0.0** (EJH manuscript deposit v4.0). Machine-specific paths from the executed
analysis were replaced by documented environment variables before public release; this does not alter
the statistical methods or deposited results.

## 1. What is NOT included (and why)

| Not included | Reason | Where to obtain |
|---|---|---|
| Oksa et al. supplementary tables | publisher copyright | supplementary files of Leukemia 2025;39:2125–2139 (doi:10.1038/s41375-025-02683-7) |
| Li et al. source data workbook | publisher copyright | Europe PMC, PMCID PMC11779914 (supplementary file of Nat Commun 2025;16:1153) |
| Oksa RNA-seq matrices | inspected only, not used for any result | GEO GSE227832 / GSE228632 |
| TARGET-ALL-P2 counts | controlled access | NCI Genomic Data Commons, project TARGET-ALL-P2 |
| EGA controlled-access datasets | not accessed or used | EGAD00001010164 (not used) |
| Institutional row-level data | individual clinical records under institutional privacy/ethics restrictions | not publicly shared; aggregated results are in `data/derived_tables/` |
| Figure 3C individual flow values | ethics permission pending (author action A8) | not deposited |

## 2. Directory map

```
environment/      Python requirements.txt, conda environment.yml, R sessionInfo()
scripts/          analysis scripts grouped by the analysis they support (see REPRODUCTION_ORDER.md)
data/derived_tables/      non-identifiable derived result tables cited in the paper
data/figure_data/         the TSV behind every figure panel (except Figure 3C individual points)
data/figure_source_data/  figure source-data workbooks, exactly as submitted to the journal
journal_tables/           Table 1 and Tables S1–S10, editable Excel, as submitted
provenance/               SHA-256 checksums and the reproduction order
```

## 3. Software environment

* Python 3.13 (numpy 2.5, pandas 3.0, scipy 1.18, statsmodels 0.14, matplotlib 3.11, openpyxl 3.1)
* R 4.6.1 (DESeq2 1.52.0, fgsea 1.38.0, msigdbr 26.1.1, data.table 1.18.4)
* Random seeds are fixed inside the scripts (20260918 / 20260919 / 20260915; bootstrap and permutation
  seeds are script-local). Re-running reproduces the reported values to the printed precision.

## 4. Reproduction order

See `provenance/REPRODUCTION_ORDER.md`. In short: acquire the published source data → run the audit and
harmonisation scripts → run the primary correspondence → run the gene-matched benchmark → run the Li95
validation → run the institutional analysis → run candidates/drugs/pathways → regenerate figure data.

## 5. Configuration

Set `EJH_PROJECT_ROOT` to a writable analysis directory and `EJH_UPSTREAM_ROOT` to legally obtained
upstream source data. Institutional scripts additionally require `EJH_INSTITUTIONAL_WORKBOOK` and
`EJH_DEIDENTIFIED_COHORT`; these restricted files are not distributed. See `scripts/README_paths.md`.

## 6. Licence

Code is released under the MIT License. Non-identifiable derived tables and figure source data are
released under CC BY 4.0. Third-party source data remain subject to their original terms and are not
redistributed.

## Archived release

The EJH manuscript reproducibility deposit is frozen as GitHub release
[`v2.0.0`](https://github.com/ABLE-Du/ETV6-RUNX1-biomarker-reanalysis/releases/tag/v2.0.0)
and archived by Zenodo:

- Version DOI: [10.5281/zenodo.22872050](https://doi.org/10.5281/zenodo.22872050)
- Concept DOI: [10.5281/zenodo.20697293](https://doi.org/10.5281/zenodo.20697293)
- Frozen code/data directory: [`ejh_manuscript_reproducibility_v2.0.0`](https://github.com/ABLE-Du/ETV6-RUNX1-biomarker-reanalysis/tree/v2.0.0/ejh_manuscript_reproducibility_v2.0.0)
