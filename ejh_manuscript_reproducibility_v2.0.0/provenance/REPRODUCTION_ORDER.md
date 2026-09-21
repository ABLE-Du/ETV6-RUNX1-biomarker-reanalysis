# Reproduction order

Run from the repository root after placing the published source data in the locations the scripts expect
(see README section 1 for provenance). Every step writes into `data/derived_tables/` or
`data/figure_data/`; figure source-data workbooks and journal tables are then regenerated.

| # | Script(s) | Produces | Paper location |
|---|---|---|---|
| 1 | `scripts/01_data_acquisition/00_fetch_public_data.py`, `00c_extract_oksa_tables.py` | raw published tables, provenance log | Methods 2.1 |
| 2 | `scripts/01_data_acquisition/_extract_li_sourcedata.py`, `_export_li95.py` | Li discovery/validation matrices, COG552 posteriors | Methods 2.1 |
| 3 | `scripts/02_audit_and_harmonisation/07_endpoint_harmonisation.py` | endpoint definitions per layer | Methods 2.1, Table 1 |
| 4 | `scripts/02_audit_and_harmonisation/08_gene_harmonisation.py` | the 219-gene shared feature space (214 exact + 5 curated renames) | Methods 2.4 |
| 5 | `scripts/03_primary_concordance/22_b1_concordance.py`, `31_b1_correction.py` | ρ = 0.794, 20,000-permutation null, 10,000-replicate bootstrap CI 0.7347–0.8394, concordance 86.8% (CI 82.19–90.87%) | Methods 2.2; Results 3.2; Figure 2A,B; Table S1 |
| 6 | `scripts/04_gene_matched_benchmark/gene_matched_sensitivity.py`, `rho_difference_bootstrap.py`, `25_b6_target_benchmark.py` | 113-gene matched ρ values, Δρ 0.448 (0.272–0.632) / 0.412 (0.242–0.591), COG552 341/211 | Methods 2.3; Results 3.2, 3.3; Figure 2C |
| 7 | `scripts/05_li95_validation/li95_analysis_v2.py` | ρ = 0.818 (CI 0.647–0.921), 95.8% concordance, 24-gene subset ρ = 0.681, Li194–Li95 25-gene ρ = 0.882, candidate q values | Methods 2.5; Results 3.3; Figure 3A; Tables S2, S3 |
| 8 | `scripts/06_institutional_cohort/local_maturation_analysis.py`, `_extract_local_flow.py`, `_dump_2020_flow.py` | CD34/CD38 Hodges–Lehmann shifts, Cliff's δ, BH q; age/WBC 2×2 tables; era CMH | Methods 2.5; Results 3.3; Figure 3C,D; Tables S6–S8 |
| 9 | `scripts/07_candidates_and_drugs/26_b5_candidate_features.py`, `32_b5_tiers_and_she.py`, `24_b7_drug_sensitivity.py` | seven pre-specified candidates; 20-agent screen with BH across 20 published P values | Methods 2.4; Results 3.4, 3.5; Figure 4, Figure 5; Table S5 |
| 10 | `scripts/08_pathways/23_b2_gsea.R`, `23_b2_pathway_convergence.py` | fgsea over the available Oksa effect universe (15,928 symbols) | Methods (Supplementary); Figure S1 |
| 11 | `scripts/09_target_de/22_DE_MRD_BALL.R` | TARGET P_BALL 191 (68/123), DESeq2 ~MRD_group and protocol-adjusted model | Methods 2.1; Table 1; Figure 2C |
| 12 | `scripts/10_figure_generation/make_figure_data.py`, `render_figures_v3_5.py` | `data/figure_data/*.tsv` and the TIFF/PDF/SVG figure set | Figures 1–5 and Figure S1 |

Notes
* Step 5 fixes the resolution floor: with 20,000 permutations the smallest attainable two-sided P value is
  5.0 × 10⁻⁵, and the paper reports `P ≤ 5.0 × 10⁻⁵` (never an equality).
* Step 8: the Hodges–Lehmann and Cliff's δ intervals are percentile bootstrap intervals (20,000 resamples),
  **not** the inversion intervals of the rank test — which is why the CD38 bootstrap CI (−80.6 to −9.0)
  excludes zero while the two-sided Mann–Whitney P = 0.057.
* Step 8 also writes only aggregated results; the row-level 2020 flow sheet is deliberately not deposited.
