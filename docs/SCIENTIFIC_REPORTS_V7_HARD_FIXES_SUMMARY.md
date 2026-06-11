# Scientific Reports v7 hard-fix summary

Output manuscript: `_ETV6-RUNX1_manuscript_v7_hard_fixes.docx`

## Hard issue 1: survival nested cross-validation

- Corrected inner penalty selection to average concordance calculated separately in each inner validation fold.
- Corrected outer evaluation to calculate concordance separately in each outer test fold before aggregation.
- Updated manuscript, Figure 5, Figure 6, Table 5, report, and machine-readable outputs.
- Corrected estimates:
  - Three-variable clinical-threshold model: mean outer-fold C-index `0.566`.
  - Clinical plus +16: mean outer-fold C-index `0.590`.
  - +16 paired increment: mean `+0.024`, median `0.000`, range `-0.292` to `+0.171`.
  - Clinical plus eligible karyotype features: mean outer-fold C-index `0.521`.

## Hard issue 2: RNA classifier permutation test

- Replaced the fixed 50-gene permutation test with a complete nested-pipeline label-permutation test.
- Each permutation reruns inner feature-number selection and outer-fold classification.
- AUC is calculated separately in each outer test fold before aggregation.
- Corrected result: mean repeated nested-CV outer-fold AUC `0.993`; complete nested-pipeline permutation `P=0.00498` using 200 permutations.
- The manuscript continues to classify this as an internal fusion-subtype signal, not prognosis.

## Hard issue 3: non-ETV6-selected contemporary registry

- Removed the 42-patient contemporary registry from the abstract, main-text results, main Figure 3, primary Table 4, and main-text limitations.
- Replaced Figure 3 with a historical-only ETV6::RUNX1-positive outcome-cohort figure.
- Retained the registry only as Supplementary Figure S6, explicitly labeled as non-ETV6::RUNX1-selected and non-validating.

## Hard issue 5: reproducibility

- Removed machine-specific paths from the corrected analysis, figure-generation, and manuscript-integration scripts.
- Added explicit `--source-root`, `--workspace-root`, and `--permutations` parameters.
- Added a one-command PowerShell reproducibility entry point.
- Generated SHA-256 input-file manifest and software-environment record.
- Added public-data source versions and commit identifiers.
- A public GitHub release and Zenodo DOI still require external publication; the manuscript retains an explicit pre-submission DOI/URL placeholder rather than making a false availability claim.

## Verification

- Full analysis completed using 200 complete nested-pipeline RNA permutations.
- Corrected figures visually inspected.
- Old Figure 3 image is absent from v7; the historical-only Figure 3 is embedded.
- Obsolete survival values and fixed-50-gene permutation wording are absent.
- Contemporary registry wording occurs only in Supplementary Figure S6.
- DOCX ZIP integrity check passed.
