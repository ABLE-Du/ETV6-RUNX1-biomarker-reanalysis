# ETV6::RUNX1 adverse-outcome biomarker reanalysis

Reproducible, multiplicity-controlled reassessment of candidate adverse-outcome
biomarkers in ETV6::RUNX1-positive pediatric B-ALL.

## Main conclusion

The integrated analysis does not establish a validated adverse-outcome
biomarker. Chromosome 16 gain remains a low-frequency external-validation
candidate, but it does not survive feature-wide FDR correction and does not
provide stable out-of-sample prognostic improvement. The robust RNA result is
an ETV6::RUNX1 fusion-subtype expression signal, not an outcome-prediction
signal.

## Corrected machine-learning results

- Three-variable clinical-threshold Cox model: mean outer-fold C-index `0.566`.
- Clinical model plus chromosome 16 gain: mean outer-fold C-index `0.590`.
- Paired +16 increment: mean `+0.024`, median `0.000`, range `-0.292` to `+0.171`.
- Primary-only RNA classifier: mean nested-CV outer-fold AUC `0.993`.
- Complete nested-pipeline RNA label-permutation test: `P=0.00498` using 200 permutations.

## Repository contents

- `scripts/`: corrected analysis, figure-generation, and manuscript-integration code.
- `results/tables/`: public-data-derived statistical and machine-learning outputs.
- `results/figures/`: public-data robustness figures in PNG, PDF, and SVG formats.
- `single_center_aggregate/`: aggregate-only single-center statistics.
- `reproducibility/`: source versions, execution instructions, and one-command runner.
- `docs/`: audit, provenance, and hard-fix summaries.

The Scientific Reports v8 manuscript-finalization changes are summarized in
[`docs/SCIENTIFIC_REPORTS_V8_IMPORTANT_POINTS_SUMMARY.md`](docs/SCIENTIFIC_REPORTS_V8_IMPORTANT_POINTS_SUMMARY.md).
The repository excludes manuscript DOCX files and individual-level
single-center data by design.

## Data privacy

This repository intentionally excludes:

- original single-center clinical workbooks;
- patient-level or pseudonymized single-center clinical tables;
- manuscript DOCX files;
- direct identifiers or protected health information.

The single-center files in `single_center_aggregate/` contain aggregate
statistics only. Public-data sample identifiers and derived statistics retain
the terms of their original sources.

## Reproduce

Install Python 3.11 and:

```powershell
pip install -r requirements.txt
.\reproducibility\run_hard_fixes.ps1 -SourceRoot "PATH_TO_PUBLIC_SOURCE_ROOT"
```

The public source root must contain the downloaded public datasets described in
[`reproducibility/DATA_SOURCES.md`](reproducibility/DATA_SOURCES.md).

## Citation and limitations

The analysis is exploratory and internally validated. It does not support
clinical risk reassignment or targeted-therapy selection. See
[`results/SCIENTIFIC_REPORTS_ROBUSTNESS_ML_REPORT.md`](results/SCIENTIFIC_REPORTS_ROBUSTNESS_ML_REPORT.md)
for the full interpretation.
