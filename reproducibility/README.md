# Reproducibility package

This package documents the corrected Scientific Reports hard-fix analysis.
The public repository reproduces the robustness analysis from downloaded public
source data. Manuscript integration additionally requires the private local
manuscript source and is therefore optional.

## Reliability corrections

- Survival nested cross-validation calculates concordance separately in every outer test fold before aggregation.
- RNA subtype classification calculates AUC separately in every outer test fold.
- The RNA label-permutation test reruns the complete nested feature-number selection and classification pipeline.
- Analysis scripts accept explicit source and workspace roots; no machine-specific path is required.
- The analysis writes SHA-256 input checksums and the software environment to `robustness_ml_analysis/`.

## Required inputs

The source root must contain the existing `public_reanalysis/` directory. The
workspace root must contain the single-center aggregate analysis outputs and the
manuscript source. Individual-level single-center data are not required for the
corrected robustness analysis.

## Run

```powershell
.\reproducibility\run_hard_fixes.ps1 -SourceRoot "PATH_TO_SOURCE_ROOT"
```

The default permutation count is 200. For a smoke test only:

```powershell
.\reproducibility\run_hard_fixes.ps1 -SourceRoot "PATH_TO_SOURCE_ROOT" -Permutations 5
```

## Primary outputs

- `robustness_ml_analysis/analysis_summary.json`
- `robustness_ml_analysis/input_file_manifest.csv`
- `robustness_ml_analysis/software_environment.json`
- `robustness_ml_analysis/tables/target_survival_nested_cv.csv`
- `robustness_ml_analysis/tables/primary_expression_nested_cv_outer_folds.csv`
- `robustness_ml_analysis/tables/primary_expression_classifier_permutation.csv`
- Corrected manuscript output when the optional `-BuildManuscript` switch is
  used with the required private local manuscript source.

## Public repository

Before submission, publish this project as a versioned GitHub release and
archive that release in Zenodo. Replace the manuscript placeholder with the
resulting repository URL and DOI.
