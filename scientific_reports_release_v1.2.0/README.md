# ETV6::RUNX1 biomarker reanalysis: Scientific Reports reproducibility release

Version 1.2.0 is the corrected reproducibility release supporting the
TARGET-derived candidate-score development disclosure and harmonized
cross-cohort applicability analyses.

Public repository:
https://github.com/ABLE-Du/ETV6-RUNX1-biomarker-reanalysis

## Included

- Full score-source and prespecification-status disclosure.
- One harmonized raw-karyotype parser for all cohorts.
- Matched complete-case C-index analysis with paired bootstrap resampling.
- Cohort baseline table, EFS event/censoring definitions, parser audit,
  parser sensitivity analyses, and full TARGET screening-universe disclosure.
- Aggregate, non-identifiable result tables.
- Publication figures and figure captions.
- Frozen statistical analysis plan and MRD endpoint harmonization rules.
- Reference-verification audit and final verified reference list.
- SHA-256 file manifest.
- Read-only release validation script.

## Excluded

- Individual-level single-center clinical records.
- Direct identifiers or re-identification keys.
- Third-party raw TARGET, NOPHO, or DepMap datasets.

Third-party public datasets must be downloaded from their original repositories
under their respective licenses and terms.

## Citation and DOI

The repository concept DOI is https://doi.org/10.5281/zenodo.20697293.
The v1.2.0 version DOI will be generated after GitHub-Zenodo archival.
Version 1.1.0 remains archived at https://doi.org/10.5281/zenodo.20699023 and
is superseded for the current revised manuscript package.

Run `python validate_release_package.py` before publishing the release.

## License

Repository-authored software and documentation are released under the MIT
License. Third-party datasets retain their original licenses.
