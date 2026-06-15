# ETV6::RUNX1 biomarker reanalysis: Scientific Reports reproducibility release

Version 1.0.0 is the reproducibility release supporting the cross-cohort
transportability and candidate-score feasibility analyses.

Public repository:
https://github.com/ABLE-Du/ETV6-RUNX1-biomarker-reanalysis

## Included

- Locked cytogenetic score definition.
- Analysis code for score assessment and bootstrap internal validation.
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

The repository URL can be cited immediately. A Zenodo DOI will be generated
after GitHub Release `v1.0.0` is published and archived. Until that action is
completed, the DOI status is `pending`.

Run `python validate_release_package.py` before publishing the release.

## License

Repository-authored software and documentation are released under the MIT
License. Third-party datasets retain their original licenses.
