# Release v1.1.0

## Purpose

This corrective release replaces the initial score analysis with a harmonized
raw-karyotype reanalysis and matched-sample incremental-performance assessment.

## Main additions

- Discloses that score construction was outcome-informed in TARGET and locked
  only before single-center assessment.
- Applies one auditable report-token definition to raw karyotypes in all
  cohorts; the variable is not labelled standard complex karyotype.
- Recomputes MRD-only and MRD-plus-score C-indices in identical complete-case
  samples using 1,000 paired bootstrap resamples.
- Reports the negative external result: no historical-cohort EFS association
  or incremental discrimination after harmonization.
- Adds component prevalence, exact confidence intervals, penalty sensitivity,
  and historical-definition agreement tables.

## Privacy

No individual-level single-center patient data are included.

## Archiving

GitHub-Zenodo integration should archive this as a new version under concept
DOI `10.5281/zenodo.20697293`. Release v1.0.0 remains immutable.
