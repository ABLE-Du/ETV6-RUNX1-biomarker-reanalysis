# Release v1.2.0

## Purpose

This corrective release expands the harmonized raw-karyotype reanalysis with
reviewer-requested reporting and reproducibility materials.

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
- Adds cohort baseline characteristics, EFS event/censoring definitions,
  parser public-example audit, parser sensitivity analyses, and the complete
  39-feature TARGET karyotype screening universe.
- Adds journal-formatted supplementary tables with rounded display estimates
  and explicit notation for non-estimable or infinite odds ratios.
- Removes the obsolete sample-size scenario table from the current submission
  package.
- Corrects MRD text parsing so values such as `<0.01` are treated as
  below-threshold measurements rather than missing.

## Privacy

No individual-level single-center patient data are included.

## Archiving

GitHub-Zenodo integration should archive this as a new version under concept
DOI `10.5281/zenodo.20697293`. Releases v1.0.0 and v1.1.0 remain immutable
but are superseded for the current revised manuscript package.
