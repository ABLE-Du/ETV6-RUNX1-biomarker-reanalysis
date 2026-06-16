# Frozen statistical analysis plan

Version: 1.1  
Freeze date: 2026-06-15  
Status: Version 1.0 was a retrospective freeze, not prospective registration.
Version 1.1 documents reviewer-driven corrections without changing score
components, weights, or threshold.

## Study objective

Assess the cross-cohort applicability and exploratory prognostic performance
of a TARGET-derived, subsequently locked diagnostic cytogenetic candidate
score in pediatric ETV6::RUNX1-positive acute lymphoblastic leukemia.

## Score definition

- Gain of chromosome 16: 2 points.
- Three or more abnormal karyotype report tokens: 1 point.
- del(6q): 1 point.
- Continuous score: 0-4 points.
- High score: total score >=1.

Score components, weights, and threshold must not be modified after this
freeze.

## Cohort roles

- TARGET: score derivation and internal feasibility assessment; never described
  as independent validation.
- Single-center historical cohort: exploratory external feasibility assessment
  for EFS.
- Single-center contemporary cohort: descriptive assessment of D19-D46 MRD
  clearance kinetics.
- NOPHO/Zenodo: external biological replication of selected CNV frequencies;
  no EFS analysis.

## Endpoints

Primary manuscript endpoint:

- Cross-cohort applicability and ascertainment heterogeneity, measured by
  score evaluability, component prevalence, and high-score prevalence.

Secondary endpoints:

- EFS association of the continuous score and high-score category.
- Cohort-specific MRD associations analyzed separately.
- Apparent and optimism-corrected C-index for MRD-only and MRD-plus-score
  penalized Cox models.

## Statistical methods

- Fisher's exact test for binary score-MRD associations.
- Penalized Cox regression for EFS associations.
- Effect estimates reported with 95% confidence intervals and exact P values.
- Five hundred bootstrap resamples for optimism correction.
- Matched complete-case analysis for each MRD-only versus MRD-plus-score
  comparison; denominators reported for every result.
- No stepwise selection, data-driven score reweighting, or pooling of different
  MRD time points.
- No interpretation of calibration curves because event counts are inadequate.

## Interpretation thresholds

- All analyses are exploratory.
- P<0.05 is not sufficient to claim clinical validation.
- TARGET results are internal derivation evidence.
- Single-center results with fewer than 20 events are considered severely
  underpowered.
- Simulated data and sample-size calculations are study-design evidence only.

## Deviations

- 2026-06-15: the historical precomputed `complex_karyotype` variable was
  replaced by the same raw-string report-token parser used for TARGET and the
  contemporary cohort. This correction was required because the old variable
  was not operationally equivalent.
- 2026-06-15: bootstrap resamples increased from 500 to 1,000 and both
  incremental models were restricted to the same complete-case sample.
- 2026-06-15: “transportability” was replaced by “applicability and
  ascertainment heterogeneity” because prevalence differences alone do not
  establish prognostic transportability.
- These corrections did not alter score components, weights, or threshold.
