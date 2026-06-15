# Frozen statistical analysis plan

Version: 1.0  
Freeze date: 2026-06-15  
Status: Frozen for the final manuscript rerun. This is a retrospective freeze,
not prospective registration.

## Study objective

Assess the feasibility and cross-cohort transportability of a predefined
diagnostic cytogenetic risk score in pediatric ETV6::RUNX1-positive acute
lymphoblastic leukemia.

## Score definition

- Gain of chromosome 16: 2 points.
- Three or more abnormal karyotype tokens: 1 point.
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

- Cross-cohort transportability, measured by score evaluability and high-score
  prevalence.

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
- Complete-case analysis for each model; denominators reported for every result.
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

Any deviation from this frozen plan must be listed, dated, justified, and
reported in the supplementary methods.
