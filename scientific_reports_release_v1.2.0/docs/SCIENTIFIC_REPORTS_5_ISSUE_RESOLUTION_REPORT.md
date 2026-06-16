# Five reviewer-identified issues: resolution report

## 1. Cohort definitions, EFS definitions, follow-up, and event composition

Resolved.

- Added `Table_3_cohort_characteristics.csv` with cohort size, diagnosis
  years, treatment protocol/era, sex, age, WBC, risk group, karyotype
  evaluability, MRD availability, EFS analyzability/events, follow-up, and
  high-score prevalence.
- Added `Supplementary_Table_EFS_definitions.csv` with time origin, time
  scale, EFS event definitions, censoring rules, event composition, and
  analysis populations.
- Updated Results and Methods to report overall EFS events, matched-model
  events, median follow-up, and cohort-specific event/censoring rules.

## 2. Parser audit and sensitivity analyses

Resolved.

- Added `Supplementary_Table_parser_audit_public_examples.csv` with
  representative public TARGET karyotype examples, normalized strings,
  detected tokens, and component calls.
- Added `Supplementary_Table_parser_audit_summary.csv` with cohort-level
  token-count summaries and classification changes under alternative parsing
  rules.
- Added `Supplementary_Table_parser_sensitivity.csv` testing two alternatives:
  excluding the defining `t(12;21)` token and counting unique reported lesions
  rather than repeated report tokens.
- Main conclusion was unchanged: historical high-score prevalence remained
  2/59 (3.4%) and no historical EFS association emerged.

## 3. TARGET screening universe disclosure

Resolved.

- Added `Supplementary_Table_TARGET_karyotype_screening_universe.csv`, the
  complete 39-feature TARGET karyotype screen.
- Added `Supplementary_Table_TARGET_screening_universe_summary.csv`, grouping
  screened features into whole-chromosome numerical, deletion, structural
  event type, structural chromosome involvement, abnormality burden, specific
  translocation, and karyotype-pattern categories.
- Updated the manuscript to state that del(6q) was retained as a biologically
  motivated, routinely reportable component despite nonsignificant TARGET EFS
  evidence.

## 4. Single-center data availability

Resolved.

- Updated Data availability to state that individual-level single-center data
  cannot be publicly released.
- Added a controlled-access route: requests should go through the
  corresponding author and require institutional review, approval of proposed
  use, ethics approval where required, and a data-use agreement.
- Clarified that requests not satisfying governance requirements are limited
  to public aggregate tables and code.

## 5. Old sample-size table

Resolved.

- Removed `Supplementary_Table_sample_size_plan.csv` from the current
  submission package because it was not used in the revised manuscript.
- Removed obsolete old-numbered result tables from the current package.
- Rebuilt the Word manuscript, figures, validation report, and SHA-256
  reproducibility manifest after removal.

## Additional data-management correction

The MRD parser now treats text values such as `<0.01` as below-threshold
measurements rather than missing. This increased contemporary D19
score-plus-MRD availability from 32 to 34 patients and did not change the
study conclusions.
