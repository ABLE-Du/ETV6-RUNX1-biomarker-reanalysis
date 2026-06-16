# MRD endpoint harmonization and non-pooling rule

## Prespecified rule

MRD measurements obtained at different treatment time points are clinically
distinct endpoints and must not be pooled into one common variable.

| Cohort | MRD time point | Threshold | Role in manuscript | Pooling rule |
|---|---|---:|---|---|
| TARGET | Day 29 | >=0.01% | Public derivation/internal assessment | Analyze separately |
| Single-center historical | Day 19 | >=0.01% | Exploratory association with historical EFS | Analyze separately |
| Single-center contemporary | Day 19 | >=0.01% | Early response endpoint | Analyze separately |
| Single-center contemporary | Day 46 | >=0.01% | Persistent MRD/delayed-clearance endpoint | Analyze separately |

## Derived clearance states

- `D19 cleared`: D19 MRD <0.01%.
- `D19 delayed / D46 cleared`: D19 MRD >=0.01% and D46 MRD <0.01%.
- `D46 persistent`: D46 MRD >=0.01%.

The contemporary D46-persistent group contains one patient. Results involving
this group are descriptive and cannot support inferential conclusions.

## Required manuscript wording

Use “MRD at the cohort-specific protocol time point” rather than implying that
TARGET Day-29 MRD, historical Day-19 MRD, and contemporary Day-46 MRD are
interchangeable measurements.
