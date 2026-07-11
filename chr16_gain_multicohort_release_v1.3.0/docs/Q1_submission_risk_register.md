# Q1 Submission Risk Register

## Current scientific position

The manuscript is written as an exploratory integrative reanalysis. Its defensible conclusion is that chromosome 16 gain is a recurrent low-frequency event in ETV6::RUNX1-positive childhood ALL, with a candidate adverse EFS association in TARGET and a RAS-MAPK transcriptional context in GSE181157. It does not claim independent outcome validation, MRD occultness, RAS-pathway mutual exclusivity, or a treatment recommendation.

## Items to resolve before submission

| Priority | Item | Why it matters | Required resolution |
| --- | --- | --- | --- |
| Critical | Freeze the current script and aggregate output tables in the public repository | The existing repository DOI predates the current expanded gseapy and manuscript build scripts | Upload this package or its reproducible analysis subset, create a tagged GitHub release, and archive that release in Zenodo before submission |
| Critical | Select a specific journal | Q1 journals have non-interchangeable limits, reference style, declaration wording, and figure-file requirements | Select the target journal before final formatting and cover-letter generation |
| High | Obtain an independent EFS-annotated ETV6::RUNX1 cohort, if feasible | TARGET is currently the only EFS discovery set, and the +16 screen-level BH q value is 0.308 | Present this work as exploratory if no independent outcome cohort can be added; do not label +16 as validated |
| High | Have a cytogeneticist verify the seven TARGET +16 reports and all GSE181157 text-derived chr16 calls | The primary TARGET exposure and GSE181157 mechanism layer rely on text-derived cytogenetic fields | Perform a blinded parser audit and append the adjudication record to the reproducibility package |
| Moderate | Confirm corresponding-author details, affiliations, author order, ethics identifiers, and funding | These are source-provided administrative fields rather than reanalysis outputs | Obtain formal author approval before upload |
| Moderate | Confirm that the target journal accepts an exploratory public-data reanalysis without independent outcome validation | Some Q1 hematology journals will desk-reject this evidence level | Use a scope-compatible target and avoid submitting to a journal that requires mechanistic experimentation or external prognostic validation |

## Recommended journal positioning

The manuscript should be positioned as an integrative genomic epidemiology/reanalysis paper, not as a definitive clinical biomarker study. Its strongest editorial features are transparent cohort-role separation, explicit multiplicity handling, cross-platform CNV recurrence, and a restrained mechanistic interpretation. Its main barrier is the lack of external EFS validation, which should be declared in the abstract, discussion, and cover letter rather than concealed.
