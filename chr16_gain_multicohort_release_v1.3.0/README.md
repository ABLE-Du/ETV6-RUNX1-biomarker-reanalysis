# ETV6::RUNX1 chromosome 16 gain multi-cohort reanalysis

This v1.3.0 release freezes public code, aggregate non-identifiable results,
figures, and manuscript source for an exploratory multi-cohort reanalysis of
chromosome 16 gain in pediatric ETV6::RUNX1-positive acute lymphoblastic leukemia.

## Evidence boundary

- TARGET is the sole EFS discovery layer. The +16 signal is exploratory: seven
  exposed patients and a karyotype-screen BH q value of 0.308.
- Oksa/NOPHO and GSE184692 support CNV recurrence/frequency only; they are not
  pooled with TARGET for EFS estimation.
- GSE181157 supports paired chr16/RAS expression context. GSE227832,
  GSE228632, and GSE87070 provide external subtype-state context only.

## Included

- Expanded gseapy, CNV-frequency, and manuscript-generation scripts.
- Aggregate non-identifiable result tables, figures, manuscript source, and
  reference-verification audit.
- SHA-256 manifest, release notes, source-data provenance, and privacy QA.

## Excluded

- Individual-level single-center records, source workbooks, identifiers, and
  cached third-party raw TARGET, GEO, and NOPHO inputs.

Source datasets must be downloaded from their original repositories under their
respective licenses. Repository-authored code and documentation use the MIT license.
