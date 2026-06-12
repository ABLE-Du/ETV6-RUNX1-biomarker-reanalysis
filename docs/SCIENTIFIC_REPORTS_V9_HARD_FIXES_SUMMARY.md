# Scientific Reports v9 hard-fix summary

## NOPHO rerun

- Reprocessed 262 public CNVkit `.cns` files in GRCh37/hg19 coordinates.
- Used the official thresholds and event rules, including retention of events
  with size `>=1 Mb`.
- Used the standard public hg19 blacklist because the exact study-specific
  modified blacklist was not publicly distributed.
- Recorded a relative-path input manifest, SHA-256 checksums, official-code
  commit, gene-coordinate source, and all workflow parameters.

Primary workflow-compatible results:

- chromosome 6 deletion: 62/262
- FOXO3 overlap: 42/262
- TNFAIP3 overlap: 26/262
- FOXO3 and TNFAIP3 co-deletion: 25/262
- FOXO3-only among nine evaluated 6q genes: 5/262
- CDKN2A-del6q co-occurrence: `P=0.0244`, `BH-FDR=0.0733`

## Statistical and manuscript corrections

- Expanded the operational multiple-testing register to include DepMap 25Q2,
  NOPHO co-occurrence, the non-duplicate karyotype-pattern screen, and the
  global multi-hit sensitivity test.
- Reported adjusted +16 `BH-FDR=0.186`.
- Reframed +16 nested cross-validation as post-selection internal sensitivity
  analysis because candidate discovery occurred before cross-validation.
- Added a sensitivity analysis excluding the source-recorded 0-month infection
  death; the favorable single-center outcome interpretation was unchanged.
- Corrected references 6 and 7 and removed unsupported published-status wording
  for earlier internal results.
- Removed the legacy NOPHO breakpoint/MCR, nominal DepMap dependency,
  virtual-knockout, and clinical-treatment framework figures from the formal
  submission package.

## Remaining author-supplied blockers

- ethics committee name, approval number, and consent/waiver statement
- funding source and grant numbers
- corresponding-author email
- CRediT author-contribution statement
- confirmation of treatment-protocol labels and follow-up cutoff date
- author review and approval of the generative-AI disclosure
