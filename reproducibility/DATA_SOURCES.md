# Public data sources and versions

- NOPHO article: https://www.nature.com/articles/s41375-025-02683-7
- NOPHO official code: https://github.com/systemsgenomics/ETV6-RUNX1_Genomics_treatment_response
- NOPHO code commit: `48607138ce0df2950d3ffc0f304246c38acf023f`
- NOPHO coordinate build: `GRCh37/hg19`
- Standard hg19 blacklist used for the workflow-compatible sensitivity rerun:
  https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/hg19-blacklist.v2.bed.gz
- Standard hg19 blacklist SHA-256:
  `1a4ba636f791936ab8952cb068f496ccbd55ec4753539547c5d3d055ed00642a`
- NOPHO Zenodo records:
  - https://zenodo.org/records/15167703
  - https://zenodo.org/records/15173882
  - https://zenodo.org/records/15174016
- DepMap official download API: https://depmap.org/portal/api/download/files
- DepMap complete releases used for cross-version sensitivity analysis: Public 25Q2 and Public 26Q1
- cBioPortal DataHub TARGET Phase II study: https://github.com/cBioPortal/datahub/tree/master/public/all_phase2_target_2018_pub
- cBioPortal DataHub commit: `78e9db506f1b0f029840a71712cb07dce1c51daf`

The exact local analysis-input files and SHA-256 checksums are generated in
`robustness_ml_analysis/input_file_manifest.csv`.

DepMap Public 24Q2 bulk objects were unavailable in the analysis environment;
the manuscript retains this limitation and does not claim exact 24Q2
falsification.

The exact NOPHO study-specific `hg19-blacklist.v2.mod.bed` file was not
distributed with the public records. The corrected GRCh37 analysis is therefore
described as official-workflow-compatible, not an exact reproduction.
