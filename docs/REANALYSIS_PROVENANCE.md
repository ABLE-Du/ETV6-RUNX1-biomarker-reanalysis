# Reanalysis provenance

Run date: 2026-06-09

## Public sources

- NOPHO article: https://www.nature.com/articles/s41375-025-02683-7
- NOPHO official code: https://github.com/systemsgenomics/ETV6-RUNX1_Genomics_treatment_response
- NOPHO code commit: `48607138ce0df2950d3ffc0f304246c38acf023f`
- NOPHO Zenodo records:
  - https://zenodo.org/records/15167703
  - https://zenodo.org/records/15173882
  - https://zenodo.org/records/15174016
- DepMap official download API: https://depmap.org/portal/api/download/files
- cBioPortal DataHub study: https://github.com/cBioPortal/datahub/tree/master/public/all_phase2_target_2018_pub
- cBioPortal DataHub commit: `78e9db506f1b0f029840a71712cb07dce1c51daf`
- Enabled bioinformatics skills: https://github.com/GPTomics/bioSkills

## Known source limitation

The manuscript-specified DepMap Public 24Q2 Figshare objects returned HTTP 403 in
this environment. DepMap sensitivity analyses therefore use official complete
25Q2 and 26Q1 bulk matrices. This limitation is explicitly retained in the
report and does not justify claiming exact 24Q2 falsification.

## Scripts

- `download_nopho_zenodo.py`
- `analyze_nopho_cnv.py`
- `download_depmap_bulk_file.py`
- `analyze_depmap_sensitivity.py`
- `analyze_target_public.py`

