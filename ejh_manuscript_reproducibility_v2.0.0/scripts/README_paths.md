# Paths used by the scripts

The public-release scripts use environment variables rather than author-machine paths:

| Variable | Purpose | Default |
|---|---|---|
| `EJH_PROJECT_ROOT` | writable analysis/output root | repository root |
| `EJH_UPSTREAM_ROOT` | legally obtained public/controlled source data | `<root>/upstream` |
| `EJH_LI_ROOT` | processed Li source-data directory | `<upstream>/li` |
| `EJH_INSTITUTIONAL_WORKBOOK` | restricted institutional source workbook | required, no default |
| `EJH_DEIDENTIFIED_COHORT` | restricted de-identified institutional table | required, no default |

The two institutional variables are intentionally unresolved in a public checkout. The corresponding
scripts document the transformation, but successful execution requires authorized local access.
