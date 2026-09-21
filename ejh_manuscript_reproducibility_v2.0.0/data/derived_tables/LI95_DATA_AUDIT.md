# Module A — Li et al. 95-patient validation cohort: data audit

**Article.** Li Z, Yang JJ, *et al.* Nature Communications 2025;16:1153.
DOI 10.1038/s41467-025-56229-7 · PMCID PMC11779914.

Full provenance (URL, bytes, SHA-256, download date) is in
`LI95_SOURCE_PROVENANCE.tsv`.

---

## 1. What was retrieved

The publisher supplementary bundle was obtained from Europe PMC
(`/PMC11779914/supplementaryFiles`). Three files were relevant:

| file | bytes | used |
|---|---|---|
| `41467_2025_56229_MOESM6_ESM.xlsx` | 2,505,150 | **yes** |
| `41467_2025_56229_MOESM4_ESM.xlsx` | 15,950 | no (legend/index) |
| `41467_2025_56229_MOESM5_ESM.pdf` | — | no (prose) |

Sheet **`Supp Figure 10`** of MOESM6 is a patient-level validation-cohort matrix:

- row 1 → 95 patient identifiers
- row 2 → C1 / C2 PAM assignment for each patient
- rows 4+ → 50 gene rows, expression values as supplied by the authors

No identifier of any kind (name, record number, date of birth, exact date) is
present in the deposited file.

---

## 2. Eligibility gate

| flag | value |
|---|---|
| `LI95_LABELS_AVAILABLE` | **TRUE** (C1 n=61, C2 n=34, recomputed) |
| `LI95_EXPRESSION_AVAILABLE` | **TRUE** (50 genes × 95 patients, 0 % missing) |
| `LI95_CLINICAL_AVAILABLE` | PARTIAL (published contingency tables: age, WBC) |
| `LI95_GENE_EFFECT_ANALYSIS_PERMITTED` | **TRUE** |
| `LI95_MOLECULAR_TRIANGULATION` | **TRUE** |

Because both patient-level labels and patient-level expression are present,
gene-level re-analysis is permitted. We used the published C1/C2 assignments
as given; we did **not** re-cluster, did **not** re-select genes, and did
**not** change any threshold.

Nothing in this module was reconstructed from a figure, a screenshot or a
digitised plot.

---

## 3. Gene universe actually available

The article states that **122 genes** were independently selected on the
validation cohort using non-negative matrix factorisation. The deposited
workbook contains **50** of them. We therefore used the 50 deposited genes and
report the size of the universe at every point a number is quoted.

Overlap with our other vectors:

| intersection | n |
|---|---|
| 50 deposited validation genes ∩ Oksa EOI coefficient table | 48 (2 genes absent from Oksa table) |
| deposited genes ∩ the frozen primary 219-gene set | 24 |
| deposited genes ∩ Oksa ∩ Li194 discovery | 24 |
| deposited genes ∩ Li194 discovery (any Oksa value) | 25 |

---

## 4. Results

### 4.1 Oksa ↔ Li95 (our re-analysis)

| gene universe | n | Spearman ρ | bootstrap 95 % CI | direction concordance | gene-identity permutation P |
|---|---|---|---|---|---|
| all deposited genes with an Oksa coefficient | 48 | **0.818** | 0.647 – 0.921 | 95.8 % (89.6 – 100) | 5.0 × 10⁻⁵ |
| restricted to the 24 genes in the frozen primary 219 | 24 | **0.681** | 0.322 – 0.895 | 91.7 % (79.2 – 100) | 6.0 × 10⁻⁴ |

Null: gene identity permuted 20,000 times (null ρ = 0.0003 ± 0.145 for the
48-gene set). Bootstrap resamples genes, so the interval describes stability
across the gene set, not patient sampling. 5.0 × 10⁻⁵ is the floor attainable
with 20,000 permutations and should be read as "no permutation exceeded the
observed statistic".

### 4.2 Li194 discovery ↔ Li95 (our re-analysis)

25 genes present in both: ρ = **0.882** (bootstrap CI 0.653 – 0.976), direction
concordance **100 %** (23/23 informative), permutation P = 5.0 × 10⁻⁵.

### 4.3 Three-way direction stability

For the 24 genes with an Oksa coefficient, a Li194 value and a Li95 value, the
sign was identical in all three contrasts for **91.7 %** (22/24).

### 4.4 Tier-1 candidates in Li95

Four of the seven frozen Tier-1 genes are present in the deposited validation
panel. All four reproduce the Oksa direction, and all four are significant in
our re-analysis after Benjamini–Hochberg correction within the 50-gene panel:

| gene | Oksa EOI coefficient | Oksa FDR | Li194 effect | **Li95 effect (C1 − C2)** | Li95 BH q | direction match |
|---|---|---|---|---|---|---|
| GATA2 | +2.329 | 0.0170 | +1.486 | **+2.939** | 2.7 × 10⁻¹⁰ | yes |
| ANPEP | +1.757 | 0.0280 | +1.196 | **+1.938** | 1.2 × 10⁻⁸ | yes |
| MPV17L | +2.506 | 0.0273 | +1.277 | **+1.998** | 1.1 × 10⁻⁶ | yes |
| SHE | +3.561 | 0.0048 | +1.323 | **+2.471** | 5.8 × 10⁻⁸ | yes |

MID1, SERPINI2 and UCK2 are **not** in the deposited 50-gene validation panel.
Their status in Li95 is `NOT_ASSESSABLE`, not `NOT_SUPPORTED`.

**CD38** (not a Tier-1 gene; pre-specified maturation marker) is present:
Li95 effect −1.501, BH q = 4.2 × 10⁻⁹ — lower in C1, i.e. lower in the adverse
state, matching the direction we observe for CD38 surface protein locally.

### 4.5 Clinical phenotype (re-tabulated from the published Source Data table)

| feature | C1 | C2 | OR (C1 vs C2) | 95 % CI | Fisher P | published (chi-squared) |
|---|---|---|---|---|---|---|
| age ≥ 5 y | 33 / 61 (54 %) | 10 / 34 (29 %) | 2.83 | 1.16 – 6.91 | 0.031 | P = 0.020; "54 % vs 29 %" |
| WBC ≥ 50 × 10⁹/L | 7 / 60 (12 %) | 11 / 34 (32 %) | 0.28 | 0.10 – 0.80 | 0.027 | "12 % vs 32 %" |

These are **Li et al.'s** validation-cohort results re-tabulated from their
deposited contingency table; the small difference in P reflects Fisher's exact
test here versus the chi-squared test in the source paper.

### 4.6 Published validation architecture (not our analysis)

> Li et al. reported: PAM assigned 61 patients (64 %) to C1 and 34 (36 %) to
> C2, and independent gene selection plus UMAP on the same 95 patients produced
> two sub-clusters consistent with the PAM classes for 87 of 95 patients (92 %).

---

## 5. COG552

The deposited COG552 source data contain **only** UMAP coordinates and PAM
posterior probabilities for 552 patients. There is no expression matrix and no
clinical variable. Re-tabulating the posteriors gives C1 n=341, C2 n=211.

Per the no-fabrication rule, COG552 is used exclusively as
`PUBLISHED_INDEPENDENT_CONTEXT`. No molecular re-analysis was attempted and no
new controlled-access acquisition was started.

---

## 6. Attribution rule for the manuscript

- Numbers in §4.1–4.4 and the Fisher tests in §4.5 are **our re-analysis of the
  published source data** → write "In our re-analysis of the deposited source
  data…".
- The 87/95 agreement and the 122-gene selection are **published results read
  from the paper** → write "Li et al. reported…".
- The two are never blurred into one another.
