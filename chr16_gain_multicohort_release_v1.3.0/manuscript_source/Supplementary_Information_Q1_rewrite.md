# Supplementary Information

## Supplementary Methods

### TARGET karyotype audit and candidate-screen universe

The TARGET source table was restricted to ETV6::RUNX1-positive patients. The diagnostic karyotype was considered evaluable when the field contained text and did not contain a documented non-assessment, failed assessment, no metaphases, or unknown result. For the primary +16 variable, whitespace was removed and the expression `(?:^|,)\\+16(?=,|/|\\[|$)` was applied. This definition identifies an explicit +16 token, avoids interpreting t(12;21) as a chromosome-16 gain, and counts each patient once even when more than one clone is recorded. The primary result is consequently a conventional-karyotype feature, not a SNP-array or sequencing-derived copy-number call.

All retained karyotype candidate features in the 141-patient/24-event EFS analysis set are supplied in Supplementary Table S1. The table retains the full screening universe, prevalence, univariable Cox estimate, and BH q value. The +16 univariable P value was 0.0088 and the corresponding BH q value was 0.3077. The q value is the appropriate result for the screen-level inference.

### Sparse-data EFS sensitivity analysis

The unpenalized Cox model is reported in the main text as the principal descriptive effect estimate because no L2 penalty was specified before inspection. Four ridge penalties were then applied to the identical 141-patient/24-event dataset. The sensitivity table reports gain16 hazard ratios, Wald confidence intervals, nominal P values, C-index, warnings, and errors. Neither an adjusted model nor a single ridge penalty was selected as confirmatory. The exposure group has seven patients and four EFS events, so all regression estimates should be read together with the Kaplan-Meier curve and event counts.

### CNV recurrence definitions

For Oksa/NOPHO CNVkit segments, broad/whole chromosome 16 gain required at least 80% of the observed chromosome 16 panel span with log2 ratio at least 0.3 and weighted chromosome 16 log2 ratio at least 0.3. For GSE184692, chromosome 16 probes were identified from GPL10150 annotation. The primary aCGH definition required mean chromosome 16 segmented log2 at least 0.20 and at least 50% of chr16 probes with segmented log2 at least 0.25. The sensitive rule used thresholds of 0.15 and 0.20, and the strict rule used thresholds of 0.25 and 0.30. All three GSE184692 rules identified four of 136 ETV6-RUNX1 samples. Exact binomial confidence intervals were calculated using the beta method.

### Expression and gene-set audit

Each expression dataset was analyzed separately. GSE181157, GSE227832, and GSE228632 RNA-seq counts were transformed to log2(counts per million + 1) and filtered at log2 counts per million at least 1 in at least 5% of samples. GSE87070 GPL570 probes were mapped to gene symbols and summarized by median expression. Gene symbols and Ensembl identifiers were mapped with the cached BioMart mappings used by the analysis script. The final GSE181157 gene-set sizes were 584 chr16 coding genes, 96 genes in the 16p/MAPK3 region, 43 RAS-MAPK core genes, 21 ERK-feedback genes, and 14-16 genes per B-cell developmental-state set.

Preranked GSEA used 1,000 phenotype permutations, fixed seed 20260710, minimum size 5, maximum size 3000, and Benjamini-Hochberg adjustment within each comparison over the eight curated gene sets. GSE181157 mutation and cytogenetic annotations were parsed from the supplied text metadata. Therefore, the data support transcriptional context, not a clinically adjudicated CNV-mutation interaction model.

### Supplementary figure legend

**Supplementary Figure S1. Sensitivity of the exploratory TARGET +16 EFS association to Cox penalization.** All estimates were calculated in the same 141-patient/24-event TARGET karyotype-EFS set. The unpenalized estimate and four L2 penalties are shown. The plot illustrates model sensitivity and is not a penalty-selection exercise.

![Supplementary Figure S1](figures/Supplementary_Figure_S1_penalty_sensitivity.png)

## Supplementary Tables

**Supplementary Table S1.** Complete TARGET diagnostic-karyotype candidate-screen universe with univariable and adjusted Cox results and BH q values. File: `Supplementary_Table_S1_TARGET_karyotype_screen.csv`.

**Supplementary Table S2.** TARGET +16 Cox penalization sensitivity analysis. File: `Supplementary_Table_S2_TARGET_penalty_sensitivity.csv`.

**Supplementary Table S3.** Platform-separated chr16-gain prevalence and exact confidence intervals. File: `Supplementary_Table_S3_chr16_prevalence.csv`.

**Supplementary Table S4.** Multi-cohort preranked GSEA output for the eight curated gene sets. File: `Supplementary_Table_S4_multicohort_prerank_GSEA.csv`.

**Supplementary Table S5.** Expression cohort metadata audit and permitted inference by dataset. File: `Supplementary_Table_S5_expression_dataset_audit.csv`.
