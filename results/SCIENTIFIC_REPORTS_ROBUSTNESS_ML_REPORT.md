# Scientific Reports robustness and machine-learning analysis

Date: 2026-06-11

## Prespecified reliability framework

- Multiplicity was reassessed within each biologically and statistically distinct hypothesis family using Benjamini-Hochberg FDR, Benjamini-Yekutieli FDR, and Holm family-wise error control.
- Prognostic machine learning was restricted to the TARGET karyotype cohort because it was the only available dataset with enough outcome events for a guarded analysis.
- All supervised performance estimates used nested cross-validation. Feature filtering, imputation, and penalty selection were performed within training folds.
- The single-center cohort (7 EFS events) and TARGET prognostic RNA subset (11 relapse-only patients) were not used to train prognostic machine-learning models.
- RNA subtype analysis was rerun using one primary diagnostic sample per patient. This evaluates ETV6::RUNX1 subtype biology, not adverse-outcome prediction.

## Main results

### Multiple testing

- Across the adverse-outcome screening families, no karyotype, CNA, prognostic RNA, paired-relapse RNA, DepMap dependency, or single-center Cox candidate survived BH-FDR at 0.05.
- The previous TARGET subtype-expression family retained 4285 BH-significant genes, but it included non-primary/repeated samples and therefore was rerun under a stricter primary-only rule.
- The primary-only expression comparison retained 1346 BH-significant genes, 485 BY-significant genes, and 140 Holm-significant genes.
- In the primary-only analysis, SLC7A11 had delta Z=0.663, nominal P=0.000907, BH-FDR=0.0229, BY-FDR=0.244, and Holm-adjusted P=1.

### Unsupervised learning

- The TARGET karyotype cohort included 141 patients and 24 EFS events.
- Hierarchical clustering selected k=4 by silhouette, but cluster separation was weak (silhouette=0.331).
- Cluster stability under repeated 80% feature subsampling was limited (median adjusted Rand index=0.114, IQR 0.114-0.454).
- The selected clusters were not reliably associated with events after accounting for the tested cluster numbers (permutation P=1; Holm P=1; log-rank P=0.942).
- Primary-only RNA did not form a significant global unsupervised separation by fusion status (label-permutation P=0.241).

### Supervised learning

- Repeated nested-CV three-variable clinical-threshold survival discrimination was C-index=0.566, calculated separately within each outer test fold.
- Adding +16 produced mean C-index=0.590; paired mean change=+0.024, median change=+0.000 (outer-fold range -0.292 to +0.171).
- Adding all eligible karyotype features produced mean C-index=0.521; paired mean change=-0.046, median change=-0.053 (outer-fold range -0.332 to +0.308).
- The primary-only RNA nearest-centroid classifier identified a sparse fusion-subtype signal with mean nested-CV outer-fold AUC=0.993 (repeat range 0.989-1.000). A complete nested feature-selection pipeline label-permutation test using 200 permutations yielded P=0.00498. Primary sample-type distribution was not detectably different by fusion status (Fisher P=0.725). These results cannot be interpreted as outcome prediction.

## Reliable conclusion

The combined multiplicity-controlled and leakage-resistant analyses do not identify a validated adverse-outcome biomarker for ETV6::RUNX1-positive pediatric B-ALL. The +16 finding remains the strongest low-frequency cytogenetic candidate, but it does not survive feature-wide multiplicity correction and does not provide stable out-of-sample prognostic improvement. Unsupervised karyotype groups are weak and unstable, and neither genome-wide CNA nor prognostic RNA screening yields an FDR-significant marker. The reproducible positive result is an ETV6::RUNX1-associated diagnostic expression program; it is a subtype signal and should not be presented as a relapse biomarker or therapeutic mechanism.

## Manuscript implication

For Scientific Reports, the manuscript should be framed as a reproducible, multi-source stress test of proposed adverse-outcome biomarkers. The defensible contribution is the clear separation between a reproducible subtype-expression signal and the absence of a validated prognostic signal after multiplicity correction, sensitivity analysis, and nested cross-validation.
