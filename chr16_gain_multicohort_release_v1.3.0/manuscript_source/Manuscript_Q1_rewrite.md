# Chromosome 16 gain in ETV6::RUNX1-positive childhood acute lymphoblastic leukemia: an integrative multi-cohort reanalysis

**Article type:** Original Article

**ChengKan Du, JiaShi Zhu, SiJian Wang, Yue Zheng, Min Liu, Hong Li**

Department of Laboratory Medicine and Department of Hematology, Shanghai Children's Hospital, School of Medicine, Shanghai Jiao Tong University, Shanghai, China

**Corresponding author:** Hong Li, MD, Department of Hematology, Shanghai Children's Hospital, School of Medicine, Shanghai Jiao Tong University, No. 24, Lane 1400, Beijing West Road, Jing'an District, Shanghai 200040, China. Email: lihonglily1978@sina.com. ORCID: 0009-0004-8324-0639.

**Running title:** Chromosome 16 gain in ETV6::RUNX1 ALL

**Keywords:** ETV6::RUNX1; childhood acute lymphoblastic leukemia; chromosome 16 gain; copy-number alteration; RAS-MAPK; event-free survival

## Abstract

ETV6::RUNX1-positive childhood acute lymphoblastic leukemia (ALL) has favorable average outcomes, but a minority of patients relapse and the contribution of uncommon secondary cytogenetic lesions remains uncertain. We performed an exploratory, multi-cohort reanalysis designed to keep prognostic discovery, copy-number recurrence, and transcriptional context analytically separate. In TARGET, 222 ETV6::RUNX1-positive patients were identified; 141 had both evaluable diagnostic karyotype and event-free survival (EFS) data. A karyotype-derived +16 call was present in 7 of 144 evaluable patients (4.9%; 95% confidence interval [CI], 2.0-9.8). EFS events occurred in 4 of 7 patients with +16 and 20 of 134 without +16; five-year EFS was 71.4% and 85.5%, respectively (log-rank P = 0.004). The univariable Cox estimate was nominally adverse (hazard ratio, 4.20; 95% CI, 1.44-12.31; P = 0.009), but the association was not retained after Benjamini-Hochberg correction across the karyotype candidate screen (q = 0.308). External CNV data supported recurrence, not outcome validation: broad or whole-chromosome 16 gain was observed in 12 of 262 Oksa/NOPHO samples and 4 of 136 GSE184692 ETV6-RUNX1 samples. The descriptive frequency across the three CNV layers was 23 of 542 (4.2%; 95% CI, 2.7-6.3). In the only expression cohort with paired chr16 and RAS annotations, GSE181157, chr16-gain B-ALL showed enrichment of chr16 dosage, 16p/MAPK3, RAS-MAPK, and ERK-feedback gene sets. These results establish chr16 gain as a recurrent low-frequency event and identify a candidate adverse EFS association and RAS-MAPK transcriptional context. Independent patient-level outcome and molecular validation are required before clinical use.

## Introduction

ETV6::RUNX1 fusion is among the most frequent initiating lesions in childhood B-cell precursor acute lymphoblastic leukemia (ALL) and is commonly associated with favorable treatment outcomes [1-4]. This average prognosis does not eliminate clinically meaningful residual risk, particularly because relapses in this subtype can occur late and may be shaped by secondary genetic changes [3,4]. The biological and clinical heterogeneity of ETV6::RUNX1-positive ALL has therefore remained an active problem for risk refinement.

Secondary copy-number alterations and cytogenetic lesions are plausible modifiers of ETV6::RUNX1 leukemia. Genome-wide studies have shown that childhood ALL contains multiple acquired genomic alterations beyond the initiating fusion [5,6]. More recently, Oksa and colleagues described genomic determinants of therapy response in ETV6::RUNX1 leukemia, reinforcing that this subtype contains clinically and molecularly relevant substructure [7]. However, uncommon whole-chromosome or broad chromosome 16 gain is difficult to study: the exposed group is small, assay platforms differ, and few public datasets jointly provide diagnostic cytogenetics, measurable residual disease (MRD), molecular annotation, and long-term event data.

The previous version of this work risked conflating several different questions. An external CNV cohort without follow-up cannot validate an EFS hazard ratio, and an expression cohort without sample-level chr16 calls cannot independently establish chr16-gain biology. We therefore restructured the analysis into non-overlapping evidence layers. TARGET was used only for EFS discovery; Oksa/NOPHO and GSE184692 were used only to assess recurrence of chr16 gain across CNV platforms; GSE181157 was used for paired expression and mutation-context analyses; and GSE227832, GSE228632, and GSE87070 were used only to describe external ETV6::RUNX1/like expression-state context. Our aim was to report what the available data support, and equally importantly, what they do not support.

## Methods

### Study design and reporting framework

This was a retrospective public-data reanalysis with a restricted local single-center audit. The report was organized according to the principles of STROBE for observational analyses and REMARK for prognostic-marker reporting [13]. The study was exploratory. Diagnostic karyotype features were evaluated in a broader candidate screen; consequently, nominal P values and Benjamini-Hochberg (BH) q values are both reported, and no statistical result was treated as practice-changing without independent outcome validation.

### TARGET ETV6::RUNX1 cohort and outcome definition

The local curated TARGET extraction contained 222 ETV6::RUNX1-positive childhood ALL patients [8,9]. Event-free survival (EFS) time was available for 219 patients, with 38 recorded events. The analysis set for karyotype-based survival estimation required an available diagnostic karyotype, valid EFS time, event status, and a parsable +16 value; it contained 141 patients and 24 events. EFS was measured from diagnosis to the first recorded qualifying event, including relapse, death, or second malignancy, with patients without a qualifying event censored at the recorded follow-up time.

The diagnostic karyotype field was treated as evaluable when it contained a non-empty result without a documented failed, uninformative, or missing assessment. The +16 exposure was defined by the patient-level regular-expression rule `(?:^|,)\\+16(?=,|/|\\[|$)` after whitespace removal. A patient was considered +16-positive when the diagnostic karyotype contained at least one +16 clone; repeated mentions across clones were not counted more than once. This conservative rule does not infer chromosome 16 gain from other cytogenetic text, and it does not convert a conventional karyotype call into an array-based CNA call.

Day-29 MRD was obtained from the curated `MRD_PERCENT_DAY_29` variable. For descriptive analyses, MRD positivity was defined as at least 0.01%. For regression sensitivity analyses, MRD was represented as log10(MRD + 0.0001). Age was dichotomized at 10 years and white blood cell count at 50 x 10^9/L only in prespecified sensitivity models; given 24 events, these multivariable estimates were interpreted as stability assessments rather than definitive adjusted effects.

### External CNV recurrence cohorts

The Oksa/NOPHO data were derived from publicly deposited CNVkit segment files associated with Oksa et al. [7]. Broad or whole-chromosome 16 gain was called using the stored broad/whole-chromosome rule: at least 80% of observed chromosome 16 panel span with log2 ratio of at least 0.3 and weighted chromosome 16 log2 ratio of at least 0.3. These files did not contain analyzable public EFS, standardized MRD, or matched RAS mutation data, and were therefore used only for prevalence estimation.

GSE184692 aCGH data were analyzed independently. ETV6-RUNX1 samples were selected from GEO subtype metadata. Chr16 probes were mapped from GPL10150 platform annotation. The primary call required a chromosome 16 mean segmented log2 value of at least 0.20 and at least 50% of chromosome 16 probes with segmented log2 value of at least 0.25. Sensitive and strict thresholds were retained as sensitivity analyses. All three thresholds produced the same number of chr16-gain calls. Because the assays, ascertainment, and call definitions differ, the cross-cohort frequency is presented as a descriptive aggregate and not as a survival meta-analysis.

### Expression and RAS-MAPK context

GSE181157, GSE227832, GSE228632, and GSE87070 were processed independently; normalized expression values were never pooled across platforms. For RNA-seq cohorts, counts were transformed to log2 counts per million and genes were retained if log2 counts per million was at least 1 in at least 5% of samples. GPL570 microarray probes in GSE87070 were mapped to symbols and collapsed by within-symbol median. Expression datasets were used at diagnosis in B-ALL reference samples as defined from their available metadata.

GSE181157 was the only expression cohort in this analysis with same-patient RNA-seq, karyotype/FISH text, and mutation text. Text-mined chr16 gain was defined by +16, trisomy 16, or chromosome-16 gain patterns in the conventional karyotype or FISH fields. RAS/MAPK mutation status was defined by a mention of KRAS, NRAS, HRAS, NF1, PTPN11, FLT3, BRAF, CBL, MAP2K1, or MAP2K2 in the mutation field. These text-derived variables were used as molecular-context annotations and not as array-level CNV calls or validated clinical assays.

Preranked gene-set enrichment analysis was run with 1,000 permutations, minimum gene-set size of 5, and a fixed seed of 20260710 [12]. Curated gene sets represented chr16 coding genes, a 16p/MAPK3-region set, RAS-MAPK core signaling, ERK-feedback genes, and B-ALL developmental-state markers. The GSE181157 chr16-gain versus no-gain comparison included 6 versus 130 B-ALL samples. Multiplicity correction was performed within each comparison across the eight curated gene sets. GSE227832, GSE228632, and GSE87070 were analyzed as external subtype-state references only because their parsed metadata lacked matched patient-level chr16-gain status.

### Statistical analysis

EFS was estimated by Kaplan-Meier methods and compared with a two-sided log-rank test. The principal regression estimate was the unpenalized univariable Cox proportional-hazards model [10]. To display dependence on sparse-data regularization, L2-penalized Cox models with penalties of 0.01, 0.05, 0.10, and 0.50 were calculated as sensitivity analyses using lifelines [11]. No penalty parameter was selected as a confirmatory primary model. Exact binomial confidence intervals were used for proportions, and Fisher exact tests compared event proportions or prevalence between cohorts. The complete diagnostic-karyotype candidate screen was controlled with BH correction. All tests were two-sided. Analyses were implemented in Python with pandas, scipy, statsmodels, lifelines, matplotlib, and gseapy.

### Ethics, data availability, and local audit

Public datasets were analyzed under the data-access conditions specified by their source repositories. The restricted single-center ETV6::RUNX1 cohort was approved by the Ethics Committee of Shanghai Children's Hospital (2015R037-F01, 2025R023-F01, and 2025R050-F01); written informed consent was obtained from guardians. The single-center cohort was not used as an independent +16 outcome validation cohort because no recorded +16-positive patient was available. Public code and aggregate derived outputs are maintained at https://github.com/ABLE-Du/ETV6-RUNX1-biomarker-reanalysis and archived at https://doi.org/10.5281/zenodo.20709667. Individual-level single-center data are restricted by privacy and institutional governance. The exact code and aggregate outputs used for this manuscript should be frozen as a tagged repository release before submission.

## Results

### A low-frequency +16 signal was observed in the TARGET EFS discovery cohort

The TARGET extraction contained 222 ETV6::RUNX1-positive patients, 219 with EFS time, and 144 with evaluable diagnostic karyotype. The karyotype-EFS analysis set comprised 141 patients with 24 events. Diagnostic +16 was detected in 7 of 144 evaluable karyotypes (4.9%; 95% CI, 2.0-9.8) and was present in 7 of 141 patients in the EFS analysis set (Table 1; Figure 1).

Four of seven +16-positive patients experienced an EFS event, compared with 20 of 134 patients without +16. The Kaplan-Meier estimate of five-year EFS was 71.4% in the +16 group and 85.5% in the no-+16 group (log-rank P = 0.004; Figure 2). The unpenalized univariable Cox model gave a hazard ratio of 4.20 (95% CI, 1.44-12.31; nominal P = 0.009). Event clustering was similarly directionally adverse by Fisher exact testing (odds ratio, 7.60; P = 0.016). However, +16 was one feature within a broader karyotype candidate screen. Its univariable BH q value was 0.308, so this result should be interpreted as a hypothesis-generating survival signal rather than a confirmed prognostic biomarker (Table 2; Supplementary Table S1).

The magnitude of the Cox estimate depended on penalization, as expected with only seven exposed patients. For the univariable model, the hazard ratio decreased from 4.20 without a penalty to 3.25 with penalty 0.10 and 1.85 with penalty 0.50 (Supplementary Figure S1 and Supplementary Table S2). This sensitivity does not negate the observed event imbalance, but it precludes selecting a single penalized estimate as a clinically stable effect size.

Day-29 MRD was available in the karyotype-evaluable subset. No +16-positive patient was MRD-positive at the 0.01% threshold, compared with 20 of 137 patients without +16 (Fisher P = 0.593). Thus, the available data do not demonstrate that +16 identifies an MRD-delayed or MRD-occult high-risk group. Multivariable models including log-transformed MRD, age, and white blood cell count remained directionally adverse, but they were based on 24 events and are presented only as sensitivity analyses (Supplementary Table S2).

### Independent CNV cohorts support recurrence of chr16 gain, not external EFS validation

External datasets showed a similar low frequency of chr16 gain across independent CNV technologies. Broad or whole-chromosome 16 gain was observed in 12 of 262 Oksa/NOPHO CNVkit samples (4.6%; 95% CI, 2.4-7.9). In GSE184692, the primary aCGH rule identified chr16 gain in 4 of 136 ETV6-RUNX1 samples (2.9%; 95% CI, 0.8-7.4). The GSE184692 estimate was statistically compatible with both TARGET and Oksa/NOPHO by Fisher exact testing (P = 0.542 and P = 0.593, respectively). Across the three CNV layers, the descriptive aggregate was 23 of 542 samples (4.2%; 95% CI, 2.7-6.3; Figure 3; Supplementary Table S3).

This concordance supports recurrence of the event and reduces the concern that the TARGET observation is solely an artifact of a single karyotype parser. It does not provide an external EFS validation because neither external CNV resource supplied a harmonized public EFS outcome suitable for analysis. Likewise, the combined frequency is descriptive and should not be interpreted as a pooled effect estimate.

### The paired expression cohort links chr16 gain to a RAS-MAPK transcriptional context

GSE181157 included 173 matched RNA-seq and metadata records, including 136 B-ALL/non-T-ALL samples. Six B-ALL samples had a text-derived chr16-gain annotation, and 47 had a RAS/MAPK mutation annotation. Five of the six chr16-gain samples also had a RAS/MAPK mutation annotation. Within the ETV6::RUNX1/like subset, one of 31 samples had chr16 gain, and that sample also carried NRAS G12D. These observations directly contradict a claim of absolute mutual exclusivity between chr16 gain and RAS/MAPK mutation in the available data.

In the B-ALL chr16-gain versus no-gain comparison, chr16 coding genes were enriched as an internal dosage positive control (normalized enrichment score [NES], 1.79; FDR q = 0.021). The 16p/MAPK3-region set (NES, 1.65; FDR q = 0.023), RAS-MAPK core set (NES, 1.69; FDR q = 0.024), and ERK-feedback set (NES, 1.60; FDR q = 0.026) were also positively enriched (Table 3; Figure 4). The RAS/MAPK-mutated B-ALL comparison independently showed strong enrichment of ERK-feedback (NES, 2.25; FDR q < 0.001) and RAS-MAPK core genes (NES, 1.93; FDR q = 0.001). These data support a transcriptional association between chr16 gain and RAS-MAPK activity, but the six exposed samples and text-derived annotations require cautious interpretation.

The three external expression cohorts broadened subtype-state context but could not test the chr16-gain mechanism directly. GSE227832 showed enrichment of a pre-BCR/pre-B program in ETV6::RUNX1/like versus other B-ALL samples (NES, 1.69; FDR q = 0.037). GSE228632 showed a directional but non-significant ERK-feedback signal (NES, 1.75; FDR q = 0.094), while GSE87070 showed B-cell developmental-state differences by single-sample enrichment but no significant RAS-MAPK or chr16-gain result. These findings are therefore contextual and not an independent replication of the GSE181157 chr16-gain comparison.

## Discussion

This multi-cohort reanalysis makes three narrow contributions. First, diagnostic +16 was associated with an adverse EFS pattern in the available TARGET karyotype subset. Second, chr16 gain occurred at a similar low frequency in two independent CNV datasets. Third, the only expression cohort with same-patient chr16 and RAS annotations linked chr16-gain B-ALL to RAS-MAPK and ERK-feedback transcriptional programs. The study therefore provides a coherent hypothesis for future validation, but not a clinically deployable marker.

The TARGET survival observation is statistically striking at the nominal level, with four events among seven +16-positive patients and a five-year EFS separation. It is nevertheless a sparse-data result. The association was not retained after BH correction across the karyotype candidate universe, and the hazard-ratio estimate changed meaningfully across reasonable L2 penalties. These facts are not peripheral qualifications; they determine the appropriate interpretation. The evidence supports a candidate adverse association, not an independently validated risk factor or a basis for treatment intensification.

The CNV analyses add a different type of support. Oksa/NOPHO and GSE184692 show that chr16 gain is not unique to the TARGET karyotype field and occurs at a frequency of approximately 3-5% across platform-separated datasets. This recurrence is important because it makes a parser artifact less likely. It does not solve the prognostic question, because outcome definitions, follow-up, and patient-level covariates were not available for harmonized survival modeling in those datasets. Keeping frequency replication separate from outcome validation avoids a common but consequential overstatement in rare-subgroup analyses.

The transcriptional data suggest a biological context rather than a causal mechanism. Chromosome 16 dosage is expected to be enriched in chr16-gain samples and should be regarded as a positive control. More informative is the concurrent enrichment of the 16p/MAPK3-region, RAS-MAPK core, and ERK-feedback sets in GSE181157. This pattern is compatible with a MAPK-active molecular context. It does not establish that MAPK3 dosage causes relapse, because exposure classification was text-derived, the chr16-gain sample count was six, and transcriptome changes can reflect co-occurring mutations or developmental state. The co-occurrence of chr16 gain and NRAS G12D in an ETV6::RUNX1 case further argues against an absolute mutual-exclusivity model.

Our analysis also corrects a potential MRD interpretation. Although no +16-positive TARGET patient was Day-29 MRD-positive, the difference from the no-+16 group was not statistically supported and was based on seven exposed patients. We therefore do not describe +16 as an MRD-occult marker. The more defensible hypothesis is that chr16 gain may identify a small biologically distinct subgroup whose clinical relevance cannot be resolved by current public data alone.

Several limitations should guide both review and future work. The primary EFS signal is based on seven exposed patients and lacks independent patient-level outcome validation. Karyotype, CNVkit, and aCGH call different biological representations of chromosome 16 gain. GSE181157 annotations were text-mined and do not replace curated array-level CNV or mutation adjudication. External expression datasets lacked matched chr16-gain annotations, and the local cohort contained no recorded +16-positive patient. Finally, the candidate-screen q value means that even the TARGET EFS association should be treated as exploratory.

The next study should be prospective and multicenter, with a harmonized diagnostic definition of whole-chromosome or broad chr16 gain, including conventional karyotype plus SNP-array or sequencing-derived CNV confirmation. At minimum, it should collect diagnostic ETV6::RUNX1 status, chr16 gain, RAS/MAPK mutations, Day-29 MRD, treatment protocol, EFS event type, and prespecified censoring rules. The current descriptive frequency of approximately 4% implies that several hundred ETV6::RUNX1-positive patients will be required to estimate a stable effect and to test whether any association is independent of MRD and treatment intensity.

In conclusion, chromosome 16 gain is a recurrent low-frequency event in ETV6::RUNX1-positive childhood ALL. TARGET provides an exploratory adverse EFS signal, while GSE181157 supports a RAS-MAPK transcriptional context. The evidence is insufficient for clinical implementation but is sufficient to justify a prospective, patient-level validation study with standardized CNV and molecular assays.

## Declarations

### Ethics approval and consent to participate

Public datasets were analyzed under their source access conditions. The single-center study was approved by the Ethics Committee of Shanghai Children's Hospital (2015R037-F01, 2025R023-F01, and 2025R050-F01). Written informed consent was obtained from guardians.

### Consent for publication

Not applicable.

### Availability of data and materials

Public code and aggregate derived outputs are maintained at https://github.com/ABLE-Du/ETV6-RUNX1-biomarker-reanalysis and archived at https://doi.org/10.5281/zenodo.20709667. TARGET data are available from their source repositories. Oksa/NOPHO CNVkit files are available through Zenodo records 10.5281/zenodo.15167703, 10.5281/zenodo.15173882, and 10.5281/zenodo.15174016. GSE184692, GSE181157, GSE227832, GSE228632, and GSE87070 are available through GEO. Individual-level single-center data are not publicly released because of privacy and institutional governance restrictions. The manuscript-specific script and aggregate output tables in this package should be deposited as a tagged public repository release before submission.

### Competing interests

The authors declare no competing interests.

### Funding

No funding to declare.

### Author contributions

ChengKan Du and JiaShi Zhu conceived the study. ChengKan Du performed the public-data reanalysis and drafted the manuscript. JiaShi Zhu curated clinical data and contributed to interpretation. SiJian Wang and Yue Zheng contributed to flow-cytometry data organization. Min Liu contributed to cytogenetic data interpretation. Hong Li supervised the study, contributed to clinical data collection, and revised the manuscript. All authors approved the final manuscript.

### Acknowledgements

We thank the patients and families who contributed to the TARGET, GEO, Oksa/NOPHO, and single-center cohorts, as well as the investigators who made the public datasets available.

## Tables

### Table 1. Cohorts, usable analysis sets, and permitted inference

| Evidence layer | Source | Usable set | Molecular variable | Permitted inference |
| --- | --- | --- | --- | --- |
| EFS discovery | TARGET | 222 ETV6::RUNX1-positive; 141 with EFS and evaluable karyotype | Diagnostic karyotype-derived +16 | Exploratory EFS association only |
| External CNV recurrence | Oksa/NOPHO | 262 ETV6::RUNX1 samples | Broad/whole chr16 gain from CNVkit segments | Prevalence/recurrence only |
| External CNV recurrence | GSE184692 | 136 ETV6-RUNX1 samples | Probe-level chr16 gain from aCGH | Prevalence/recurrence only |
| Paired molecular context | GSE181157 | 136 B-ALL samples; 31 ETV6::RUNX1/like | RNA-seq, chr16 karyotype/FISH text, mutation text | Transcriptional/RAS context only |
| External expression state | GSE227832, GSE228632, GSE87070 | 301, 60, and 574 diagnostic B-ALL/BCP samples, respectively | Expression subtype labels | External expression-state context only |
| Local restricted audit | Shanghai Children's Hospital | 115 ETV6::RUNX1/TEL-AML1 patients | Recorded karyotype +16 | Not used for +16 outcome validation because no +16-positive case was recorded |

### Table 2. TARGET +16 EFS discovery results and interpretation boundaries

| Analysis | Result | Interpretation |
| --- | --- | --- |
| Diagnostic +16 frequency | 7/144 (4.9%; 95% CI, 2.0-9.8) | Low-frequency karyotype event |
| EFS analysis set | 141 patients; 24 events | Available karyotype-EFS subset |
| EFS events | 4/7 with +16 versus 20/134 without +16 | Sparse exposed group |
| Five-year EFS | 71.4% with +16 versus 85.5% without +16 | Kaplan-Meier estimate |
| Log-rank test | P = 0.004 | Nominal comparison |
| Univariable Cox model | HR, 4.20; 95% CI, 1.44-12.31; P = 0.009 | Nominal, unadjusted estimate |
| Karyotype-screen multiplicity | BH q = 0.308 | Does not support a confirmed screen-positive biomarker |
| Day-29 MRD >= 0.01% | 0/7 with +16 versus 20/137 without +16; Fisher P = 0.593 | No evidence of MRD enrichment |

### Table 3. GSE181157 B-ALL chr16-gain expression analysis

| Curated gene set | Comparison | NES | FDR q value | Interpretation |
| --- | --- | ---: | ---: | --- |
| chr16 coding genes | chr16 gain versus no gain (6 versus 130) | 1.79 | 0.021 | Expected dosage positive control |
| 16p/MAPK3 region | chr16 gain versus no gain (6 versus 130) | 1.65 | 0.023 | Positive enrichment |
| RAS-MAPK core | chr16 gain versus no gain (6 versus 130) | 1.69 | 0.024 | Positive enrichment |
| ERK feedback | chr16 gain versus no gain (6 versus 130) | 1.60 | 0.026 | Positive enrichment |
| RAS-MAPK core | RAS/MAPK-mutated versus wild-type (47 versus 89) | 1.93 | 0.001 | Positive mutation-associated control |
| ERK feedback | RAS/MAPK-mutated versus wild-type (47 versus 89) | 2.25 | <0.001 | Positive mutation-associated control |

## Figure legends

**Figure 1. Integrated public-data design and evidentiary boundaries.** TARGET was used for exploratory EFS discovery. Oksa/NOPHO and GSE184692 were used for external CNV recurrence only. GSE181157 was the sole expression cohort with same-patient chr16 and RAS annotations and was used for molecular context. GSE227832, GSE228632, and GSE87070 were used for external ETV6::RUNX1/like expression-state context. No outcome-free cohort was pooled into survival modeling.

![Figure 1](figures/Figure_1_study_design.png)

**Figure 2. Diagnostic +16 and EFS in the TARGET karyotype-EFS subset.** The analysis included 141 ETV6::RUNX1-positive patients with evaluable diagnostic karyotype, EFS time, and event status. Shaded bands show 95% confidence intervals. At-risk counts are shown in 20-month intervals. The log-rank P value is nominal and the survival association remained exploratory after the karyotype candidate-screen multiplicity adjustment.

![Figure 2](figures/Figure_2_TARGET_EFS.png)

**Figure 3. Recurrent low-frequency chr16 gain across independent CNV platforms.** Points and horizontal lines show prevalence and exact 95% confidence intervals. The dashed line shows the descriptive frequency across the three non-overlapping CNV layers, 23/542 (4.2%). It is not a pooled EFS or RAS/MAPK estimate.

![Figure 3](figures/Figure_3_external_CNV_recurrence.png)

**Figure 4. Expression-context analysis across public cohorts.** Circle color indicates the direction of the normalized enrichment score (NES; orange, positive; blue, negative); bold values have within-comparison FDR q < 0.05. GSE181157 provides same-patient chr16/RAS expression context. The external expression datasets lack matched chr16-gain calls and are displayed only to describe ETV6::RUNX1/like subtype-state context.

![Figure 4](figures/Figure_4_expression_context.png)

## References

1. Shurtleff SA, Buijs A, Behm FG, et al. TEL/AML1 fusion resulting from a cryptic t(12;21) is the most common genetic lesion in pediatric ALL and defines a subgroup of patients with an excellent prognosis. Leukemia. 1995;9:1985-1989. PMID: 8609706.

2. Borkhardt A, Cazzaniga G, Viehmann S, et al. Incidence and clinical relevance of TEL/AML1 fusion genes in children with acute lymphoblastic leukemia enrolled in the German and Italian multicenter therapy trials. Blood. 1997;90:571-577. doi:10.1182/blood.v90.2.571.

3. Loh ML, Goldwasser MA, Silverman LB, et al. Prospective analysis of TEL/AML1-positive patients treated on Dana-Farber Cancer Institute Consortium Protocol 95-01. Blood. 2006;107:4508-4513. doi:10.1182/blood-2005-08-3451.

4. Forestier E, Heyman M, Andersen MK, et al. Outcome of ETV6/RUNX1-positive childhood acute lymphoblastic leukaemia in the NOPHO-ALL-1992 protocol: frequent late relapses but good overall survival. Br J Haematol. 2008;140:665-672. doi:10.1111/j.1365-2141.2008.06980.x.

5. Mullighan CG, Goorha S, Radtke I, et al. Genome-wide analysis of genetic alterations in acute lymphoblastic leukaemia. Nature. 2007;446:758-764. doi:10.1038/nature05690.

6. Bokemeyer A, Eckert C, Meyr F, et al. Copy number genome alterations are associated with treatment response and outcome in relapsed childhood ETV6/RUNX1-positive acute lymphoblastic leukemia. Haematologica. 2014;99:706-714. doi:10.3324/haematol.2012.072470.

7. Oksa L, et al. Genomic determinants of therapy response in ETV6::RUNX1 leukemia. Leukemia. 2025;39:2125-2139. doi:10.1038/s41375-025-02683-7.

8. Cerami E, Gao J, Dogrusoz U, et al. The cBio cancer genomics portal: an open platform for exploring multidimensional cancer genomics data. Cancer Discov. 2012;2:401-404. doi:10.1158/2159-8290.CD-12-0095.

9. Gao J, Aksoy BA, Dogrusoz U, et al. Integrative analysis of complex cancer genomics and clinical profiles using the cBioPortal. Sci Signal. 2013;6:pl1. doi:10.1126/scisignal.2004088.

10. Cox DR. Regression models and life-tables. J R Stat Soc Series B Stat Methodol. 1972;34:187-220. doi:10.1111/j.2517-6161.1972.tb00899.x.

11. Davidson-Pilon C. lifelines: survival analysis in Python. J Open Source Softw. 2019;4:1317. doi:10.21105/joss.01317.

12. Subramanian A, Tamayo P, Mootha VK, et al. Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles. Proc Natl Acad Sci U S A. 2005;102:15545-15550. doi:10.1073/pnas.0506580102.

13. McShane LM, Altman DG, Sauerbrei W, et al. Reporting recommendations for tumor marker prognostic studies (REMARK). J Natl Cancer Inst. 2005;97:1180-1184. doi:10.1093/jnci/dji237.
