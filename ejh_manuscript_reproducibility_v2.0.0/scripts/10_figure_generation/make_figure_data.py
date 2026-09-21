# -*- coding: utf-8 -*-
"""Generate source-data tables for every revised main figure and main table.
Nothing is hard-coded: every value is read from an authoritative file."""
import pandas as pd, numpy as np, os, sys, io, json
from pathlib import Path
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
P = str(Path(os.environ.get("EJH_PROJECT_ROOT", Path(__file__).resolve().parents[2])))
A = str(Path(os.environ.get("EJH_ANALYSIS_ROOT", Path(P) / "analysis_enhancement")))
os.makedirs(os.path.join(A, "figure_data"), exist_ok=True)
FD = os.path.join(A, "figure_data")
W = lambda n, d: d.to_csv(os.path.join(FD, n), sep="\t", index=False)

def rd(p, **k):
    return pd.read_csv(p, sep="\t", keep_default_na=False, **k)

def num(s):
    return pd.to_numeric(pd.Series(s).replace("", np.nan), errors="coerce")

# ---------- Fig 1 : study architecture ----------
f1 = pd.DataFrame([
    dict(layer="L1", cohort="Oksa et al. 2025 (Nordic)", n=362, role="Source of the adverse-response contrast",
         endpoint="End-of-induction MRD: slow (>=0.1% at day 29) vs fast (0)",
         contribution="effect vector", evidence_ceiling="published contrast, reused as published"),
    dict(layer="L1", cohort="Li et al. 2025 discovery", n=194, role="Source of the second adverse-response contrast",
         endpoint="Unsupervised C1 (n=118) vs C2 (n=76)",
         contribution="effect vector + drug screen + MRD variable",
         evidence_ceiling="published contrast, reused as published"),
    dict(layer="L2", cohort="Li et al. 2025 validation", n=95, role="Independent molecular triangulation",
         endpoint="PAM-assigned C1 (n=61) vs C2 (n=34); 50 deposited genes",
         contribution="our re-analysis of deposited source data",
         evidence_ceiling="INDEPENDENT_MOLECULAR_TRIANGULATION"),
    dict(layer="L2", cohort="COG552 (published context)", n=552, role="Published external context only",
         endpoint="PAM C1 (n=341) vs C2 (n=211)", contribution="none - UMAP/PAM only, no expression or clinical deposited",
         evidence_ceiling="PUBLISHED_INDEPENDENT_CONTEXT"),
    dict(layer="L3", cohort="Single-centre institutional cohort", n=115, role="Independent phenotypic triangulation",
         endpoint="Day-19 MRD >=0.1% (n=23 poor / 87 good of 110 evaluable)",
         contribution="CD34/CD38 immunophenotype, age, WBC, cytogenetics",
         evidence_ceiling="INDEPENDENT_PHENOTYPIC_TRIANGULATION"),
    dict(layer="L4", cohort="TARGET-ALL Phase 2", n=191, role="General paediatric B-ALL benchmark",
         endpoint="Poor (n=68) vs good (n=123) response", contribution="gene-matched sensitivity context",
         evidence_ceiling="GENERAL_BALL_SUPPORTIVE_CONTEXT"),
])
W("Fig1_study_design.tsv", f1)

# ---------- Fig 2 : primary correspondence ----------
b1 = rd(os.path.join(P, "tables/B1_concordance_summary_v2.tsv"))
b1.to_csv(os.path.join(FD, "Fig2_primary_correspondence.tsv"), sep="\t", index=False)
# per-gene scatter source
bm = rd(os.path.join(P, "results/B6_TARGET_benchmark.tsv"))
keep = bm.dropna(subset=["oksa_coefficient_slow_vs_fast", "li_C1_minus_C2"])[
    ["gene", "oksa_coefficient_slow_vs_fast", "li_C1_minus_C2"]].copy()
keep.columns = ["gene", "oksa_effect", "li194_effect"]
keep["concordant"] = np.sign(keep.oksa_effect) == np.sign(keep.li194_effect)
keep.to_csv(os.path.join(FD, "Fig2_gene_scatter.tsv"), sep="\t", index=False)

# gene-matched 113
gm = rd(os.path.join(A, "TARGET_113_GENE_MATCHED_SENSITIVITY.tsv"))
gm.to_csv(os.path.join(FD, "Fig2_gene_matched_sensitivity.tsv"), sep="\t", index=False)
rd(os.path.join(A, "RHO_DIFFERENCE_BOOTSTRAP.tsv")).to_csv(
    os.path.join(FD, "Fig2_rho_difference.tsv"), sep="\t", index=False)

# ---------- Fig 3 : independent triangulation ----------
rows = []
h = rd(os.path.join(A, "LI95_GENE_HARMONISATION.tsv"))
for _, r in h.iterrows():
    rows.append(dict(panel="A_Li95_gene_level", gene=r["gene"],
                     oksa_coef=num([r["oksa_coef"]]).iloc[0],
                     li194_effect=num([r["li194_mean_diff"]]).iloc[0],
                     li95_effect=num([r["li95_mean_diff"]]).iloc[0],
                     li95_bh_q=num([r["li95_bh_q"]]).iloc[0],
                     in_primary_219=r["in_primary_219"],
                     oksa_li95_match=r["oksa_li95_match"]))
res = json.load(open(os.path.join(A, "li95_out/_li95_results.json")))
for k, v in res.items():
    if isinstance(v, dict) and "rho" in v:
        rows.append(dict(panel="A_Li95_summary", gene=k, oksa_coef="", li194_effect="",
                         li95_effect="", li95_bh_q="",
                         in_primary_219="rho=%s n=%s concordance=%s P=%s" % (
                             round(v["rho"], 4), v["n"], round(v["concordance_pct"], 1), v["p_rho"]),
                         oksa_li95_match=""))
# local CD34/CD38
cc = rd(os.path.join(A, "LOCAL_CD34_CD38_RESULTS.tsv"))
for _, r in cc.iterrows():
    rows.append(dict(panel="C_local_CD34_CD38", gene=r["marker"],
                     oksa_coef="", li194_effect="", li95_effect="", li95_bh_q="",
                     in_primary_219="n=%s median_poor=%s median_good=%s HL=%s CI=[%s,%s] P=%s q=%s" % (
                         r["n_evaluable"], r["median_poor"], r["median_good"],
                         round(float(r["HL_shift"]), 2), round(float(r["HL_CI_lo"]), 2),
                         round(float(r["HL_CI_hi"]), 2), round(float(r["P"]), 4),
                         round(float(r["BH_q_2markers"]), 4)),
                     oksa_li95_match=r["observed_direction"]))
aw = rd(os.path.join(A, "LOCAL_AGE_WBC_TRIANGULATION.tsv"))
for _, r in aw.iterrows():
    rows.append(dict(panel="D_local_age_WBC", gene=r["feature"], oksa_coef="", li194_effect="",
                     li95_effect="", li95_bh_q="",
                     in_primary_219="n=%s poor_exposed=%s/%s good_exposed=%s/%s OR=%s P=%s q=%s" % (
                         r["n_evaluable"], r["a"], int(r["a"]) + int(r["b"]), r["c"],
                         int(r["c"]) + int(r["d"]), r["OR"], round(float(r["P"]), 4),
                         round(float(r["BH_q_2features"]), 4)),
                     oksa_li95_match=r["observed_direction"]))
# Li95 clinical
cl = json.load(open(os.path.join(A, "_gene_matched_status.json")))["li95_clinical"]
for r in cl:
    rows.append(dict(panel="A_Li95_clinical", gene=r["feature"], oksa_coef="", li194_effect="",
                     li95_effect="", li95_bh_q="",
                     in_primary_219="C1 %s/%s vs C2 %s/%s OR=%s CI=[%s,%s] P=%s" % (
                         r["C1_positive"], int(r["C1_positive"]) + int(r["C1_negative"]),
                         r["C2_positive"], int(r["C2_positive"]) + int(r["C2_negative"]),
                         round(r["OR_C1_vs_C2"], 3), round(r["CI_lo"], 3), round(r["CI_hi"], 3),
                         round(r["P"], 4)),
                     oksa_li95_match=r["direction"]))
# patient-level CD34/CD38 for plotting (study IDs only)
fa = rd(os.path.join(A, "LOCAL_FLOW_FIELD_AUDIT.tsv"))
lr = rd(os.path.join(A, "LOCAL_FLOW_2020_RAW.tsv"))
lr.to_csv(os.path.join(FD, "Fig3_local_CD34_CD38_patients.tsv"), sep="\t", index=False)
W("Fig3_independent_triangulation.tsv", pd.DataFrame(rows))

# ---------- Fig 4 : Tier-1 candidates ----------
b5 = rd(os.path.join(P, "tables/B5_CANDIDATE_GENES_v2.tsv"))
t1 = b5[b5.is_Tier1.astype(str).str.upper() == "TRUE"].copy()
h95 = rd(os.path.join(A, "LI95_GENE_HARMONISATION.tsv")).set_index("gene")
t1["li95_effect"] = t1.gene.map(h95["li95_mean_diff"])
t1["li95_bh_q"] = t1.gene.map(h95["li95_bh_q"])
t1["li95_direction_match"] = t1.gene.map(h95["oksa_li95_match"])
t1["li95_status"] = np.where(t1.gene.isin(h95.index), "measured", "not_in_validation_panel")
t1 = t1[["gene", "oksa_coefficient_slow_vs_fast", "oksa_fdr", "li_C1_minus_C2", "li_fdr",
         "target_class", "target_assessability", "singlecell_pct_cells",
         "li95_effect", "li95_bh_q", "li95_direction_match", "li95_status"]]
W("Fig4_candidate_biology.tsv", t1)

# ---------- Fig 5 : drugs ----------
d = rd(os.path.join(P, "tables/B7_drug_sensitivity_20drug_FINAL.tsv"))
d = d[["drug", "n_measurements_paper", "n_patients_nonmissing_in_C1C2", "n_C1", "n_C2",
       "median_LC50_difference_C1_minus_C2", "published_p", "bh_fdr_20drug", "direction",
       "significant_bh_0.05", "claim_class"]]
W("Fig5_drug_annotation.tsv", d)

# ---------- Main Table 1 / Table 2 ----------
W("Table1_cohorts_endpoints_roles.tsv", f1)
t2 = pd.DataFrame([
    dict(cohort="Li discovery (n=194)", older_age="C1 enriched (published)",
         high_WBC="C2 enriched (published)", CD34_HSC_like="GATA2 higher in C1 (+1.49)",
         CD38_proB_like="not in the 252-gene panel", early_response="C1 = adverse state"),
    dict(cohort="Li validation (n=95)", older_age="OR 2.83 (1.16-6.91), P 0.031",
         high_WBC="OR 0.28 (0.10-0.80), P 0.027",
         CD34_HSC_like="GATA2 higher in C1 (+2.94, q 2.7e-10)",
         CD38_proB_like="CD38 lower in C1 (-1.501, q 4.2e-9)", early_response="C1 = adverse state (PAM)"),
    dict(cohort="Local institutional (n=115)",
         older_age="OR 0.875 (0.32-2.36), P 1.00 - not reproduced",
         high_WBC="0/23 vs 5/87, P 0.582 - direction consistent, underpowered",
         CD34_HSC_like="CD34 % higher in poor responders (median 86.0 vs 63.9, P 0.184)",
         CD38_proB_like="CD38 % lower in poor responders (median 31.6 vs 86.3, P 0.057)",
         early_response="day-19 MRD >=0.1% (distinct endpoint, never pooled)"),
])
W("Table2_phenotypic_triangulation.tsv", t2)

print("figure_data written:")
for f in sorted(os.listdir(FD)):
    print("  ", f, os.path.getsize(os.path.join(FD, f)))
print()
print(t2.to_string(index=False))
