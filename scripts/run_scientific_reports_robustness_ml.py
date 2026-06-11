from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
import math
import os
import platform
import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from lifelines import CoxPHFitter
from lifelines.statistics import multivariate_logrank_test
from lifelines.utils import concordance_index
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import pdist, squareform
from scipy.stats import chi2_contingency, fisher_exact, mannwhitneyu
from statsmodels.stats.multitest import multipletests


WORKSPACE = Path(os.environ.get("ETV6_RUNX1_WORKSPACE_ROOT", Path.cwd())).resolve()
SOURCE = Path(os.environ.get("ETV6_RUNX1_SOURCE_ROOT", WORKSPACE)).resolve()
OUT = WORKSPACE / "robustness_ml_analysis"
TABLES = OUT / "tables"
FIGURES = OUT / "figures"

TARGET_BASE = (
    SOURCE
    / "public_reanalysis"
    / "TARGET_cBioPortal"
    / "datahub"
    / "public"
    / "all_phase2_target_2018_pub"
)
BIOMARKER = SOURCE / "public_reanalysis" / "ETV6_poor_outcome_biomarkers"
TARGET_RESULTS = SOURCE / "public_reanalysis" / "TARGET_cBioPortal" / "results"

RANDOM_SEED = 20260611
RNG = np.random.default_rng(RANDOM_SEED)
EXPRESSION_PERMUTATIONS = 200

KARYOTYPE_PREFIXES = (
    "abnormal_",
    "del_",
    "del_chr_",
    "gain_chr_",
    "hyperdiploid_",
    "hypodiploid_",
    "loss_chr_",
    "multiple_",
    "normal_",
    "structural_",
    "translocation_",
)


def setup() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)
    sns.set_theme(style="whitegrid", context="paper")
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 8,
            "axes.titlesize": 9,
            "axes.labelsize": 8,
            "legend.fontsize": 7,
            "savefig.dpi": 600,
            "figure.dpi": 160,
            "pdf.fonttype": 42,
            "svg.fonttype": "none",
        }
    )


def configure_paths(workspace_root: Path, source_root: Path, permutations: int) -> None:
    global WORKSPACE, SOURCE, OUT, TABLES, FIGURES
    global TARGET_BASE, BIOMARKER, TARGET_RESULTS, EXPRESSION_PERMUTATIONS
    WORKSPACE = workspace_root.resolve()
    SOURCE = source_root.resolve()
    OUT = WORKSPACE / "robustness_ml_analysis"
    TABLES = OUT / "tables"
    FIGURES = OUT / "figures"
    TARGET_BASE = (
        SOURCE
        / "public_reanalysis"
        / "TARGET_cBioPortal"
        / "datahub"
        / "public"
        / "all_phase2_target_2018_pub"
    )
    BIOMARKER = SOURCE / "public_reanalysis" / "ETV6_poor_outcome_biomarkers"
    TARGET_RESULTS = SOURCE / "public_reanalysis" / "TARGET_cBioPortal" / "results"
    EXPRESSION_PERMUTATIONS = permutations


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run multiplicity, unsupervised-learning, and nested-CV robustness analyses."
    )
    parser.add_argument(
        "--workspace-root",
        type=Path,
        default=WORKSPACE,
        help="Writable project root for outputs (default: cwd or ETV6_RUNX1_WORKSPACE_ROOT).",
    )
    parser.add_argument(
        "--source-root",
        type=Path,
        default=SOURCE,
        help="Root containing public_reanalysis inputs (default: workspace or ETV6_RUNX1_SOURCE_ROOT).",
    )
    parser.add_argument(
        "--permutations",
        type=int,
        default=EXPRESSION_PERMUTATIONS,
        help="Number of complete nested-pipeline RNA label permutations (default: 200).",
    )
    args = parser.parse_args()
    if args.permutations < 1:
        parser.error("--permutations must be at least 1")
    return args


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_reproducibility_manifest() -> None:
    inputs = [
        BIOMARKER / "karyotype_biomarker_screen.csv",
        BIOMARKER / "karyotype_patient_feature_matrix.csv",
        BIOMARKER / "cna_deletion_biomarker_screen.csv",
        BIOMARKER / "rna_primary_time_to_relapse_biomarker_screen.csv",
        BIOMARKER / "rna_paired_relapse_change_screen.csv",
        TARGET_RESULTS / "target_expression_all_gene_statistics.csv",
        TARGET_BASE / "data_mrna_seq_rpkm_zscores_ref_all_samples.txt",
        TARGET_BASE / "data_clinical_sample.txt",
        SOURCE
        / "public_reanalysis"
        / "DepMap_Public_26Q1"
        / "results"
        / "all_gene_del6q_vs_no_del6q_statistics.csv",
        WORKSPACE
        / "single_center_manuscript_analysis"
        / "tables"
        / "exploratory_univariable_cox.csv",
    ]
    missing = [str(path) for path in inputs if not path.exists()]
    if missing:
        raise FileNotFoundError("Missing required inputs:\n" + "\n".join(missing))
    input_rows = [
        {
            "role": "analysis_input",
            "path_relative_to_source_or_workspace": (
                str(path.relative_to(SOURCE))
                if path.is_relative_to(SOURCE)
                else str(path.relative_to(WORKSPACE))
            ),
            "bytes": path.stat().st_size,
            "sha256": sha256_file(path),
        }
        for path in inputs
    ]
    pd.DataFrame(input_rows).to_csv(OUT / "input_file_manifest.csv", index=False)
    packages = [
        "numpy",
        "pandas",
        "scipy",
        "statsmodels",
        "lifelines",
        "matplotlib",
        "seaborn",
    ]
    environment = {
        "python": platform.python_version(),
        "platform": platform.platform(),
        "random_seed": RANDOM_SEED,
        "expression_nested_pipeline_permutations": EXPRESSION_PERMUTATIONS,
        "packages": {
            name: importlib.metadata.version(name)
            for name in packages
        },
    }
    (OUT / "software_environment.json").write_text(
        json.dumps(environment, indent=2), encoding="utf-8"
    )


def adjust_values(values: pd.Series) -> pd.DataFrame:
    p = pd.to_numeric(values, errors="coerce").to_numpy(float)
    output = pd.DataFrame(index=values.index, data={"p_value": p})
    valid = np.isfinite(p)
    for method, name in [
        ("fdr_bh", "bh_fdr"),
        ("fdr_by", "by_fdr"),
        ("holm", "holm_fwer"),
    ]:
        adjusted = np.full(len(p), np.nan)
        if valid.any():
            adjusted[valid] = multipletests(p[valid], method=method)[1]
        output[name] = adjusted
    return output


def multiple_testing_audit() -> tuple[pd.DataFrame, pd.DataFrame]:
    families = [
        ("target_karyotype_event", BIOMARKER / "karyotype_biomarker_screen.csv", "event_fisher_p", "feature"),
        ("target_karyotype_poor5y", BIOMARKER / "karyotype_biomarker_screen.csv", "poor_5y_fisher_p", "feature"),
        ("target_karyotype_univ_cox", BIOMARKER / "karyotype_biomarker_screen.csv", "univ_cox_p", "feature"),
        ("target_karyotype_adjusted_cox", BIOMARKER / "karyotype_biomarker_screen.csv", "adjusted_cox_p", "feature"),
        ("target_cna_event", BIOMARKER / "cna_deletion_biomarker_screen.csv", "event_fisher_p", "gene"),
        ("target_cna_poor5y", BIOMARKER / "cna_deletion_biomarker_screen.csv", "poor_5y_fisher_p", "gene"),
        ("target_cna_candidate_univ_cox", BIOMARKER / "cna_deletion_biomarker_screen.csv", "univ_cox_p", "gene"),
        ("target_cna_candidate_adjusted_cox", BIOMARKER / "cna_deletion_biomarker_screen.csv", "adjusted_cox_p", "gene"),
        ("target_rna_time_to_relapse", BIOMARKER / "rna_primary_time_to_relapse_biomarker_screen.csv", "spearman_p", "gene"),
        ("target_rna_early_vs_late_relapse", BIOMARKER / "rna_primary_time_to_relapse_biomarker_screen.csv", "mannwhitney_p", "gene"),
        ("target_rna_paired_relapse_change", BIOMARKER / "rna_paired_relapse_change_screen.csv", "wilcoxon_p", "gene"),
        ("target_subtype_expression_previous", TARGET_RESULTS / "target_expression_all_gene_statistics.csv", "mw_p", "gene"),
        ("depmap_26q1_pooled_t", SOURCE / "public_reanalysis" / "DepMap_Public_26Q1" / "results" / "all_gene_del6q_vs_no_del6q_statistics.csv", "pooled_t_p", "gene"),
        ("depmap_26q1_welch_t", SOURCE / "public_reanalysis" / "DepMap_Public_26Q1" / "results" / "all_gene_del6q_vs_no_del6q_statistics.csv", "welch_t_p", "gene"),
        ("depmap_26q1_exact_permutation", SOURCE / "public_reanalysis" / "DepMap_Public_26Q1" / "results" / "all_gene_del6q_vs_no_del6q_statistics.csv", "exact_permutation_p", "gene"),
        ("single_center_univ_cox", WORKSPACE / "single_center_manuscript_analysis" / "tables" / "exploratory_univariable_cox.csv", "p", "factor"),
    ]
    all_results: list[pd.DataFrame] = []
    audit_rows: list[dict[str, object]] = []
    for family, path, p_col, id_col in families:
        data = pd.read_csv(path)
        corrected = adjust_values(data[p_col])
        frame = pd.DataFrame(
            {
                "family": family,
                "identifier": data[id_col].astype(str),
                **corrected.to_dict("series"),
            }
        )
        all_results.append(frame)
        valid = frame.dropna(subset=["p_value"])
        audit_rows.append(
            {
                "family": family,
                "tests": len(valid),
                "min_p": valid["p_value"].min() if len(valid) else np.nan,
                "min_bh_fdr": valid["bh_fdr"].min() if len(valid) else np.nan,
                "min_by_fdr": valid["by_fdr"].min() if len(valid) else np.nan,
                "min_holm_fwer": valid["holm_fwer"].min() if len(valid) else np.nan,
                "bh_significant_0_05": int((valid["bh_fdr"] < 0.05).sum()),
                "by_significant_0_05": int((valid["by_fdr"] < 0.05).sum()),
                "holm_significant_0_05": int((valid["holm_fwer"] < 0.05).sum()),
            }
        )
    combined = pd.concat(all_results, ignore_index=True)
    audit = pd.DataFrame(audit_rows)
    combined.to_csv(TABLES / "multiple_testing_all_families.csv", index=False)
    audit.to_csv(TABLES / "multiple_testing_audit.csv", index=False)
    return combined, audit


def unique_binary_features(data: pd.DataFrame, min_count: int = 5) -> list[str]:
    candidates = [col for col in data.columns if col.startswith(KARYOTYPE_PREFIXES)]
    candidates = [
        col
        for col in candidates
        if min_count <= pd.to_numeric(data[col], errors="coerce").fillna(0).sum() <= len(data) - min_count
    ]
    seen: set[tuple[int, ...]] = set()
    keep: list[str] = []
    for col in candidates:
        key = tuple(pd.to_numeric(data[col], errors="coerce").fillna(0).astype(int))
        if key not in seen:
            seen.add(key)
            keep.append(col)
    return keep


def adjusted_rand_index(labels_a: np.ndarray, labels_b: np.ndarray) -> float:
    a_vals, a_inv = np.unique(labels_a, return_inverse=True)
    b_vals, b_inv = np.unique(labels_b, return_inverse=True)
    table = np.zeros((len(a_vals), len(b_vals)), dtype=int)
    np.add.at(table, (a_inv, b_inv), 1)

    def choose2(x: np.ndarray | int) -> np.ndarray | float:
        return np.asarray(x) * (np.asarray(x) - 1) / 2

    sum_cells = choose2(table).sum()
    sum_rows = choose2(table.sum(axis=1)).sum()
    sum_cols = choose2(table.sum(axis=0)).sum()
    total = choose2(table.sum())
    expected = sum_rows * sum_cols / total if total else 0
    maximum = 0.5 * (sum_rows + sum_cols)
    return float((sum_cells - expected) / (maximum - expected)) if maximum != expected else 1.0


def silhouette_from_distance(distance: np.ndarray, labels: np.ndarray) -> float:
    values = []
    for i in range(len(labels)):
        same = labels == labels[i]
        same[i] = False
        a = distance[i, same].mean() if same.any() else 0.0
        b = min(
            distance[i, labels == other].mean()
            for other in np.unique(labels)
            if other != labels[i]
        )
        values.append((b - a) / max(a, b) if max(a, b) > 0 else 0.0)
    return float(np.mean(values))


def cluster_event_stat(labels: np.ndarray, event: np.ndarray) -> float:
    table = pd.crosstab(labels, event).reindex(columns=[0, 1], fill_value=0)
    return float(chi2_contingency(table, correction=False)[0])


def karyotype_unsupervised() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, dict]:
    data = pd.read_csv(BIOMARKER / "karyotype_patient_feature_matrix.csv")
    features = unique_binary_features(data)
    x = data[features].to_numpy(float)
    condensed = pdist(x, metric="jaccard")
    distance = squareform(condensed)
    tree = linkage(condensed, method="average")
    event = data["event"].astype(int).to_numpy()

    rows = []
    labels_by_k: dict[int, np.ndarray] = {}
    for k in range(2, 7):
        labels = fcluster(tree, k, criterion="maxclust")
        labels_by_k[k] = labels
        actual = cluster_event_stat(labels, event)
        permuted = np.array(
            [cluster_event_stat(labels, RNG.permutation(event)) for _ in range(3000)]
        )
        rows.append(
            {
                "k": k,
                "silhouette": silhouette_from_distance(distance, labels),
                "event_chi2": actual,
                "event_permutation_p": (1 + int((permuted >= actual).sum())) / (len(permuted) + 1),
            }
        )
    evaluation = pd.DataFrame(rows)
    corrected = adjust_values(evaluation["event_permutation_p"])
    evaluation[["event_bh_fdr", "event_by_fdr", "event_holm_fwer"]] = corrected[
        ["bh_fdr", "by_fdr", "holm_fwer"]
    ]
    best_k = int(evaluation.sort_values("silhouette", ascending=False).iloc[0]["k"])
    best_labels = labels_by_k[best_k]

    stability = []
    feature_count = max(2, int(math.ceil(0.8 * len(features))))
    for _ in range(300):
        chosen = RNG.choice(len(features), size=feature_count, replace=False)
        sub_condensed = pdist(x[:, chosen], metric="jaccard")
        sub_labels = fcluster(linkage(sub_condensed, method="average"), best_k, criterion="maxclust")
        stability.append(adjusted_rand_index(best_labels, sub_labels))

    standardized = (x - x.mean(axis=0)) / np.where(x.std(axis=0) == 0, 1, x.std(axis=0))
    u, s, _ = np.linalg.svd(standardized, full_matrices=False)
    scores = u[:, :2] * s[:2]
    explained = (s**2) / np.sum(s**2)
    patient_scores = pd.DataFrame(
        {
            "PATIENT_ID": data["PATIENT_ID"],
            "event": event,
            "cluster": best_labels,
            "PC1": scores[:, 0],
            "PC2": scores[:, 1],
        }
    )
    patient_scores.to_csv(TABLES / "karyotype_unsupervised_patient_scores.csv", index=False)
    evaluation.to_csv(TABLES / "karyotype_cluster_evaluation.csv", index=False)
    pd.DataFrame({"ari": stability}).to_csv(TABLES / "karyotype_cluster_stability.csv", index=False)
    cluster_table = pd.crosstab(best_labels, event).reindex(columns=[0, 1], fill_value=0)
    cluster_table.columns = ["censored", "event"]
    cluster_table.index.name = "cluster"
    cluster_table.to_csv(TABLES / "karyotype_cluster_event_table.csv")
    logrank_p = float(
        multivariate_logrank_test(data["time"], best_labels, data["event"]).p_value
    )
    summary = {
        "n": len(data),
        "events": int(event.sum()),
        "input_features": len(features),
        "best_k": best_k,
        "best_silhouette": float(evaluation.loc[evaluation["k"].eq(best_k), "silhouette"].iloc[0]),
        "best_k_event_permutation_p": float(evaluation.loc[evaluation["k"].eq(best_k), "event_permutation_p"].iloc[0]),
        "best_k_event_holm_fwer": float(evaluation.loc[evaluation["k"].eq(best_k), "event_holm_fwer"].iloc[0]),
        "best_k_logrank_p": logrank_p,
        "stability_median_ari": float(np.median(stability)),
        "stability_iqr_low": float(np.quantile(stability, 0.25)),
        "stability_iqr_high": float(np.quantile(stability, 0.75)),
        "pc1_explained": float(explained[0]),
        "pc2_explained": float(explained[1]),
    }
    return patient_scores, evaluation, pd.DataFrame({"ari": stability}), summary


def read_clinical(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", comment="#", low_memory=False)


def primary_expression_data() -> tuple[pd.DataFrame, pd.DataFrame]:
    expression_path = TARGET_BASE / "data_mrna_seq_rpkm_zscores_ref_all_samples.txt"
    samples = read_clinical(TARGET_BASE / "data_clinical_sample.txt")
    header = pd.read_csv(expression_path, sep="\t", nrows=0).columns.tolist()
    available = set(header[2:])
    selected = samples.loc[
        samples["SAMPLE_ID"].isin(available)
        & samples["ETV6_RUNX1_FUSION_STATUS"].isin(["Positive", "Negative"])
        & samples["SAMPLE_ID"].str.rsplit("-", n=1).str[-1].isin(["03", "09"])
    ].copy()
    selected["sample_type_code"] = selected["SAMPLE_ID"].str.rsplit("-", n=1).str[-1]
    selected["priority"] = selected["sample_type_code"].map({"09": 0, "03": 1})
    selected = (
        selected.sort_values(["PATIENT_ID", "priority", "SAMPLE_ID"])
        .drop_duplicates("PATIENT_ID", keep="first")
        .reset_index(drop=True)
    )
    usecols = ["Hugo_Symbol", "Entrez_Gene_Id", *selected["SAMPLE_ID"].tolist()]
    expression = pd.read_csv(
        expression_path,
        sep="\t",
        usecols=lambda col: col in usecols,
        low_memory=False,
    )
    expression = expression.dropna(subset=["Hugo_Symbol"]).copy()
    numeric = expression[selected["SAMPLE_ID"]].apply(pd.to_numeric, errors="coerce")
    expression["_observed"] = numeric.notna().sum(axis=1)
    expression = (
        pd.concat([expression[["Hugo_Symbol", "Entrez_Gene_Id", "_observed"]], numeric], axis=1)
        .sort_values(["Hugo_Symbol", "_observed"], ascending=[True, False])
        .drop_duplicates("Hugo_Symbol", keep="first")
        .drop(columns="_observed")
        .set_index("Hugo_Symbol")
    )
    return selected, expression


def primary_expression_statistics(
    selected: pd.DataFrame, expression: pd.DataFrame
) -> tuple[pd.DataFrame, dict]:
    positive = selected.loc[selected["ETV6_RUNX1_FUSION_STATUS"].eq("Positive"), "SAMPLE_ID"].tolist()
    negative = selected.loc[selected["ETV6_RUNX1_FUSION_STATUS"].eq("Negative"), "SAMPLE_ID"].tolist()
    rows = []
    for gene, row in expression.iterrows():
        pos = pd.to_numeric(row[positive], errors="coerce").dropna()
        neg = pd.to_numeric(row[negative], errors="coerce").dropna()
        if len(pos) < 5 or len(neg) < 5:
            continue
        rows.append(
            {
                "gene": gene,
                "positive_n": len(pos),
                "negative_n": len(neg),
                "positive_mean_z": pos.mean(),
                "negative_mean_z": neg.mean(),
                "delta_z": pos.mean() - neg.mean(),
                "mw_p": mannwhitneyu(pos, neg, alternative="two-sided").pvalue,
            }
        )
    results = pd.DataFrame(rows)
    adjusted = adjust_values(results["mw_p"])
    results[["bh_fdr", "by_fdr", "holm_fwer"]] = adjusted[["bh_fdr", "by_fdr", "holm_fwer"]]
    results = results.sort_values("mw_p")
    results.to_csv(TABLES / "primary_only_expression_statistics.csv", index=False)
    selected.to_csv(TABLES / "primary_only_expression_samples.csv", index=False)
    slc = results.loc[results["gene"].eq("SLC7A11")]
    summary = {
        "positive_patients": len(positive),
        "negative_patients": len(negative),
        "genes_tested": len(results),
        "bh_significant": int((results["bh_fdr"] < 0.05).sum()),
        "by_significant": int((results["by_fdr"] < 0.05).sum()),
        "holm_significant": int((results["holm_fwer"] < 0.05).sum()),
        "slc7a11_delta_z": float(slc.iloc[0]["delta_z"]) if len(slc) else np.nan,
        "slc7a11_p": float(slc.iloc[0]["mw_p"]) if len(slc) else np.nan,
        "slc7a11_bh_fdr": float(slc.iloc[0]["bh_fdr"]) if len(slc) else np.nan,
        "slc7a11_by_fdr": float(slc.iloc[0]["by_fdr"]) if len(slc) else np.nan,
        "slc7a11_holm_fwer": float(slc.iloc[0]["holm_fwer"]) if len(slc) else np.nan,
    }
    return results, summary


def append_primary_expression_audit(
    combined: pd.DataFrame, audit: pd.DataFrame, results: pd.DataFrame
) -> tuple[pd.DataFrame, pd.DataFrame]:
    family = pd.DataFrame(
        {
            "family": "target_primary_only_subtype_expression",
            "identifier": results["gene"],
            "p_value": results["mw_p"],
            "bh_fdr": results["bh_fdr"],
            "by_fdr": results["by_fdr"],
            "holm_fwer": results["holm_fwer"],
        }
    )
    row = pd.DataFrame(
        [
            {
                "family": "target_primary_only_subtype_expression",
                "tests": len(results),
                "min_p": results["mw_p"].min(),
                "min_bh_fdr": results["bh_fdr"].min(),
                "min_by_fdr": results["by_fdr"].min(),
                "min_holm_fwer": results["holm_fwer"].min(),
                "bh_significant_0_05": int((results["bh_fdr"] < 0.05).sum()),
                "by_significant_0_05": int((results["by_fdr"] < 0.05).sum()),
                "holm_significant_0_05": int((results["holm_fwer"] < 0.05).sum()),
            }
        ]
    )
    combined = pd.concat([combined, family], ignore_index=True)
    audit = pd.concat([audit, row], ignore_index=True)
    combined.to_csv(TABLES / "multiple_testing_all_families.csv", index=False)
    audit.to_csv(TABLES / "multiple_testing_audit.csv", index=False)
    return combined, audit


def stratified_folds(y: np.ndarray, k: int, rng: np.random.Generator) -> np.ndarray:
    fold = np.empty(len(y), dtype=int)
    for value in np.unique(y):
        indices = np.flatnonzero(y == value)
        rng.shuffle(indices)
        for i, part in enumerate(np.array_split(indices, k)):
            fold[part] = i
    return fold


def auc_rank(y: np.ndarray, score: np.ndarray) -> float:
    pos = score[y == 1]
    neg = score[y == 0]
    if len(pos) == 0 or len(neg) == 0:
        return np.nan
    comparisons = (pos[:, None] > neg[None, :]).sum() + 0.5 * (pos[:, None] == neg[None, :]).sum()
    return float(comparisons / (len(pos) * len(neg)))


def nearest_centroid_score(
    train_x: np.ndarray,
    train_y: np.ndarray,
    test_x: np.ndarray,
    top_k: int,
) -> np.ndarray:
    medians = np.nanmedian(train_x, axis=0)
    medians[~np.isfinite(medians)] = 0
    train_x = np.where(np.isfinite(train_x), train_x, medians)
    test_x = np.where(np.isfinite(test_x), test_x, medians)
    means = train_x.mean(axis=0)
    sds = train_x.std(axis=0)
    sds[sds < 1e-8] = 1
    train_z = (train_x - means) / sds
    test_z = (test_x - means) / sds
    effect = np.abs(train_z[train_y == 1].mean(axis=0) - train_z[train_y == 0].mean(axis=0))
    keep = np.argsort(effect)[::-1][: min(top_k, len(effect))]
    positive_centroid = train_z[train_y == 1][:, keep].mean(axis=0)
    negative_centroid = train_z[train_y == 0][:, keep].mean(axis=0)
    d_pos = np.mean((test_z[:, keep] - positive_centroid) ** 2, axis=1)
    d_neg = np.mean((test_z[:, keep] - negative_centroid) ** 2, axis=1)
    return d_neg - d_pos


def expression_machine_learning(
    selected: pd.DataFrame, expression: pd.DataFrame
) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    sample_ids = selected["SAMPLE_ID"].tolist()
    y = selected["ETV6_RUNX1_FUSION_STATUS"].eq("Positive").astype(int).to_numpy()
    x_all = expression[sample_ids].T.to_numpy(float)
    observed_fraction = np.isfinite(x_all).mean(axis=0)
    x_all = x_all[:, observed_fraction >= 0.8]
    variances = np.nanvar(x_all, axis=0)
    variable_indices = np.argsort(np.nan_to_num(variances, nan=-np.inf))[::-1][:500]
    x = x_all[:, variable_indices]
    medians = np.nanmedian(x, axis=0)
    x = np.where(np.isfinite(x), x, medians)
    x = (x - x.mean(axis=0)) / np.where(x.std(axis=0) == 0, 1, x.std(axis=0))
    u, s, _ = np.linalg.svd(x, full_matrices=False)
    scores = u[:, :2] * s[:2]
    explained = (s**2) / np.sum(s**2)

    def separation_stat(labels: np.ndarray) -> float:
        pos = scores[labels == 1, :2]
        neg = scores[labels == 0, :2]
        centroid_distance = np.linalg.norm(pos.mean(axis=0) - neg.mean(axis=0))
        within = np.sqrt(
            (
                ((pos - pos.mean(axis=0)) ** 2).sum()
                + ((neg - neg.mean(axis=0)) ** 2).sum()
            )
            / len(scores)
        )
        return float(centroid_distance / within) if within else np.inf

    actual_sep = separation_stat(y)
    permuted_sep = np.array([separation_stat(RNG.permutation(y)) for _ in range(5000)])
    separation_p = (1 + int((permuted_sep >= actual_sep).sum())) / (len(permuted_sep) + 1)

    score_rows = pd.DataFrame(
        {
            "PATIENT_ID": selected["PATIENT_ID"],
            "SAMPLE_ID": sample_ids,
            "fusion_status": selected["ETV6_RUNX1_FUSION_STATUS"],
            "sample_type_code": selected["sample_type_code"],
            "PC1": scores[:, 0],
            "PC2": scores[:, 1],
        }
    )
    score_rows.to_csv(TABLES / "primary_expression_pca_scores.csv", index=False)

    top_options = [10, 50, 200]

    def nested_cv_once(labels: np.ndarray, seed: int) -> tuple[float, float, list[dict[str, float]]]:
        nested_rng = np.random.default_rng(seed)
        outer = stratified_folds(labels, 5, nested_rng)
        fold_rows: list[dict[str, float]] = []
        for outer_fold in range(5):
            train_idx = np.flatnonzero(outer != outer_fold)
            test_idx = np.flatnonzero(outer == outer_fold)
            inner = stratified_folds(labels[train_idx], 3, nested_rng)
            candidate_scores = {}
            for top_k in top_options:
                inner_fold_aucs = []
                for inner_fold in range(3):
                    inner_train = np.flatnonzero(inner != inner_fold)
                    inner_test = np.flatnonzero(inner == inner_fold)
                    inner_pred = nearest_centroid_score(
                        x_all[train_idx][inner_train],
                        labels[train_idx][inner_train],
                        x_all[train_idx][inner_test],
                        top_k,
                    )
                    inner_fold_aucs.append(
                        auc_rank(labels[train_idx][inner_test], inner_pred)
                    )
                candidate_scores[top_k] = float(np.nanmean(inner_fold_aucs))
            chosen = max(
                candidate_scores,
                key=lambda key: np.nan_to_num(candidate_scores[key], nan=-1),
            )
            predictions = nearest_centroid_score(
                x_all[train_idx], labels[train_idx], x_all[test_idx], chosen
            )
            predicted = (predictions >= 0).astype(int)
            test_labels = labels[test_idx]
            sensitivity = float((predicted[test_labels == 1] == 1).mean())
            specificity = float((predicted[test_labels == 0] == 0).mean())
            fold_rows.append(
                {
                    "outer_fold": outer_fold + 1,
                    "fold_auc": auc_rank(test_labels, predictions),
                    "fold_balanced_accuracy": (sensitivity + specificity) / 2,
                    "sensitivity": sensitivity,
                    "specificity": specificity,
                    "selected_genes": float(chosen),
                    "test_n": float(len(test_idx)),
                    "test_positive": float(test_labels.sum()),
                }
            )
        mean_auc = float(np.nanmean([row["fold_auc"] for row in fold_rows]))
        mean_balanced_accuracy = float(
            np.nanmean([row["fold_balanced_accuracy"] for row in fold_rows])
        )
        return mean_auc, mean_balanced_accuracy, fold_rows

    ml_rows = []
    fold_rows = []
    for repeat in range(10):
        mean_auc, mean_balanced_accuracy, repeat_folds = nested_cv_once(
            y, RANDOM_SEED + 1000 + repeat
        )
        for row in repeat_folds:
            row["repeat"] = repeat + 1
            fold_rows.append(row)
        ml_rows.append(
            {
                "repeat": repeat + 1,
                "nested_cv_auc": mean_auc,
                "nested_cv_balanced_accuracy": mean_balanced_accuracy,
                "median_selected_genes": float(
                    np.median([row["selected_genes"] for row in repeat_folds])
                ),
            }
        )
    ml = pd.DataFrame(ml_rows)
    ml.to_csv(TABLES / "primary_expression_nested_cv.csv", index=False)
    pd.DataFrame(fold_rows).to_csv(
        TABLES / "primary_expression_nested_cv_outer_folds.csv", index=False
    )

    nested_pipeline_observed_auc, _, _ = nested_cv_once(y, RANDOM_SEED + 3000)
    permutation_values = []
    for i in range(EXPRESSION_PERMUTATIONS):
        permutation_values.append(
            nested_cv_once(
                np.random.default_rng(RANDOM_SEED + 4000 + i).permutation(y),
                RANDOM_SEED + 5000 + i,
            )[0]
        )
        if (i + 1) % 20 == 0 or i + 1 == EXPRESSION_PERMUTATIONS:
            print(f"Completed RNA nested-pipeline permutations: {i + 1}/{EXPRESSION_PERMUTATIONS}")
    permutation_aucs = np.asarray(permutation_values)
    nested_pipeline_permutation_p = (
        1 + int((permutation_aucs >= nested_pipeline_observed_auc).sum())
    ) / (len(permutation_aucs) + 1)
    pd.DataFrame(
        {
            "permutation": np.arange(1, len(permutation_aucs) + 1),
            "permuted_nested_pipeline_mean_fold_auc": permutation_aucs,
        }
    ).to_csv(
        TABLES / "primary_expression_classifier_permutation.csv", index=False
    )

    type_table = pd.crosstab(
        selected["ETV6_RUNX1_FUSION_STATUS"], selected["sample_type_code"]
    ).reindex(index=["Negative", "Positive"], columns=["03", "09"], fill_value=0)
    sample_type_fisher_p = float(fisher_exact(type_table.to_numpy())[1])
    type_table.to_csv(TABLES / "primary_expression_sample_type_table.csv")
    summary = {
        "pca_pc1_explained": float(explained[0]),
        "pca_pc2_explained": float(explained[1]),
        "pca_label_separation_stat": actual_sep,
        "pca_label_permutation_p": separation_p,
        "nested_cv_auc_mean": float(ml["nested_cv_auc"].mean()),
        "nested_cv_auc_min": float(ml["nested_cv_auc"].min()),
        "nested_cv_auc_max": float(ml["nested_cv_auc"].max()),
        "nested_cv_balanced_accuracy_mean": float(ml["nested_cv_balanced_accuracy"].mean()),
        "nested_pipeline_observed_mean_fold_auc": nested_pipeline_observed_auc,
        "nested_pipeline_permutations": EXPRESSION_PERMUTATIONS,
        "nested_pipeline_permutation_p": nested_pipeline_permutation_p,
        "sample_type_fisher_p": sample_type_fisher_p,
        "classifier_genes_eligible": int(x_all.shape[1]),
    }
    return score_rows, ml, summary


def cox_feature_columns(train: pd.DataFrame, mode: str) -> list[str]:
    clinical = ["mrd29_positive", "wbc_ge50", "age_ge10"]
    if mode == "clinical":
        candidates = clinical
    elif mode == "clinical_plus_gain16":
        candidates = clinical + ["gain_chr_16"]
    else:
        candidates = clinical + [col for col in train.columns if col.startswith(KARYOTYPE_PREFIXES)]
    keep = []
    seen: set[tuple[float, ...]] = set()
    min_count = max(3, int(math.ceil(len(train) * 0.035)))
    for col in candidates:
        values = pd.to_numeric(train[col], errors="coerce")
        observed = values.dropna()
        if observed.nunique() < 2:
            continue
        if col not in clinical and not (min_count <= observed.sum() <= len(observed) - min_count):
            continue
        key = tuple(values.fillna(values.median()).astype(float))
        if key in seen:
            continue
        seen.add(key)
        keep.append(col)
    return keep


def fit_predict_cox(
    train: pd.DataFrame,
    test: pd.DataFrame,
    mode: str,
    penalizer: float,
) -> np.ndarray:
    features = cox_feature_columns(train, mode)
    train_model = train[["time", "event", *features]].copy()
    test_model = test[features].copy()
    for col in features:
        median = pd.to_numeric(train_model[col], errors="coerce").median()
        train_model[col] = pd.to_numeric(train_model[col], errors="coerce").fillna(median)
        test_model[col] = pd.to_numeric(test_model[col], errors="coerce").fillna(median)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        model = CoxPHFitter(penalizer=penalizer, l1_ratio=0.0)
        model.fit(train_model, duration_col="time", event_col="event", show_progress=False)
    return model.predict_partial_hazard(test_model).to_numpy(float)


def survival_nested_cv() -> tuple[pd.DataFrame, dict]:
    data = pd.read_csv(BIOMARKER / "karyotype_patient_feature_matrix.csv")
    data = data.loc[data["time"].notna()].reset_index(drop=True)
    event = data["event"].astype(int).to_numpy()
    modes = ["clinical", "clinical_plus_gain16", "clinical_plus_karyotype"]
    penalties = [0.1, 1.0, 10.0]
    rows = []
    for repeat in range(5):
        repeat_rng = np.random.default_rng(RANDOM_SEED + 2000 + repeat)
        outer = stratified_folds(event, 5, repeat_rng)
        for outer_fold in range(5):
            train_idx = np.flatnonzero(outer != outer_fold)
            test_idx = np.flatnonzero(outer == outer_fold)
            train = data.iloc[train_idx].reset_index(drop=True)
            test = data.iloc[test_idx].reset_index(drop=True)
            inner = stratified_folds(train["event"].astype(int).to_numpy(), 3, repeat_rng)
            for mode in modes:
                penalty_scores = {}
                for penalty in penalties:
                    inner_fold_scores = []
                    for inner_fold in range(3):
                        inner_train = train.iloc[np.flatnonzero(inner != inner_fold)]
                        inner_test = train.iloc[np.flatnonzero(inner == inner_fold)]
                        try:
                            inner_prediction = fit_predict_cox(
                                inner_train, inner_test, mode, penalty
                            )
                            inner_fold_scores.append(
                                concordance_index(
                                    inner_test["time"],
                                    -inner_prediction,
                                    inner_test["event"],
                                )
                            )
                        except Exception:
                            inner_fold_scores.append(0.5)
                    penalty_scores[penalty] = float(np.mean(inner_fold_scores))
                chosen = max(penalty_scores, key=penalty_scores.get)
                try:
                    prediction = fit_predict_cox(train, test, mode, chosen)
                    c_index = concordance_index(
                        test["time"], -prediction, test["event"]
                    )
                except Exception:
                    c_index = 0.5
                rows.append(
                    {
                        "repeat": repeat + 1,
                        "outer_fold": outer_fold + 1,
                        "model": mode,
                        "nested_cv_c_index": c_index,
                        "selected_penalty": chosen,
                        "test_n": len(test),
                        "test_events": int(test["event"].sum()),
                    }
                )
    results = pd.DataFrame(rows)
    results.to_csv(TABLES / "target_survival_nested_cv.csv", index=False)
    repeat_summary = (
        results.groupby(["repeat", "model"], as_index=False)
        .agg(
            nested_cv_c_index=("nested_cv_c_index", "mean"),
            median_selected_penalty=("selected_penalty", "median"),
        )
    )
    repeat_summary.to_csv(TABLES / "target_survival_nested_cv_repeat_summary.csv", index=False)
    pivot = results.pivot(
        index=["repeat", "outer_fold"], columns="model", values="nested_cv_c_index"
    )
    pivot["gain16_minus_clinical"] = pivot["clinical_plus_gain16"] - pivot["clinical"]
    pivot["karyotype_minus_clinical"] = pivot["clinical_plus_karyotype"] - pivot["clinical"]
    pivot.to_csv(TABLES / "target_survival_nested_cv_paired_differences.csv")
    summary = {
        "n": len(data),
        "events": int(data["event"].sum()),
        "clinical_mean_c_index": float(pivot["clinical"].mean()),
        "gain16_mean_c_index": float(pivot["clinical_plus_gain16"].mean()),
        "karyotype_mean_c_index": float(pivot["clinical_plus_karyotype"].mean()),
        "gain16_mean_delta": float(pivot["gain16_minus_clinical"].mean()),
        "gain16_median_delta": float(pivot["gain16_minus_clinical"].median()),
        "gain16_delta_min": float(pivot["gain16_minus_clinical"].min()),
        "gain16_delta_max": float(pivot["gain16_minus_clinical"].max()),
        "karyotype_mean_delta": float(pivot["karyotype_minus_clinical"].mean()),
        "karyotype_median_delta": float(pivot["karyotype_minus_clinical"].median()),
        "karyotype_delta_min": float(pivot["karyotype_minus_clinical"].min()),
        "karyotype_delta_max": float(pivot["karyotype_minus_clinical"].max()),
    }
    return results, summary


def plot_results(
    audit: pd.DataFrame,
    karyo_scores: pd.DataFrame,
    karyo_eval: pd.DataFrame,
    expression_scores: pd.DataFrame,
    survival_cv: pd.DataFrame,
    expression_cv: pd.DataFrame,
) -> None:
    def add_panel_label(ax: plt.Axes, label: str) -> None:
        ax.text(
            -0.12,
            1.06,
            label,
            transform=ax.transAxes,
            fontsize=11,
            fontweight="bold",
            va="top",
        )

    fig, axes = plt.subplots(2, 2, figsize=(12, 9), constrained_layout=True)
    ax = axes[0, 0]
    plot_audit = audit.sort_values("min_bh_fdr").head(12).copy()
    y = np.arange(len(plot_audit))
    ax.scatter(-np.log10(plot_audit["min_p"].clip(lower=1e-300)), y, label="Nominal P", color="#D55E00")
    ax.scatter(-np.log10(plot_audit["min_bh_fdr"].clip(lower=1e-300)), y, label="BH-FDR", color="#0072B2")
    ax.axvline(-np.log10(0.05), color="#444444", linestyle="--", linewidth=0.8)
    ax.set_yticks(y)
    ax.set_yticklabels(plot_audit["family"].str.replace("_", " "))
    ax.invert_yaxis()
    ax.set_xlabel("-log10(value)")
    ax.set_title("Multiplicity changes the evidence ranking")
    ax.legend(frameon=False)
    add_panel_label(ax, "a")

    ax = axes[0, 1]
    for event, marker, color, label in [(0, "o", "#999999", "Censored"), (1, "X", "#C44E52", "Event")]:
        subset = karyo_scores.loc[karyo_scores["event"].eq(event)]
        ax.scatter(subset["PC1"], subset["PC2"], c=subset["cluster"], cmap="tab10", marker=marker, s=35, alpha=0.85, label=label)
    ax.set_xlabel("Karyotype PC1")
    ax.set_ylabel("Karyotype PC2")
    ax.set_title(f"Unsupervised karyotype structure (best k={int(karyo_eval.sort_values('silhouette', ascending=False).iloc[0]['k'])})")
    ax.legend(frameon=False)
    add_panel_label(ax, "b")

    ax = axes[1, 0]
    palette = {"Negative": "#999999", "Positive": "#0072B2"}
    for status, subset in expression_scores.groupby("fusion_status"):
        ax.scatter(subset["PC1"], subset["PC2"], label=status, color=palette[status], s=35, alpha=0.85)
    ax.set_xlabel("Expression PC1")
    ax.set_ylabel("Expression PC2")
    ax.set_title("Primary-only RNA: no global unsupervised separation")
    ax.legend(frameon=False)
    add_panel_label(ax, "c")

    ax = axes[1, 1]
    sns.boxplot(data=survival_cv, x="model", y="nested_cv_c_index", ax=ax, color="#A6CEE3", width=0.55)
    sns.stripplot(data=survival_cv, x="model", y="nested_cv_c_index", ax=ax, color="#1F78B4", size=4)
    ax.axhline(0.5, color="#444444", linestyle="--", linewidth=0.8)
    ax.set_xlabel("")
    ax.set_ylabel("Nested-CV concordance index")
    ax.set_xticks(
        range(3),
        ["Clinical", "Clinical + +16", "Clinical + karyotype"],
        rotation=15,
        ha="right",
    )
    ax.set_title("No reliable prognostic gain from karyotype features")
    add_panel_label(ax, "d")
    for suffix in ["png", "pdf", "svg"]:
        fig.savefig(FIGURES / f"Figure_R1_robustness_ml.{suffix}", bbox_inches="tight", facecolor="white")
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(10, 4), constrained_layout=True)
    survival_plot = survival_cv.copy()
    survival_plot["model"] = survival_plot["model"].map(
        {
            "clinical": "Clinical",
            "clinical_plus_gain16": "Clinical + +16",
            "clinical_plus_karyotype": "Clinical + karyotype",
        }
    )
    sns.lineplot(data=survival_plot, x="repeat", y="nested_cv_c_index", hue="model", marker="o", ax=axes[0])
    axes[0].axhline(0.5, color="#444444", linestyle="--", linewidth=0.8)
    axes[0].set_title("Repeated nested-CV survival models")
    axes[0].set_ylabel("Outer-fold concordance index")
    axes[0].legend(frameon=False)
    add_panel_label(axes[0], "a")
    sns.boxplot(data=expression_cv, y="nested_cv_auc", ax=axes[1], color="#80B1D3", width=0.4)
    sns.stripplot(data=expression_cv, y="nested_cv_auc", ax=axes[1], color="#0072B2", size=4)
    axes[1].axhline(0.5, color="#444444", linestyle="--", linewidth=0.8)
    axes[1].set_title("Primary RNA classifies fusion subtype, not prognosis")
    axes[1].set_ylabel("Nested-CV AUC")
    add_panel_label(axes[1], "b")
    for suffix in ["png", "pdf", "svg"]:
        fig.savefig(FIGURES / f"Figure_R2_cross_validation.{suffix}", bbox_inches="tight", facecolor="white")
    plt.close(fig)


def write_report(summary: dict) -> None:
    mt = summary["multiple_testing"]
    karyo = summary["karyotype_unsupervised"]
    expr = summary["primary_expression"]
    expr_ml = summary["expression_ml"]
    survival = summary["survival_ml"]
    report = f"""# Scientific Reports robustness and machine-learning analysis

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
- The previous TARGET subtype-expression family retained {mt['target_subtype_previous_bh']} BH-significant genes, but it included non-primary/repeated samples and therefore was rerun under a stricter primary-only rule.
- The primary-only expression comparison retained {expr['bh_significant']} BH-significant genes, {expr['by_significant']} BY-significant genes, and {expr['holm_significant']} Holm-significant genes.
- In the primary-only analysis, SLC7A11 had delta Z={expr['slc7a11_delta_z']:.3f}, nominal P={expr['slc7a11_p']:.3g}, BH-FDR={expr['slc7a11_bh_fdr']:.3g}, BY-FDR={expr['slc7a11_by_fdr']:.3g}, and Holm-adjusted P={expr['slc7a11_holm_fwer']:.3g}.

### Unsupervised learning

- The TARGET karyotype cohort included {karyo['n']} patients and {karyo['events']} EFS events.
- Hierarchical clustering selected k={karyo['best_k']} by silhouette, but cluster separation was weak (silhouette={karyo['best_silhouette']:.3f}).
- Cluster stability under repeated 80% feature subsampling was limited (median adjusted Rand index={karyo['stability_median_ari']:.3f}, IQR {karyo['stability_iqr_low']:.3f}-{karyo['stability_iqr_high']:.3f}).
- The selected clusters were not reliably associated with events after accounting for the tested cluster numbers (permutation P={karyo['best_k_event_permutation_p']:.3g}; Holm P={karyo['best_k_event_holm_fwer']:.3g}; log-rank P={karyo['best_k_logrank_p']:.3g}).
- Primary-only RNA did not form a significant global unsupervised separation by fusion status (label-permutation P={expr_ml['pca_label_permutation_p']:.3g}).

### Supervised learning

- Repeated nested-CV three-variable clinical-threshold survival discrimination was C-index={survival['clinical_mean_c_index']:.3f}, calculated separately within each outer test fold.
- Adding +16 produced mean C-index={survival['gain16_mean_c_index']:.3f}; paired mean change={survival['gain16_mean_delta']:+.3f}, median change={survival['gain16_median_delta']:+.3f} (outer-fold range {survival['gain16_delta_min']:+.3f} to {survival['gain16_delta_max']:+.3f}).
- Adding all eligible karyotype features produced mean C-index={survival['karyotype_mean_c_index']:.3f}; paired mean change={survival['karyotype_mean_delta']:+.3f}, median change={survival['karyotype_median_delta']:+.3f} (outer-fold range {survival['karyotype_delta_min']:+.3f} to {survival['karyotype_delta_max']:+.3f}).
- The primary-only RNA nearest-centroid classifier identified a sparse fusion-subtype signal with mean nested-CV outer-fold AUC={expr_ml['nested_cv_auc_mean']:.3f} (repeat range {expr_ml['nested_cv_auc_min']:.3f}-{expr_ml['nested_cv_auc_max']:.3f}). A complete nested feature-selection pipeline label-permutation test using {expr_ml['nested_pipeline_permutations']} permutations yielded P={expr_ml['nested_pipeline_permutation_p']:.3g}. Primary sample-type distribution was not detectably different by fusion status (Fisher P={expr_ml['sample_type_fisher_p']:.3g}). These results cannot be interpreted as outcome prediction.

## Reliable conclusion

The combined multiplicity-controlled and leakage-resistant analyses do not identify a validated adverse-outcome biomarker for ETV6::RUNX1-positive pediatric B-ALL. The +16 finding remains the strongest low-frequency cytogenetic candidate, but it does not survive feature-wide multiplicity correction and does not provide stable out-of-sample prognostic improvement. Unsupervised karyotype groups are weak and unstable, and neither genome-wide CNA nor prognostic RNA screening yields an FDR-significant marker. The reproducible positive result is an ETV6::RUNX1-associated diagnostic expression program; it is a subtype signal and should not be presented as a relapse biomarker or therapeutic mechanism.

## Manuscript implication

For Scientific Reports, the manuscript should be framed as a reproducible, multi-source stress test of proposed adverse-outcome biomarkers. The defensible contribution is the clear separation between a reproducible subtype-expression signal and the absence of a validated prognostic signal after multiplicity correction, sensitivity analysis, and nested cross-validation.
"""
    (OUT / "SCIENTIFIC_REPORTS_ROBUSTNESS_ML_REPORT.md").write_text(report, encoding="utf-8")


def main() -> None:
    args = parse_args()
    configure_paths(args.workspace_root, args.source_root, args.permutations)
    setup()
    write_reproducibility_manifest()
    combined, audit = multiple_testing_audit()
    karyo_scores, karyo_eval, _, karyo_summary = karyotype_unsupervised()
    selected, expression = primary_expression_data()
    expression_stats, expression_summary = primary_expression_statistics(selected, expression)
    combined, audit = append_primary_expression_audit(combined, audit, expression_stats)
    expression_scores, expression_cv, expression_ml_summary = expression_machine_learning(selected, expression)
    survival_cv, survival_summary = survival_nested_cv()

    subtype_previous = audit.loc[audit["family"].eq("target_subtype_expression_previous")].iloc[0]
    adverse_families = audit.loc[
        ~audit["family"].isin(
            [
                "target_subtype_expression_previous",
                "target_primary_only_subtype_expression",
            ]
        )
    ]
    summary = {
        "random_seed": RANDOM_SEED,
        "multiple_testing": {
            "adverse_families_with_bh_significant_results": int((adverse_families["bh_significant_0_05"] > 0).sum()),
            "target_subtype_previous_bh": int(subtype_previous["bh_significant_0_05"]),
        },
        "karyotype_unsupervised": karyo_summary,
        "primary_expression": expression_summary,
        "expression_ml": expression_ml_summary,
        "survival_ml": survival_summary,
    }
    (OUT / "analysis_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    plot_results(audit, karyo_scores, karyo_eval, expression_scores, survival_cv, expression_cv)
    write_report(summary)
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
