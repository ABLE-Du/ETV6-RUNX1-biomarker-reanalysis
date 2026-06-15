from __future__ import annotations

import argparse
from pathlib import Path
import math
import re
import warnings

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from lifelines import CoxPHFitter
from lifelines.statistics import proportional_hazard_test
from lifelines.utils import concordance_index
from scipy.stats import beta, contingency, fisher_exact


ROOT = Path(__file__).resolve().parents[1]
WORKSPACE = ROOT.parent
INPUT = WORKSPACE / "genetic_risk_score" / "inputs"
SCREEN = WORKSPACE / "efs_factor_reanalysis" / "inputs" / "karyotype_biomarker_screen.csv"
TABLE = ROOT / "tables"
FIGURE = ROOT / "figures"
TABLE.mkdir(exist_ok=True)
FIGURE.mkdir(exist_ok=True)

SEED = 20260615
N_BOOT = 1000
RIDGE_EFS = 0.05
RIDGE_INCREMENTAL = 0.5
MRD_OFFSET = 0.0001

STRUCTURAL_TOKEN = re.compile(r"(?:DEL|ADD|T|INV|DUP|DER|DIC|IDIC|I|R)\(", re.I)
NUMERICAL_TOKEN = re.compile(r"(?:^|[,/])(?:\+|-)(?:\d+|X|Y|MAR)(?=[,/\[(]|$)", re.I)
FAILED_KARYOTYPE = re.compile(
    r"未做|未查|无分裂|未见分裂|无分裂相|未发现可供分析|失败|不详|未知|NOT\s*DONE|FAILED|UNKNOWN|NAN|NONE",
    re.I,
)


def normalize_karyotype(value: object) -> str:
    if pd.isna(value):
        return ""
    text = str(value).upper()
    replacements = {
        "（": "(",
        "）": ")",
        "，": ",",
        "；": ";",
        "：": ":",
        "＋": "+",
        "－": "-",
        "［": "[",
        "］": "]",
        " ": "",
    }
    for source, target in replacements.items():
        text = text.replace(source, target)
    return text


def parse_karyotype(value: object) -> dict[str, int | float]:
    text = normalize_karyotype(value)
    evaluable = int(bool(text) and not FAILED_KARYOTYPE.search(text))
    structural_n = len(STRUCTURAL_TOKEN.findall(text))
    numerical_n = len(NUMERICAL_TOKEN.findall(text))
    abnormal_n = structural_n + numerical_n
    parsed = {
        "karyotype_evaluable": evaluable,
        "structural_token_n": structural_n,
        "numerical_token_n": numerical_n,
        "abnormal_report_token_n": abnormal_n,
        "gain16": int(bool(re.search(r"(?:^|[,/])\+16(?=[,/\[(]|$)", text))),
        "abnormal_report_tokens_ge3": int(abnormal_n >= 3),
        "del6q": int(bool(re.search(r"DEL\(6\)\([^)]*Q|DEL\(6Q|DEL6Q", text))),
    }
    if not evaluable:
        for key in [
            "structural_token_n",
            "numerical_token_n",
            "abnormal_report_token_n",
            "gain16",
            "abnormal_report_tokens_ge3",
            "del6q",
        ]:
            parsed[key] = np.nan
    return parsed


def add_score(data: pd.DataFrame, karyotype_col: str) -> pd.DataFrame:
    parsed = pd.DataFrame([parse_karyotype(value) for value in data[karyotype_col]], index=data.index)
    out = pd.concat([data.drop(columns=parsed.columns, errors="ignore").copy(), parsed], axis=1)
    out["candidate_score"] = (
        2 * out["gain16"] + out["abnormal_report_tokens_ge3"] + out["del6q"]
    )
    out["candidate_score_high"] = (out["candidate_score"] >= 1).astype(float)
    out.loc[out["karyotype_evaluable"] == 0, ["candidate_score", "candidate_score_high"]] = np.nan
    return out


def contains_etv6(row: pd.Series) -> bool:
    columns = [
        "融合基因",
        "其它融合基因",
        "融合基因(RNAseq报告)",
        "分子分型(RNAseq报告)",
        "FISH",
        "FISH其它",
    ]
    text = " ".join(str(row.get(column, "")) for column in columns)
    return bool(re.search(r"ETV6|TEL.?AML|12;21", text, re.I))


def assert_historical_alignment(raw: pd.DataFrame, analytic: pd.DataFrame) -> None:
    if len(raw) != len(analytic):
        raise ValueError("Historical raw and analytic tables have different row counts.")
    checks = [
        (raw["来源表"].fillna("").astype(str), analytic["来源表"].fillna("").astype(str), "source"),
        (
            pd.to_datetime(raw["诊断日期"], errors="coerce").dt.strftime("%Y-%m-%d").fillna(""),
            pd.to_datetime(analytic["diagnosis_date"], errors="coerce").dt.strftime("%Y-%m-%d").fillna(""),
            "diagnosis_date",
        ),
        (raw["性别"].fillna("").astype(str), analytic["性别"].fillna("").astype(str), "sex"),
    ]
    for left, right, label in checks:
        if not left.reset_index(drop=True).equals(right.reset_index(drop=True)):
            raise ValueError(f"Historical row alignment failed for {label}.")
    for raw_col, analytic_col in [("年龄", "age_years"), ("WBC", "wbc")]:
        left = pd.to_numeric(raw[raw_col], errors="coerce")
        right = pd.to_numeric(analytic[analytic_col], errors="coerce")
        if not np.isclose(left, right, equal_nan=True, rtol=1e-6, atol=1e-6).all():
            raise ValueError(f"Historical row alignment failed for {raw_col}.")


def load_cohorts() -> dict[str, pd.DataFrame]:
    target = pd.read_csv(INPUT / "target_etv6_positive_patient_derived_variables.csv")
    target = add_score(target, "KARYOTYPE")
    target["time"] = pd.to_numeric(target["time"], errors="coerce")
    target["event"] = pd.to_numeric(target["event"], errors="coerce")
    target["mrd"] = pd.to_numeric(target["MRD_PERCENT_DAY_29"], errors="coerce")
    target["mrd_delayed"] = (target["mrd"] >= 0.01).astype(float)
    target.loc[target["mrd"].isna(), "mrd_delayed"] = np.nan

    workbook = INPUT / "ALL2020-2025_single_center_cleaned_pseudonymized.xlsx"
    raw_historical = pd.read_excel(workbook, sheet_name="结局数据_合并")
    raw_historical = raw_historical.loc[raw_historical["TEL/AML"].eq("阳性")].copy().reset_index(drop=True)
    historical = pd.read_csv(INPUT / "historical_outcome_analytic_cohort_pseudonymized.csv")
    historical = historical.reset_index(drop=True)
    assert_historical_alignment(raw_historical, historical)
    historical = historical.rename(
        columns={
            "gain16": "old_gain16",
            "del6q": "old_del6q",
            "complex_karyotype": "old_complex_karyotype",
        }
    )
    historical["raw_karyotype"] = raw_historical["染色体"]
    historical = add_score(historical, "raw_karyotype")
    historical["time"] = pd.to_numeric(historical["efs_months"], errors="coerce")
    historical["event"] = pd.to_numeric(historical["efs_event"], errors="coerce")
    historical["mrd"] = pd.to_numeric(historical["d19_mrd"], errors="coerce")
    historical["mrd_delayed"] = (historical["mrd"] >= 0.01).astype(float)
    historical.loc[historical["mrd"].isna(), "mrd_delayed"] = np.nan

    registry = pd.read_excel(workbook, sheet_name="登记数据_合并")
    contemporary = registry.loc[registry.apply(contains_etv6, axis=1)].copy()
    contemporary = add_score(contemporary, "核型")
    contemporary["mrd"] = pd.to_numeric(contemporary["D19MRD(%)"], errors="coerce")
    contemporary["mrd_delayed"] = (contemporary["mrd"] >= 0.01).astype(float)
    contemporary.loc[contemporary["mrd"].isna(), "mrd_delayed"] = np.nan
    contemporary["d46_mrd"] = pd.to_numeric(contemporary["D46MRD(%)"], errors="coerce")
    contemporary["d46_persistent"] = (contemporary["d46_mrd"] >= 0.01).astype(float)
    contemporary.loc[contemporary["d46_mrd"].isna(), "d46_persistent"] = np.nan

    return {
        "TARGET": target,
        "Single-center historical": historical,
        "Single-center contemporary": contemporary,
    }


def exact_binomial_ci(x: int, n: int) -> tuple[float, float]:
    if n == 0:
        return np.nan, np.nan
    low = 0.0 if x == 0 else beta.ppf(0.025, x, n - x + 1)
    high = 1.0 if x == n else beta.ppf(0.975, x + 1, n - x)
    return float(low), float(high)


def build_prevalence_table(cohorts: dict[str, pd.DataFrame]) -> pd.DataFrame:
    rows = []
    variables = ["gain16", "abnormal_report_tokens_ge3", "del6q", "candidate_score_high"]
    target = cohorts["TARGET"]
    for cohort, data in cohorts.items():
        for variable in variables:
            values = data[variable].dropna().astype(int)
            positive = int(values.sum())
            n = len(values)
            low, high = exact_binomial_ci(positive, n)
            row = {
                "cohort": cohort,
                "variable": variable,
                "positive_n": positive,
                "evaluable_n": n,
                "prevalence": positive / n if n else np.nan,
                "ci_low": low,
                "ci_high": high,
                "fisher_p_vs_TARGET": np.nan,
            }
            if cohort != "TARGET":
                target_values = target[variable].dropna().astype(int)
                table = [
                    [int(target_values.sum()), int(len(target_values) - target_values.sum())],
                    [positive, n - positive],
                ]
                row["fisher_p_vs_TARGET"] = fisher_exact(table).pvalue
            rows.append(row)
    return pd.DataFrame(rows)


def build_historical_agreement(historical: pd.DataFrame) -> pd.DataFrame:
    rows = []
    old_to_new = [
        ("old_gain16", "gain16"),
        ("old_complex_karyotype", "abnormal_report_tokens_ge3"),
        ("old_del6q", "del6q"),
    ]
    for old, new in old_to_new:
        use = historical[[old, new]].dropna().astype(int)
        rows.append(
            {
                "old_variable": old,
                "harmonized_variable": new,
                "n_compared": len(use),
                "old_positive_n": int(use[old].sum()),
                "harmonized_positive_n": int(use[new].sum()),
                "agreement_n": int((use[old] == use[new]).sum()),
                "percent_agreement": float((use[old] == use[new]).mean()),
            }
        )
    old_score = (
        2 * pd.to_numeric(historical["old_gain16"], errors="coerce")
        + pd.to_numeric(historical["old_complex_karyotype"], errors="coerce")
        + pd.to_numeric(historical["old_del6q"], errors="coerce")
    )
    old_high = (old_score >= 1).astype(float)
    old_high.loc[historical["karyotype_class"].eq("Uninformative")] = np.nan
    use = pd.DataFrame({"old_high": old_high, "new_high": historical["candidate_score_high"]}).dropna().astype(int)
    rows.append(
        {
            "old_variable": "old_score_high",
            "harmonized_variable": "candidate_score_high",
            "n_compared": len(use),
            "old_positive_n": int(use["old_high"].sum()),
            "harmonized_positive_n": int(use["new_high"].sum()),
            "agreement_n": int((use["old_high"] == use["new_high"]).sum()),
            "percent_agreement": float((use["old_high"] == use["new_high"]).mean()),
        }
    )
    return pd.DataFrame(rows)


def fit_cox(data: pd.DataFrame, variables: list[str], penalizer: float) -> tuple[CoxPHFitter, pd.DataFrame]:
    use = data[["time", "event", *variables]].dropna().copy()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        model = CoxPHFitter(penalizer=penalizer)
        model.fit(use, duration_col="time", event_col="event")
    return model, use


def cox_row(data: pd.DataFrame, variable: str, cohort: str, penalizer: float) -> dict[str, object]:
    use = data[["time", "event", variable]].dropna()
    row: dict[str, object] = {
        "cohort": cohort,
        "variable": variable,
        "penalizer": penalizer,
        "n": len(use),
        "events": int(use["event"].sum()),
    }
    if len(use) < 20 or use["event"].sum() < 5 or use[variable].nunique() < 2:
        row["error"] = "insufficient data"
        return row
    try:
        model, use = fit_cox(use, [variable], penalizer)
        summary = model.summary.loc[variable]
        row.update(
            {
                "hr": math.exp(summary["coef"]),
                "ci_low": math.exp(summary["coef lower 95%"]),
                "ci_high": math.exp(summary["coef upper 95%"]),
                "p": summary["p"],
                "c_index": model.concordance_index_,
            }
        )
        if cohort == "TARGET":
            row["ph_test_p"] = proportional_hazard_test(model, use, time_transform="rank").summary.loc[
                variable, "p"
            ]
    except Exception as exc:
        row["error"] = str(exc)
    return row


def build_efs_tables(cohorts: dict[str, pd.DataFrame]) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows = []
    sensitivity = []
    for cohort in ["TARGET", "Single-center historical"]:
        data = cohorts[cohort]
        for variable in ["candidate_score", "candidate_score_high"]:
            rows.append(cox_row(data, variable, cohort, RIDGE_EFS))
            for penalty in [0.0, 0.01, 0.05, 0.1, 0.5]:
                sensitivity.append(cox_row(data, variable, cohort, penalty))
    return pd.DataFrame(rows), pd.DataFrame(sensitivity)


def binary_association(data: pd.DataFrame, exposure: str, outcome: str, cohort: str) -> dict[str, object]:
    use = data[[exposure, outcome]].dropna().astype(int)
    table = pd.crosstab(use[exposure], use[outcome]).reindex(index=[0, 1], columns=[0, 1], fill_value=0)
    conditional = contingency.odds_ratio(table.to_numpy(), kind="conditional")
    interval = conditional.confidence_interval(confidence_level=0.95)
    return {
        "cohort": cohort,
        "outcome": outcome,
        "exposure": exposure,
        "n": len(use),
        "outcomes": int(use[outcome].sum()),
        "exposed_n": int(use[exposure].sum()),
        "exposed_outcomes": int(((use[exposure] == 1) & (use[outcome] == 1)).sum()),
        "conditional_odds_ratio": conditional.statistic,
        "exact_ci_low": interval.low,
        "exact_ci_high": interval.high,
        "fisher_p": fisher_exact(table.to_numpy()).pvalue,
    }


def build_mrd_table(cohorts: dict[str, pd.DataFrame]) -> pd.DataFrame:
    rows = []
    for cohort, outcome in [
        ("TARGET", "mrd_delayed"),
        ("Single-center historical", "mrd_delayed"),
        ("Single-center contemporary", "mrd_delayed"),
        ("Single-center contemporary", "d46_persistent"),
    ]:
        for exposure in ["candidate_score_high", "gain16", "abnormal_report_tokens_ge3", "del6q"]:
            rows.append(binary_association(cohorts[cohort], exposure, outcome, cohort))
    return pd.DataFrame(rows)


def fit_c_index(train: pd.DataFrame, test: pd.DataFrame, variables: list[str]) -> tuple[float, float, CoxPHFitter]:
    model, _ = fit_cox(train, variables, RIDGE_INCREMENTAL)
    train_risk = model.predict_partial_hazard(train[variables]).to_numpy()
    test_risk = model.predict_partial_hazard(test[variables]).to_numpy()
    return (
        concordance_index(train["time"], -train_risk, train["event"]),
        concordance_index(test["time"], -test_risk, test["event"]),
        model,
    )


def matched_incremental_analysis(data: pd.DataFrame, cohort: str) -> tuple[list[dict[str, object]], dict[str, object]]:
    use = data[["time", "event", "mrd", "candidate_score"]].dropna().copy().reset_index(drop=True)
    use["log_mrd"] = np.log10(use["mrd"] + MRD_OFFSET)
    variables = {"MRD only": ["log_mrd"], "MRD + candidate score": ["log_mrd", "candidate_score"]}
    fitted: dict[str, tuple[float, float, CoxPHFitter]] = {
        name: fit_c_index(use, use, model_vars) for name, model_vars in variables.items()
    }
    rows = []
    for name, model_vars in variables.items():
        model = fitted[name][2]
        row = {
            "cohort": cohort,
            "model": name,
            "n": len(use),
            "events": int(use["event"].sum()),
            "penalizer": RIDGE_INCREMENTAL,
            "apparent_c_index": fitted[name][0],
            "candidate_score_hr": np.nan,
            "candidate_score_ci_low": np.nan,
            "candidate_score_ci_high": np.nan,
            "candidate_score_p": np.nan,
        }
        if "candidate_score" in model_vars:
            summary = model.summary.loc["candidate_score"]
            row.update(
                {
                    "candidate_score_hr": math.exp(summary["coef"]),
                    "candidate_score_ci_low": math.exp(summary["coef lower 95%"]),
                    "candidate_score_ci_high": math.exp(summary["coef upper 95%"]),
                    "candidate_score_p": summary["p"],
                }
            )
        rows.append(row)

    rng = np.random.default_rng(SEED)
    optimism = {name: [] for name in variables}
    test_delta = []
    failed = 0
    for _ in range(N_BOOT):
        sample = use.iloc[rng.integers(0, len(use), len(use))].reset_index(drop=True)
        if sample["event"].sum() < 3 or sample["log_mrd"].nunique() < 2 or sample["candidate_score"].nunique() < 2:
            failed += 1
            continue
        try:
            boot = {}
            for name, model_vars in variables.items():
                boot[name] = fit_c_index(sample, use, model_vars)
                optimism[name].append(boot[name][0] - boot[name][1])
            test_delta.append(boot["MRD + candidate score"][1] - boot["MRD only"][1])
        except Exception:
            failed += 1

    for row in rows:
        values = np.asarray(optimism[row["model"]])
        row["bootstrap_requested"] = N_BOOT
        row["bootstrap_valid"] = len(values)
        row["bootstrap_failed"] = failed
        row["mean_optimism"] = float(np.mean(values))
        row["optimism_corrected_c_index"] = row["apparent_c_index"] - row["mean_optimism"]

    by_model = {row["model"]: row for row in rows}
    delta = {
        "cohort": cohort,
        "n": len(use),
        "events": int(use["event"].sum()),
        "apparent_delta_c": (
            by_model["MRD + candidate score"]["apparent_c_index"]
            - by_model["MRD only"]["apparent_c_index"]
        ),
        "optimism_corrected_delta_c": (
            by_model["MRD + candidate score"]["optimism_corrected_c_index"]
            - by_model["MRD only"]["optimism_corrected_c_index"]
        ),
        "paired_bootstrap_test_delta_median": float(np.median(test_delta)),
        "paired_bootstrap_test_delta_ci_low": float(np.quantile(test_delta, 0.025)),
        "paired_bootstrap_test_delta_ci_high": float(np.quantile(test_delta, 0.975)),
        "bootstrap_valid": len(test_delta),
        "interpretation": "conditional internal assessment of a fixed post hoc score",
    }
    return rows, delta


def build_incremental_tables(cohorts: dict[str, pd.DataFrame]) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows = []
    deltas = []
    for cohort in ["TARGET", "Single-center historical"]:
        cohort_rows, delta = matched_incremental_analysis(cohorts[cohort], cohort)
        rows.extend(cohort_rows)
        deltas.append(delta)
    return pd.DataFrame(rows), pd.DataFrame(deltas)


def build_source_disclosure() -> pd.DataFrame:
    screen = pd.read_csv(SCREEN).set_index("feature")
    specifications = [
        (
            "gain16",
            "Whole-chromosome +16 in the diagnostic karyotype report",
            2,
            "gain_chr_16",
            "Selected after TARGET EFS screening; assigned two points heuristically because it was the strongest screened karyotype signal and was rare.",
        ),
        (
            "abnormal_report_tokens_ge3",
            "At least three regex-detected structural or numerical abnormality tokens in the complete reported karyotype string",
            1,
            "abnormal_tokens_ge3",
            "Selected after TARGET EFS screening as a report-burden candidate; this is not a standard complex-karyotype classification.",
        ),
        (
            "del6q",
            "Reported deletion involving the long arm of chromosome 6",
            1,
            "del_6q",
            "Selected after TARGET review because of biological interest and availability in routine reports; TARGET EFS evidence was not significant.",
        ),
    ]
    rows = []
    for variable, definition, points, feature, rationale in specifications:
        result = screen.loc[feature]
        rows.append(
            {
                "component": variable,
                "operational_definition": definition,
                "points": points,
                "selection_source": "Outcome-informed TARGET exploratory screen",
                "target_screen_n": int(result["univ_n"]),
                "target_events": int(result["univ_events"]),
                "target_univariable_hr": result["univ_hr"],
                "target_univariable_p": result["univ_cox_p"],
                "target_univariable_bh_fdr": result["univ_cox_fdr"],
                "weight_rationale": rationale,
                "prespecification_status": "Locked only before single-center assessment; not prospectively prespecified before TARGET analysis",
            }
        )
    rows.append(
        {
            "component": "candidate_score_high",
            "operational_definition": "Total score >=1; any component is sufficient",
            "points": np.nan,
            "selection_source": "Heuristic threshold selected after TARGET review",
            "weight_rationale": "Chosen for sensitivity and feasibility; the two-point +16 weight does not affect this binary classification.",
            "prespecification_status": "Locked only before single-center assessment; not prospectively prespecified before TARGET analysis",
        }
    )
    return pd.DataFrame(rows)


def build_variable_dictionary() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "variable": "karyotype_evaluable",
                "uniform_definition": "Non-empty reported karyotype without an explicit failed, unknown, or no-metaphase statement",
                "handling": "Unevaluable reports are missing, never score zero",
            },
            {
                "variable": "abnormal_report_token_n",
                "uniform_definition": "Count of all regex-detected structural and numerical abnormality tokens in the complete reported string",
                "handling": "Includes reported t(12;21) and repeated tokens across clones; report-burden measure, not unique-lesion count",
            },
            {
                "variable": "gain16",
                "uniform_definition": "Whole-chromosome +16 token",
                "handling": "Two candidate-score points",
            },
            {
                "variable": "abnormal_report_tokens_ge3",
                "uniform_definition": "abnormal_report_token_n >=3",
                "handling": "One candidate-score point; not labelled complex karyotype",
            },
            {
                "variable": "del6q",
                "uniform_definition": "Deletion token involving chromosome 6q",
                "handling": "One candidate-score point",
            },
            {
                "variable": "candidate_score",
                "uniform_definition": "2*gain16 + abnormal_report_tokens_ge3 + del6q",
                "handling": "Range 0-4; fixed after TARGET-informed construction",
            },
            {
                "variable": "candidate_score_high",
                "uniform_definition": "candidate_score >=1",
                "handling": "Any component qualifies; exploratory binary threshold",
            },
            {
                "variable": "log_mrd",
                "uniform_definition": "log10(cohort-specific MRD percentage + 0.0001)",
                "handling": "TARGET Day 29 and historical Day 19 analyzed separately, never pooled",
            },
        ]
    )


def build_figures(prevalence: pd.DataFrame, incremental: pd.DataFrame, deltas: pd.DataFrame) -> None:
    selected = prevalence[prevalence["variable"] == "candidate_score_high"].copy()
    selected["percent"] = 100 * selected["prevalence"]
    selected["low"] = 100 * selected["ci_low"]
    selected["high"] = 100 * selected["ci_high"]
    colors = ["#0072B2", "#D55E00", "#009E73"]
    fig, ax = plt.subplots(figsize=(8.5, 5.2))
    x = np.arange(len(selected))
    ax.bar(x, selected["percent"], color=colors, width=0.65)
    ax.errorbar(
        x,
        selected["percent"],
        yerr=[selected["percent"] - selected["low"], selected["high"] - selected["percent"]],
        fmt="none",
        ecolor="black",
        capsize=5,
        linewidth=1.2,
    )
    for i, row in enumerate(selected.itertuples()):
        ax.text(i, row.percent + 4, f"{row.positive_n}/{row.evaluable_n}", ha="center", fontsize=10)
    ax.set_xticks(x, selected["cohort"])
    ax.set_ylabel("High candidate-score prevalence (%)")
    ax.set_title("Cross-cohort applicability and ascertainment heterogeneity")
    ax.set_ylim(0, max(70, selected["high"].max() + 10))
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(FIGURE / "Figure_2_score_applicability.png", dpi=300)
    fig.savefig(FIGURE / "Figure_2_score_applicability.pdf")
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.8), gridspec_kw={"width_ratios": [1.5, 1]})
    cohorts = ["TARGET", "Single-center historical"]
    x = np.arange(len(cohorts))
    width = 0.19
    offsets = [-1.5 * width, -0.5 * width, 0.5 * width, 1.5 * width]
    series = [
        ("MRD apparent", "MRD only", "apparent_c_index", "#56B4E9"),
        ("MRD corrected", "MRD only", "optimism_corrected_c_index", "#0072B2"),
        ("MRD + score apparent", "MRD + candidate score", "apparent_c_index", "#F0C36E"),
        ("MRD + score corrected", "MRD + candidate score", "optimism_corrected_c_index", "#E69F00"),
    ]
    for offset, (label, model, value, color) in zip(offsets, series):
        values = [
            incremental.loc[(incremental["cohort"] == cohort) & (incremental["model"] == model), value].iloc[0]
            for cohort in cohorts
        ]
        axes[0].bar(x + offset, values, width=width, color=color, label=label)
    axes[0].axhline(0.5, color="black", linestyle="--", linewidth=1)
    axes[0].set_xticks(
        x,
        [
            f"TARGET\nn=141, events=24",
            f"Single-center historical\nn=56, events=5",
        ],
    )
    axes[0].set_ylabel("Harrell C-index")
    axes[0].set_ylim(0.45, 0.78)
    axes[0].set_title("Matched complete-case models")
    axes[0].legend(frameon=False, fontsize=8, ncol=2)
    axes[0].spines[["top", "right"]].set_visible(False)

    y = np.arange(len(deltas))
    error = np.vstack(
        [
            deltas["paired_bootstrap_test_delta_median"] - deltas["paired_bootstrap_test_delta_ci_low"],
            deltas["paired_bootstrap_test_delta_ci_high"] - deltas["paired_bootstrap_test_delta_median"],
        ]
    )
    axes[1].errorbar(
        deltas["paired_bootstrap_test_delta_median"],
        y,
        xerr=error,
        fmt="o",
        color="#0072B2",
        capsize=4,
    )
    axes[1].axvline(0, color="black", linestyle="--", linewidth=1)
    axes[1].set_yticks(y, [f"{row.cohort}\nn={row.n}, events={row.events}" for row in deltas.itertuples()])
    axes[1].set_xlabel("Paired bootstrap test ΔC\n(MRD + score minus MRD)")
    axes[1].set_title("Uncertainty in incremental discrimination")
    axes[1].spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(FIGURE / "Figure_3_matched_cindex.png", dpi=300)
    fig.savefig(FIGURE / "Figure_3_matched_cindex.pdf")
    plt.close(fig)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run harmonized karyotype scoring and matched-sample C-index analyses."
    )
    parser.add_argument(
        "--input-root",
        type=Path,
        default=INPUT,
        help="Directory containing the TARGET CSV, historical CSV, and restricted single-center workbook.",
    )
    parser.add_argument(
        "--screen-file",
        type=Path,
        default=SCREEN,
        help="TARGET exploratory karyotype-screen CSV used for score-source disclosure.",
    )
    parser.add_argument(
        "--output-root",
        type=Path,
        default=ROOT,
        help="Output directory; tables/ and figures/ subdirectories will be created.",
    )
    return parser.parse_args()


def main() -> None:
    global INPUT, SCREEN, ROOT, TABLE, FIGURE
    args = parse_args()
    INPUT = args.input_root.resolve()
    SCREEN = args.screen_file.resolve()
    ROOT = args.output_root.resolve()
    TABLE = ROOT / "tables"
    FIGURE = ROOT / "figures"
    TABLE.mkdir(parents=True, exist_ok=True)
    FIGURE.mkdir(parents=True, exist_ok=True)

    cohorts = load_cohorts()
    prevalence = build_prevalence_table(cohorts)
    agreement = build_historical_agreement(cohorts["Single-center historical"])
    efs, sensitivity = build_efs_tables(cohorts)
    mrd = build_mrd_table(cohorts)
    incremental, deltas = build_incremental_tables(cohorts)
    source = build_source_disclosure()
    definitions = build_variable_dictionary()

    source.to_csv(TABLE / "Table_1_score_source_disclosure.csv", index=False)
    definitions.to_csv(TABLE / "Table_2_uniform_variable_definitions.csv", index=False)
    prevalence.to_csv(TABLE / "Table_3_component_prevalence.csv", index=False)
    efs.to_csv(TABLE / "Table_4_harmonized_efs_associations.csv", index=False)
    incremental.to_csv(TABLE / "Table_5_matched_incremental_models.csv", index=False)
    deltas.to_csv(TABLE / "Table_6_paired_cindex_delta.csv", index=False)
    agreement.to_csv(TABLE / "Supplementary_Table_historical_definition_agreement.csv", index=False)
    sensitivity.to_csv(TABLE / "Supplementary_Table_penalty_sensitivity.csv", index=False)
    mrd.to_csv(TABLE / "Supplementary_Table_harmonized_MRD_associations.csv", index=False)
    build_figures(prevalence, incremental, deltas)

    audit = [
        "# Harmonized reanalysis audit",
        "",
        "- One parser was applied to raw TARGET, historical single-center, and contemporary single-center karyotype strings.",
        "- Historical raw-to-analytic row alignment passed exact checks for source, diagnosis date, sex, age, and WBC.",
        "- The report-token burden variable is explicitly not a standard complex-karyotype classification.",
        f"- Matched complete-case C-index models used ridge penalizer {RIDGE_INCREMENTAL} and {N_BOOT} paired bootstrap resamples.",
        "- TARGET bootstrap results are conditional internal assessment of a fixed, post hoc candidate score; feature selection was not repeated.",
        "",
    ]
    for row in deltas.itertuples():
        audit.append(
            f"- {row.cohort}: matched n={row.n}, events={row.events}; apparent delta C={row.apparent_delta_c:.3f}; "
            f"optimism-corrected delta C={row.optimism_corrected_delta_c:.3f}; paired bootstrap test delta "
            f"median={row.paired_bootstrap_test_delta_median:.3f} "
            f"(95% interval {row.paired_bootstrap_test_delta_ci_low:.3f} to {row.paired_bootstrap_test_delta_ci_high:.3f})."
        )
    (ROOT / "HARMONIZED_REANALYSIS_AUDIT.md").write_text("\n".join(audit), encoding="utf-8")


if __name__ == "__main__":
    main()
