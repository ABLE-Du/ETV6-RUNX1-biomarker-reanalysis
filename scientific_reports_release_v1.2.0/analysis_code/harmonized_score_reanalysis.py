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
STRUCTURAL_DETAIL = re.compile(r"(?:DEL|ADD|T|INV|DUP|DER|DIC|IDIC|I|R)\([^,\[/]+", re.I)
NUMERICAL_DETAIL = re.compile(r"(?:^|[,/])((?:\+|-)(?:\d+|X|Y|MAR))(?=[,/\[(]|$)", re.I)
T12_21_TOKEN = re.compile(r"T\(12;21\*?\)", re.I)
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


def parse_mrd_value(value: object) -> float:
    if pd.isna(value):
        return np.nan
    if isinstance(value, (int, float, np.integer, np.floating)):
        return float(value)
    text = str(value).strip().replace("%", "")
    if not text:
        return np.nan
    if text.startswith("<"):
        try:
            return float(text[1:]) / 2
        except ValueError:
            return 0.0
    match = re.search(r"[-+]?\d+(?:\.\d+)?", text)
    return float(match.group()) if match else np.nan


def is_t12_21_token(text: str, start: int) -> bool:
    return bool(T12_21_TOKEN.match(text[start:]))


def parse_karyotype(value: object) -> dict[str, int | float]:
    text = normalize_karyotype(value)
    evaluable = int(bool(text) and not FAILED_KARYOTYPE.search(text))
    structural_matches = list(STRUCTURAL_TOKEN.finditer(text))
    structural_n = len(structural_matches)
    numerical_n = len(NUMERICAL_TOKEN.findall(text))
    abnormal_n = structural_n + numerical_n
    structural_excluding_t12_n = sum(
        1 for match in structural_matches if not is_t12_21_token(text, match.start())
    )
    abnormal_excluding_t12_n = structural_excluding_t12_n + numerical_n
    unique_lesions = {
        match.group(0).upper()
        for match in STRUCTURAL_DETAIL.finditer(text)
    } | {
        match.group(1).upper()
        for match in NUMERICAL_DETAIL.finditer(text)
    }
    unique_lesion_n = len(unique_lesions)
    parsed = {
        "karyotype_evaluable": evaluable,
        "structural_token_n": structural_n,
        "numerical_token_n": numerical_n,
        "abnormal_report_token_n": abnormal_n,
        "abnormal_report_token_n_excluding_t12_21": abnormal_excluding_t12_n,
        "abnormal_unique_lesion_n": unique_lesion_n,
        "gain16": int(bool(re.search(r"(?:^|[,/])\+16(?=[,/\[(]|$)", text))),
        "abnormal_report_tokens_ge3": int(abnormal_n >= 3),
        "abnormal_report_tokens_ge3_excluding_t12_21": int(abnormal_excluding_t12_n >= 3),
        "abnormal_unique_lesions_ge3": int(unique_lesion_n >= 3),
        "del6q": int(bool(re.search(r"DEL\(6\)\([^)]*Q|DEL\(6Q|DEL6Q", text))),
    }
    if not evaluable:
        for key in [
            "structural_token_n",
            "numerical_token_n",
            "abnormal_report_token_n",
            "abnormal_report_token_n_excluding_t12_21",
            "abnormal_unique_lesion_n",
            "gain16",
            "abnormal_report_tokens_ge3",
            "abnormal_report_tokens_ge3_excluding_t12_21",
            "abnormal_unique_lesions_ge3",
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
    target["mrd"] = target["MRD_PERCENT_DAY_29"].map(parse_mrd_value)
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
    historical["mrd"] = historical["d19_mrd"].map(parse_mrd_value)
    historical["mrd_delayed"] = (historical["mrd"] >= 0.01).astype(float)
    historical.loc[historical["mrd"].isna(), "mrd_delayed"] = np.nan

    registry = pd.read_excel(workbook, sheet_name="登记数据_合并")
    contemporary = registry.loc[registry.apply(contains_etv6, axis=1)].copy()
    contemporary = add_score(contemporary, "核型")
    contemporary["mrd"] = contemporary["D19MRD(%)"].map(parse_mrd_value)
    contemporary["mrd_delayed"] = (contemporary["mrd"] >= 0.01).astype(float)
    contemporary.loc[contemporary["mrd"].isna(), "mrd_delayed"] = np.nan
    contemporary["d46_mrd"] = contemporary["D46MRD(%)"].map(parse_mrd_value)
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


def fmt_n_pct(count: int, denominator: int) -> str:
    if denominator == 0:
        return "NA"
    return f"{count}/{denominator} ({100 * count / denominator:.1f}%)"


def fmt_median_iqr(values: pd.Series) -> str:
    use = pd.to_numeric(values, errors="coerce").dropna()
    if use.empty:
        return "Not available"
    q1, median, q3 = use.quantile([0.25, 0.5, 0.75])
    return f"{median:.1f} ({q1:.1f}-{q3:.1f})"


def fmt_counts(values: pd.Series, denominator: int | None = None, max_levels: int = 6) -> str:
    use = values.fillna("Missing").astype(str)
    if use.empty:
        return "Not available"
    denominator = denominator or len(use)
    counts = use.value_counts()
    parts = [f"{level}: {count} ({100 * count / denominator:.1f}%)" for level, count in counts.iloc[:max_levels].items()]
    if len(counts) > max_levels:
        parts.append(f"Other levels: {int(counts.iloc[max_levels:].sum())} ({100 * counts.iloc[max_levels:].sum() / denominator:.1f}%)")
    return "; ".join(parts)


def diagnosis_year_range(values: pd.Series) -> str:
    dates = pd.to_datetime(values, errors="coerce").dropna()
    if dates.empty:
        return "Not available"
    years = dates.dt.year
    return f"{int(years.min())}-{int(years.max())}"


def build_cohort_characteristics(cohorts: dict[str, pd.DataFrame]) -> pd.DataFrame:
    target = cohorts["TARGET"]
    historical = cohorts["Single-center historical"]
    contemporary = cohorts["Single-center contemporary"]

    risk_contemporary = contemporary["最终危险度_标准"].fillna(contemporary["最终危险度"]) if "最终危险度_标准" in contemporary else pd.Series(dtype=object)
    rows = [
        {
            "characteristic": "Patients, n",
            "TARGET": str(len(target)),
            "Single-center historical": str(len(historical)),
            "Single-center contemporary": str(len(contemporary)),
        },
        {
            "characteristic": "Diagnosis years",
            "TARGET": "Not available in analysis file",
            "Single-center historical": diagnosis_year_range(historical["diagnosis_date"]),
            "Single-center contemporary": diagnosis_year_range(contemporary["诊断日期"]),
        },
        {
            "characteristic": "Treatment protocol or era",
            "TARGET": fmt_counts(target["PROTOCOL"], len(target)),
            "Single-center historical": fmt_counts(historical["来源表"], len(historical)),
            "Single-center contemporary": fmt_counts(contemporary["来源表"], len(contemporary)),
        },
        {
            "characteristic": "Male sex",
            "TARGET": fmt_n_pct(int(target["SEX"].eq("Male").sum()), int(target["SEX"].notna().sum())),
            "Single-center historical": fmt_n_pct(int(historical["性别"].eq("男").sum()), int(historical["性别"].notna().sum())),
            "Single-center contemporary": fmt_n_pct(int(contemporary["性别"].eq("男").sum()), int(contemporary["性别"].notna().sum())),
        },
        {
            "characteristic": "Age at diagnosis, years, median (IQR)",
            "TARGET": fmt_median_iqr(target["AGE"]),
            "Single-center historical": fmt_median_iqr(historical["age_years"]),
            "Single-center contemporary": fmt_median_iqr(contemporary["诊断年龄_年"]),
        },
        {
            "characteristic": "WBC, x10^9/L, median (IQR)",
            "TARGET": fmt_median_iqr(target["WBC"]),
            "Single-center historical": fmt_median_iqr(historical["wbc"]),
            "Single-center contemporary": fmt_median_iqr(contemporary["WBC（×109）"] if "WBC（×109）" in contemporary else pd.Series(dtype=float)),
        },
        {
            "characteristic": "Clinical risk group",
            "TARGET": "Not available in analysis file",
            "Single-center historical": fmt_counts(historical["risk_group"], len(historical)),
            "Single-center contemporary": fmt_counts(risk_contemporary, len(contemporary)),
        },
        {
            "characteristic": "Evaluable karyotype",
            "TARGET": fmt_n_pct(int(target["karyotype_evaluable"].sum()), len(target)),
            "Single-center historical": fmt_n_pct(int(historical["karyotype_evaluable"].sum()), len(historical)),
            "Single-center contemporary": fmt_n_pct(int(contemporary["karyotype_evaluable"].sum()), len(contemporary)),
        },
        {
            "characteristic": "MRD available at primary cohort time point",
            "TARGET": fmt_n_pct(int(target["mrd"].notna().sum()), len(target)),
            "Single-center historical": fmt_n_pct(int(historical["mrd"].notna().sum()), len(historical)),
            "Single-center contemporary": fmt_n_pct(int(contemporary["mrd"].notna().sum()), len(contemporary)),
        },
        {
            "characteristic": "D46 MRD available",
            "TARGET": "Not applicable",
            "Single-center historical": "Not applicable",
            "Single-center contemporary": fmt_n_pct(int(contemporary["d46_mrd"].notna().sum()), len(contemporary)),
        },
        {
            "characteristic": "EFS analyzable patients/events",
            "TARGET": f"{int(target[['time', 'event']].dropna().shape[0])}/{int(target['event'].fillna(0).sum())}",
            "Single-center historical": f"{int(historical[['time', 'event']].dropna().shape[0])}/{int(historical['event'].fillna(0).sum())}",
            "Single-center contemporary": "Not analyzed",
        },
        {
            "characteristic": "EFS follow-up, months, median (IQR)",
            "TARGET": fmt_median_iqr(target["time"] / 30.4375),
            "Single-center historical": fmt_median_iqr(historical["time"]),
            "Single-center contemporary": "Not analyzed",
        },
        {
            "characteristic": "High candidate score among evaluable karyotypes",
            "TARGET": fmt_n_pct(int(target["candidate_score_high"].fillna(0).sum()), int(target["candidate_score_high"].notna().sum())),
            "Single-center historical": fmt_n_pct(int(historical["candidate_score_high"].fillna(0).sum()), int(historical["candidate_score_high"].notna().sum())),
            "Single-center contemporary": fmt_n_pct(int(contemporary["candidate_score_high"].fillna(0).sum()), int(contemporary["candidate_score_high"].notna().sum())),
        },
    ]
    return pd.DataFrame(rows)


def build_efs_definition_table(cohorts: dict[str, pd.DataFrame]) -> pd.DataFrame:
    target = cohorts["TARGET"].copy()
    historical = cohorts["Single-center historical"].copy()
    target_events = (
        target.loc[target["time"].notna(), "FIRST_EVENT"]
        .fillna("Censored/no first event recorded")
        .replace({"Censored": "Censored/no first event recorded"})
    )
    historical_events = (
        historical.get("EFS事件", pd.Series(dtype=object))
        .fillna("Missing")
        .replace(
            {
                "存活": "Censored/event-free alive",
                "复发": "Relapse",
                "骨髓复发": "Bone marrow relapse",
                "感染死亡": "Infection-related death",
                "骨髓复发、睾丸复发": "Bone marrow and testicular relapse",
            }
        )
    )
    return pd.DataFrame(
        [
            {
                "cohort": "TARGET",
                "time_origin": "Diagnosis/enrolment as encoded in TARGET time-to-event variable",
                "time_scale": "Days in source analysis; converted to months only for descriptive follow-up",
                "efs_event_definition": "FIRST_EVENT of relapse, death, or second malignant neoplasm counted as EFS event",
                "censoring_rule": "Patients without a qualifying first event were censored at the recorded follow-up/time variable",
                "event_composition": fmt_counts(target_events, int(target["time"].notna().sum()), max_levels=8),
                "analysis_population": f"{int(target[['time', 'event']].dropna().shape[0])} with EFS time/event; {int(target['event'].fillna(0).sum())} events overall; matched MRD-score model n=141/events=24",
            },
            {
                "cohort": "Single-center historical",
                "time_origin": "Diagnosis date",
                "time_scale": "Months",
                "efs_event_definition": "Source EFS event variable; non-survival entries counted as EFS events",
                "censoring_rule": "Patients recorded as alive/event-free were censored at last follow-up or calculated EFS time",
                "event_composition": fmt_counts(historical_events, len(historical_events), max_levels=8),
                "analysis_population": f"{int(historical[['time', 'event']].dropna().shape[0])} with EFS time/event; {int(historical['event'].fillna(0).sum())} events overall; matched MRD-score model n=56/events=5",
            },
            {
                "cohort": "Single-center contemporary",
                "time_origin": "Not applicable",
                "time_scale": "Not applicable",
                "efs_event_definition": "EFS not analyzed in this cohort",
                "censoring_rule": "Not applicable",
                "event_composition": "Not applicable",
                "analysis_population": "Contemporary cohort used for descriptive D19 and D46 MRD analyses only",
            },
        ]
    )


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


def token_audit(value: object) -> dict[str, object]:
    text = normalize_karyotype(value)
    structural_tokens = [match.group(0).upper() for match in STRUCTURAL_DETAIL.finditer(text)]
    numerical_tokens = [match.group(1).upper() for match in NUMERICAL_DETAIL.finditer(text)]
    t12_tokens = [
        match.group(0).upper()
        for match in STRUCTURAL_DETAIL.finditer(text)
        if T12_21_TOKEN.match(match.group(0))
    ]
    return {
        "normalized_karyotype": text,
        "structural_tokens": "; ".join(structural_tokens),
        "numerical_tokens": "; ".join(numerical_tokens),
        "t12_21_tokens": "; ".join(t12_tokens),
        "unique_abnormal_lesions": "; ".join(sorted(set(structural_tokens + numerical_tokens))),
    }


def build_parser_audit_tables(cohorts: dict[str, pd.DataFrame]) -> tuple[pd.DataFrame, pd.DataFrame]:
    target = cohorts["TARGET"].copy()
    categories = [
        ("normal_or_no_detected_abnormality", target["abnormal_report_token_n"].fillna(0).eq(0)),
        ("report_token_burden_ge3", target["abnormal_report_tokens_ge3"].eq(1)),
        ("gain16_component", target["gain16"].eq(1)),
        ("del6q_component", target["del6q"].eq(1)),
        ("contains_t12_21_token", target["KARYOTYPE"].map(lambda x: bool(T12_21_TOKEN.search(normalize_karyotype(x))))),
    ]
    example_rows = []
    seen_karyotypes: set[str] = set()
    for category, mask in categories:
        candidates = target.loc[mask & target["karyotype_evaluable"].eq(1)].copy()
        if candidates.empty:
            continue
        row = candidates.iloc[0]
        if str(row["KARYOTYPE"]) in seen_karyotypes and len(candidates) > 1:
            row = candidates.iloc[1]
        seen_karyotypes.add(str(row["KARYOTYPE"]))
        audit = token_audit(row["KARYOTYPE"])
        example_rows.append(
            {
                "example_category": category,
                "source": "TARGET public karyotype string; no patient identifier included",
                "raw_karyotype": row["KARYOTYPE"],
                **audit,
                "abnormal_report_token_n": row["abnormal_report_token_n"],
                "abnormal_report_token_n_excluding_t12_21": row["abnormal_report_token_n_excluding_t12_21"],
                "abnormal_unique_lesion_n": row["abnormal_unique_lesion_n"],
                "gain16": row["gain16"],
                "abnormal_report_tokens_ge3": row["abnormal_report_tokens_ge3"],
                "del6q": row["del6q"],
                "candidate_score": row["candidate_score"],
            }
        )

    summary_rows = []
    for cohort, data in cohorts.items():
        use = data.loc[data["karyotype_evaluable"].eq(1)].copy()
        n = len(use)
        t12_any = data.apply(
            lambda row: bool(T12_21_TOKEN.search(normalize_karyotype(row.get("KARYOTYPE", row.get("raw_karyotype", row.get("核型", "")))))),
            axis=1,
        )
        score_excluding_t12 = (
            2 * use["gain16"] + use["abnormal_report_tokens_ge3_excluding_t12_21"] + use["del6q"]
        )
        score_unique = 2 * use["gain16"] + use["abnormal_unique_lesions_ge3"] + use["del6q"]
        summary_rows.append(
            {
                "cohort": cohort,
                "evaluable_karyotype_n": n,
                "abnormal_report_token_n_median_iqr": fmt_median_iqr(use["abnormal_report_token_n"]),
                "abnormal_report_token_n_excluding_t12_21_median_iqr": fmt_median_iqr(use["abnormal_report_token_n_excluding_t12_21"]),
                "abnormal_unique_lesion_n_median_iqr": fmt_median_iqr(use["abnormal_unique_lesion_n"]),
                "contains_t12_21_token_n": int(t12_any.loc[use.index].sum()),
                "primary_high_score_n": int(use["candidate_score_high"].sum()),
                "high_score_excluding_t12_21_n": int((score_excluding_t12 >= 1).sum()),
                "high_score_unique_lesions_n": int((score_unique >= 1).sum()),
                "changed_high_score_excluding_t12_21_n": int(((score_excluding_t12 >= 1).astype(int) != use["candidate_score_high"].astype(int)).sum()),
                "changed_high_score_unique_lesions_n": int(((score_unique >= 1).astype(int) != use["candidate_score_high"].astype(int)).sum()),
            }
        )
    return pd.DataFrame(example_rows), pd.DataFrame(summary_rows)


def score_for_burden(data: pd.DataFrame, burden_col: str) -> pd.Series:
    return 2 * data["gain16"] + data[burden_col] + data["del6q"]


def build_parser_sensitivity(cohorts: dict[str, pd.DataFrame]) -> pd.DataFrame:
    definitions = [
        ("primary_report_tokens", "abnormal_report_tokens_ge3"),
        ("exclude_t12_21_report_tokens", "abnormal_report_tokens_ge3_excluding_t12_21"),
        ("unique_reported_lesions", "abnormal_unique_lesions_ge3"),
    ]
    rows = []
    for definition, burden_col in definitions:
        for cohort, data in cohorts.items():
            use = data.loc[data["karyotype_evaluable"].eq(1)].copy()
            score = score_for_burden(use, burden_col)
            high = score >= 1
            row = {
                "sensitivity_definition": definition,
                "cohort": cohort,
                "burden_variable": burden_col,
                "evaluable_n": len(use),
                "high_score_n": int(high.sum()),
                "high_score_prevalence": float(high.mean()) if len(use) else np.nan,
                "efs_n": np.nan,
                "efs_events": np.nan,
                "continuous_score_hr": np.nan,
                "continuous_score_ci_low": np.nan,
                "continuous_score_ci_high": np.nan,
                "continuous_score_p": np.nan,
            }
            if cohort in {"TARGET", "Single-center historical"}:
                temp = use[["time", "event"]].copy()
                temp["sensitivity_score"] = score
                fitted = cox_row(temp, "sensitivity_score", cohort, RIDGE_EFS)
                row.update(
                    {
                        "efs_n": fitted.get("n"),
                        "efs_events": fitted.get("events"),
                        "continuous_score_hr": fitted.get("hr"),
                        "continuous_score_ci_low": fitted.get("ci_low"),
                        "continuous_score_ci_high": fitted.get("ci_high"),
                        "continuous_score_p": fitted.get("p"),
                    }
                )
            rows.append(row)
    return pd.DataFrame(rows)


def screen_category(feature: str) -> str:
    if feature.startswith("gain_chr") or feature.startswith("loss_chr"):
        return "whole_chromosome_numerical"
    if feature.startswith("del_") or feature.startswith("del_chr"):
        return "deletion"
    if feature.startswith("structural_type"):
        return "structural_event_type"
    if feature.startswith("structural_chr"):
        return "structural_chromosome_involvement"
    if feature in {"structural_ge3", "abnormal_tokens_ge3"}:
        return "abnormality_burden"
    if feature.startswith("translocation"):
        return "specific_translocation"
    if feature in {"hypodiploid_clone", "normal_clone_present", "normal_karyotype_only", "multiple_clones"}:
        return "karyotype_pattern"
    return "other"


def build_screening_universe_tables() -> tuple[pd.DataFrame, pd.DataFrame]:
    screen = pd.read_csv(SCREEN)
    screen.insert(1, "feature_category", screen["feature"].map(screen_category))
    summary = (
        screen.groupby("feature_category")
        .agg(
            screened_feature_n=("feature", "count"),
            min_univariable_p=("univ_cox_p", "min"),
            min_univariable_bh_fdr=("univ_cox_fdr", "min"),
        )
        .reset_index()
        .sort_values(["screened_feature_n", "feature_category"], ascending=[False, True])
    )
    return screen, summary


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
                "variable": "abnormal_report_token_n_excluding_t12_21",
                "uniform_definition": "Same report-token count after excluding the defining t(12;21) token",
                "handling": "Sensitivity analysis for whether the defining fusion should contribute to report burden",
            },
            {
                "variable": "abnormal_unique_lesion_n",
                "uniform_definition": "Number of unique regex-detected structural/numerical lesion strings after duplicate removal",
                "handling": "Sensitivity analysis for repeated clone-level reporting",
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
    cohort_characteristics = build_cohort_characteristics(cohorts)
    efs_definitions = build_efs_definition_table(cohorts)
    prevalence = build_prevalence_table(cohorts)
    agreement = build_historical_agreement(cohorts["Single-center historical"])
    efs, sensitivity = build_efs_tables(cohorts)
    mrd = build_mrd_table(cohorts)
    incremental, deltas = build_incremental_tables(cohorts)
    source = build_source_disclosure()
    definitions = build_variable_dictionary()
    parser_examples, parser_summary = build_parser_audit_tables(cohorts)
    parser_sensitivity = build_parser_sensitivity(cohorts)
    screen_full, screen_summary = build_screening_universe_tables()

    source.to_csv(TABLE / "Table_1_score_source_disclosure.csv", index=False)
    definitions.to_csv(TABLE / "Table_2_uniform_variable_definitions.csv", index=False)
    cohort_characteristics.to_csv(TABLE / "Table_3_cohort_characteristics.csv", index=False)
    prevalence.to_csv(TABLE / "Table_4_component_prevalence.csv", index=False)
    efs.to_csv(TABLE / "Table_5_harmonized_efs_associations.csv", index=False)
    incremental.to_csv(TABLE / "Table_6_matched_incremental_models.csv", index=False)
    deltas.to_csv(TABLE / "Table_7_paired_cindex_delta.csv", index=False)
    efs_definitions.to_csv(TABLE / "Supplementary_Table_EFS_definitions.csv", index=False)
    parser_examples.to_csv(TABLE / "Supplementary_Table_parser_audit_public_examples.csv", index=False)
    parser_summary.to_csv(TABLE / "Supplementary_Table_parser_audit_summary.csv", index=False)
    parser_sensitivity.to_csv(TABLE / "Supplementary_Table_parser_sensitivity.csv", index=False)
    screen_full.to_csv(TABLE / "Supplementary_Table_TARGET_karyotype_screening_universe.csv", index=False)
    screen_summary.to_csv(TABLE / "Supplementary_Table_TARGET_screening_universe_summary.csv", index=False)
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
