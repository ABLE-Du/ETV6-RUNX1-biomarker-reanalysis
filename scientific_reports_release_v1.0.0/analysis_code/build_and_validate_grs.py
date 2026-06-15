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
from scipy.stats import fisher_exact, mannwhitneyu
from statsmodels.stats.proportion import proportion_effectsize
from statsmodels.stats.power import NormalIndPower


ROOT = Path(__file__).resolve().parent
INP = ROOT / "inputs"
OUT = ROOT / "results"
OUT.mkdir(exist_ok=True)

STRUCTURAL = re.compile(r"(?:del|add|t|inv|dup|der|dic|idic|i|r)\(", re.I)
ABNORMAL = re.compile(
    r"(?:del|add|t|inv|dup|der|dic|idic|i|r)\(|(?:^|,)[+-](?:\d+|X|Y|mar)",
    re.I,
)


def parse_karyotype(value: object) -> dict:
    text = "" if pd.isna(value) else str(value).replace(" ", "")
    failed = bool(re.search(r"未做|未查|无分裂|未见分裂|失败|不详|未知|nan|none", text, re.I))
    evaluable = int(bool(text) and not failed)
    abnormal_n = len(ABNORMAL.findall(text))
    return {
        "karyotype_evaluable": evaluable,
        "gain16": int(bool(re.search(r"(?:^|,)\+16(?=,|/|\[|$)", text, re.I))),
        "high_abnormal_burden": int(abnormal_n >= 3),
        "del6q": int(bool(re.search(r"del\(6\)\([^)]*q", text, re.I))),
        "abnormal_token_n": abnormal_n,
    }


def add_score(df: pd.DataFrame, karyotype_col: str) -> pd.DataFrame:
    parsed = pd.DataFrame([parse_karyotype(x) for x in df[karyotype_col]], index=df.index)
    out = pd.concat([df.drop(columns=parsed.columns, errors="ignore").copy(), parsed], axis=1)
    out["grs"] = 2 * out["gain16"] + out["high_abnormal_burden"] + out["del6q"]
    out["grs_high"] = (out["grs"] >= 1).astype(int)
    for col in ["gain16", "high_abnormal_burden", "del6q", "grs", "grs_high"]:
        out.loc[out["karyotype_evaluable"] == 0, col] = np.nan
    return out


def binary_summary(df: pd.DataFrame, exposure: str, outcome: str, cohort: str) -> dict:
    use = df[[exposure, outcome]].dropna()
    table = pd.crosstab(use[exposure].astype(int), use[outcome].astype(int)).reindex(
        index=[0, 1], columns=[0, 1], fill_value=0
    )
    odds, p = fisher_exact(table.to_numpy())
    sensitivity = table.loc[1, 1] / table[1].sum() if table[1].sum() else np.nan
    specificity = table.loc[0, 0] / table[0].sum() if table[0].sum() else np.nan
    return {
        "cohort": cohort,
        "exposure": exposure,
        "outcome": outcome,
        "n": len(use),
        "outcomes": int(use[outcome].sum()),
        "exposed_n": int(use[exposure].sum()),
        "exposed_outcomes": int(((use[exposure] == 1) & (use[outcome] == 1)).sum()),
        "unexposed_outcomes": int(((use[exposure] == 0) & (use[outcome] == 1)).sum()),
        "odds_ratio": odds,
        "fisher_p": p,
        "auc": (sensitivity + specificity) / 2 if np.isfinite(sensitivity) and np.isfinite(specificity) else np.nan,
    }


def cox_summary(df: pd.DataFrame, exposure: str, cohort: str) -> dict:
    use = df[["time", "event", exposure]].dropna()
    result = {"cohort": cohort, "exposure": exposure, "n": len(use), "events": int(use.event.sum())}
    if len(use) < 20 or use.event.sum() < 5 or use[exposure].nunique() < 2:
        return result
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            model = CoxPHFitter(penalizer=0.05)
            model.fit(use, "time", "event")
        row = model.summary.loc[exposure]
        result.update({
            "hr": math.exp(row["coef"]),
            "ci_low": math.exp(row["coef lower 95%"]),
            "ci_high": math.exp(row["coef upper 95%"]),
            "p": row["p"],
            "c_index": model.concordance_index_,
        })
    except Exception as exc:
        result["error"] = str(exc)
    return result


def compare_mrd_models(df: pd.DataFrame, mrd_col: str, cohort: str) -> pd.DataFrame:
    use = df[["time", "event", mrd_col, "grs"]].dropna().copy()
    use["log_mrd"] = np.log10(pd.to_numeric(use[mrd_col], errors="coerce") + 0.0001)
    rows = []
    for variables in [["log_mrd"], ["log_mrd", "grs"]]:
        data = use[["time", "event", *variables]].dropna()
        row = {"cohort": cohort, "model": "+".join(variables), "n": len(data), "events": int(data.event.sum())}
        if len(data) >= 20 and data.event.sum() >= 5:
            try:
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    model = CoxPHFitter(penalizer=0.5).fit(data, "time", "event")
                row["c_index"] = model.concordance_index_
                for variable in variables:
                    row[f"{variable}_hr"] = math.exp(model.params_[variable])
                    row[f"{variable}_p"] = model.summary.loc[variable, "p"]
            except Exception as exc:
                row["error"] = str(exc)
        rows.append(row)
    return pd.DataFrame(rows)


def target_analysis() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    d = pd.read_csv(INP / "target_etv6_positive_patient_derived_variables.csv")
    d = d[d["KARYOTYPE"].notna()].copy()
    d = add_score(d, "KARYOTYPE")
    d["mrd_delayed"] = (pd.to_numeric(d["MRD_PERCENT_DAY_29"], errors="coerce") >= 0.01).astype(float)
    d.loc[d["MRD_PERCENT_DAY_29"].isna(), "mrd_delayed"] = np.nan
    d["time"] = pd.to_numeric(d["time"], errors="coerce")
    d["event"] = pd.to_numeric(d["event"], errors="coerce")
    binary = pd.DataFrame([
        binary_summary(d, "gain16", "mrd_delayed", "TARGET"),
        binary_summary(d, "high_abnormal_burden", "mrd_delayed", "TARGET"),
        binary_summary(d, "del6q", "mrd_delayed", "TARGET"),
        binary_summary(d, "grs_high", "mrd_delayed", "TARGET"),
    ])
    cox = pd.DataFrame([
        cox_summary(d, "gain16", "TARGET"),
        cox_summary(d, "high_abnormal_burden", "TARGET"),
        cox_summary(d, "del6q", "TARGET"),
        cox_summary(d, "grs", "TARGET"),
        cox_summary(d, "grs_high", "TARGET"),
    ])
    return d, binary, cox


def contains_etv6(row: pd.Series, columns: list[str]) -> bool:
    text = " ".join(str(row.get(c, "")) for c in columns)
    return bool(re.search(r"ETV6|TEL.?AML|12;21", text, re.I))


def single_center_analysis() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    workbook = INP / "ALL2020-2025_single_center_cleaned_pseudonymized.xlsx"
    registry = pd.read_excel(workbook, sheet_name=2)
    fusion_cols = ["融合基因", "其它融合基因", "融合基因(RNAseq报告)", "分子分型(RNAseq报告)", "FISH", "FISH其它"]
    registry["etv6_positive"] = registry.apply(lambda r: contains_etv6(r, fusion_cols), axis=1)
    reg = registry.loc[registry["etv6_positive"]].copy()
    reg = add_score(reg, "核型")
    reg["d19"] = pd.to_numeric(reg["D19MRD(%)"], errors="coerce")
    reg["d46"] = pd.to_numeric(reg["D46MRD(%)"], errors="coerce")
    reg["d19_delayed"] = (reg["d19"] >= 0.01).astype(float)
    reg.loc[reg["d19"].isna(), "d19_delayed"] = np.nan
    reg["d46_persistent"] = (reg["d46"] >= 0.01).astype(float)
    reg.loc[reg["d46"].isna(), "d46_persistent"] = np.nan
    reg["clearance_class"] = np.where(
        reg["d46_persistent"] == 1, 2,
        np.where(reg["d19_delayed"] == 1, 1, np.where(reg["d19_delayed"] == 0, 0, np.nan)),
    )

    hist = pd.read_csv(INP / "historical_outcome_analytic_cohort_pseudonymized.csv")
    hist["karyotype_evaluable"] = (~hist["karyotype_class"].eq("Uninformative")).astype(int)
    hist["gain16"] = pd.to_numeric(hist["gain16"], errors="coerce")
    hist["del6q"] = pd.to_numeric(hist["del6q"], errors="coerce")
    hist["high_abnormal_burden"] = pd.to_numeric(hist["complex_karyotype"], errors="coerce")
    hist["grs"] = 2 * hist["gain16"] + hist["high_abnormal_burden"] + hist["del6q"]
    hist["grs_high"] = (hist["grs"] >= 1).astype(float)
    for col in ["gain16", "del6q", "high_abnormal_burden", "grs", "grs_high"]:
        hist.loc[hist["karyotype_evaluable"] == 0, col] = np.nan
    hist["time"] = pd.to_numeric(hist["efs_months"], errors="coerce")
    hist["event"] = pd.to_numeric(hist["efs_event"], errors="coerce")
    hist["d19"] = pd.to_numeric(hist["d19_mrd"], errors="coerce")
    hist["d19_delayed"] = (hist["d19"] >= 0.01).astype(float)
    hist.loc[hist["d19"].isna(), "d19_delayed"] = np.nan

    binary_rows = []
    for cohort, data, outcome in [
        ("single_center_contemporary", reg, "d19_delayed"),
        ("single_center_contemporary", reg, "d46_persistent"),
        ("single_center_historical", hist, "d19_delayed"),
    ]:
        for exposure in ["gain16", "high_abnormal_burden", "del6q", "grs_high"]:
            binary_rows.append(binary_summary(data, exposure, outcome, cohort))
    binary = pd.DataFrame(binary_rows)
    cox = pd.DataFrame([
        cox_summary(hist, "grs", "single_center_historical"),
        cox_summary(hist, "grs_high", "single_center_historical"),
    ])
    return reg, hist, binary, cox


def sample_size_table(score_prevalence: float, baseline_mrd_delay: float, event_rate: float) -> pd.DataFrame:
    rows = []
    power = NormalIndPower()
    for prevalence in sorted(set([0.05, 0.10, 0.25, round(score_prevalence, 3)])):
        for odds_ratio in [1.5, 2.0, 3.0]:
            p0 = baseline_mrd_delay
            p1 = odds_ratio * p0 / (1 - p0 + odds_ratio * p0)
            effect = abs(proportion_effectsize(p1, p0))
            ratio = prevalence / (1 - prevalence)
            n0 = power.solve_power(effect_size=effect, power=0.8, alpha=0.05, ratio=ratio)
            rows.append({
                "endpoint": "MRD_clearance_delay", "assumed_effect": f"OR={odds_ratio}",
                "high_score_prevalence": prevalence, "baseline_risk": p0,
                "total_n_for_80pct_power": math.ceil(n0 * (1 + ratio)),
            })
        for hr in [1.5, 2.0, 3.0]:
            events = (1.959964 + 0.841621) ** 2 / (prevalence * (1 - prevalence) * math.log(hr) ** 2)
            rows.append({
                "endpoint": "EFS", "assumed_effect": f"HR={hr}",
                "high_score_prevalence": prevalence, "baseline_risk": event_rate,
                "required_events": math.ceil(events),
                "total_n_for_80pct_power": math.ceil(events / event_rate),
            })
    return pd.DataFrame(rows)


def main() -> None:
    target, target_binary, target_cox = target_analysis()
    reg, hist, sc_binary, sc_cox = single_center_analysis()
    binary = pd.concat([target_binary, sc_binary], ignore_index=True)
    cox = pd.concat([target_cox, sc_cox], ignore_index=True)
    prevalence = float(target["grs_high"].mean())
    baseline_delay = float(reg["d19_delayed"].mean()) if reg["d19_delayed"].notna().any() else 0.5
    event_rate = float(hist["event"].sum() / hist["event"].notna().sum())
    sizes = sample_size_table(prevalence, baseline_delay, event_rate)
    model_comparison = pd.concat([
        compare_mrd_models(target, "MRD_PERCENT_DAY_29", "TARGET"),
        compare_mrd_models(hist, "d19", "single_center_historical"),
    ], ignore_index=True)

    target.to_csv(OUT / "target_grs_patient_level.csv", index=False)
    reg.to_csv(OUT / "single_center_contemporary_grs.csv", index=False)
    hist.to_csv(OUT / "single_center_historical_grs.csv", index=False)
    binary.to_csv(OUT / "grs_mrd_associations.csv", index=False)
    cox.to_csv(OUT / "grs_efs_associations.csv", index=False)
    sizes.to_csv(OUT / "multicenter_sample_size_plan.csv", index=False)
    model_comparison.to_csv(OUT / "mrd_vs_mrd_plus_grs_models.csv", index=False)

    scorecard = pd.DataFrame([
        {"component": "Gain of chromosome 16 (+16)", "points": 2, "rationale": "Strongest TARGET karyotype EFS signal; rare"},
        {"component": "Three or more abnormal karyotype tokens", "points": 1, "rationale": "Captures cytogenetic complexity/burden"},
        {"component": "del(6q)", "points": 1, "rationale": "Biologically replicated CNV; exploratory prognosis"},
    ])
    scorecard.to_csv(OUT / "locked_grs_scorecard.csv", index=False)

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.8))
    cohorts = [("TARGET", target), ("SC contemporary", reg), ("SC historical", hist)]
    axes[0].bar([x[0] for x in cohorts], [100 * x[1]["grs_high"].mean() for x in cohorts], color=["#0072B2", "#009E73", "#D55E00"])
    axes[0].set(ylabel="High genetic-risk score (%)", title="Score transportability")
    reg_plot = reg.dropna(subset=["clearance_class"])
    labels = ["D19 cleared", "D19 delayed / D46 cleared", "D46 persistent"]
    medians = [reg_plot.loc[reg_plot.clearance_class == i, "grs"].mean() for i in range(3)]
    counts = [int((reg_plot.clearance_class == i).sum()) for i in range(3)]
    bars = axes[1].bar(labels, medians, color=["#56B4E9", "#E69F00", "#CC79A7"])
    for bar, count in zip(bars, counts):
        axes[1].text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.03, f"n={count}", ha="center", fontsize=9)
    axes[1].set(ylabel="Mean genetic-risk score", title="Exploratory score across MRD clearance states")
    axes[1].tick_params(axis="x", rotation=15)
    fig.tight_layout()
    fig.savefig(OUT / "Figure_GRS_validation.png", dpi=300)
    fig.savefig(OUT / "Figure_GRS_validation.pdf")
    plt.close(fig)

    caption = (
        "Figure. Transportability and exploratory MRD-clearance distribution of the predefined "
        "cytogenetic risk score. Left: high-score prevalence differs substantially between TARGET "
        "and the single-center cohorts. Right: mean score across contemporary MRD-clearance states; "
        "the D46-persistent group contains only one patient and must not be interpreted as validation."
    )
    (OUT / "Figure_GRS_validation_caption.txt").write_text(caption, encoding="utf-8")

    target_grs = binary[(binary.cohort == "TARGET") & (binary.exposure == "grs_high") & (binary.outcome == "mrd_delayed")].iloc[0]
    sc_grs = binary[(binary.cohort == "single_center_contemporary") & (binary.exposure == "grs_high") & (binary.outcome == "d19_delayed")].iloc[0]
    report = f"""# 预先定义组合遗传风险评分：构建与可行性验证

## 锁定评分

- `+16`：2分
- `≥3个核型异常事件`：1分
- `del(6q)`：1分
- 总分范围0–4分；总分≥1定义遗传高风险。

评分仅依赖诊断期常规核型，权重在查看单中心MRD/结局关联前锁定。

## 数据可迁移性

- TARGET有核型记录患者：{len(target)}例；可评价{int(target.grs_high.notna().sum())}例；遗传高风险比例{target.grs_high.mean():.1%}。
- 单中心当代ETV6::RUNX1阳性登记患者：{len(reg)}例；可评价{int(reg.grs_high.notna().sum())}例；遗传高风险比例{reg.grs_high.mean():.1%}。
- 单中心历史ETV6::RUNX1队列：{len(hist)}例；可评价{int(hist.grs_high.notna().sum())}例；遗传高风险比例{hist.grs_high.mean():.1%}。

## MRD清除关联

- TARGET：高评分与Day-29 MRD延迟的OR={target_grs.odds_ratio:.2f}，Fisher P={target_grs.fisher_p:.3g}。
- 单中心当代队列：高评分与D19 MRD延迟的OR={sc_grs.odds_ratio:.2f}，Fisher P={sc_grs.fisher_p:.3g}。
- 单中心D46持续阳性病例极少，不能可靠检验评分与D46清除失败的关联。

## EFS与增量价值

- TARGET中评分每增加1分的EFS HR见`grs_efs_associations.csv`；该队列同时参与候选构建，因此属于内部衍生证据。
- 单中心历史队列的方向可用于外部可行性评价，但高评分患者和事件极少，不能作为正式验证。
- `mrd_vs_mrd_plus_grs_models.csv`比较MRD单独模型与MRD+评分模型；表观C-index变化不得解释为临床有效性，需多中心外部验证。

## 解释

该评分目前属于预先定义、可实施但尚未验证的候选评分。若TARGET与单中心关联不一致，应报告为“运输失败/证据不足”，不能修改权重追求显著性。建议多中心研究以D19-D46 MRD清除延迟为主要终点，以EFS为长期次要终点。
"""
    (OUT / "GENETIC_RISK_SCORE_REPORT.md").write_text(report, encoding="utf-8")


if __name__ == "__main__":
    main()
