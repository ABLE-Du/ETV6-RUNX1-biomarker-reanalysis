from pathlib import Path
import warnings

import numpy as np
import pandas as pd
from lifelines import CoxPHFitter
from lifelines.utils import concordance_index


ROOT = Path(__file__).resolve().parent
INP = ROOT / "results"
OUT = ROOT / "results"
RNG = np.random.default_rng(20260615)
N_BOOT = 500


def fit_and_predict(train: pd.DataFrame, test: pd.DataFrame, variables: list[str]):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        model = CoxPHFitter(penalizer=0.5)
        model.fit(train[["time", "event", *variables]], "time", "event")
    train_risk = model.predict_partial_hazard(train[variables]).to_numpy()
    test_risk = model.predict_partial_hazard(test[variables]).to_numpy()
    train_c = concordance_index(train["time"], -train_risk, train["event"])
    test_c = concordance_index(test["time"], -test_risk, test["event"])
    return train_c, test_c


def validate(data: pd.DataFrame, variables: list[str], cohort: str) -> dict:
    use = data[["time", "event", *variables]].dropna().reset_index(drop=True)
    apparent, _ = fit_and_predict(use, use, variables)
    optimism = []
    test_values = []
    failed = 0
    for _ in range(N_BOOT):
        sample = use.iloc[RNG.integers(0, len(use), len(use))].copy().reset_index(drop=True)
        if sample["event"].sum() < 3 or any(sample[v].nunique() < 2 for v in variables):
            failed += 1
            continue
        try:
            app_boot, test_original = fit_and_predict(sample, use, variables)
            optimism.append(app_boot - test_original)
            test_values.append(test_original)
        except Exception:
            failed += 1
    optimism = np.asarray(optimism)
    corrected = apparent - np.mean(optimism) if len(optimism) else np.nan
    return {
        "cohort": cohort,
        "model": "+".join(variables),
        "n": len(use),
        "events": int(use["event"].sum()),
        "events_per_parameter": use["event"].sum() / len(variables),
        "bootstrap_requested": N_BOOT,
        "bootstrap_valid": len(optimism),
        "bootstrap_failed": failed,
        "apparent_c_index": apparent,
        "mean_optimism": float(np.mean(optimism)) if len(optimism) else np.nan,
        "optimism_corrected_c_index": corrected,
        "bootstrap_test_c_low": float(np.quantile(test_values, 0.025)) if len(test_values) else np.nan,
        "bootstrap_test_c_high": float(np.quantile(test_values, 0.975)) if len(test_values) else np.nan,
        "reliability_flag": (
            "severely_underpowered"
            if use["event"].sum() < 20
            else "exploratory"
            if use["event"].sum() < 100
            else "adequate_for_internal_validation"
        ),
    }


def main() -> None:
    target = pd.read_csv(INP / "target_grs_patient_level.csv")
    target["time"] = pd.to_numeric(target["time"], errors="coerce")
    target["event"] = pd.to_numeric(target["event"], errors="coerce")
    target["log_mrd"] = np.log10(pd.to_numeric(target["MRD_PERCENT_DAY_29"], errors="coerce") + 0.0001)

    historical = pd.read_csv(INP / "single_center_historical_grs.csv")
    historical["time"] = pd.to_numeric(historical["time"], errors="coerce")
    historical["event"] = pd.to_numeric(historical["event"], errors="coerce")
    historical["log_mrd"] = np.log10(pd.to_numeric(historical["d19"], errors="coerce") + 0.0001)

    rows = []
    for cohort, data in [("TARGET", target), ("single_center_historical", historical)]:
        rows.append(validate(data, ["log_mrd"], cohort))
        rows.append(validate(data, ["log_mrd", "grs"], cohort))
    result = pd.DataFrame(rows)
    result.to_csv(OUT / "bootstrap_internal_validation.csv", index=False)

    report = [
        "# Bootstrap internal validation",
        "",
        "Five hundred bootstrap resamples were requested. Penalized Cox models used the same fixed predictors in every resample.",
        "The optimism-corrected C-index is reported as an internal-validation estimate, not external validation.",
        "",
        "|Cohort|Model|n/events|Apparent C|Mean optimism|Corrected C|Reliability|",
        "|---|---|---:|---:|---:|---:|---|",
    ]
    for r in result.itertuples():
        report.append(
            f"|{r.cohort}|{r.model}|{r.n}/{r.events}|{r.apparent_c_index:.3f}|"
            f"{r.mean_optimism:.3f}|{r.optimism_corrected_c_index:.3f}|{r.reliability_flag}|"
        )
    report += [
        "",
        "Calibration curves and calibration slopes were not interpreted because both datasets contain far fewer than 100 EFS events, and the historical cohort contains only five complete-case events.",
    ]
    (OUT / "BOOTSTRAP_INTERNAL_VALIDATION_REPORT.md").write_text("\n".join(report), encoding="utf-8")


if __name__ == "__main__":
    main()
