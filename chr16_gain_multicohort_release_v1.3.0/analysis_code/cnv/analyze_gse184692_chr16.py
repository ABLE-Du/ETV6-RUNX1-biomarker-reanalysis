from __future__ import annotations

import csv
import gzip
import json
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact
from statsmodels.stats.proportion import proportion_confint


ROOT = Path(__file__).resolve().parent
PROJECT = ROOT.parent

SERIES_MATRIX = ROOT / "geo_matrix" / "GSE184692_series_matrix.txt.gz"
SEGMENTED = ROOT / "geo_supp" / "GSE184692_Segmented.txt.gz"
GPL = ROOT / "geo_platforms" / "GPL10150_family.soft.gz"
TARGET = PROJECT / "genetic_risk_score" / "results" / "target_grs_patient_level.csv"
OKSA = PROJECT / "database_expansion_search_20260623" / "oksa_chr16_gain_screen.csv"


def read_gse184692_samples() -> pd.DataFrame:
    metadata: dict[str, list[list[str]]] = {}
    with gzip.open(SERIES_MATRIX, "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("!series_matrix_table_begin"):
                break
            if line.startswith("!Sample_"):
                parts = next(csv.reader([line.rstrip("\n")], delimiter="\t"))
                key = parts[0].replace("!Sample_", "")
                vals = [v.strip('"') for v in parts[1:]]
                metadata.setdefault(key, []).append(vals)

    sample_titles = metadata.get("title", [[]])[0]
    geo = metadata.get("geo_accession", [[]])[0]
    chars = metadata.get("characteristics_ch1", [])
    rows = []
    for idx, title in enumerate(sample_titles):
        row = {"sample_title": title, "geo_accession": geo[idx] if idx < len(geo) else None}
        for char_row in chars:
            if idx >= len(char_row):
                continue
            value = char_row[idx]
            if ":" in value:
                key, val = value.split(":", 1)
                row[key.strip().lower().replace(" ", "_")] = val.strip()
        rows.append(row)
    return pd.DataFrame(rows)


def read_gpl_chr16_probe_ids() -> set[str]:
    chr16_ids: set[str] = set()
    table_started = False
    header: list[str] | None = None
    with gzip.open(GPL, "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("!platform_table_begin"):
                table_started = True
                continue
            if line.startswith("!platform_table_end"):
                break
            if not table_started:
                continue
            parts = next(csv.reader([line.rstrip("\n")], delimiter="\t"))
            if header is None:
                header = parts
                continue
            if not parts or len(parts) < len(header):
                continue
            rec = dict(zip(header, parts))
            if rec.get("CONTROL_TYPE") not in {"FALSE", "false", ""}:
                continue
            chrom = rec.get("CHROMOSOME", "")
            if chrom in {"chr16", "16"}:
                probe_id = rec.get("ID")
                if probe_id:
                    chr16_ids.add(probe_id)
    return chr16_ids


def summarize_chr16_segmented(chr16_ids: set[str], sample_meta: pd.DataFrame) -> pd.DataFrame:
    sample_cols = sample_meta["sample_title"].tolist()
    accum = {
        s: {
            "chr16_probe_n": 0,
            "sum_log2": 0.0,
            "n_ge_0_15": 0,
            "n_ge_0_20": 0,
            "n_ge_0_25": 0,
            "n_ge_0_30": 0,
            "n_le_minus_0_30": 0,
        }
        for s in sample_cols
    }

    chunks = pd.read_csv(SEGMENTED, sep="\t", compression="gzip", chunksize=10000)
    for chunk in chunks:
        chunk = chunk[chunk["ID_REF"].isin(chr16_ids)]
        if chunk.empty:
            continue
        keep = ["ID_REF", *[c for c in sample_cols if c in chunk.columns]]
        chunk = chunk[keep]
        for col in keep[1:]:
            vals = pd.to_numeric(chunk[col], errors="coerce").dropna().to_numpy()
            if vals.size == 0:
                continue
            a = accum[col]
            a["chr16_probe_n"] += int(vals.size)
            a["sum_log2"] += float(vals.sum())
            a["n_ge_0_15"] += int((vals >= 0.15).sum())
            a["n_ge_0_20"] += int((vals >= 0.20).sum())
            a["n_ge_0_25"] += int((vals >= 0.25).sum())
            a["n_ge_0_30"] += int((vals >= 0.30).sum())
            a["n_le_minus_0_30"] += int((vals <= -0.30).sum())

    rows = []
    for sample, a in accum.items():
        n = a["chr16_probe_n"]
        row = {"sample_title": sample, **a}
        row["chr16_mean_log2"] = a["sum_log2"] / n if n else np.nan
        for key in ["n_ge_0_15", "n_ge_0_20", "n_ge_0_25", "n_ge_0_30", "n_le_minus_0_30"]:
            row[key.replace("n_", "frac_")] = a[key] / n if n else np.nan
        row["chr16_gain_sensitive"] = bool(row["chr16_mean_log2"] >= 0.15 and row["frac_ge_0_20"] >= 0.50)
        row["chr16_gain_primary"] = bool(row["chr16_mean_log2"] >= 0.20 and row["frac_ge_0_25"] >= 0.50)
        row["chr16_gain_strict"] = bool(row["chr16_mean_log2"] >= 0.25 and row["frac_ge_0_30"] >= 0.50)
        rows.append(row)

    out = pd.DataFrame(rows)
    return sample_meta.merge(out, on="sample_title", how="left")


def ci(k: int, n: int) -> tuple[float, float]:
    if n == 0:
        return math.nan, math.nan
    lo, hi = proportion_confint(k, n, alpha=0.05, method="beta")
    return float(lo), float(hi)


def prevalence_row(dataset: str, definition: str, k: int, n: int) -> dict:
    lo, hi = ci(k, n)
    return {
        "dataset": dataset,
        "definition": definition,
        "chr16_gain": k,
        "total": n,
        "prevalence_percent": 100 * k / n if n else np.nan,
        "ci95_low_percent": 100 * lo if n else np.nan,
        "ci95_high_percent": 100 * hi if n else np.nan,
        "formatted": f"{k}/{n} ({100*k/n:.1f}%; 95% CI {100*lo:.1f}-{100*hi:.1f})" if n else "NA",
    }


def main() -> None:
    (ROOT / "tables").mkdir(exist_ok=True)
    (ROOT / "figures").mkdir(exist_ok=True)

    samples = read_gse184692_samples()
    samples.to_csv(ROOT / "tables" / "GSE184692_sample_metadata.csv", index=False, encoding="utf-8-sig")

    chr16_ids = read_gpl_chr16_probe_ids()
    (ROOT / "tables" / "GPL10150_chr16_probe_ids.txt").write_text("\n".join(sorted(chr16_ids)), encoding="utf-8")

    gse = summarize_chr16_segmented(chr16_ids, samples)
    gse.to_csv(ROOT / "tables" / "GSE184692_chr16_probe_summary.csv", index=False, encoding="utf-8-sig")

    subtype_summary = (
        gse.groupby("all_subtype", dropna=False)
        .agg(
            n=("sample_title", "size"),
            chr16_gain_sensitive=("chr16_gain_sensitive", "sum"),
            chr16_gain_primary=("chr16_gain_primary", "sum"),
            chr16_gain_strict=("chr16_gain_strict", "sum"),
            median_chr16_mean_log2=("chr16_mean_log2", "median"),
        )
        .reset_index()
        .sort_values("n", ascending=False)
    )
    subtype_summary.to_csv(ROOT / "tables" / "GSE184692_chr16_gain_by_subtype.csv", index=False, encoding="utf-8-sig")

    etv6 = gse[gse["all_subtype"].eq("ETV6-RUNX1")].copy()
    target = pd.read_csv(TARGET)
    target_k = target[(target["karyotype_evaluable"] == 1) & target["gain16"].isin([0, 1])]
    oksa = pd.read_csv(OKSA)

    rows = [
        prevalence_row(
            "TARGET diagnostic karyotype",
            "curated karyotype +16/gain16",
            int(target_k["gain16"].sum()),
            int(len(target_k)),
        ),
        prevalence_row(
            "Oksa/NOPHO CNVkit",
            "broad/whole-chromosome chr16 gain in CNVkit segments",
            int(oksa["broad_chr16_gain_rule"].sum()),
            int(len(oksa)),
        ),
        prevalence_row(
            "GSE184692 aCGH ETV6-RUNX1",
            "primary rule: chr16 mean log2 >=0.20 and >=50% chr16 probes >=0.25",
            int(etv6["chr16_gain_primary"].sum()),
            int(len(etv6)),
        ),
        prevalence_row(
            "GSE184692 aCGH ETV6-RUNX1",
            "sensitive rule: chr16 mean log2 >=0.15 and >=50% chr16 probes >=0.20",
            int(etv6["chr16_gain_sensitive"].sum()),
            int(len(etv6)),
        ),
        prevalence_row(
            "GSE184692 aCGH ETV6-RUNX1",
            "strict rule: chr16 mean log2 >=0.25 and >=50% chr16 probes >=0.30",
            int(etv6["chr16_gain_strict"].sum()),
            int(len(etv6)),
        ),
    ]
    rows.append(
        prevalence_row(
            "Combined prevalence only",
            "TARGET + Oksa/NOPHO + GSE184692 primary rule; not pooled for EFS/RAS",
            int(target_k["gain16"].sum()) + int(oksa["broad_chr16_gain_rule"].sum()) + int(etv6["chr16_gain_primary"].sum()),
            int(len(target_k)) + int(len(oksa)) + int(len(etv6)),
        )
    )
    prevalence = pd.DataFrame(rows)

    primary = prevalence.iloc[2]
    comparison_rows = []
    for idx in [0, 1]:
        ref = prevalence.iloc[idx]
        odds, p = fisher_exact(
            [
                [int(primary["chr16_gain"]), int(primary["total"] - primary["chr16_gain"])],
                [int(ref["chr16_gain"]), int(ref["total"] - ref["chr16_gain"])],
            ],
            alternative="two-sided",
        )
        comparison_rows.append(
            {
                "comparison": f"GSE184692 primary vs {ref['dataset']}",
                "odds_ratio": odds,
                "fisher_two_sided_p": p,
            }
        )
    comp = pd.DataFrame(comparison_rows)
    prevalence.to_csv(ROOT / "tables" / "expanded_chr16_gain_prevalence_comparison.csv", index=False, encoding="utf-8-sig")
    comp.to_csv(ROOT / "tables" / "expanded_chr16_gain_prevalence_fisher.csv", index=False, encoding="utf-8-sig")

    plot_df = prevalence.iloc[[0, 1, 2]].copy()
    fig, ax = plt.subplots(figsize=(6.5, 4.0))
    x = np.arange(len(plot_df))
    y = plot_df["prevalence_percent"].to_numpy()
    yerr = np.vstack(
        [
            y - plot_df["ci95_low_percent"].to_numpy(),
            plot_df["ci95_high_percent"].to_numpy() - y,
        ]
    )
    colors = ["#009E73", "#CC79A7", "#0072B2"]
    ax.bar(x, y, yerr=yerr, capsize=4, color=colors, edgecolor="white")
    ax.set_xticks(x)
    ax.set_xticklabels(["TARGET\nkaryotype", "Oksa/NOPHO\nCNVkit", "GSE184692\naCGH"], rotation=0)
    ax.set_ylabel("chr16 gain prevalence (%)")
    ax.set_title("External validation of chr16 gain prevalence")
    ax.grid(axis="y", alpha=0.2)
    for i, row in plot_df.iterrows():
        loc = list(plot_df.index).index(i)
        ax.text(loc, row["ci95_high_percent"] + 0.8, f"{int(row['chr16_gain'])}/{int(row['total'])}", ha="center", fontsize=9)
    fig.tight_layout()
    fig.savefig(ROOT / "figures" / "Figure_expanded_chr16_gain_prevalence.png", dpi=300)
    fig.savefig(ROOT / "figures" / "Figure_expanded_chr16_gain_prevalence.pdf")
    plt.close(fig)

    metrics = {
        "gse184692_total_samples": int(len(gse)),
        "gse184692_etv6_runx1_samples": int(len(etv6)),
        "gse184692_chr16_probe_count": int(len(chr16_ids)),
        "gse184692_etv6_chr16_gain_primary": int(etv6["chr16_gain_primary"].sum()),
        "gse184692_etv6_chr16_gain_sensitive": int(etv6["chr16_gain_sensitive"].sum()),
        "gse184692_etv6_chr16_gain_strict": int(etv6["chr16_gain_strict"].sum()),
        "subtype_counts": samples["all_subtype"].value_counts(dropna=False).to_dict(),
    }
    (ROOT / "GSE184692_chr16_metrics.json").write_text(json.dumps(metrics, indent=2, ensure_ascii=False), encoding="utf-8")


if __name__ == "__main__":
    main()
