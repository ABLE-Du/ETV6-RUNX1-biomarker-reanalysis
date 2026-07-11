from __future__ import annotations

import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
import numpy as np
import pandas as pd
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
from statsmodels.stats.proportion import proportion_confint


ROOT = Path(__file__).resolve().parent
PROJECT = ROOT.parent
FIGURES = ROOT / "figures"
TABLES = ROOT / "tables"
FIGURES.mkdir(parents=True, exist_ok=True)
TABLES.mkdir(parents=True, exist_ok=True)

TARGET = PROJECT / "genetic_risk_score" / "results" / "target_grs_patient_level.csv"
PENALTY = PROJECT / "bjh_v10_framework_reanalysis_20260623" / "target_gain16_cox_penalty_sensitivity.csv"
PREVALENCE = PROJECT / "chr16_ras_gseapy_expanded_20260710" / "tables" / "chr16_gain_public_prevalence_layer.csv"
GSEA = PROJECT / "chr16_ras_gseapy_expanded_20260710" / "tables" / "gseapy_prerank_summary_all_datasets.csv"


COLORS = {
    "blue": "#0072B2",
    "orange": "#D55E00",
    "green": "#009E73",
    "purple": "#CC79A7",
    "yellow": "#E69F00",
    "gray": "#6C757D",
    "light_gray": "#F4F6F8",
    "dark": "#202124",
}


def style() -> None:
    plt.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 8.5,
            "axes.labelsize": 9,
            "axes.titlesize": 10,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "legend.fontsize": 8,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "axes.spines.top": False,
            "axes.spines.right": False,
        }
    )


def save(fig: plt.Figure, stem: str) -> None:
    fig.savefig(FIGURES / f"{stem}.png", dpi=600, bbox_inches="tight", facecolor="white")
    fig.savefig(FIGURES / f"{stem}.pdf", bbox_inches="tight", facecolor="white")
    plt.close(fig)


def panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(-0.11, 1.06, label, transform=ax.transAxes, fontweight="bold", fontsize=11, va="top")


def figure_1_study_design() -> None:
    fig, ax = plt.subplots(figsize=(7.2, 3.5))
    ax.set_axis_off()
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    boxes = [
        (0.02, 0.54, 0.215, 0.32, COLORS["blue"], "TARGET discovery", "222 ETV6::RUNX1+\n141 with EFS + karyotype\n7 +16; 24 events\nEFS signal discovery"),
        (0.275, 0.54, 0.215, 0.32, COLORS["green"], "External CNV recurrence", "Oksa/NOPHO: 12/262\nGSE184692: 4/136\nPlatform-separated\nfrequency replication"),
        (0.53, 0.54, 0.215, 0.32, COLORS["purple"], "Paired molecular context", "GSE181157 B-ALL\n6 chr16-gain; 47 RAS-mut\nRNA-seq + text-mined\nchr16/RAS annotations"),
        (0.785, 0.54, 0.195, 0.32, COLORS["yellow"], "External expression state", "GSE227832, GSE228632\nand GSE87070\n1,059 expression samples\nNo patient-level +16/EFS"),
    ]
    for x, y, width, height, color, title, body in boxes:
        patch = FancyBboxPatch(
            (x, y), width, height, boxstyle="round,pad=0.012,rounding_size=0.016",
            facecolor="white", edgecolor=color, linewidth=1.5,
        )
        ax.add_patch(patch)
        ax.add_patch(
            FancyBboxPatch(
                (x, y + height - 0.075), width, 0.075, boxstyle="round,pad=0.012,rounding_size=0.016",
                facecolor=color, edgecolor=color, linewidth=0,
            )
        )
        ax.text(x + width / 2, y + height - 0.037, title, ha="center", va="center", color="white", fontweight="bold", fontsize=7.3)
        ax.text(x + width / 2, y + 0.115, body, ha="center", va="center", color=COLORS["dark"], fontsize=7.6, linespacing=1.35)
    for x1, x2 in [(0.235, 0.275), (0.49, 0.53), (0.745, 0.785)]:
        ax.add_patch(FancyArrowPatch((x1, 0.70), (x2, 0.70), arrowstyle="->", mutation_scale=11, linewidth=1, color=COLORS["gray"]))

    note = FancyBboxPatch((0.08, 0.12), 0.84, 0.21, boxstyle="round,pad=0.014,rounding_size=0.015", facecolor=COLORS["light_gray"], edgecolor="#D0D7DE", linewidth=0.8)
    ax.add_patch(note)
    ax.text(0.5, 0.255, "Evidence rule", ha="center", va="center", fontweight="bold", fontsize=8.5, color=COLORS["dark"])
    ax.text(
        0.5,
        0.175,
        "EFS was estimated only in TARGET. CNV cohorts were not pooled for survival, and expression cohorts were not treated as independent chr16-gain validation.",
        ha="center", va="center", fontsize=8, color=COLORS["dark"], wrap=True,
    )
    ax.set_title("Integrated public-data design and evidentiary boundaries", loc="left", pad=10, fontweight="bold")
    save(fig, "Figure_1_study_design")


def risk_counts(df: pd.DataFrame, times: list[int]) -> list[int]:
    return [int((df["time"] >= t).sum()) for t in times]


def figure_2_km() -> None:
    raw = pd.read_csv(TARGET)
    data = raw[(raw["karyotype_evaluable"] == 1) & raw["gain16"].isin([0, 1])].copy()
    data["time"] = pd.to_numeric(data["time"], errors="coerce")
    data["event"] = pd.to_numeric(data["event"], errors="coerce")
    data = data.dropna(subset=["time", "event", "gain16"])
    # TARGET's curated EFS variable is stored in days; display it in months.
    data["time"] = data["time"] / 30.4375
    data["gain16"] = data["gain16"].astype(int)
    no_gain = data[data["gain16"] == 0]
    gain = data[data["gain16"] == 1]
    test = logrank_test(gain["time"], no_gain["time"], event_observed_A=gain["event"], event_observed_B=no_gain["event"])
    km_no = KaplanMeierFitter(label="No +16")
    km_gain = KaplanMeierFitter(label="+16")
    km_no.fit(no_gain["time"], no_gain["event"])
    km_gain.fit(gain["time"], gain["event"])

    fig = plt.figure(figsize=(6.9, 4.7))
    ax = fig.add_axes([0.12, 0.27, 0.84, 0.65])
    ax_risk = fig.add_axes([0.12, 0.06, 0.84, 0.14])
    km_no.plot_survival_function(ax=ax, ci_show=True, ci_alpha=0.15, color=COLORS["blue"], linewidth=2.1)
    km_gain.plot_survival_function(ax=ax, ci_show=True, ci_alpha=0.15, color=COLORS["orange"], linewidth=2.1)
    ax.set_xlim(0, 145)
    ax.set_ylim(0.35, 1.02)
    ax.set_xlabel("Months from diagnosis")
    ax.set_ylabel("Event-free survival probability")
    ax.grid(axis="y", color="#E6E6E6", linewidth=0.7)
    ax.legend(frameon=False, loc="lower left")
    ax.text(
        0.985,
        0.05,
        f"+16: 4/7 events\nNo +16: 20/134 events\nLog-rank P = {test.p_value:.3f}",
        transform=ax.transAxes, ha="right", va="bottom", fontsize=8,
        bbox={"boxstyle": "round,pad=0.35", "facecolor": "white", "edgecolor": "#C8CDD3"},
    )
    ax.set_title("TARGET EFS discovery subset", loc="left", pad=8, fontweight="bold")
    panel_label(ax, "A")

    times = [0, 20, 40, 60, 80, 100, 120, 140]
    ax_risk.set_xlim(0, 145)
    ax_risk.set_ylim(-0.6, 2.0)
    ax_risk.axis("off")
    ax_risk.text(-10, 1.5, "At risk", ha="right", va="center", fontweight="bold", fontsize=8)
    ax_risk.text(-10, 0.8, "No +16", ha="right", va="center", color=COLORS["blue"], fontsize=8)
    ax_risk.text(-10, 0.1, "+16", ha="right", va="center", color=COLORS["orange"], fontsize=8)
    for t, n0, n1 in zip(times, risk_counts(no_gain, times), risk_counts(gain, times)):
        ax_risk.text(t, 0.8, str(n0), ha="center", va="center", fontsize=8)
        ax_risk.text(t, 0.1, str(n1), ha="center", va="center", fontsize=8)
        ax_risk.text(t, 1.55, str(t), ha="center", va="center", fontsize=7, color=COLORS["gray"])
    save(fig, "Figure_2_TARGET_EFS")

    five_year = []
    for name, km, grp in [("No +16", km_no, no_gain), ("+16", km_gain, gain)]:
        five_year.append({"group": name, "n": len(grp), "events": int(grp["event"].sum()), "efs_60_months": float(km.predict(60))})
    pd.DataFrame(five_year).to_csv(TABLES / "Figure_2_TARGET_EFS_data.csv", index=False)


def figure_3_prevalence() -> None:
    dat = pd.read_csv(PREVALENCE)
    dat = dat[dat["dataset"] != "Combined prevalence only"].copy()
    dat = dat.drop_duplicates(subset=["dataset"], keep="first").reset_index(drop=True)
    dat["label"] = ["TARGET\n(karyotype)", "Oksa/NOPHO\n(CNVkit)", "GSE184692\n(aCGH)"]
    dat["p"] = dat["chr16_gain"] / dat["total"]
    lo, hi = proportion_confint(dat["chr16_gain"], dat["total"], alpha=0.05, method="beta")
    dat["lo"] = lo
    dat["hi"] = hi

    fig, ax = plt.subplots(figsize=(6.6, 3.3))
    y = np.arange(len(dat))[::-1]
    for i, (_, row) in enumerate(dat.iterrows()):
        ypos = y[i]
        ax.plot([row["lo"] * 100, row["hi"] * 100], [ypos, ypos], color=COLORS["gray"], lw=1.4, zorder=1)
        ax.scatter(row["p"] * 100, ypos, color=[COLORS["blue"], COLORS["green"], COLORS["purple"]][i], s=55, zorder=2, edgecolor="white", linewidth=0.6)
        ax.text(10.35, ypos, f"{int(row['chr16_gain'])}/{int(row['total'])} ({row['p']*100:.1f}%)", va="center", fontsize=8)
    ax.axvline(4.24, color=COLORS["dark"], ls="--", lw=1.0)
    ax.set_yticks(y, dat["label"])
    ax.set_xlim(0, 13)
    ax.set_ylim(-0.25, 2.65)
    ax.set_xlabel("chr16 gain frequency (%) with exact 95% CI")
    ax.grid(axis="x", color="#E6E6E6", linewidth=0.7)
    ax.set_title("Recurrent low-frequency chr16 gain across independent CNV platforms", loc="left", pad=8, fontweight="bold")
    panel_label(ax, "A")
    ax.text(0.01, 0.02, "Dashed line: descriptive combined frequency, 4.2% (23/542); not pooled for EFS or RAS/MAPK testing.", fontsize=7.0, color=COLORS["gray"], transform=ax.transAxes, va="bottom")
    save(fig, "Figure_3_external_CNV_recurrence")
    dat.to_csv(TABLES / "Figure_3_external_CNV_recurrence_data.csv", index=False)


def short_comparison(value: str) -> str:
    mapping = {
        "B_ALL_chr16_gain_vs_no_gain": "GSE181157\nchr16 gain",
        "B_ALL_RASmut_vs_RASwt": "GSE181157\nRAS/MAPK-mut",
        "GSE227832_ETV6_RUNX1_or_like_vs_other_BALL": "GSE227832\nETV6::RUNX1/like",
        "GSE228632_ETV6_RUNX1_or_like_vs_other_BALL": "GSE228632\nETV6::RUNX1/like",
        "GSE87070_ETV6_RUNX1_vs_other_BCP": "GSE87070\nETV6::RUNX1",
    }
    return mapping[value]


def figure_4_molecular_context() -> None:
    dat = pd.read_csv(GSEA)
    selected_terms = ["CHR16_DOSAGE_EXPRESSED", "CHR16_16P_MAPK3_REGION", "RAS_MAPK_CORE", "MAPK_ERK_FEEDBACK"]
    term_labels = ["chr16 coding\n(positive control)", "16p/MAPK3 region", "RAS-MAPK core", "ERK feedback"]
    comparisons = [
        "B_ALL_chr16_gain_vs_no_gain",
        "B_ALL_RASmut_vs_RASwt",
        "GSE227832_ETV6_RUNX1_or_like_vs_other_BALL",
        "GSE228632_ETV6_RUNX1_or_like_vs_other_BALL",
        "GSE87070_ETV6_RUNX1_vs_other_BCP",
    ]
    subset = dat[dat["comparison"].isin(comparisons) & dat["Term"].isin(selected_terms)].copy()
    subset["Term"] = pd.Categorical(subset["Term"], categories=selected_terms, ordered=True)
    subset["comparison"] = pd.Categorical(subset["comparison"], categories=comparisons, ordered=True)
    pivot = subset.pivot(index="Term", columns="comparison", values="NES").reindex(index=selected_terms, columns=comparisons)
    fdr = subset.pivot(index="Term", columns="comparison", values="FDR q-val").reindex(index=selected_terms, columns=comparisons)

    fig = plt.figure(figsize=(7.1, 4.5))
    ax = fig.add_axes([0.12, 0.2, 0.64, 0.7])
    for yi, term in enumerate(selected_terms[::-1]):
        for xi, comp in enumerate(comparisons):
            nes = pivot.loc[term, comp]
            q = fdr.loc[term, comp]
            if pd.isna(nes):
                continue
            color = COLORS["orange"] if nes >= 0 else COLORS["blue"]
            size = 30 + min(-math.log10(max(float(q), 1e-4)), 4.0) * 38
            ax.scatter(xi, yi, s=size, color=color, alpha=0.78, edgecolor=COLORS["dark"], linewidth=0.35)
            ax.text(xi, yi, f"{nes:.2f}", ha="center", va="center", fontsize=6.7, color="white" if abs(nes) > 1.45 else COLORS["dark"], fontweight="bold" if q < 0.05 else "normal")
    ax.axvline(1.5, color="#C7CDD4", lw=1, ls="--")
    ax.set_xticks(range(len(comparisons)), [short_comparison(x) for x in comparisons], rotation=0, ha="center", fontsize=7.4)
    ax.set_yticks(range(len(selected_terms)), term_labels[::-1])
    ax.set_xlim(-0.55, len(comparisons) - 0.45)
    ax.set_ylim(-0.65, len(selected_terms) - 0.35)
    ax.set_title("Preranked pathway enrichment across expression cohorts", loc="left", pad=8, fontweight="bold")
    panel_label(ax, "A")
    ax.text(0.5, -0.32, "Sample sizes: 6 vs 130; 47 vs 89; 49 vs 252; 35 vs 25; 172 vs 402. Circle color indicates NES direction; bold values have FDR < 0.05.", transform=ax.transAxes, ha="center", fontsize=6.8, color=COLORS["gray"])

    ax2 = fig.add_axes([0.80, 0.2, 0.18, 0.7])
    ax2.axis("off")
    ax2.text(0.0, 0.96, "Interpretation", fontweight="bold", fontsize=9, va="top")
    ax2.text(0.0, 0.79, "Direct same-patient\nchr16/RAS evidence", fontsize=7.7, fontweight="bold", va="top")
    ax2.text(0.0, 0.65, "GSE181157 showed enrichment of chr16 dosage, 16p/MAPK3, RAS-MAPK and ERK feedback in chr16-gain B-ALL.", fontsize=7.3, va="top", wrap=True)
    ax2.text(0.0, 0.41, "External state context", fontsize=7.7, fontweight="bold", va="top")
    ax2.text(0.0, 0.27, "The three external datasets do not carry matched chr16-gain annotations. They support subtype-state context only, not independent mechanistic validation.", fontsize=7.3, va="top", wrap=True)
    save(fig, "Figure_4_expression_context")
    subset.to_csv(TABLES / "Figure_4_expression_context_data.csv", index=False)


def figure_s1_penalty_sensitivity() -> None:
    dat = pd.read_csv(PENALTY)
    dat = dat[dat["variables"] == "gain16"].copy().sort_values("penalizer")
    fig, ax = plt.subplots(figsize=(5.9, 3.1))
    y = np.arange(len(dat))[::-1]
    for i, (_, row) in enumerate(dat.iterrows()):
        ypos = y[i]
        ax.plot([row["gain16_ci_low"], row["gain16_ci_high"]], [ypos, ypos], color=COLORS["gray"], lw=1.3)
        ax.scatter(row["gain16_hr"], ypos, s=43, color=COLORS["purple"], edgecolor="white", linewidth=0.5, zorder=3)
        ax.text(13.2, ypos, f"{row['gain16_hr']:.2f} ({row['gain16_ci_low']:.2f}-{row['gain16_ci_high']:.2f})", va="center", fontsize=7.5)
    ax.axvline(1, color=COLORS["dark"], ls="--", lw=1)
    ax.set_xscale("log")
    ax.set_xlim(0.55, 20)
    ax.set_yticks(y, ["Unpenalized" if x == 0 else f"Ridge penalty = {x:g}" for x in dat["penalizer"]])
    ax.set_xlabel("Hazard ratio for +16 (log scale)")
    ax.set_title("Sensitivity of the exploratory TARGET association to Cox penalization", loc="left", pad=8, fontweight="bold")
    ax.text(0.0, -0.27, "All estimates use the same 141-patient/24-event subset. No penalty value was pre-specified; results are sensitivity analyses rather than confirmatory model selection.", transform=ax.transAxes, fontsize=7.2, color=COLORS["gray"])
    save(fig, "Supplementary_Figure_S1_penalty_sensitivity")
    dat.to_csv(TABLES / "Supplementary_Figure_S1_penalty_sensitivity_data.csv", index=False)


def main() -> None:
    style()
    figure_1_study_design()
    figure_2_km()
    figure_3_prevalence()
    figure_4_molecular_context()
    figure_s1_penalty_sensitivity()


if __name__ == "__main__":
    main()
