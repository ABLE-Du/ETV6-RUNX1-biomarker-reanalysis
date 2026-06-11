from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from lifelines import KaplanMeierFitter
from lifelines.statistics import multivariate_logrank_test


COLORS = {
    "blue": "#2C7FB8",
    "teal": "#238B8D",
    "orange": "#F28E2B",
    "red": "#C44E52",
    "purple": "#7A5195",
    "lightgray": "#D9D9D9",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate the historical-only single-center outcome figure."
    )
    parser.add_argument("--workspace-root", type=Path, default=Path.cwd())
    return parser.parse_args()


def panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(
        -0.13,
        1.08,
        label,
        transform=ax.transAxes,
        fontsize=13,
        fontweight="bold",
        va="top",
    )


def main() -> None:
    root = parse_args().workspace_root.resolve()
    input_path = (
        root
        / "single_center_manuscript_analysis"
        / "tables"
        / "historical_outcome_analytic_cohort_pseudonymized.csv"
    )
    output_dir = root / "single_center_manuscript_analysis" / "figures"
    output_dir.mkdir(parents=True, exist_ok=True)
    outcome = pd.read_csv(input_path)

    sns.set_theme(style="whitegrid", context="paper")
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 8,
            "axes.titlesize": 9,
            "axes.labelsize": 8,
            "legend.fontsize": 7,
            "savefig.dpi": 600,
            "pdf.fonttype": 42,
            "svg.fonttype": "none",
        }
    )

    fig = plt.figure(figsize=(13, 9), constrained_layout=True)
    gs = fig.add_gridspec(2, 2, height_ratios=[1, 1.08])
    axes = [fig.add_subplot(gs[row, col]) for row in range(2) for col in range(2)]
    ax_a, ax_b, ax_c, ax_d = axes

    era = outcome["来源表"].astype(str).map(
        lambda value: "ALL2015" if "2015" in value else "ALL2005/2009"
    )
    composition = (
        pd.DataFrame({"cohort": era, "event": outcome["efs_event"].astype(int)})
        .groupby("cohort", as_index=False)
        .agg(n=("event", "size"), events=("event", "sum"))
        .set_index("cohort")
        .reindex(["ALL2005/2009", "ALL2015"])
        .reset_index()
    )
    bars = ax_a.barh(
        composition["cohort"],
        composition["n"],
        color=[COLORS["blue"], COLORS["teal"]],
    )
    for bar, row in zip(bars, composition.itertuples()):
        ax_a.text(
            row.n + 0.8,
            bar.get_y() + bar.get_height() / 2,
            f"n={row.n}; events={row.events}",
            va="center",
        )
    ax_a.set_xlim(0, max(composition["n"]) + 9)
    ax_a.set_xlabel("Patients")
    ax_a.set_title("Historical ETV6::RUNX1-positive outcome cohort")
    sns.despine(ax=ax_a)
    panel_label(ax_a, "a")

    for event, time, label, color in [
        ("efs_event", "efs_months", "EFS", COLORS["red"]),
        ("os_event", "os_months", "OS", COLORS["blue"]),
    ]:
        subset = outcome[[event, time]].dropna()
        km = KaplanMeierFitter().fit(subset[time], subset[event], label=label)
        km.plot_survival_function(ax=ax_b, ci_show=True, color=color, linewidth=1.8)
    ax_b.set_xlim(0, 210)
    ax_b.set_ylim(0.65, 1.01)
    ax_b.set_xlabel("Months from diagnosis")
    ax_b.set_ylabel("Survival probability")
    ax_b.set_title(
        f"Overall outcomes (n={len(outcome)}; "
        f"{int(outcome.efs_event.sum())} EFS events; {int(outcome.os_event.sum())} deaths)"
    )
    ax_b.legend(frameon=False, loc="lower left")
    sns.despine(ax=ax_b)
    panel_label(ax_b, "b")

    risk_colors = {"LR": COLORS["blue"], "IR": COLORS["orange"], "HR": COLORS["red"]}
    risk_data = outcome[outcome["risk_group"].isin(risk_colors)].copy()
    for group in ["LR", "IR", "HR"]:
        subset = risk_data[risk_data["risk_group"].eq(group)]
        km = KaplanMeierFitter().fit(
            subset["efs_months"],
            subset["efs_event"],
            label=f"{group} (n={len(subset)}, events={int(subset.efs_event.sum())})",
        )
        km.plot_survival_function(
            ax=ax_c, ci_show=False, color=risk_colors[group], linewidth=1.8
        )
    p_value = multivariate_logrank_test(
        risk_data["efs_months"], risk_data["risk_group"], risk_data["efs_event"]
    ).p_value
    ax_c.text(
        0.98,
        0.96,
        f"Global log-rank P={p_value:.3f}",
        transform=ax_c.transAxes,
        ha="right",
        va="top",
    )
    ax_c.set_xlim(0, 210)
    ax_c.set_ylim(0.4, 1.01)
    ax_c.set_xlabel("Months from diagnosis")
    ax_c.set_ylabel("EFS probability")
    ax_c.set_title("EFS by source-recorded risk group")
    ax_c.legend(frameon=False, loc="lower left")
    sns.despine(ax=ax_c)
    panel_label(ax_c, "c")

    events = outcome[outcome["efs_event"].eq(1)].sort_values("efs_months").copy()
    event_text = events["EFS事件"].fillna("").astype(str)
    events["event_class"] = np.where(
        event_text.str.contains("复发"),
        "Relapse",
        np.where(event_text.str.contains("感染"), "Infection-related death", "Other event"),
    )
    event_colors = {
        "Relapse": COLORS["red"],
        "Infection-related death": COLORS["orange"],
        "Other event": COLORS["purple"],
    }
    y = np.arange(len(events))
    ax_d.hlines(y, 0, events["efs_months"], color=COLORS["lightgray"], linewidth=3)
    for index, row in enumerate(events.itertuples()):
        ax_d.scatter(
            row.efs_months,
            index,
            s=45,
            color=event_colors[row.event_class],
            edgecolor="white",
            linewidth=0.5,
            zorder=3,
        )
        status = "death" if row.os_event else "alive/censored"
        ax_d.text(
            row.efs_months + 3,
            index,
            f"{row.event_class}; {status}",
            va="center",
            fontsize=7,
        )
    ax_d.set_yticks(y)
    ax_d.set_yticklabels(events["study_id"])
    ax_d.set_xlabel("Months to first EFS event")
    ax_d.set_title("Patient-level event spectrum")
    ax_d.set_xlim(0, max(events["efs_months"].max() + 55, 105))
    sns.despine(ax=ax_d)
    panel_label(ax_d, "d")

    fig.suptitle(
        "Single-center ETV6::RUNX1-positive historical outcome cohort",
        fontsize=12,
        fontweight="bold",
    )
    for suffix in ["png", "pdf", "svg"]:
        fig.savefig(
            output_dir / f"Figure_SC1_historical_only_outcomes.{suffix}",
            bbox_inches="tight",
            facecolor="white",
        )
    plt.close(fig)


if __name__ == "__main__":
    main()
