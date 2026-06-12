from __future__ import annotations

import argparse
import gzip
import hashlib
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests


DEFAULT_SOURCE_ROOT = Path("public_data/NOPHO_CNV/raw")
DEFAULT_BLACKLIST = Path("reproducibility/nopho_grch37/input/hg19-blacklist.v2.bed.gz")
DEFAULT_OUTPUT = Path("reproducibility/nopho_grch37/results")
OFFICIAL_CODE_COMMIT = "48607138ce0df2950d3ffc0f304246c38acf023f"

# Ensembl GRCh37 REST lookup/symbol/homo_sapiens, retrieved 2026-06-12.
# Coordinates are recorded as returned by Ensembl. One-base endpoint differences do
# not affect overlap calls for the >=1-Mb events evaluated here.
GENES = {
    "FOXO3": ("6", 108881038, 109005977),
    "TNFAIP3": ("6", 138188351, 138204449),
    "HACE1": ("6", 105175968, 105307794),
    "GRIK2": ("6", 101846664, 102517958),
    "EPHA7": ("6", 93949738, 94129265),
    "MYB": ("6", 135502453, 135540311),
    "PRDM1": ("6", 106534195, 106557814),
    "BACH2": ("6", 90636248, 91006627),
    "LATS1": ("6", 149979289, 150039392),
    "CDKN2A": ("9", 21967751, 21995300),
    "PAX5": ("9", 36833272, 37034103),
    "IKZF1": ("7", 50343720, 50472799),
}

SENSITIVITY_THRESHOLDS = (-0.3, -0.4, -0.5)
SENSITIVITY_MIN_SIZES = (0, 500_000, 1_000_000)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="GRCh37-compatible reanalysis of the public NOPHO CNVkit segment files."
    )
    parser.add_argument("--source-root", type=Path, default=DEFAULT_SOURCE_ROOT)
    parser.add_argument("--blacklist", type=Path, default=DEFAULT_BLACKLIST)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_segments(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path, sep="\t")
    frame["chromosome"] = frame["chromosome"].astype(str).str.replace("chr", "", regex=False)
    frame = frame.loc[frame["chromosome"].isin([str(i) for i in range(1, 23)])].copy()
    for column in ("start", "end", "log2"):
        frame[column] = pd.to_numeric(frame[column])
    return frame.sort_values(["chromosome", "start", "end"]).reset_index(drop=True)


def load_blacklist(path: Path) -> pd.DataFrame:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as handle:
        frame = pd.read_csv(
            handle,
            sep="\t",
            header=None,
            usecols=[0, 1, 2],
            names=["chromosome", "start", "end"],
        )
    frame["chromosome"] = frame["chromosome"].astype(str).str.replace("chr", "", regex=False)
    return frame


def call_and_merge(
    segments: pd.DataFrame,
    blacklist: pd.DataFrame,
    del_threshold: float = -0.4,
    amp_threshold: float = 0.3,
    merge_gap: int = 10_000,
    min_size: int = 1_000_000,
) -> pd.DataFrame:
    work = segments.copy()
    work["type"] = np.where(
        work["log2"] <= del_threshold,
        "DEL",
        np.where(work["log2"] >= amp_threshold, "AMP", None),
    )

    events: list[dict[str, object]] = []
    current: dict[str, object] | None = None
    for row in work.itertuples(index=False):
        if row.type is None:
            current = None
            continue
        item = {
            "chromosome": str(row.chromosome),
            "start": int(row.start),
            "end": int(row.end),
            "type": row.type,
            "min_log2": float(row.log2),
            "max_log2": float(row.log2),
        }
        if (
            current is not None
            and current["chromosome"] == item["chromosome"]
            and current["type"] == item["type"]
            and int(item["start"]) - int(current["end"]) < merge_gap
        ):
            current["end"] = max(int(current["end"]), int(item["end"]))
            current["min_log2"] = min(float(current["min_log2"]), float(item["min_log2"]))
            current["max_log2"] = max(float(current["max_log2"]), float(item["max_log2"]))
        else:
            events.append(item)
            current = events[-1]

    columns = ["chromosome", "start", "end", "type", "min_log2", "max_log2", "length"]
    if not events:
        return pd.DataFrame(columns=columns)
    output = pd.DataFrame(events)
    output["length"] = output["end"] - output["start"]
    # The official R code removes events <1 Mb, so exactly 1-Mb events are retained.
    output = output.loc[output["length"] >= min_size].reset_index(drop=True)

    keep: list[bool] = []
    for event in output.itertuples(index=False):
        regions = blacklist.loc[blacklist["chromosome"] == str(event.chromosome)]
        fractions = [
            max(0, min(event.end, region.end) - max(event.start, region.start)) / event.length
            for region in regions.itertuples(index=False)
        ]
        # Mirrors bedtools intersect -f 0.5 -wa -v used by the official workflow.
        keep.append(max(fractions, default=0) < 0.5)
    return output.loc[keep].reset_index(drop=True)


def overlaps(events: pd.DataFrame, gene: str, event_type: str = "DEL") -> bool:
    chromosome, start, end = GENES[gene]
    if events.empty:
        return False
    matches = events.loc[
        (events["chromosome"] == chromosome)
        & (events["type"] == event_type)
        & (events["start"] < end)
        & (events["end"] > start)
    ]
    return not matches.empty


def raw_segment_overlaps(segments: pd.DataFrame, gene: str, threshold: float = -0.4) -> bool:
    chromosome, start, end = GENES[gene]
    matches = segments.loc[
        (segments["chromosome"] == chromosome)
        & (segments["log2"] <= threshold)
        & (segments["start"] < end)
        & (segments["end"] > start)
    ]
    return not matches.empty


def summarize_matrix(matrix: pd.DataFrame, method: str) -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "method": method,
                "gene": gene,
                "n_deleted": int(matrix[gene].sum()),
                "n_total": len(matrix),
                "percent": 100 * matrix[gene].mean(),
            }
            for gene in GENES
        ]
    )


def plot_nopho_summary(
    official: pd.DataFrame, sensitivity_summary: pd.DataFrame, output_path: Path
) -> None:
    genes_6q = ["BACH2", "EPHA7", "GRIK2", "HACE1", "PRDM1", "FOXO3", "MYB", "TNFAIP3", "LATS1"]
    counts = official[genes_6q].sum().sort_values(ascending=False)
    threshold_view = sensitivity_summary.loc[sensitivity_summary["min_size"] == 1_000_000]

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.6), constrained_layout=True)
    axes[0].bar(counts.index, counts.values, color="#3C78A8")
    axes[0].set_ylabel("Cases with overlapping deletion")
    axes[0].set_title("a  GRCh37 gene-overlap counts")
    axes[0].tick_params(axis="x", rotation=55)

    for gene, color in [("del6q", "#2F4B7C"), ("FOXO3", "#D45087"), ("TNFAIP3", "#FFA600")]:
        axes[1].plot(
            threshold_view["del_threshold"],
            threshold_view[gene],
            marker="o",
            linewidth=2,
            label=gene,
            color=color,
        )
    axes[1].invert_xaxis()
    axes[1].set_xlabel("Deletion log2 threshold")
    axes[1].set_ylabel("Cases at >=1-Mb size filter")
    axes[1].set_title("b  Threshold sensitivity")
    axes[1].legend(frameon=False)
    fig.suptitle("NOPHO CNV sensitivity reanalysis (GRCh37/hg19)", fontsize=13)
    fig.savefig(output_path, dpi=300)
    plt.close(fig)


def main() -> None:
    args = parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    files = sorted(args.source_root.glob("*/*.cns"))
    if len(files) != 262:
        raise RuntimeError(f"Expected 262 CNS files, found {len(files)}")
    if not args.blacklist.exists():
        raise FileNotFoundError(args.blacklist)

    blacklist = load_blacklist(args.blacklist)
    input_manifest = pd.DataFrame(
        [
            {
                "sample": path.name.replace(".tumor.merged.cns", ""),
                "relative_path": path.relative_to(args.source_root).as_posix(),
                "size_bytes": path.stat().st_size,
                "sha256": sha256(path),
            }
            for path in files
        ]
    )
    input_manifest.to_csv(args.output / "input_cns_manifest.csv", index=False)

    official_rows: list[dict[str, object]] = []
    raw_rows: list[dict[str, object]] = []
    event_tables: list[pd.DataFrame] = []
    sensitivity_rows: list[dict[str, object]] = []
    for path in files:
        sample = path.name.replace(".tumor.merged.cns", "")
        segments = load_segments(path)
        events = call_and_merge(segments, blacklist)
        if not events.empty:
            table = events.copy()
            table.insert(0, "sample", sample)
            event_tables.append(table)

        official: dict[str, object] = {
            "sample": sample,
            "del6q": bool(((events["chromosome"] == "6") & (events["type"] == "DEL")).any()),
        }
        raw: dict[str, object] = {"sample": sample}
        for gene in GENES:
            official[gene] = overlaps(events, gene)
            raw[gene] = raw_segment_overlaps(segments, gene)
        official_rows.append(official)
        raw_rows.append(raw)

        for threshold in SENSITIVITY_THRESHOLDS:
            for min_size in SENSITIVITY_MIN_SIZES:
                sensitivity_events = call_and_merge(
                    segments,
                    blacklist,
                    del_threshold=threshold,
                    min_size=min_size,
                )
                sensitivity_rows.append(
                    {
                        "sample": sample,
                        "del_threshold": threshold,
                        "min_size": min_size,
                        "del6q": bool(
                            (
                                (sensitivity_events["chromosome"] == "6")
                                & (sensitivity_events["type"] == "DEL")
                            ).any()
                        ),
                        **{gene: overlaps(sensitivity_events, gene) for gene in GENES},
                    }
                )

    official = pd.DataFrame(official_rows).sort_values("sample")
    raw = pd.DataFrame(raw_rows).sort_values("sample")
    events = pd.concat(event_tables, ignore_index=True)
    sensitivity = pd.DataFrame(sensitivity_rows)
    official.to_csv(args.output / "sample_gene_deletion_matrix_grch37_ge1Mb.csv", index=False)
    raw.to_csv(args.output / "sample_gene_deletion_matrix_grch37_raw_segments.csv", index=False)
    events.to_csv(args.output / "merged_filtered_events_grch37_ge1Mb.csv", index=False)
    sensitivity.to_csv(args.output / "threshold_size_sensitivity_matrix_grch37.csv", index=False)

    frequencies = pd.concat(
        [
            summarize_matrix(official, "grch37_threshold_-0.4_ge1Mb"),
            summarize_matrix(raw, "grch37_raw_segment_threshold_-0.4"),
        ],
        ignore_index=True,
    )
    frequencies.to_csv(args.output / "gene_deletion_frequencies_grch37.csv", index=False)

    foxo = official["FOXO3"]
    tnf = official["TNFAIP3"]
    del6q = official["del6q"]
    evaluated_6q = ["BACH2", "EPHA7", "GRIK2", "HACE1", "PRDM1", "FOXO3", "MYB", "TNFAIP3", "LATS1"]
    patterns = {
        "total_samples": len(official),
        "del6q_grch37_ge1Mb": int(del6q.sum()),
        "FOXO3_deleted_grch37_ge1Mb": int(foxo.sum()),
        "TNFAIP3_deleted_grch37_ge1Mb": int(tnf.sum()),
        "FOXO3_and_TNFAIP3": int((foxo & tnf).sum()),
        "FOXO3_not_TNFAIP3": int((foxo & ~tnf).sum()),
        "TNFAIP3_not_FOXO3": int((tnf & ~foxo).sum()),
        "FOXO3_only_among_9_evaluated_6q_genes": int(
            (foxo & ~official[[gene for gene in evaluated_6q if gene != "FOXO3"]].any(axis=1)).sum()
        ),
        "del6q_involving_FOXO3": int((del6q & foxo).sum()),
    }
    (args.output / "patterns_grch37.json").write_text(json.dumps(patterns, indent=2), encoding="utf-8")

    cooccurrence_rows = []
    for gene in ("CDKN2A", "PAX5", "IKZF1"):
        table = pd.crosstab(official["del6q"], official[gene]).reindex(
            index=[False, True], columns=[False, True], fill_value=0
        )
        odds_ratio, p_value = fisher_exact(table.to_numpy())
        cooccurrence_rows.append(
            {
                "gene": gene,
                "deleted_n": int(official[gene].sum()),
                "deleted_percent": 100 * official[gene].mean(),
                "table_no_del6q_no_gene": int(table.loc[False, False]),
                "table_no_del6q_gene": int(table.loc[False, True]),
                "table_del6q_no_gene": int(table.loc[True, False]),
                "table_del6q_gene": int(table.loc[True, True]),
                "fisher_odds_ratio": odds_ratio,
                "fisher_p": p_value,
            }
        )
    cooccurrence = pd.DataFrame(cooccurrence_rows)
    cooccurrence["bh_fdr"] = multipletests(cooccurrence["fisher_p"], method="fdr_bh")[1]
    cooccurrence.to_csv(args.output / "del6q_cooccurrence_fisher_grch37.csv", index=False)

    sensitivity_summary = (
        sensitivity.groupby(["del_threshold", "min_size"])[["del6q", *GENES]]
        .sum()
        .reset_index()
    )
    sensitivity_summary.to_csv(args.output / "threshold_size_sensitivity_summary_grch37.csv", index=False)
    plot_nopho_summary(official, sensitivity_summary, args.output / "nopho_grch37_sensitivity.png")

    provenance = {
        "analysis_label": "official-workflow-compatible GRCh37/hg19 reanalysis",
        "coordinate_build": "GRCh37/hg19",
        "source_root": "user-supplied NOPHO raw root; files listed by relative_path in input_cns_manifest.csv",
        "input_cns_files": len(files),
        "input_manifest": "input_cns_manifest.csv",
        "blacklist": args.blacklist.name,
        "blacklist_sha256": sha256(args.blacklist),
        "blacklist_source": "https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/hg19-blacklist.v2.bed.gz",
        "official_code_commit": OFFICIAL_CODE_COMMIT,
        "official_parameters": {
            "deletion_log2_lte": -0.4,
            "gain_log2_gte": 0.3,
            "merge_same_type_gap_lt_bp": 10_000,
            "minimum_event_size_bp": 1_000_000,
            "minimum_event_size_inclusive": True,
            "discard_if_single_blacklist_interval_overlap_fraction_gte": 0.5,
            "chromosomes": "autosomes 1-22",
        },
        "gene_coordinates": {
            "source": "Ensembl GRCh37 REST API",
            "retrieved": "2026-06-12",
            "genes": GENES,
        },
    }
    (args.output / "provenance_grch37.json").write_text(
        json.dumps(provenance, indent=2), encoding="utf-8"
    )

    print(json.dumps(patterns, indent=2))
    print("\nGRCh37 gene frequencies:")
    print(frequencies.loc[frequencies["method"] == "grch37_threshold_-0.4_ge1Mb"].to_string(index=False))
    print("\nCo-occurrence with BH-FDR:")
    print(cooccurrence.to_string(index=False))


if __name__ == "__main__":
    main()
