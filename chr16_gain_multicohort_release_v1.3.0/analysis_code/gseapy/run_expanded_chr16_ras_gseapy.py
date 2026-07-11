from __future__ import annotations

import csv
import gzip
import json
import math
import re
import tarfile
from io import BytesIO
from pathlib import Path

import gseapy as gp
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats


ROOT = Path(__file__).resolve().parent
PROJECT = ROOT.parent
OLD_GSEAPY = PROJECT / "chr16_ras_gseapy_closure_20260710"
GEO = PROJECT / "expanded_database_validation_20260630" / "geo_supp"
EXPANDED_VALIDATION = PROJECT / "expanded_database_validation_20260630"

TABLES = ROOT / "tables"
FIGURES = ROOT / "figures"
GSEA_OUT = ROOT / "gseapy_prerank"
GENESETS = ROOT / "gene_sets"
for d in [TABLES, FIGURES, GSEA_OUT, GENESETS]:
    d.mkdir(parents=True, exist_ok=True)


GSE181 = GEO / "GSE181157"
GSE181_RAW_TAR = GSE181 / "GSE181157_RAW.tar"
GSE181_META_CSV = GSE181 / "GSE181157_SampleMetadata.csv"
GSE227832_COUNTS = GEO / "GSE227832_RNAseq_read_counts.txt.gz"
GSE227832_MATRICES = sorted(GEO.glob("GSE227832-GPL*_series_matrix.txt.gz"))
GSE228632_COUNTS = GEO / "GSE228632_RNAseq_read_counts.txt.gz"
GSE228632_MATRIX = GEO / "GSE228632_series_matrix.txt.gz"
GSE87070_MATRIX = GEO / "GSE87070_series_matrix.txt.gz"
GPL570_ANNOT = PROJECT / "CARS_reviewer_response_external_validation_20260629" / "geo_download" / "GPL570.annot.gz"


RAS_PATHWAY_GENES = [
    "KRAS", "NRAS", "HRAS", "NF1", "PTPN11", "BRAF", "RAF1", "ARAF",
    "MAP2K1", "MAP2K2", "MAPK1", "MAPK3", "SOS1", "SOS2", "GRB2", "SHC1",
    "RASA1", "RASA2", "RASGRP1", "RASGRP2", "RASGRP3", "CBL", "CBLB",
    "DUSP4", "DUSP5", "DUSP6", "SPRY1", "SPRY2", "SPRY4",
    "ELK1", "ETS1", "ETS2", "FOS", "FOSB", "JUN", "JUNB", "JUND",
    "ETV4", "ETV5", "EGR1", "EGR2", "EGR3", "CCND1",
]

MAPK_ERK_FEEDBACK_GENES = [
    "DUSP4", "DUSP5", "DUSP6", "DUSP7", "SPRY1", "SPRY2", "SPRY4",
    "ETV4", "ETV5", "FOS", "FOSB", "JUN", "JUNB", "JUND",
    "EGR1", "EGR2", "EGR3", "CCND1", "MYC", "IER3", "BTG2",
]

B_EARLY_PROGENITOR_GENES = [
    "CD34", "KIT", "FLT3", "IL7R", "DNTT", "RAG1", "RAG2", "VPREB1",
    "VPREB3", "IGLL1", "MME", "PAX5", "EBF1", "TCF3", "IKZF1", "LEF1",
]

B_PRE_BCR_PRE_B_GENES = [
    "CD19", "CD79A", "CD79B", "BLNK", "VPREB1", "VPREB3", "IGLL1",
    "PAX5", "EBF1", "TCF3", "RAG1", "RAG2", "MME", "BCL11A", "POU2AF1",
]

B_MATURE_B_CELL_GENES = [
    "MS4A1", "CD22", "CD37", "CD40", "CD74", "BANK1", "BLK", "SPIB",
    "FCRL1", "FCRL5", "CR2", "CD83", "TNFRSF13B", "TNFRSF13C", "AICDA",
]

B_STEM_CELL_LIKE_GENES = [
    "PROM1", "CD34", "KIT", "FLT3", "MECOM", "GATA2", "LMO2", "HHEX",
    "MEIS1", "HOXA9", "ERG", "SOX4", "MYB", "BMI1", "MPL",
]

FOCUS_GENESETS = [
    "CHR16_DOSAGE_EXPRESSED",
    "CHR16_16P_MAPK3_REGION",
    "RAS_MAPK_CORE",
    "MAPK_ERK_FEEDBACK",
    "B_ALL_EARLY_PROGENITOR",
    "B_ALL_PRE_BCR_PRE_B",
    "B_ALL_MATURE_B",
    "B_ALL_STEM_CELL_LIKE",
]

TERM_LABELS = {
    "CHR16_DOSAGE_EXPRESSED": "chr16 coding",
    "CHR16_16P_MAPK3_REGION": "16p/MAPK3 region",
    "RAS_MAPK_CORE": "RAS-MAPK core",
    "MAPK_ERK_FEEDBACK": "ERK feedback",
    "B_ALL_EARLY_PROGENITOR": "B early progenitor",
    "B_ALL_PRE_BCR_PRE_B": "pre-BCR/pre-B",
    "B_ALL_MATURE_B": "mature B",
    "B_ALL_STEM_CELL_LIKE": "stem-like",
}


def clean_text(value: object) -> str:
    if pd.isna(value):
        return ""
    text = str(value).strip()
    if len(text) >= 2 and text[0] == text[-1] == '"':
        text = text[1:-1]
    return text.strip()


def safe_name(text: str) -> str:
    text = text.lower().strip()
    text = re.sub(r"[^a-z0-9]+", "_", text)
    text = re.sub(r"_+", "_", text)
    return text.strip("_")


def split_geo_line(line: str) -> list[str]:
    return next(csv.reader([line.rstrip("\n\r")], delimiter="\t", quotechar='"'))


def has_chr16_gain(text: object) -> bool:
    value = clean_text(text)
    if not value:
        return False
    patterns = [
        r"\+16(?!\d)",
        r"trisomy\s*16",
        r"gain\s*(?:of\s*)?(?:chromosome\s*)?16",
        r"chr(?:omosome)?\s*16\s*gain",
    ]
    return any(re.search(p, value, flags=re.I) for p in patterns)


def has_ras_pathway_mutation(text: object) -> bool:
    value = clean_text(text)
    if not value:
        return False
    genes = ["KRAS", "NRAS", "NF1", "PTPN11", "FLT3", "BRAF", "CBL", "MAP2K1", "MAP2K2", "HRAS"]
    return any(re.search(rf"\b{re.escape(g)}\b", value, flags=re.I) for g in genes)


def parse_series_matrix_metadata(paths: list[Path], dataset: str) -> pd.DataFrame:
    rows = []
    for path in paths:
        sample_fields: dict[str, list[str]] = {}
        characteristic_lines: list[list[str]] = []
        with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
            for line in handle:
                if line.startswith("!series_matrix_table_begin"):
                    break
                if not line.startswith("!Sample_"):
                    continue
                parts = split_geo_line(line)
                key = parts[0].lstrip("!")
                values = [clean_text(v) for v in parts[1:]]
                if key == "Sample_characteristics_ch1":
                    characteristic_lines.append(values)
                else:
                    sample_fields[key] = values
        if not sample_fields:
            continue
        n = max(len(v) for v in sample_fields.values())
        for i in range(n):
            row = {"dataset": dataset, "series_matrix_file": path.name}
            for key, values in sample_fields.items():
                row[safe_name(key.replace("Sample_", ""))] = values[i] if i < len(values) else ""
            for values in characteristic_lines:
                if i >= len(values):
                    continue
                val = values[i]
                if ":" in val:
                    k, v = val.split(":", 1)
                    k = safe_name(k)
                    v = v.strip()
                    if k in row and row[k] and row[k] != v:
                        row[k] = f"{row[k]}; {v}"
                    else:
                        row[k] = v
                else:
                    key = f"characteristic_{len([x for x in row if x.startswith('characteristic_')]) + 1}"
                    row[key] = val
            title = row.get("title", "")
            bracket = re.search(r"\[([^\]]+)\]", title)
            row["sample_id"] = bracket.group(1).strip() if bracket else (title or row.get("geo_accession", ""))
            row["geo_accession"] = row.get("geo_accession", "")
            rows.append(row)
    meta = pd.DataFrame(rows)
    if meta.empty:
        return meta
    meta = meta.drop_duplicates(subset=["sample_id", "geo_accession", "series_matrix_file"])
    return meta


def read_gse181_metadata() -> pd.DataFrame:
    meta = pd.read_csv(GSE181_META_CSV, header=1, encoding="utf-8-sig")
    meta.columns = [re.sub(r"\s+", " ", str(c)).strip() for c in meta.columns]
    meta = meta[meta["DFCI ID"].notna()].copy()
    meta["dataset"] = "GSE181157"
    meta["sample_id"] = meta["DFCI ID"].astype(str).str.strip()
    meta["predicted_subtype"] = meta["Predicted subtype by RNA-seq"].astype(str).str.strip()
    meta["conventional_subtype"] = meta["Conventional methods subtype"].astype(str).str.strip()
    meta["diagnosis"] = meta["Diagnosis"].astype(str).str.strip()
    meta["fusion_etv6_runx1"] = meta["Fusions"].astype(str).str.contains("ETV6-RUNX1", case=False, na=False)
    meta["etv6_runx1_or_like"] = (
        meta["fusion_etv6_runx1"]
        | meta["predicted_subtype"].str.contains(r"ETV6[-:]?RUNX1|ETV6::RUNX1", case=False, regex=True, na=False)
    )
    meta["ras_pathway_mutation"] = meta["Mutations"].apply(has_ras_pathway_mutation)
    meta["karyotype_chr16_gain"] = meta["Karyotype"].apply(has_chr16_gain)
    meta["fish_chr16_gain"] = meta["Conventional cytogenetics (FISH)"].apply(has_chr16_gain)
    meta["chr16_gain_text"] = meta["karyotype_chr16_gain"] | meta["fish_chr16_gain"]
    meta["is_t_all"] = (
        meta["predicted_subtype"].str.contains("T-ALL", case=False, na=False)
        | meta["conventional_subtype"].str.contains("T-ALL", case=False, na=False)
        | meta["diagnosis"].str.contains("Pre-T", case=False, na=False)
    )
    meta["is_b_all"] = ~meta["is_t_all"]
    meta["is_diagnosis"] = True
    meta["include_b_all_reference"] = meta["is_b_all"]
    return meta


def read_gse181_counts() -> pd.DataFrame:
    series = []
    with tarfile.open(GSE181_RAW_TAR, "r") as tf:
        members = [m for m in tf.getmembers() if m.name.endswith(".counts.txt.gz")]
        for member in members:
            sample = re.sub(r"^GSM\d+_", "", member.name)
            sample = re.sub(r"\.counts\.txt\.gz$", "", sample)
            raw = tf.extractfile(member).read()
            text = gzip.decompress(raw)
            s = pd.read_csv(BytesIO(text), sep="\t", header=None, names=["gene", sample])
            s["gene"] = s["gene"].astype(str).str.replace(r"\.\d+$", "", regex=True)
            s[sample] = pd.to_numeric(s[sample], errors="coerce").fillna(0.0).astype(float)
            series.append(s.set_index("gene")[sample])
    counts = pd.concat(series, axis=1).groupby(level=0).sum()
    counts.index.name = "gene"
    return counts


def read_count_matrix(path: Path) -> pd.DataFrame:
    counts = pd.read_csv(path, sep="\t", compression="gzip")
    counts = counts.rename(columns={counts.columns[0]: "gene"})
    counts["gene"] = counts["gene"].astype(str).str.replace(r"\.\d+$", "", regex=True)
    counts = counts.set_index("gene")
    counts = counts.apply(pd.to_numeric, errors="coerce").fillna(0.0)
    counts = counts.groupby(level=0).sum()
    return counts


def log2_cpm(counts: pd.DataFrame) -> tuple[pd.DataFrame, pd.Series]:
    lib_sizes = counts.sum(axis=0)
    cpm = counts.div(lib_sizes.replace(0, np.nan), axis=1) * 1_000_000
    return np.log2(cpm.fillna(0) + 1.0), lib_sizes


def standardize_gse227832_metadata() -> pd.DataFrame:
    meta = parse_series_matrix_metadata(GSE227832_MATRICES, "GSE227832")
    if meta.empty:
        return meta
    pred = meta.get("revised_predicted_genotype", pd.Series("", index=meta.index)).fillna("")
    genotype = meta.get("genotype", pd.Series("", index=meta.index)).fillna("")
    meta["predicted_subtype"] = np.where(pred.astype(str).str.len() > 0, pred, genotype)
    disease = meta.get("disease", pd.Series("", index=meta.index)).astype(str)
    state = meta.get("disease_state", pd.Series("", index=meta.index)).astype(str)
    meta["is_b_all"] = disease.str.contains("B acute lymphoblastic leukemia", case=False, na=False)
    meta["is_t_all"] = disease.str.contains("T acute lymphoblastic leukemia", case=False, na=False) | meta[
        "predicted_subtype"
    ].astype(str).str.contains("T-ALL", case=False, na=False)
    meta["is_diagnosis"] = state.str.contains("diagnosis", case=False, na=False) | meta["title"].astype(str).str.contains(
        "diagnosis", case=False, na=False
    )
    meta["etv6_runx1_or_like"] = meta["predicted_subtype"].astype(str).str.contains(
        r"ETV6::RUNX1|ETV6-RUNX1|ETV6.*RUNX1", case=False, regex=True, na=False
    )
    meta["include_b_all_reference"] = meta["is_b_all"] & ~meta["is_t_all"] & meta["is_diagnosis"]
    meta["chr16_gain_text"] = False
    meta["ras_pathway_mutation"] = False
    return meta


def standardize_gse228632_metadata() -> pd.DataFrame:
    meta = parse_series_matrix_metadata([GSE228632_MATRIX], "GSE228632")
    if meta.empty:
        return meta
    pred = meta.get("predicted_subtype", pd.Series("", index=meta.index)).fillna("")
    meta["predicted_subtype"] = pred
    disease = meta.get("disease", pd.Series("", index=meta.index)).astype(str)
    state = meta.get("disease_state", pd.Series("", index=meta.index)).astype(str)
    meta["is_b_all"] = disease.str.contains("B acute lymphoblastic leukemia", case=False, na=False)
    meta["is_t_all"] = disease.str.contains("T acute lymphoblastic leukemia", case=False, na=False) | pred.astype(
        str
    ).str.contains("T-ALL", case=False, na=False)
    meta["is_diagnosis"] = state.str.contains("diagnosis", case=False, na=False) | meta["title"].astype(str).str.contains(
        "diagnosis", case=False, na=False
    )
    meta["etv6_runx1_or_like"] = pred.astype(str).str.contains(
        r"ETV6::RUNX1|ETV6-RUNX1|ETV6.*RUNX1", case=False, regex=True, na=False
    )
    meta["include_b_all_reference"] = meta["is_b_all"] & ~meta["is_t_all"] & meta["is_diagnosis"]
    meta["chr16_gain_text"] = False
    meta["ras_pathway_mutation"] = False
    return meta


def standardize_gse87070_metadata() -> pd.DataFrame:
    meta = parse_series_matrix_metadata([GSE87070_MATRIX], "GSE87070")
    if meta.empty:
        return meta
    meta["sample_id"] = meta["geo_accession"]
    meta["predicted_subtype"] = meta.get("all_subtype", pd.Series("", index=meta.index)).fillna("")
    meta["is_b_all"] = meta["predicted_subtype"].isin(["ETV6-RUNX1", "Other BCP"])
    meta["is_t_all"] = meta["predicted_subtype"].str.contains("T", case=False, na=False)
    meta["is_diagnosis"] = True
    meta["etv6_runx1_or_like"] = meta["predicted_subtype"].eq("ETV6-RUNX1")
    meta["include_b_all_reference"] = meta["predicted_subtype"].isin(["ETV6-RUNX1", "Other BCP"])
    meta["chr16_gain_text"] = False
    meta["ras_pathway_mutation"] = False
    return meta


def read_gpl570_annotation() -> pd.DataFrame:
    header_line = None
    with gzip.open(GPL570_ANNOT, "rt", encoding="utf-8", errors="replace") as handle:
        for i, line in enumerate(handle):
            if line.startswith("ID\t"):
                header_line = i
                break
    if header_line is None:
        raise RuntimeError(f"Cannot find GPL570 annotation header in {GPL570_ANNOT}")
    annot = pd.read_csv(GPL570_ANNOT, sep="\t", compression="gzip", skiprows=header_line, dtype=str)
    annot = annot[["ID", "Gene symbol"]].copy()
    annot["gene_symbol"] = (
        annot["Gene symbol"].fillna("").astype(str).str.split(r"\s*///\s*", regex=True).str[0].str.strip()
    )
    annot = annot[(annot["ID"].notna()) & (annot["gene_symbol"].str.len() > 0)]
    annot = annot[~annot["gene_symbol"].isin(["---", "nan"])]
    annot = annot.drop_duplicates(subset=["ID"])
    return annot[["ID", "gene_symbol"]]


def series_matrix_table_info(path: Path) -> tuple[int, int | None]:
    data_rows: int | None = None
    header_line: int | None = None
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for i, line in enumerate(handle):
            if line.startswith("!Sample_data_row_count"):
                parts = split_geo_line(line)
                counts = [clean_text(x) for x in parts[1:] if clean_text(x).isdigit()]
                if counts:
                    data_rows = int(counts[0])
            if line.startswith("!series_matrix_table_begin"):
                header_line = i + 1
                break
    if header_line is None:
        raise RuntimeError(f"Cannot find expression table in {path}")
    return header_line, data_rows


def read_gse87070_symbol_expression() -> pd.DataFrame:
    header_line, data_rows = series_matrix_table_info(GSE87070_MATRIX)
    expr = pd.read_csv(
        GSE87070_MATRIX,
        sep="\t",
        compression="gzip",
        skiprows=header_line,
        nrows=data_rows,
    )
    expr = expr.rename(columns={expr.columns[0]: "ID"})
    annot = read_gpl570_annotation()
    expr = expr.merge(annot, on="ID", how="inner")
    numeric_cols = [c for c in expr.columns if c not in {"ID", "gene_symbol"}]
    expr[numeric_cols] = expr[numeric_cols].apply(pd.to_numeric, errors="coerce")
    symbol_expr = expr.groupby("gene_symbol")[numeric_cols].median()
    symbol_expr = symbol_expr.replace([np.inf, -np.inf], np.nan).dropna(how="all")
    symbol_expr.index.name = "gene"
    return symbol_expr


def read_cached_mapping() -> pd.DataFrame:
    path = OLD_GSEAPY / "tables" / "biomart_symbol_mapping.csv"
    mapping = pd.read_csv(path)
    mapping["ensembl_gene_id"] = mapping["ensembl_gene_id"].astype(str).str.replace(r"\.\d+$", "", regex=True)
    return mapping


def read_cached_chr16() -> pd.DataFrame:
    path = OLD_GSEAPY / "tables" / "biomart_chr16_genes.csv"
    chr16 = pd.read_csv(path)
    chr16["ensembl_gene_id"] = chr16["ensembl_gene_id"].astype(str).str.replace(r"\.\d+$", "", regex=True)
    return chr16


def symbol_set_to_ensembl(mapping: pd.DataFrame, symbols: list[str], available: set[str]) -> list[str]:
    sub = mapping[mapping["external_gene_name"].isin(symbols)].copy()
    coding = sub[sub["gene_biotype"].eq("protein_coding")]
    if not coding.empty:
        sub = coding
    genes = set(sub["ensembl_gene_id"].astype(str)) & available
    return sorted(genes)


def write_gmt(gene_sets: dict[str, list[str]], path: Path) -> None:
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for term, genes in gene_sets.items():
            handle.write("\t".join([term, "custom_curated"] + sorted(set(map(str, genes)))) + "\n")


def build_rna_gene_sets(expr: pd.DataFrame) -> dict[str, list[str]]:
    mapping = read_cached_mapping()
    chr16 = read_cached_chr16()
    available = set(expr.index.astype(str))
    expressed = set(expr.index[(expr >= 1.0).sum(axis=1) >= max(5, int(math.ceil(0.05 * expr.shape[1])) )].astype(str))
    chr16 = chr16[chr16["gene_biotype"].eq("protein_coding")].copy()
    chr16["start_position"] = pd.to_numeric(chr16["start_position"], errors="coerce")
    chr16["end_position"] = pd.to_numeric(chr16["end_position"], errors="coerce")
    chr16_all = sorted(set(chr16["ensembl_gene_id"].astype(str)) & available & expressed)
    chr16_region = chr16[
        chr16["start_position"].between(25_000_000, 35_000_000)
        | chr16["end_position"].between(25_000_000, 35_000_000)
    ]
    gene_sets = {
        "CHR16_DOSAGE_EXPRESSED": chr16_all,
        "CHR16_16P_MAPK3_REGION": sorted(set(chr16_region["ensembl_gene_id"].astype(str)) & available & expressed),
        "RAS_MAPK_CORE": symbol_set_to_ensembl(mapping, RAS_PATHWAY_GENES, available),
        "MAPK_ERK_FEEDBACK": symbol_set_to_ensembl(mapping, MAPK_ERK_FEEDBACK_GENES, available),
        "B_ALL_EARLY_PROGENITOR": symbol_set_to_ensembl(mapping, B_EARLY_PROGENITOR_GENES, available),
        "B_ALL_PRE_BCR_PRE_B": symbol_set_to_ensembl(mapping, B_PRE_BCR_PRE_B_GENES, available),
        "B_ALL_MATURE_B": symbol_set_to_ensembl(mapping, B_MATURE_B_CELL_GENES, available),
        "B_ALL_STEM_CELL_LIKE": symbol_set_to_ensembl(mapping, B_STEM_CELL_LIKE_GENES, available),
    }
    return {k: v for k, v in gene_sets.items() if len(v) >= 5}


def build_symbol_gene_sets(expr: pd.DataFrame) -> dict[str, list[str]]:
    chr16 = read_cached_chr16()
    chr16 = chr16[chr16["gene_biotype"].eq("protein_coding")].dropna(subset=["external_gene_name"]).copy()
    chr16["start_position"] = pd.to_numeric(chr16["start_position"], errors="coerce")
    chr16["end_position"] = pd.to_numeric(chr16["end_position"], errors="coerce")
    available = set(expr.index.astype(str))
    expressed = set(expr.index[expr.notna().sum(axis=1) >= max(5, int(math.ceil(0.05 * expr.shape[1])) )].astype(str))
    chr16_symbols = sorted(set(chr16["external_gene_name"].astype(str)) & available & expressed)
    chr16_region = chr16[
        chr16["start_position"].between(25_000_000, 35_000_000)
        | chr16["end_position"].between(25_000_000, 35_000_000)
    ]
    gene_sets = {
        "CHR16_DOSAGE_EXPRESSED": chr16_symbols,
        "CHR16_16P_MAPK3_REGION": sorted(set(chr16_region["external_gene_name"].astype(str)) & available & expressed),
        "RAS_MAPK_CORE": sorted(set(RAS_PATHWAY_GENES) & available),
        "MAPK_ERK_FEEDBACK": sorted(set(MAPK_ERK_FEEDBACK_GENES) & available),
        "B_ALL_EARLY_PROGENITOR": sorted(set(B_EARLY_PROGENITOR_GENES) & available),
        "B_ALL_PRE_BCR_PRE_B": sorted(set(B_PRE_BCR_PRE_B_GENES) & available),
        "B_ALL_MATURE_B": sorted(set(B_MATURE_B_CELL_GENES) & available),
        "B_ALL_STEM_CELL_LIKE": sorted(set(B_STEM_CELL_LIKE_GENES) & available),
    }
    return {k: v for k, v in gene_sets.items() if len(v) >= 5}


def filter_expr(expr: pd.DataFrame, gene_sets: dict[str, list[str]]) -> pd.DataFrame:
    expr = expr.copy()
    expr.index = expr.index.astype(str)
    keep = set(expr.index[expr.var(axis=1, skipna=True) > 1e-8])
    for genes in gene_sets.values():
        keep.update(genes)
    out = expr.loc[sorted(keep & set(expr.index))]
    return out.replace([np.inf, -np.inf], np.nan).dropna(how="all")


def filter_rna_expr(expr: pd.DataFrame, gene_sets: dict[str, list[str]]) -> pd.DataFrame:
    expr = expr.copy()
    expr.index = expr.index.astype(str)
    min_samples = max(5, int(math.ceil(0.05 * expr.shape[1])))
    expressed = (expr >= 1.0).sum(axis=1) >= min_samples
    keep = set(expr.index[expressed])
    for genes in gene_sets.values():
        keep.update(genes)
    out = expr.loc[sorted(keep & set(expr.index))]
    return out.replace([np.inf, -np.inf], np.nan).dropna(how="all")


def rank_genes(expr: pd.DataFrame, pos: list[str], neg: list[str]) -> pd.DataFrame:
    x = expr[pos].to_numpy(dtype=float)
    y = expr[neg].to_numpy(dtype=float)
    nx, ny = x.shape[1], y.shape[1]
    mean_x = np.nanmean(x, axis=1)
    mean_y = np.nanmean(y, axis=1)
    var_x = np.nanvar(x, axis=1, ddof=1)
    var_y = np.nanvar(y, axis=1, ddof=1)
    denom = np.sqrt(var_x / nx + var_y / ny)
    stat = np.divide(mean_x - mean_y, denom, out=np.zeros_like(denom), where=denom > 0)
    rank = pd.DataFrame(
        {
            "gene": expr.index.astype(str),
            "rank_stat": stat,
            "mean_positive": mean_x,
            "mean_negative": mean_y,
            "mean_difference": mean_x - mean_y,
        }
    )
    rank = rank.replace([np.inf, -np.inf], np.nan).dropna(subset=["rank_stat"])
    return rank.sort_values("rank_stat", ascending=False)


def bh_fdr(pvals: list[float]) -> list[float]:
    p = np.asarray([np.nan if pd.isna(x) else float(x) for x in pvals], dtype=float)
    q = np.full_like(p, np.nan, dtype=float)
    mask = ~np.isnan(p)
    if not mask.any():
        return q.tolist()
    vals = p[mask]
    order = np.argsort(vals)
    ranked = vals[order]
    n = len(ranked)
    adj = ranked * n / (np.arange(n) + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    adj = np.clip(adj, 0, 1)
    out = np.empty(n)
    out[order] = adj
    q[mask] = out
    return q.tolist()


def cliffs_delta(x: np.ndarray, y: np.ndarray) -> float:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    x = x[~np.isnan(x)]
    y = y[~np.isnan(y)]
    if len(x) == 0 or len(y) == 0:
        return np.nan
    diff = np.sign(x[:, None] - y[None, :]).sum()
    return float(diff / (len(x) * len(y)))


def get_comparisons(dataset: str, meta: pd.DataFrame, expr: pd.DataFrame) -> list[dict[str, object]]:
    m = meta.drop_duplicates(subset=["sample_id"]).set_index("sample_id")
    m = m.loc[[s for s in expr.columns if s in m.index]].copy()
    comparisons: list[dict[str, object]] = []
    if dataset == "GSE181157":
        configs = [
            ("B_ALL_chr16_gain_vs_no_gain", m["is_b_all"], m["chr16_gain_text"], "chr16 gain", "no chr16 gain"),
            ("B_ALL_RASmut_vs_RASwt", m["is_b_all"], m["ras_pathway_mutation"], "RAS/MAPK mutated", "RAS/MAPK wild-type"),
            (
                "ETV6_like_RASmut_vs_RASwt",
                m["etv6_runx1_or_like"],
                m["ras_pathway_mutation"],
                "ETV6/RAS mutated",
                "ETV6/RAS wild-type",
            ),
        ]
    elif dataset == "GSE87070":
        configs = [
            (
                "GSE87070_ETV6_RUNX1_vs_other_BCP",
                m["include_b_all_reference"],
                m["etv6_runx1_or_like"],
                "ETV6-RUNX1",
                "Other BCP",
            )
        ]
    else:
        configs = [
            (
                f"{dataset}_ETV6_RUNX1_or_like_vs_other_BALL",
                m["include_b_all_reference"],
                m["etv6_runx1_or_like"],
                "ETV6::RUNX1/like",
                "other B-ALL",
            )
        ]
    for name, subset, group, pos_label, neg_label in configs:
        available = list(m.index[subset.fillna(False)])
        available = [s for s in available if s in expr.columns]
        pos = [s for s in available if bool(group.loc[s])]
        neg = [s for s in available if not bool(group.loc[s])]
        comparisons.append(
            {
                "dataset": dataset,
                "comparison": name,
                "positive_label": pos_label,
                "negative_label": neg_label,
                "positive": pos,
                "negative": neg,
            }
        )
    return comparisons


def run_prerank_for_dataset(
    dataset: str,
    expr: pd.DataFrame,
    meta: pd.DataFrame,
    gene_sets: dict[str, list[str]],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows = []
    rank_rows = []
    for cfg in get_comparisons(dataset, meta, expr):
        pos = cfg["positive"]
        neg = cfg["negative"]
        comp = str(cfg["comparison"])
        if len(pos) < 2 or len(neg) < 2:
            continue
        rank = rank_genes(expr, pos, neg)
        rank.insert(0, "dataset", dataset)
        rank.insert(1, "comparison", comp)
        rank.to_csv(TABLES / f"ranking_{comp}.csv", index=False, encoding="utf-8-sig")
        rank_rows.append(rank.head(500))
        outdir = GSEA_OUT / comp
        outdir.mkdir(parents=True, exist_ok=True)
        pre = gp.prerank(
            rnk=rank[["gene", "rank_stat"]],
            gene_sets=gene_sets,
            outdir=str(outdir),
            min_size=5,
            max_size=3000,
            permutation_num=1000,
            threads=4,
            no_plot=True,
            seed=20260710,
            verbose=False,
        )
        res = pre.res2d.copy()
        res.insert(0, "dataset", dataset)
        res.insert(1, "comparison", comp)
        res.insert(2, "positive_label", cfg["positive_label"])
        res.insert(3, "negative_label", cfg["negative_label"])
        res.insert(4, "n_positive", len(pos))
        res.insert(5, "n_negative", len(neg))
        rows.append(res)
    result = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()
    ranks = pd.concat(rank_rows, ignore_index=True) if rank_rows else pd.DataFrame()
    return result, ranks


def run_ssgsea(expr: pd.DataFrame, gene_sets: dict[str, list[str]]) -> pd.DataFrame:
    data = expr.copy()
    data.index.name = "gene_name"
    data = data.reset_index()
    ss = gp.ssgsea(
        data=data,
        gene_sets=gene_sets,
        outdir=None,
        sample_norm_method="rank",
        correl_norm_type="rank",
        min_size=5,
        max_size=3000,
        threads=4,
        no_plot=True,
        seed=20260710,
        verbose=False,
    )
    res = ss.res2d.copy()
    if {"Name", "Term", "NES"}.issubset(res.columns):
        scores = res.pivot(index="Name", columns="Term", values="NES")
    elif {"Name", "Term", "ES"}.issubset(res.columns):
        scores = res.pivot(index="Name", columns="Term", values="ES")
    else:
        scores = res
    scores = scores.apply(pd.to_numeric, errors="coerce")
    scores.index.name = "sample_id"
    return scores


def compare_ssgsea_for_dataset(dataset: str, scores: pd.DataFrame, meta: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for cfg in get_comparisons(dataset, meta, scores.T):
        pos = [s for s in cfg["positive"] if s in scores.index]
        neg = [s for s in cfg["negative"] if s in scores.index]
        comp = str(cfg["comparison"])
        for term in scores.columns:
            x = scores.loc[pos, term].astype(float).dropna().to_numpy()
            y = scores.loc[neg, term].astype(float).dropna().to_numpy()
            p = np.nan
            if len(x) >= 2 and len(y) >= 2:
                p = stats.mannwhitneyu(x, y, alternative="two-sided").pvalue
            rows.append(
                {
                    "dataset": dataset,
                    "comparison": comp,
                    "positive_label": cfg["positive_label"],
                    "negative_label": cfg["negative_label"],
                    "gene_set": term,
                    "n_positive": len(x),
                    "n_negative": len(y),
                    "median_positive": np.nanmedian(x) if len(x) else np.nan,
                    "median_negative": np.nanmedian(y) if len(y) else np.nan,
                    "median_difference": (np.nanmedian(x) - np.nanmedian(y)) if len(x) and len(y) else np.nan,
                    "cliffs_delta": cliffs_delta(x, y),
                    "mannwhitney_p": p,
                }
            )
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    out["bh_fdr_within_comparison"] = np.nan
    for _, idx in out.groupby("comparison").groups.items():
        out.loc[idx, "bh_fdr_within_comparison"] = bh_fdr(out.loc[idx, "mannwhitney_p"].tolist())
    return out


def dataset_audit_row(dataset: str, expr: pd.DataFrame, meta: pd.DataFrame, evidence_tier: str, data_type: str) -> dict[str, object]:
    m = meta.drop_duplicates(subset=["sample_id"]).set_index("sample_id")
    matched = [s for s in expr.columns if s in m.index]
    mm = m.loc[matched].copy()
    return {
        "dataset": dataset,
        "data_type": data_type,
        "evidence_tier": evidence_tier,
        "expression_samples": expr.shape[1],
        "matched_metadata_samples": len(matched),
        "features_after_gene_level_processing": expr.shape[0],
        "diagnostic_b_all_reference_n": int(mm.get("include_b_all_reference", pd.Series(False, index=mm.index)).sum()),
        "etv6_runx1_or_like_n": int(mm.get("etv6_runx1_or_like", pd.Series(False, index=mm.index)).sum()),
        "chr16_gain_text_n": int(mm.get("chr16_gain_text", pd.Series(False, index=mm.index)).sum()),
        "ras_mapk_mutation_n": int(mm.get("ras_pathway_mutation", pd.Series(False, index=mm.index)).sum()),
    }


def save_fig(fig: plt.Figure, stem: str) -> None:
    fig.savefig(FIGURES / f"{stem}.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIGURES / f"{stem}.pdf", bbox_inches="tight")
    plt.close(fig)


def plot_prerank_dotplot(prerank: pd.DataFrame) -> None:
    if prerank.empty:
        return
    df = prerank.copy()
    term_col = "Term" if "Term" in df.columns else "gene_set"
    nes_col = "NES"
    fdr_col = "FDR q-val" if "FDR q-val" in df.columns else "fdr"
    df = df[df[term_col].isin(FOCUS_GENESETS)].copy()
    if df.empty:
        return
    df["term_label"] = df[term_col].map(TERM_LABELS).fillna(df[term_col])
    df["fdr"] = pd.to_numeric(df[fdr_col], errors="coerce")
    df["nes"] = pd.to_numeric(df[nes_col], errors="coerce")
    df["minus_log10_fdr"] = -np.log10(df["fdr"].clip(lower=1e-4))
    order_terms = [TERM_LABELS[t] for t in FOCUS_GENESETS if TERM_LABELS[t] in set(df["term_label"])]
    order_comps = list(dict.fromkeys(df["comparison"].tolist()))
    if order_terms:
        df["term_label"] = pd.Categorical(df["term_label"], categories=order_terms, ordered=True)
    if order_comps:
        df["comparison"] = pd.Categorical(df["comparison"], categories=order_comps, ordered=True)
    sns.set_theme(style="whitegrid", context="paper", font_scale=0.9)
    fig, ax = plt.subplots(figsize=(12.5, 5.8))
    cmap = sns.diverging_palette(240, 10, as_cmap=True)
    sns.scatterplot(
        data=df,
        x="comparison",
        y="term_label",
        hue="nes",
        size="minus_log10_fdr",
        sizes=(35, 260),
        palette=cmap,
        hue_norm=(-2.5, 2.5),
        edgecolor="black",
        linewidth=0.35,
        ax=ax,
    )
    ax.axhline(1.5, color="none")
    ax.set_xlabel("")
    ax.set_ylabel("")
    for label in ax.get_xticklabels():
        label.set_rotation(35)
        label.set_ha("right")
    ax.set_title("Preranked gseapy results across public expression cohorts")
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles, labels, bbox_to_anchor=(1.02, 1), loc="upper left", frameon=True, title="NES / -log10(FDR)")
    save_fig(fig, "Figure_1_multicohort_prerank_NES")


def plot_ssgsea_heatmap(ss_tests: pd.DataFrame) -> None:
    if ss_tests.empty:
        return
    df = ss_tests[ss_tests["gene_set"].isin(FOCUS_GENESETS)].copy()
    df["term_label"] = df["gene_set"].map(TERM_LABELS).fillna(df["gene_set"])
    mat = df.pivot(index="term_label", columns="comparison", values="median_difference")
    if mat.empty:
        return
    order = [TERM_LABELS[t] for t in FOCUS_GENESETS if TERM_LABELS[t] in mat.index]
    mat = mat.loc[order]
    sns.set_theme(style="white", context="paper", font_scale=0.9)
    fig, ax = plt.subplots(figsize=(12.8, 5.8))
    sns.heatmap(
        mat,
        cmap=sns.diverging_palette(240, 10, as_cmap=True),
        center=0,
        linewidths=0.5,
        linecolor="white",
        annot=True,
        fmt=".2f",
        cbar_kws={"label": "ssGSEA NES median difference"},
        ax=ax,
    )
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_title("ssGSEA directionality across direct and external expression cohorts")
    for label in ax.get_xticklabels():
        label.set_rotation(35)
        label.set_ha("right")
    save_fig(fig, "Figure_2_multicohort_ssGSEA_heatmap")


def plot_evidence_map(audit: pd.DataFrame, prevalence: pd.DataFrame) -> None:
    rows = []
    for _, r in audit.iterrows():
        rows.append(
            {
                "dataset": r["dataset"],
                "samples": r["matched_metadata_samples"],
                "evidence": r["evidence_tier"],
                "type": r["data_type"],
            }
        )
    if not prevalence.empty:
        prev_plot = prevalence[~prevalence["dataset"].astype(str).str.contains("Combined prevalence", case=False, na=False)]
        for _, r in prev_plot.iterrows():
            rows.append(
                {
                    "dataset": str(r["dataset"]).replace(" diagnostic karyotype", ""),
                    "samples": int(r["total"]),
                    "evidence": "external/public CNV prevalence",
                    "type": "CNV/karyotype",
                }
            )
    df = pd.DataFrame(rows)
    if df.empty:
        return
    tier_order = [
        "same-patient chr16/RAS/expression",
        "external expression subtype reference",
        "external/public CNV prevalence",
    ]
    df["evidence"] = pd.Categorical(df["evidence"], categories=tier_order, ordered=True)
    sns.set_theme(style="whitegrid", context="paper", font_scale=1.0)
    fig, ax = plt.subplots(figsize=(10.5, 4.8))
    sns.barplot(data=df.sort_values("evidence"), x="dataset", y="samples", hue="evidence", dodge=False, ax=ax)
    ax.set_ylabel("Samples contributing to evidence layer")
    ax.set_xlabel("")
    ax.set_title("Expanded evidence map: direct closure versus external context")
    for label in ax.get_xticklabels():
        label.set_rotation(30)
        label.set_ha("right")
    ax.legend(title="")
    save_fig(fig, "Figure_3_expanded_evidence_map")


def load_cnv_prevalence() -> pd.DataFrame:
    path = EXPANDED_VALIDATION / "tables" / "expanded_chr16_gain_prevalence_comparison.csv"
    if not path.exists():
        return pd.DataFrame()
    return pd.read_csv(path)


def markdown_table(df: pd.DataFrame) -> str:
    if df.empty:
        return ""
    out = df.copy()
    out = out.replace([np.inf, -np.inf], np.nan)
    out = out.fillna("")
    cols = [str(c) for c in out.columns]

    def cell(value: object) -> str:
        text = str(value)
        text = text.replace("|", "\\|")
        text = text.replace("\n", "<br>")
        return text

    lines = [
        "| " + " | ".join(cols) + " |",
        "| " + " | ".join(["---"] * len(cols)) + " |",
    ]
    for _, row in out.iterrows():
        lines.append("| " + " | ".join(cell(row[c]) for c in out.columns) + " |")
    return "\n".join(lines)


def write_report(
    audit: pd.DataFrame,
    prerank: pd.DataFrame,
    ss_tests: pd.DataFrame,
    prevalence: pd.DataFrame,
    gene_set_sizes: pd.DataFrame,
) -> None:
    term_col = "Term" if "Term" in prerank.columns else "gene_set"
    fdr_col = "FDR q-val" if "FDR q-val" in prerank.columns else "fdr"
    lines = [
        "# Expanded gseapy chr16/RAS-MAPK/B-ALL closure",
        "",
        "## Scope",
        "",
        "This analysis expands the expression evidence base beyond GSE181157 by adding GSE227832, GSE228632, and GSE87070. GSE181157 remains the only available public cohort in this workspace with same-patient RNA-seq, text-mined chr16 gain, and RAS/MAPK mutation fields. External expression datasets are therefore used as subtype/state context, not as direct chr16 gain or RAS mutation validation.",
        "",
        "## Dataset audit",
        "",
        markdown_table(audit),
        "",
        "## Gene-set audit",
        "",
        markdown_table(gene_set_sizes),
        "",
        "## chr16 gain public prevalence layer",
        "",
    ]
    if prevalence.empty:
        lines.append("No precomputed CNV prevalence table was found.")
    else:
        prev_cols = [c for c in ["dataset", "definition", "chr16_gain", "total", "formatted"] if c in prevalence.columns]
        lines.append(markdown_table(prevalence[prev_cols]))
    lines += ["", "## Preranked gseapy summary", ""]
    if prerank.empty:
        lines.append("No prerank results were generated.")
    else:
        keep = prerank[prerank[term_col].isin(FOCUS_GENESETS)].copy()
        keep["NES"] = pd.to_numeric(keep["NES"], errors="coerce")
        keep[fdr_col] = pd.to_numeric(keep[fdr_col], errors="coerce")
        keep = keep.sort_values(["comparison", fdr_col, "NES"], ascending=[True, True, False])
        cols = ["dataset", "comparison", "positive_label", "negative_label", "n_positive", "n_negative", term_col, "NES", "NOM p-val", fdr_col]
        cols = [c for c in cols if c in keep.columns]
        lines.append(markdown_table(keep[cols]))
    lines += ["", "## ssGSEA group-test summary", ""]
    if ss_tests.empty:
        lines.append("No ssGSEA group tests were generated.")
    else:
        keep = ss_tests[ss_tests["gene_set"].isin(FOCUS_GENESETS)].copy()
        keep = keep.sort_values(["comparison", "bh_fdr_within_comparison", "gene_set"])
        cols = [
            "dataset",
            "comparison",
            "positive_label",
            "negative_label",
            "gene_set",
            "n_positive",
            "n_negative",
            "median_difference",
            "cliffs_delta",
            "mannwhitney_p",
            "bh_fdr_within_comparison",
        ]
        lines.append(markdown_table(keep[cols]))
    lines += [
        "",
        "## Interpretation boundary",
        "",
        "- The direct chr16 gain/RAS-MAPK expression loop is still supported primarily by GSE181157, where chr16-gain B-ALL shows positive prerank enrichment for chr16 dosage and RAS-MAPK gene sets.",
        "- GSE227832, GSE228632, and GSE87070 expand the B-ALL expression-state context. They test whether ETV6::RUNX1/like cases sit in a reproducible B-lineage/RAS-MAPK transcriptional background, but they do not independently validate chr16 gain as a patient-level exposure.",
        "- The CNV prevalence layer (TARGET/Oksa-NOPHO/GSE184692) supports that chr16 gain is a recurrent low-frequency event, but it is not pooled with expression cohorts for survival or mutation interaction modeling.",
        "",
    ]
    (ROOT / "expanded_gseapy_closure_report.md").write_text("\n".join(lines), encoding="utf-8")


def append_analysis_log(audit: pd.DataFrame, prerank: pd.DataFrame, ss_tests: pd.DataFrame) -> None:
    log = PROJECT / "analysis_log.md"
    term_col = "Term" if "Term" in prerank.columns else "gene_set"
    fdr_col = "FDR q-val" if "FDR q-val" in prerank.columns else "fdr"
    direct_line = ""
    if not prerank.empty and term_col in prerank.columns:
        direct = prerank[
            (prerank["comparison"].eq("B_ALL_chr16_gain_vs_no_gain"))
            & (prerank[term_col].isin(["CHR16_DOSAGE_EXPRESSED", "RAS_MAPK_CORE", "MAPK_ERK_FEEDBACK"]))
        ].copy()
        if not direct.empty:
            direct[fdr_col] = pd.to_numeric(direct[fdr_col], errors="coerce")
            direct["NES"] = pd.to_numeric(direct["NES"], errors="coerce")
            direct_line = "; ".join(
                f"{r[term_col]} NES={r['NES']:.2f}, FDR={r[fdr_col]:.3g}" for _, r in direct.iterrows()
            )
    entry = [
        "",
        "## 2026-07-10 expanded gseapy chr16/RAS-MAPK/B-ALL closure",
        "",
        f"- Datasets processed: {', '.join(audit['dataset'].astype(str).tolist())}.",
        f"- Matched metadata samples by dataset: {audit[['dataset', 'matched_metadata_samples']].to_dict(orient='records')}.",
        f"- Direct GSE181157 chr16-gain prerank highlights: {direct_line or 'not available'}.",
        f"- Generated tables in `{TABLES}` and figures in `{FIGURES}`.",
        "- Evidence boundary: external expression cohorts expand biological context; same-patient chr16 gain/RAS validation remains limited to GSE181157, with CNV prevalence supported separately by TARGET/Oksa-NOPHO/GSE184692 summaries.",
        "",
    ]
    with log.open("a", encoding="utf-8", newline="\n") as handle:
        handle.write("\n".join(entry))


def process_rnaseq_dataset(dataset: str, counts: pd.DataFrame, meta: pd.DataFrame, evidence_tier: str) -> dict[str, object]:
    common = [c for c in counts.columns if c in set(meta["sample_id"].astype(str))]
    counts = counts[common]
    expr, lib_sizes = log2_cpm(counts)
    gene_sets = build_rna_gene_sets(expr)
    expr = filter_rna_expr(expr, gene_sets)
    write_gmt(gene_sets, GENESETS / f"{dataset}_custom_gene_sets.gmt")
    scores = run_ssgsea(expr, gene_sets)
    scores.to_csv(TABLES / f"{dataset}_ssgsea_scores.csv", encoding="utf-8-sig")
    ss_tests = compare_ssgsea_for_dataset(dataset, scores, meta)
    prerank, ranks = run_prerank_for_dataset(dataset, expr, meta, gene_sets)
    meta.to_csv(TABLES / f"{dataset}_harmonized_metadata.csv", index=False, encoding="utf-8-sig")
    pd.DataFrame({"sample_id": lib_sizes.index, "library_size": lib_sizes.values}).to_csv(
        TABLES / f"{dataset}_library_sizes.csv", index=False, encoding="utf-8-sig"
    )
    gs_sizes = pd.DataFrame(
        [{"dataset": dataset, "gene_set": k, "n_genes": len(v), "id_space": "Ensembl"} for k, v in gene_sets.items()]
    )
    return {
        "audit": dataset_audit_row(dataset, expr, meta, evidence_tier, "RNA-seq counts/log2CPM"),
        "prerank": prerank,
        "ranks": ranks,
        "ss_tests": ss_tests,
        "gene_set_sizes": gs_sizes,
    }


def process_microarray_dataset(dataset: str, expr: pd.DataFrame, meta: pd.DataFrame, evidence_tier: str) -> dict[str, object]:
    common = [c for c in expr.columns if c in set(meta["sample_id"].astype(str))]
    expr = expr[common]
    gene_sets = build_symbol_gene_sets(expr)
    expr = filter_expr(expr, gene_sets)
    write_gmt(gene_sets, GENESETS / f"{dataset}_custom_gene_sets.gmt")
    scores = run_ssgsea(expr, gene_sets)
    scores.to_csv(TABLES / f"{dataset}_ssgsea_scores.csv", encoding="utf-8-sig")
    ss_tests = compare_ssgsea_for_dataset(dataset, scores, meta)
    prerank, ranks = run_prerank_for_dataset(dataset, expr, meta, gene_sets)
    meta.to_csv(TABLES / f"{dataset}_harmonized_metadata.csv", index=False, encoding="utf-8-sig")
    gs_sizes = pd.DataFrame(
        [{"dataset": dataset, "gene_set": k, "n_genes": len(v), "id_space": "HGNC symbol"} for k, v in gene_sets.items()]
    )
    return {
        "audit": dataset_audit_row(dataset, expr, meta, evidence_tier, "GPL570 microarray/symbol median"),
        "prerank": prerank,
        "ranks": ranks,
        "ss_tests": ss_tests,
        "gene_set_sizes": gs_sizes,
    }


def main() -> None:
    outputs: list[dict[str, object]] = []

    gse181_meta = read_gse181_metadata()
    gse181_counts = read_gse181_counts()
    outputs.append(
        process_rnaseq_dataset(
            "GSE181157",
            gse181_counts,
            gse181_meta,
            "same-patient chr16/RAS/expression",
        )
    )

    gse227832_meta = standardize_gse227832_metadata()
    gse227832_counts = read_count_matrix(GSE227832_COUNTS)
    outputs.append(
        process_rnaseq_dataset(
            "GSE227832",
            gse227832_counts,
            gse227832_meta,
            "external expression subtype reference",
        )
    )

    gse228632_meta = standardize_gse228632_metadata()
    gse228632_counts = read_count_matrix(GSE228632_COUNTS)
    outputs.append(
        process_rnaseq_dataset(
            "GSE228632",
            gse228632_counts,
            gse228632_meta,
            "external expression subtype reference",
        )
    )

    gse87070_meta = standardize_gse87070_metadata()
    gse87070_expr = read_gse87070_symbol_expression()
    outputs.append(
        process_microarray_dataset(
            "GSE87070",
            gse87070_expr,
            gse87070_meta,
            "external expression subtype reference",
        )
    )

    audit = pd.DataFrame([o["audit"] for o in outputs])
    prerank = pd.concat([o["prerank"] for o in outputs if not o["prerank"].empty], ignore_index=True)
    ss_tests = pd.concat([o["ss_tests"] for o in outputs if not o["ss_tests"].empty], ignore_index=True)
    ranks = pd.concat([o["ranks"] for o in outputs if not o["ranks"].empty], ignore_index=True)
    gene_set_sizes = pd.concat([o["gene_set_sizes"] for o in outputs], ignore_index=True)
    prevalence = load_cnv_prevalence()

    audit.to_csv(TABLES / "dataset_audit.csv", index=False, encoding="utf-8-sig")
    prerank.to_csv(TABLES / "gseapy_prerank_summary_all_datasets.csv", index=False, encoding="utf-8-sig")
    ss_tests.to_csv(TABLES / "ssgsea_group_tests_all_datasets.csv", index=False, encoding="utf-8-sig")
    ranks.to_csv(TABLES / "top500_ranked_genes_all_comparisons.csv", index=False, encoding="utf-8-sig")
    gene_set_sizes.to_csv(TABLES / "gene_set_audit_all_datasets.csv", index=False, encoding="utf-8-sig")
    if not prevalence.empty:
        prevalence.to_csv(TABLES / "chr16_gain_public_prevalence_layer.csv", index=False, encoding="utf-8-sig")

    harmonized_meta = pd.concat(
        [
            gse181_meta,
            gse227832_meta,
            gse228632_meta,
            gse87070_meta,
        ],
        ignore_index=True,
        sort=False,
    )
    keep_cols = [
        "dataset",
        "sample_id",
        "geo_accession",
        "title",
        "disease",
        "disease_state",
        "predicted_subtype",
        "include_b_all_reference",
        "is_b_all",
        "is_t_all",
        "is_diagnosis",
        "etv6_runx1_or_like",
        "chr16_gain_text",
        "ras_pathway_mutation",
    ]
    harmonized_meta[[c for c in keep_cols if c in harmonized_meta.columns]].to_csv(
        TABLES / "harmonized_sample_metadata_all_expression_datasets.csv",
        index=False,
        encoding="utf-8-sig",
    )

    plot_prerank_dotplot(prerank)
    plot_ssgsea_heatmap(ss_tests)
    plot_evidence_map(audit, prevalence)
    write_report(audit, prerank, ss_tests, prevalence, gene_set_sizes)
    append_analysis_log(audit, prerank, ss_tests)

    manifest = {
        "analysis_date": "2026-07-10",
        "python_gseapy": gp.__version__,
        "datasets": audit.to_dict(orient="records"),
        "outputs": {
            "report": str(ROOT / "expanded_gseapy_closure_report.md"),
            "tables": str(TABLES),
            "figures": str(FIGURES),
        },
        "evidence_boundary": (
            "Only GSE181157 has same-patient RNA-seq plus chr16/RAS annotations in this workspace; "
            "GSE227832, GSE228632, and GSE87070 are external expression-state references."
        ),
    }
    (ROOT / "analysis_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
