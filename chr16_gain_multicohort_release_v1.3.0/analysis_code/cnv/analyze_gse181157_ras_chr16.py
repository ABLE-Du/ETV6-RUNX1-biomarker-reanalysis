from __future__ import annotations

import json
import re
from pathlib import Path

import pandas as pd
from scipy.stats import fisher_exact


ROOT = Path(__file__).resolve().parent
XLSX = ROOT / "geo_supp" / "GSE181157" / "GSE181157_SampleMetadata.xlsx"

RAS_PATTERN = re.compile(r"\b(KRAS|NRAS|PTPN11|NF1|FLT3|BRAF|CBL|MAP2K1|MAP2K2|HRAS)\b", re.I)


def has_chr16_gain(text: object) -> bool:
    if pd.isna(text):
        return False
    s = str(text)
    # Avoid matching unrelated chromosome positions such as t(5;16) without a plus sign.
    return bool(re.search(r"(?<!\d)\+16(?!\d)", s))


def has_ras_pathway_mutation(text: object) -> bool:
    if pd.isna(text):
        return False
    s = str(text).strip()
    if not s or s.lower() in {"no mutation", "none", "nan"}:
        return False
    return bool(RAS_PATTERN.search(s))


def main() -> None:
    (ROOT / "tables").mkdir(exist_ok=True)

    raw = pd.read_excel(XLSX, sheet_name="S1 Patient characteristics", header=1)
    raw = raw[raw["ID"].apply(lambda x: str(x).isdigit())].copy()
    raw["ID"] = raw["ID"].astype(int)
    raw.columns = [str(c).strip() for c in raw.columns]

    raw["rna_predicted_etv6_runx1_like"] = raw["Predicted subtype by RNA-seq"].eq("ETV6-RUNX1_or_like")
    raw["conventional_etv6_runx1"] = raw["Conventional methods subtype"].eq("ETV6-RUNX1")
    raw["fusion_etv6_runx1"] = raw["Fusions"].astype(str).str.contains("ETV6-RUNX1", case=False, regex=False)
    raw["ras_pathway_mutation"] = raw["Mutations"].apply(has_ras_pathway_mutation)
    raw["karyotype_chr16_gain"] = raw["Karyotype"].apply(has_chr16_gain)
    raw["fish_chr16_gain"] = raw["Conventional cytogenetics  (FISH)"].apply(has_chr16_gain)
    raw["any_chr16_gain_text"] = raw["karyotype_chr16_gain"] | raw["fish_chr16_gain"]

    cols = [
        "ID",
        "DFCI ID",
        "Age at Dx (years)",
        "Highest WBC at Dx prior to systemic treatment",
        "NCI Risk",
        "Initial Risk",
        "Final Risk",
        "TP1 MRD",
        "TP2 MRD",
        "Conventional cytogenetics  (FISH)",
        "Karyotype",
        "Conventional methods subtype",
        "Fusions",
        "Mutations",
        "Predicted subtype by RNA-seq",
        "rna_predicted_etv6_runx1_like",
        "fusion_etv6_runx1",
        "ras_pathway_mutation",
        "any_chr16_gain_text",
    ]
    raw[cols].to_csv(ROOT / "tables" / "GSE181157_patient_level_flags.csv", index=False, encoding="utf-8-sig")

    subsets = {
        "all_173": raw,
        "rna_predicted_ETV6_RUNX1_or_like": raw[raw["rna_predicted_etv6_runx1_like"]],
        "fusion_ETV6_RUNX1": raw[raw["fusion_etv6_runx1"]],
        "conventional_ETV6_RUNX1": raw[raw["conventional_etv6_runx1"]],
    }
    rows = []
    for label, df in subsets.items():
        rows.append(
            {
                "subset": label,
                "n": len(df),
                "ras_pathway_mutation_n": int(df["ras_pathway_mutation"].sum()),
                "chr16_gain_text_n": int(df["any_chr16_gain_text"].sum()),
                "chr16_gain_and_ras_n": int((df["any_chr16_gain_text"] & df["ras_pathway_mutation"]).sum()),
                "chr16_gain_without_ras_n": int((df["any_chr16_gain_text"] & ~df["ras_pathway_mutation"]).sum()),
                "ras_without_chr16_gain_n": int((~df["any_chr16_gain_text"] & df["ras_pathway_mutation"]).sum()),
            }
        )
    summary = pd.DataFrame(rows)
    summary.to_csv(ROOT / "tables" / "GSE181157_ras_chr16_summary.csv", index=False, encoding="utf-8-sig")

    etv6 = raw[raw["rna_predicted_etv6_runx1_like"]].copy()
    a = int((etv6["any_chr16_gain_text"] & etv6["ras_pathway_mutation"]).sum())
    b = int((etv6["any_chr16_gain_text"] & ~etv6["ras_pathway_mutation"]).sum())
    c = int((~etv6["any_chr16_gain_text"] & etv6["ras_pathway_mutation"]).sum())
    d = int((~etv6["any_chr16_gain_text"] & ~etv6["ras_pathway_mutation"]).sum())
    odds, p_two = fisher_exact([[a, b], [c, d]], alternative="two-sided")
    _, p_less = fisher_exact([[a, b], [c, d]], alternative="less")
    contingency = pd.DataFrame(
        [
            {
                "subset": "rna_predicted_ETV6_RUNX1_or_like",
                "table_layout": "[[chr16_gain_RASmut,chr16_gain_RASwt],[noChr16gain_RASmut,noChr16gain_RASwt]]",
                "chr16_gain_ras_mut": a,
                "chr16_gain_ras_wt": b,
                "no_chr16_gain_ras_mut": c,
                "no_chr16_gain_ras_wt": d,
                "odds_ratio": odds,
                "fisher_two_sided_p": p_two,
                "fisher_less_p_mutual_exclusivity": p_less,
            }
        ]
    )
    contingency.to_csv(ROOT / "tables" / "GSE181157_etv6_chr16_ras_contingency.csv", index=False, encoding="utf-8-sig")

    examples = etv6[etv6["any_chr16_gain_text"] | etv6["ras_pathway_mutation"]][cols]
    examples.to_csv(ROOT / "tables" / "GSE181157_ETV6_like_chr16_or_ras_examples.csv", index=False, encoding="utf-8-sig")

    metrics = {
        "total_patients": int(len(raw)),
        "rna_predicted_etv6_runx1_or_like": int(len(etv6)),
        "etv6_like_ras_mutated": int(etv6["ras_pathway_mutation"].sum()),
        "etv6_like_chr16_gain_text": int(etv6["any_chr16_gain_text"].sum()),
        "etv6_like_chr16_gain_and_ras": int((etv6["any_chr16_gain_text"] & etv6["ras_pathway_mutation"]).sum()),
        "note": "chr16 gain is text-mined from conventional karyotype/FISH fields, not array-level CNV calls.",
    }
    (ROOT / "GSE181157_ras_chr16_metrics.json").write_text(json.dumps(metrics, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
