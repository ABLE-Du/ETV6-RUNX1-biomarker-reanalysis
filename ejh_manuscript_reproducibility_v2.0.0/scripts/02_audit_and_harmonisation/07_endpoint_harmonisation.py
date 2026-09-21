#!/usr/bin/env python -u
# -*- coding: utf-8 -*-
"""07_endpoint_harmonisation.py -> tables/07_ENDPOINT_HARMONISATION.tsv

Round-1 deliverable 7: the endpoint table that makes MRD pooling impossible.

Every cohort measures "poor early response" at a different timepoint and threshold.
The project rule is that these are never merged, but a rule in a document is easy to
violate by accident.  So this script turns the rule into data plus an assertion:
each row is one (cohort, timepoint, threshold) triple, and the script refuses to
finish if two rows collide on the same (timepoint, threshold) pair, because that
would be the shape that invites a silent merge.

It also records the Li day-numbering discrepancy discovered in `04_li_audit.py`.
"""
from __future__ import annotations

import csv
import os
import json
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OKSA = ROOT / "intermediate" / "oksa"
LI = ROOT / "intermediate" / "li"
TABLES = ROOT / "tables"
REPORTS = ROOT / "reports"
UPSTREAM = Path(os.environ.get("EJH_UPSTREAM_ROOT", Path(__file__).resolve().parents[2] / "upstream"))

COLUMNS = ["cohort", "endpoint_id", "timepoint_day", "threshold", "poor_definition",
           "good_definition", "intermediate_definition", "n_total", "n_evaluable",
           "n_poor", "comparability_class", "may_be_pooled_with", "source",
           "caveat"]


def say(*a):
    print(*a, flush=True)


def main() -> int:
    TABLES.mkdir(parents=True, exist_ok=True)
    s1 = pd.read_csv(OKSA / "SuppTable1_response_groups.tsv", sep="\t",
                     dtype=str, keep_default_na=False).drop_duplicates("Patient ID")
    mrd = pd.read_csv(LI / "li_mrd.tsv", sep="\t", dtype=str, keep_default_na=False)
    local = pd.read_csv(UPSTREAM / "local_cohort" / "02_LOCAL_DEIDENTIFIED_COHORT.tsv",
                        sep="\t", dtype=str, keep_default_na=False)
    pball = pd.read_csv(UPSTREAM / "analysis" / "revision_01" / "cohorts" / "P_BALL.tsv",
                        sep="\t", dtype=str, keep_default_na=False)

    # ---- Oksa ---------------------------------------------------------------
    o_mid = s1["Mid-induction response"].str.strip().value_counts().to_dict()
    o_eoi = s1["EOI response"].str.strip().value_counts().to_dict()
    o_eoc = s1["EOC MRD status"].str.strip().value_counts().to_dict()

    # ---- Li -----------------------------------------------------------------
    li_d19 = int((~mrd["D19MRD"].str.strip().isin(["", "NA", "nan"])).sum())
    li_d46 = int((~mrd["D46MRD"].str.strip().isin(["", "NA", "nan"])).sum())

    # ---- local --------------------------------------------------------------
    ev = local["D19_MRD_ge_0p1"].astype(str).str.upper().isin(["TRUE", "FALSE"])
    l_poor = int(local.loc[ev, "D19_MRD_ge_0p1"].astype(str).str.upper().eq("TRUE").sum())

    # ---- TARGET -------------------------------------------------------------
    t_poor = int(pball["MRD_group"].astype(str).str.strip().eq("Poor").sum())

    rows = [
        dict(cohort="OKSA", endpoint_id="OKSA_d15", timepoint_day=15,
             threshold="10%",
             poor_definition="slow: MRD >= 10%",
             good_definition="fast: MRD < 10%", intermediate_definition="not defined",
             n_total=len(s1), n_evaluable=int(sum(o_mid.get(k, 0) for k in ("fast", "slow"))),
             n_poor=int(o_mid.get("slow", 0)),
             comparability_class="non_poolable",
             may_be_pooled_with="",
             source="Oksa Supp Table 1 (Mid-induction response)",
             caveat="day-15 label is 'unknown' for all 8 publicly linkable expression patients"),
        dict(cohort="OKSA", endpoint_id="OKSA_d29_EOI", timepoint_day=29,
             threshold="0.1%",
             poor_definition="slow: MRD >= 0.1%",
             good_definition="fast: MRD = 0 (undetectable/negative)",
             intermediate_definition="intermediate: 0 < MRD < 0.1%",
             n_total=len(s1), n_evaluable=int(sum(o_eoi.get(k, 0)
                                                  for k in ("fast", "intermediate", "slow"))),
             n_poor=int(o_eoi.get("slow", 0)),
             comparability_class="non_poolable", may_be_pooled_with="",
             source="Oksa Supp Table 1 (EOI response)",
             caveat="three-level ordinal; the primary contrast is slow vs fast"),
        dict(cohort="OKSA", endpoint_id="OKSA_d79_EOC", timepoint_day=79,
             threshold="detectable",
             poor_definition="slow: MRD > 0%", good_definition="fast: MRD = 0%",
             intermediate_definition="not defined",
             n_total=len(s1),
             n_evaluable=int(sum(o_eoc.get(k, 0) for k in ("positive", "negative"))),
             n_poor=int(o_eoc.get("positive", 0)),
             comparability_class="non_poolable", may_be_pooled_with="",
             source="Oksa Supp Table 1 (EOC MRD status)",
             caveat="different timepoint and qualitative threshold"),
        dict(cohort="LI", endpoint_id="LI_early_induction", timepoint_day="19_or_15_UNRESOLVED",
             threshold="0.01% (negative <0.01%)",
             poor_definition="MRD positive (>= 0.01%)",
             good_definition="MRD negative (< 0.01%)", intermediate_definition="not defined",
             n_total=194, n_evaluable=li_d19, n_poor="",
             comparability_class="non_poolable", may_be_pooled_with="",
             source="Li Supp workbook, column D19MRD",
             caveat="NAME COLLISION RISK: the supplementary column is D19MRD while the main text says day 15. "
                    "Never equate with the institutional day-19 endpoint."),
        dict(cohort="LI", endpoint_id="LI_EOI_day46", timepoint_day="46_or_42_UNRESOLVED",
             threshold="0.01% (negative <0.01%)",
             poor_definition="MRD positive (>= 0.01%)",
             good_definition="MRD negative (< 0.01%)", intermediate_definition="not defined",
             n_total=194, n_evaluable=li_d46, n_poor="",
             comparability_class="non_poolable", may_be_pooled_with="",
             source="Li Supp workbook, column D46MRD",
             caveat="column name says day 46; main text says end of induction day 42"),
        dict(cohort="LOCAL", endpoint_id="LOCAL_d19", timepoint_day=19,
             threshold="0.1%",
             poor_definition="day-19 MRD >= 0.1%",
             good_definition="day-19 MRD < 0.1%", intermediate_definition="not defined",
             n_total=len(local), n_evaluable=int(ev.sum()), n_poor=l_poor,
             comparability_class="non_poolable", may_be_pooled_with="",
             source="local_cohort/02_LOCAL_DEIDENTIFIED_COHORT.tsv",
             caveat="<0.1% is a reporting threshold, not documented MRD negativity"),
        dict(cohort="TARGET", endpoint_id="TARGET_d29", timepoint_day=29,
             threshold="0.01%",
             poor_definition="day-29 MRD >= 0.01%",
             good_definition="day-29 MRD < 0.01%", intermediate_definition="not defined",
             n_total=len(pball), n_evaluable=len(pball), n_poor=t_poor,
             comparability_class="non_poolable", may_be_pooled_with="",
             source="analysis/revision_01/cohorts/P_BALL.tsv",
             caveat="general B-ALL, not ETV6-specific; supportive role only"),
    ]

    # ---- the assertion that makes silent pooling impossible -----------------
    seen: dict[tuple, str] = {}
    collisions = []
    for r in rows:
        key = (r["timepoint_day"], r["threshold"])
        if key in seen:
            collisions.append((seen[key], r["endpoint_id"], key))
        seen[key] = r["endpoint_id"]

    for r in rows:
        assert set(r.keys()) == set(COLUMNS), f"column mismatch in {r['endpoint_id']}"

    out = TABLES / "07_ENDPOINT_HARMONISATION.tsv"
    with out.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLUMNS, delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    say(f"[write] {out.relative_to(ROOT)}  ({len(rows)} endpoints x {len(COLUMNS)} columns)")

    nunique = len({(r["timepoint_day"], r["threshold"]) for r in rows})
    say(f"[assert] distinct (timepoint, threshold) pairs = {nunique} of {len(rows)}")
    if collisions:
        for a, b, k in collisions:
            say(f"  [COLLISION] {a} and {b} share {k}")

    md = f"""# 07 -- Endpoint harmonisation

Machine-readable form: `tables/07_ENDPOINT_HARMONISATION.tsv`.
Every row is one (cohort, timepoint, threshold) triple. **No two rows may be aggregated.**

| cohort | endpoint | day | threshold | poor definition | n total | n evaluable | n poor |
|---|---|---|---|---|---|---|---|
"""
    for r in rows:
        md += (f"| {r['cohort']} | `{r['endpoint_id']}` | {r['timepoint_day']} | "
               f"{r['threshold']} | {r['poor_definition']} | {r['n_total']} | "
               f"{r['n_evaluable']} | {r['n_poor']} |\n")

    md += f"""
Distinct (timepoint, threshold) pairs: **{nunique} of {len(rows)}**
{"(collisions: " + str(collisions) + ")" if collisions else "(no collisions: every endpoint is uniquely keyed, so a merge cannot happen by accident)"}

## Notes that each row carries

| endpoint | caveat |
|---|---|
"""
    for r in rows:
        md += f"| `{r['endpoint_id']}` | {r['caveat']} |\n"

    md += """
## The one collision risk worth naming

The institutional cohort's endpoint is **day 19 at 0.1 %**. The Li supplementary
workbook's early-MRD column is literally named **`D19MRD`**. Two different things with the
same day number in their names is the most likely single error in this project, and the
`caveat` column on `LI_early_induction` exists specifically to stop it. Li's column is
recorded as `19_or_15_UNRESOLVED` because the paper's own text calls the same measurement
"day 15" while the column header says D19 — the audit refuses to resolve that, and any
downstream table must carry the note forward.

## What the assertion does not cover

The assertion only catches two endpoints that share a (day, threshold) key. It cannot catch
a *wrong* aggregation that renames one endpoint first. That is why each Round-2 script that
consumes this table must also carry the endpoint id into its output rows, so a merged row
would be visibly missing a valid endpoint id.
"""

    (REPORTS / "07_ENDPOINT_HARMONISATION.md").write_text(md, encoding="utf-8")
    say(f"[write] reports/07_ENDPOINT_HARMONISATION.md")

    (TABLES / "_ENDPOINT_HARMONISATION.json").write_text(json.dumps({
        "n_endpoints": len(rows), "distinct_keys": nunique, "collisions": collisions,
        "endpoints": {r["endpoint_id"]: {"day": r["timepoint_day"],
                                         "threshold": r["threshold"],
                                         "n_evaluable": r["n_evaluable"],
                                         "n_poor": r["n_poor"]} for r in rows},
    }, indent=2, ensure_ascii=False), encoding="utf-8")

    return 0 if not collisions else 1


if __name__ == "__main__":
    sys.exit(main())
