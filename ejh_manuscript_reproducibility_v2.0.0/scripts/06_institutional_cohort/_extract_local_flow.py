# -*- coding: utf-8 -*-
"""
Module B - local immunophenotype extraction.
READ ONLY on the source workbook. Only study_id (LOCAL-xxx) is written out;
no name / record number / ID / date of birth / exact dates are exported.

Step 1: verify row-order mapping between the workbook '2020' sheet and study_id
        by cross-checking three independent fields (age, WBC, D19 MRD).
Step 2: extract CD34 / CD38 as reported.
"""
import openpyxl, pandas as pd, sys, io, re, os
from pathlib import Path
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

SRC = Path(os.environ["EJH_INSTITUTIONAL_WORKBOOK"])
COH = Path(os.environ["EJH_DEIDENTIFIED_COHORT"])

coh = pd.read_csv(COH, sep="\t", keep_default_na=False, dtype=str)
sub = coh[coh.source_sheet == "2020"].reset_index(drop=True)
print("cohort 2020 n =", len(sub))

wb = openpyxl.load_workbook(SRC, read_only=True, data_only=True)
ws = wb["2020"]
rows = list(ws.iter_rows(max_row=33, values_only=True))
hdr = {i: str(c).strip() for i, c in enumerate(rows[0]) if c is not None}
inv = {v: k for k, v in hdr.items()}

def col(name):
    return inv.get(name)

print("col idx: age", col("诊断年龄"), "WBC", col("WBC（×109）"),
      "D19MRD", col("D19MRD(%)"), "CD34", col("CD34"), "CD38", col("CD38"))

def num(v):
    if v is None:
        return None
    s = str(v).strip()
    if s in ("", "/", "-", "—", "无", "NA"):
        return None
    try:
        return float(s)
    except Exception:
        return None

# --- step 1: verify row order on three independent fields ---
print("\n=== ROW-ORDER VERIFICATION (workbook r_k  vs  LOCAL-xxx #k) ===")
ok = 0
mismatch = []
for k in range(1, 33):
    r = rows[k]
    rec = sub.iloc[k - 1]
    a_w, w_w, m_w = num(r[col("诊断年龄")]), num(r[col("WBC（×109）")]), num(r[col("D19MRD(%)")])
    a_c = float(rec["age_at_diagnosis_years_num"]) if rec["age_at_diagnosis_years_num"] else None
    w_c = float(rec["WBC_num"]) if rec["WBC_num"] else None
    m_c = float(rec["D19_MRD_continuous"]) if rec["D19_MRD_continuous"] else None
    same = (a_w == a_c or (a_w is None and a_c is None)) and \
           (w_w == w_c or (w_w is None and w_c is None)) and \
           (m_w == m_c or (m_w is None and m_c is None))
    ok += same
    if not same:
        mismatch.append((rec["study_id"], a_w, a_c, w_w, w_c, m_w, m_c))
print("rows matching on all three fields: %d / 32" % ok)
for m in mismatch:
    print("  MISMATCH", m)

# --- step 2: extract ---
print("\n=== CD34 / CD38 EXTRACTION ===")
out = []
for k in range(1, 33):
    r = rows[k]
    rec = sub.iloc[k - 1]
    c34_raw = r[col("CD34")]
    c38_raw = r[col("CD38")]
    out.append({
        "study_id": rec["study_id"],
        "treatment_era": rec["treatment_era"],
        "CD34_raw": "" if c34_raw is None else str(c34_raw).strip(),
        "CD38_raw": "" if c38_raw is None else str(c38_raw).strip(),
        "D19_MRD_percent": rec["D19_MRD_continuous"],
        "D19_MRD_ge_0p1": rec["D19_MRD_ge_0p1"],
        "age_years": rec["age_at_diagnosis_years_num"],
        "WBC": rec["WBC_num"],
        "sex": rec["sex"],
    })
df = pd.DataFrame(out)
n34 = (df["CD34_raw"].isin(["", "/", "-", "无"])).sum()
n38 = (df["CD38_raw"].isin(["", "/", "-", "无"])).sum()
print("CD34 missing:", n34, " available:", 32 - n34)
print("CD38 missing:", n38, " available:", 32 - n38)
print(df.to_string())
df.to_csv("LOCAL_FLOW_2020_RAW.tsv", sep="\t", index=False)
print("\nwritten LOCAL_FLOW_2020_RAW.tsv")
wb.close()
