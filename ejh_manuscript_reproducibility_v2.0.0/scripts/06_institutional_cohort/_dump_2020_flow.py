# -*- coding: utf-8 -*-
"""Dump per-antigen flow values from the 2020 sheet (only sheet with antigen columns).
READ ONLY. Only study_id (LOCAL-xxx) is written; no identifiers."""
import openpyxl, sys, io, json, os
from pathlib import Path
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

SRC = Path(os.environ["EJH_INSTITUTIONAL_WORKBOOK"])

wb = openpyxl.load_workbook(SRC, read_only=True, data_only=True)
ws = wb["2020"]
rows = list(ws.iter_rows(max_row=33, values_only=True))
hdr = rows[0]
IDX = {i: str(c) for i, c in enumerate(hdr) if c is not None}
cols = [34, 37, 40, 41, 42, 43, 45, 46, 47, 58, 63, 65, 66, 68, 70, 72, 75, 13]
print("HEADER:")
for i in cols:
    print("  col%-4d %s" % (i, IDX.get(i)))
print()
print("ROWS (patient rows start at row index 1):")
for r in range(1, len(rows)):
    row = rows[r]
    if all(v is None for v in row):
        continue
    vals = {IDX.get(i): row[i] for i in cols if i < len(row)}
    print(" r%-3d" % r, {k: (str(v)[:26] if v is not None else "") for k, v in vals.items()})
wb.close()
