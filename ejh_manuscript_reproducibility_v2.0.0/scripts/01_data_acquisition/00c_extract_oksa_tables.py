#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""00c_extract_oksa_tables.py -- parse the Oksa 2025 supplementary workbook once.

The workbook (data/oksa_supp/Oksa2025_supplementary_tables.xlsx, 7.3 MB, 25 sheets)
is expensive to parse, so it is converted to TSV here under intermediate/oksa/ and
every downstream audit script reads the TSVs instead.

Header detection is deliberately explicit rather than assumed: the sheets do NOT
share a header row (Table 1 has a 5-line preamble, Table 20 has 3, Table 22 has
4).  For each sheet the header is taken to be the first row within the first 8
rows whose non-empty cell count is maximal and strictly greater than the row
above it.

Three small sheets defeat that heuristic -- they carry a merged two-line header
or an alignment row above the real header (Table 9, 12, 15).  For those an
explicit override is used, and the fact that an override was needed is recorded
in the QC file, so no sheet is silently mis-parsed.

Every decision is written to intermediate/oksa/_SHEET_QC.tsv together with the
raw first 8 rows, so the detection is auditable by eye.

Outputs
  intermediate/oksa/_SHEET_QC.tsv          header-row decision + raw preamble
  intermediate/oksa/<table>.tsv            one per sheet (tab separated)
"""
from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import openpyxl

ROOT = Path(__file__).resolve().parents[1]
XLSX = ROOT / "data" / "oksa_supp" / "Oksa2025_supplementary_tables.xlsx"
OUT = ROOT / "intermediate" / "oksa"
QC = OUT / "_SHEET_QC.tsv"
MAX_PREAMBLE_SCAN = 8

# Manual overrides.  Each entry states the 1-based header row; optional
# `columns` replaces the parsed header entirely (used when the source has a
# merged two-line header that cannot be reconstructed automatically), optional
# `rename` fixes individual cells, optional `drop_first_col_values` removes
# alignment/label rows that would otherwise be read as data.
HEADER_OVERRIDE: dict[str, dict] = {
    "SuppTable9_Musica_mutsignature": {
        "row": 5,                                  # "Case ID | GE0310 | ..."
        "rename": {0: "Signature"},
        "drop_first_col_values": {"Signature"},    # r6 is a merged column label
    },
    "SuppTable12_CellCycle_stats": {
        "row": 6,
        "columns": ["gene_set_score", "mid_induction_p", "mid_induction_cor",
                    "eoi_p", "eoi_cor"],
    },
    "SuppTable15_IGK_stats": {
        "row": 5,                                  # "Responder groups tested | p-value"
    },
    # ---- multi-block wide sheets -------------------------------------------
    # These carry a second header row above the field names, holding the label of
    # each contrast / block.  Without composing the two rows the contrast labels
    # are lost and the `.1 / .2 / .3` suffixes pandas invents are meaningless --
    # which would silently make rung B1 read the wrong contrast.
    "SuppTable22_DE_analyses": {
        "compose": {"block_row": 4, "field_row": 5, "n_cols": 22},
    },
    "SuppTable7_panel_CNV": {
        "compose": {"block_row": 4, "field_row": 5, "n_cols": 12},
    },
    # ---- sheets whose layout cannot be reconstructed safely ---------------
    # A banner row plus two unlabelled side-by-side blocks: the block boundary is
    # not recoverable from the sheet, so these are extracted but flagged unusable
    # rather than guessed at.  Neither is needed by any planned analysis.
    "SuppTable11_Mut_signature_stats": {"row": 5, "manual_review": True},
    "SuppTable16_CNV_CRISPR_summary": {"row": 5, "manual_review": True},
}


def say(*a):
    print(*a, flush=True)


def cell(v) -> str:
    if v is None:
        return ""
    s = str(v)
    return s.replace("\t", " ").replace("\n", " ").strip()


def detect_header(rows: list[tuple]) -> int:
    """Return the 0-based index of the header row."""
    counts = [sum(1 for v in r if v is not None and str(v).strip() != "")
              for r in rows[:MAX_PREAMBLE_SCAN]]
    if not counts:
        return 0
    best = max(counts)
    for i, c in enumerate(counts):
        if c == best and c >= 3:
            prev = counts[i - 1] if i > 0 else -1
            if c > prev:
                return i
    return counts.index(best) if best >= 3 else 0


def main() -> int:
    if not XLSX.exists():
        say(f"[FAIL] missing {XLSX}")
        return 1
    OUT.mkdir(parents=True, exist_ok=True)
    wb = openpyxl.load_workbook(XLSX, read_only=True, data_only=True)

    qc_rows = []
    written = []
    for name in wb.sheetnames:
        ws = wb[name]
        it = ws.iter_rows(values_only=True)
        head_raw = []
        for i, r in enumerate(it):
            head_raw.append(r)
            if i >= MAX_PREAMBLE_SCAN - 1:
                break

        ov = HEADER_OVERRIDE.get(name)
        manual_review = bool((ov or {}).get("manual_review"))

        if ov and "compose" in ov:
            # two-row header: a block label row above a field-name row
            cp = ov["compose"]
            br, fr, nc = cp["block_row"], cp["field_row"], cp["n_cols"]
            ws2 = wb[name]
            it2 = ws2.iter_rows(values_only=True)
            blk_raw = fld_raw = None
            for i, r in enumerate(it2):
                if i == br - 1:
                    blk_raw = [cell(v) for v in r[:nc]]
                if i == fr - 1:
                    fld_raw = [cell(v) for v in r[:nc]]
                    break
            blk_raw = (blk_raw or []) + [""] * nc
            fld_raw = (fld_raw or []) + [""] * nc
            header = []
            cur = ""
            for j in range(nc):
                if blk_raw[j]:
                    cur = blk_raw[j]
                f = fld_raw[j]
                header.append(f"{cur} | {f}" if (cur and f) else (f or cur))
            while header and header[-1] == "":
                header.pop()
            ncol = len(header)
            method = "composed_two_row"
            h = fr - 1
            # continue reading data from the field row onwards
            ws3 = wb[name]
            it3 = ws3.iter_rows(values_only=True)
            for _ in range(fr):
                next(it3, None)
            data_iter = it3
        else:
            if ov:
                h = ov["row"] - 1
                method = "manual_override"
            else:
                h = detect_header(head_raw)
                method = "auto"

            ws2 = wb[name]
            it2 = ws2.iter_rows(values_only=True)
            for _ in range(h):
                next(it2, None)
            raw_header = [cell(v) for v in next(it2, ())]
            while raw_header and raw_header[-1] == "":
                raw_header.pop()
            ncol = len(raw_header)

            if ov and "columns" in ov:
                header = list(ov["columns"])
            else:
                header = list(raw_header)
            if ov and "rename" in ov:
                for idx, newname in ov["rename"].items():
                    if idx < len(header):
                        header[idx] = newname
            data_iter = it2

        drop = (ov or {}).get("drop_first_col_values", set())

        fname = OUT / f"{name}.tsv"
        n = 0
        with fname.open("w", encoding="utf-8", newline="") as fh:
            w = csv.writer(fh, delimiter="\t", lineterminator="\n")
            w.writerow(header)
            for r in data_iter:
                vals = [cell(v) for v in r[:ncol]]
                while len(vals) < ncol:
                    vals.append("")
                if not any(vals):
                    continue
                if vals[0] in drop:
                    continue
                w.writerow(vals)
                n += 1
        written.append((name, ncol, n, fname.name))

        qc_rows.append({
            "sheet": name,
            "header_detection": method,
            "manual_review_required": manual_review,
            "detected_header_row_1based": h + 1,
            "n_header_cols": ncol,
            "n_data_rows": n,
            "header_preview": " | ".join(header[:12]),
            "raw_first_8_rows": " ~~ ".join(
                " / ".join(cell(v) for v in r[:6]) for r in head_raw),
        })
        flag = "  [MANUAL REVIEW]" if manual_review else ""
        say(f"[ok] {name:<34} {method:<17} header_row={h+1:<3} cols={ncol:<3} rows={n}{flag}")

    with QC.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(qc_rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(qc_rows)
    say(f"[write] {QC.relative_to(ROOT)}  ({len(qc_rows)} sheets)")

    manifest = {
        "source": str(XLSX.relative_to(ROOT)),
        "n_sheets": len(written),
        "tables": [{"sheet": s, "cols": c, "rows": r, "file": f} for s, c, r, f in written],
    }
    (OUT / "_EXTRACT_MANIFEST.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8")

    expected = {"SuppTable1_response_groups": 362, "SuppTable_2_fusion_clas": 359,
                "SuppTable22_DE_analyses": 19000}
    fails = []
    got = {s: r for s, _, r, _ in written}
    for sheet, lo in expected.items():
        if got.get(sheet, 0) < lo * 0.95:
            fails.append(f"{sheet}: rows={got.get(sheet)} expected>={lo}")
    for f in fails:
        say(f"  [FAIL] {f}")
    say(f"[done] sheet_row_checks={'PASS' if not fails else 'FAIL'}")
    return 0 if not fails else 1


if __name__ == "__main__":
    sys.exit(main())
