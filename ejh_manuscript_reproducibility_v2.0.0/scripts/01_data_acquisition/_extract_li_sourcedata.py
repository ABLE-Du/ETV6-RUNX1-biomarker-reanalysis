"""Extract all Li et al. official Source Data needed for the EJH enhancement.

Source file: 41467_2025_56229_MOESM6_ESM.xlsx (Nature Communications Source Data)
Retrieved from Europe PMC PMC11779914 supplementaryFiles.

Outputs (all verbatim from source; nothing invented or imputed):
  li95_out/LI194_matrix.tsv          194 discovery patients x genes (with Subtype/Age/Sex/WBC rows)
  li95_out/LI194_labels.tsv
  li95_out/LI95_validation_age_table.tsv   (Figure 1e / 1f cross-tabs)
  li95_out/COG552_pam_calls.tsv      (Supp Figure 4, published context)
"""
import openpyxl, os, json

SRC = 'li95_raw/41467_2025_56229_MOESM6_ESM.xlsx'
os.makedirs('li95_out', exist_ok=True)


def sheet_rows(name):
    wb = openpyxl.load_workbook(SRC, read_only=True)
    ws = wb[name]
    rows = [list(r) for r in ws.iter_rows(values_only=True)]
    wb.close()
    return rows


# ---------------- discovery cohort 194 x genes ----------------
rows = sheet_rows('Supp Figure 1')
ids = [str(c) for c in rows[1][1:] if c is not None]
meta = {}
gene_rows = []
for r in rows[2:]:
    key = r[0]
    if key is None:
        continue
    key = str(key)
    vals = [('' if c is None else str(c)) for c in r[1:1 + len(ids)]]
    if key in ('Subtype', 'Age', 'Sex', 'WBC'):
        meta[key] = vals
    else:
        gene_rows.append((key, vals))
print('discovery: n_patients', len(ids), 'n_gene_rows', len(gene_rows), 'meta', list(meta))

with open('li95_out/LI194_matrix.tsv', 'w', encoding='utf-8') as f:
    f.write('gene\t' + '\t'.join(ids) + '\n')
    for g, vals in gene_rows:
        f.write(g + '\t' + '\t'.join(vals) + '\n')

with open('li95_out/LI194_labels.tsv', 'w', encoding='utf-8') as f:
    f.write('patient_id\tSubtype\tAge\tSex\tWBC\n')
    for i, pid in enumerate(ids):
        f.write('%s\t%s\t%s\t%s\t%s\n' % (pid,
                meta.get('Subtype', [''] * len(ids))[i],
                meta.get('Age', [''] * len(ids))[i],
                meta.get('Sex', [''] * len(ids))[i],
                meta.get('WBC', [''] * len(ids))[i]))

# ---------------- validation cohort age / WBC tables (Figure 1e, 1f) ----------
rows = sheet_rows('Figure 1')
hdr = rows[1]


def block(col0, ncols):
    out = []
    for r in rows[2:]:
        vals = r[col0:col0 + ncols]
        if vals[0] is None:
            continue
        out.append([('' if c is None else str(c)) for c in vals])
    return out


# Figure 1e: cols 16-18 (Age group, C1, C2); Figure 1f: cols 17? locate by header
idx_e = [i for i, h in enumerate(hdr) if h == 'Age group']
print('Age group header at cols', idx_e)
tables = {}
for c in idx_e:
    lab = hdr[c + 1] if c + 1 < len(hdr) else '?'
    lab2 = hdr[c + 2] if c + 2 < len(hdr) else '?'
    out = []
    for r in rows[2:]:
        if r[c] is None:
            continue
        out.append([str(r[c]), ('' if r[c + 1] is None else str(r[c + 1])),
                    ('' if r[c + 2] is None else str(r[c + 2]))])
    tables[(lab, lab2)] = out

with open('li95_out/LI95_validation_age_WBC_table.tsv', 'w', encoding='utf-8') as f:
    f.write('table_c1_col\ttable_c2_col\tlevel\tC1\tC2\n')
    for (a, b), out in tables.items():
        for lv, c1, c2 in out:
            f.write(f'{a}\t{b}\t{lv}\t{c1}\t{c2}\n')
for k, v in tables.items():
    print('  block', k, v)

# ---------------- COG 552 PAM calls (Supp Figure 4) ----------------
rows = sheet_rows('Supp Figure 4')
cog = []
for r in rows[2:]:
    if r[0] is None:
        continue
    cog.append([('' if c is None else str(c)) for c in r[:5]])
with open('li95_out/COG552_pam_calls.tsv', 'w', encoding='utf-8') as f:
    f.write('USI\tUMAP_1\tUMAP_2\tPAM_C1\tPAM_C2\n')
    for c in cog:
        f.write('\t'.join(c) + '\n')
n_c1 = sum(1 for c in cog if c[3] and float(c[3]) > 0.5)
print('COG552 rows', len(cog), 'n with PAM_C1>0.5', n_c1)

print('DONE')
