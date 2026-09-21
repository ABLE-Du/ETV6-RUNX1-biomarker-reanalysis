"""Export the Li et al. 95-patient validation-cohort patient-level data
from the official Nature Communications Source Data file.

Source: 41467_2025_56229_MOESM6_ESM.xlsx, sheet 'Supp Figure 10'
Downloaded from Europe PMC (PMC11779914) supplementaryFiles.

No values are invented, rescaled or imputed.
"""
import openpyxl, hashlib, os, json

SRC = 'li95_raw/41467_2025_56229_MOESM6_ESM.xlsx'
SHEET = 'Supp Figure 10'

wb = openpyxl.load_workbook(SRC, read_only=True)
ws = wb[SHEET]
rows = list(ws.iter_rows(values_only=True))
wb.close()

hdr = rows[1]
sub = rows[2]
ids = [str(c) for c in hdr[3:] if c is not None]
labs = [str(c) for c in sub[3:] if c is not None]
assert len(ids) == len(labs) == 95, (len(ids), len(labs))

genes, dirn, rank, mat = [], [], [], []
for r in rows[3:]:
    if r[2] is None:
        continue
    genes.append(str(r[2])); dirn.append(str(r[0])); rank.append(r[1])
    vals = []
    for c in r[3:3 + len(ids)]:
        vals.append('' if c is None or str(c).strip() == '' or str(c) == 'NA' else str(c))
    mat.append(vals)

os.makedirs('li95_out', exist_ok=True)

# long format: gene x patient
with open('li95_out/LI95_expression_long.tsv', 'w', encoding='utf-8') as f:
    f.write('gene\tpublished_direction_in_validation\trank\tpatient_id\tvalidation_subtype\texpression\n')
    for g, d, rk, vals in zip(genes, dirn, rank, mat):
        for pid, lb, v in zip(ids, labs, vals):
            f.write(f'{g}\t{d}\t{rk}\t{pid}\t{lb}\t{v}\n')

# wide format
with open('li95_out/LI95_expression_matrix.tsv', 'w', encoding='utf-8') as f:
    f.write('gene\t' + '\t'.join(ids) + '\n')
    for g, vals in zip(genes, mat):
        f.write(g + '\t' + '\t'.join(vals) + '\n')

with open('li95_out/LI95_patient_labels.tsv', 'w', encoding='utf-8') as f:
    f.write('patient_id\tvalidation_subtype\tid_format\n')
    for pid, lb in zip(ids, labs):
        fmt = 'sj_prefixed' if pid.startswith('SJ') else 'numeric'
        f.write(f'{pid}\t{lb}\t{fmt}\n')

# missingness per gene
miss = {}
for g, vals in zip(genes, mat):
    miss[g] = sum(1 for v in vals if v == '')
with open('li95_out/LI95_gene_missingness.tsv', 'w', encoding='utf-8') as f:
    f.write('gene\tpublished_direction\tn_missing_of_95\tpct_missing\n')
    for g in genes:
        f.write(f'{g}\t{dirn[genes.index(g)]}\t{miss[g]}\t{100*miss[g]/95:.1f}\n')

with open('li95_out/_sha256.txt', 'w', encoding='utf-8') as f:
    h = hashlib.sha256(open(SRC, 'rb').read()).hexdigest()
    f.write(f'{SRC}\t{h}\t{os.path.getsize(SRC)}\n')

print('n_genes', len(genes), 'n_patients', len(ids))
from collections import Counter
print('labels', Counter(labs))
print('id formats', Counter(['SJ' if p.startswith('SJ') else 'num' for p in ids]))
print('missing per gene (max)', max(miss.values()))
print('genes:', genes)
