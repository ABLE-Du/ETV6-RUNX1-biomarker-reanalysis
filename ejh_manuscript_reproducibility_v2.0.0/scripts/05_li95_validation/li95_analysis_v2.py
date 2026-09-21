"""Module A (final): Li95 independent-validation analysis, full 50-gene coverage.

Oksa EOI slow-vs-fast coefficients are read for EVERY gene in the official Oksa
Supplementary Table 22 (19,588 rows / 15,928 unique symbols), so the Oksa-vs-Li95
direction test uses all 50 genes with patient-level validation expression, not
only the 24 that happen to sit in the 219-gene primary universe.

Nothing is invented, rescaled across platforms, or imputed.
  - Oksa  : published EOI slow-vs-fast coefficient (Oksa Supplementary Table 22)
  - Li194 : C1-minus-C2 mean difference, recomputed from official Source Data
  - Li95  : C1-minus-C2 mean difference, recomputed from official Source Data
Only RANKS and SIGNS are compared across platforms; magnitudes are never pooled.
"""
import numpy as np, pandas as pd
from scipy import stats
import os, json

RNG = np.random.default_rng(20260919)
OUT = 'li95_out'

# ---------------------------------------------------------------- Oksa DE table
ok = pd.read_csv('../intermediate/oksa/SuppTable22_DE_analyses.tsv', sep='\t',
                 keep_default_na=False, dtype=str)
cols = list(ok.columns)
eoi_coef = [c for c in cols if c.startswith('EOI slow VS fast')][0]
eoi_p = [c for c in cols if c.startswith('EOI slow VS fast')][1]
eoi_fdr = [c for c in cols if c.startswith('EOI slow VS fast')][2]
print('Oksa EOI columns:', eoi_coef, '|', eoi_p, '|', eoi_fdr)

okx = ok[['Gene symbol', eoi_coef, eoi_p, eoi_fdr]].copy()
okx.columns = ['gene', 'oksa_coef', 'oksa_p', 'oksa_fdr']
for c in ['oksa_coef', 'oksa_p', 'oksa_fdr']:
    okx[c] = pd.to_numeric(okx[c], errors='coerce')
okx = okx.dropna(subset=['oksa_coef'])
okx['gene_u'] = okx['gene'].str.upper()
okx = okx.drop_duplicates('gene_u').set_index('gene_u')
print('Oksa genes with an EOI coefficient:', len(okx))

# ---------------------------------------------------------------- Li95
e95 = pd.read_csv(f'{OUT}/LI95_expression_long.tsv', sep='\t', keep_default_na=False)
e95 = e95[e95.expression != '']
e95['expr'] = pd.to_numeric(e95['expression'], errors='coerce')
lab95 = e95[['patient_id', 'validation_subtype']].drop_duplicates().set_index('patient_id')
sub95 = lab95['validation_subtype']
wide = e95.pivot_table(index='gene', columns='patient_id', values='expr', aggfunc='first')
c1 = [c for c in wide.columns if sub95[c] == 'C1']
c2 = [c for c in wide.columns if sub95[c] == 'C2']
print('Li95: n_C1', len(c1), 'n_C2', len(c2), 'n_genes', wide.shape[0])

d = pd.DataFrame(index=wide.index)
d['li95_mean_diff'] = wide[c1].mean(axis=1) - wide[c2].mean(axis=1)
d['li95_median_diff'] = wide[c1].median(axis=1) - wide[c2].median(axis=1)
mw = [stats.mannwhitneyu(wide.loc[g, c1], wide.loc[g, c2], alternative='two-sided')
      for g in wide.index]
d['li95_mwu_p'] = [x.pvalue for x in mw]
d['li95_rank_biserial'] = [2 * x.statistic / (len(c1) * len(c2)) - 1 for x in mw]
d['published_direction'] = e95.groupby('gene')['published_direction_in_validation'].first()
d['gene_u'] = [g.upper() for g in d.index]

# ---------------------------------------------------------------- Li194 (official Source Data)
m = pd.read_csv(f'{OUT}/LI194_matrix.tsv', sep='\t', keep_default_na=False, index_col=0)
lb = pd.read_csv(f'{OUT}/LI194_labels.tsv', sep='\t', keep_default_na=False).set_index('patient_id')
s = lb.loc[m.columns, 'Subtype']
c1d = [c for c in m.columns if s[c] == 'C1']; c2d = [c for c in m.columns if s[c] == 'C2']
M = m.apply(pd.to_numeric, errors='coerce')
li194 = pd.DataFrame({'li194_mean_diff': M[c1d].mean(axis=1) - M[c2d].mean(axis=1)})
li194['gene_u'] = [g.upper() for g in li194.index]
li194 = li194.drop_duplicates('gene_u').set_index('gene_u')

# ---------------------------------------------------------------- primary 219 membership
prim = pd.read_csv('../tables/08_OKSA_LI_gene_overlap.tsv', sep='\t', keep_default_na=False)
prim_u = set(prim.iloc[:, 1].astype(str).str.upper())
print('primary 219 symbols:', len(prim_u))

# ---------------------------------------------------------------- harmonisation
h = d.copy()
h = h.join(okx[['oksa_coef', 'oksa_p', 'oksa_fdr']], on='gene_u')
h = h.join(li194[['li194_mean_diff']], on='gene_u')
h['in_primary_219'] = h['gene_u'].isin(prim_u)
h = h.reset_index().rename(columns={'index': 'gene', 'level_0': 'gene'})
h = h[['gene', 'published_direction', 'in_primary_219', 'oksa_coef', 'oksa_p', 'oksa_fdr',
       'li194_mean_diff', 'li95_mean_diff', 'li95_median_diff', 'li95_rank_biserial', 'li95_mwu_p']]
p = h['li95_mwu_p'].to_numpy(float)
order = np.argsort(p); r = p[order]; mm = len(p)
q = np.minimum.accumulate((r * mm / (np.arange(mm) + 1))[::-1])[::-1]
qq = np.empty(mm); qq[order] = np.clip(q, 0, 1)
h['li95_bh_q'] = qq
h['oksa_direction'] = np.where(h.oksa_coef > 0, 'up_in_slow', 'down_in_slow')
h['li95_direction'] = np.where(h.li95_mean_diff > 0, 'up_in_C1', 'down_in_C1')
h['oksa_li95_match'] = np.sign(h.oksa_coef) == np.sign(h.li95_mean_diff)
h['li194_li95_match'] = np.where(h.li194_mean_diff.notna(),
                                 np.sign(h.li194_mean_diff) == np.sign(h.li95_mean_diff), None)
h = h.sort_values('gene')
h.to_csv(f'{OUT}/LI95_GENE_HARMONISATION.tsv', sep='\t', index=False)


def concordance_test(x, y, label, nperm=20000, nboot=10000):
    mask = np.isfinite(x) & np.isfinite(y)
    x = x[mask]; y = y[mask]; n = len(x)
    rho = stats.spearmanr(x, y).statistic
    conc = 100 * np.mean(np.sign(x) == np.sign(y))
    nr = np.empty(nperm); nc = np.empty(nperm)
    for i in range(nperm):
        pp = RNG.permutation(n)
        nr[i] = stats.spearmanr(x, y[pp]).statistic
        nc[i] = 100 * np.mean(np.sign(x) == np.sign(y[pp]))
    p_rho = (np.sum(np.abs(nr) >= abs(rho)) + 1) / (nperm + 1)
    p_conc = (np.sum(nc >= conc) + 1) / (nperm + 1)
    br = np.empty(nboot); bc = np.empty(nboot)
    for i in range(nboot):
        idx = RNG.integers(0, n, n)
        br[i] = stats.spearmanr(x[idx], y[idx]).statistic
        bc[i] = 100 * np.mean(np.sign(x[idx]) == np.sign(y[idx]))
    res = dict(label=label, n=int(n), rho=float(rho), concordance_pct=float(conc),
               null_rho_mean=float(nr.mean()), null_rho_sd=float(nr.std(ddof=1)),
               null_conc_mean=float(nc.mean()), null_conc_sd=float(nc.std(ddof=1)),
               p_rho=float(p_rho), p_concordance=float(p_conc),
               rho_ci=[float(np.percentile(br, 2.5)), float(np.percentile(br, 97.5))],
               conc_ci=[float(np.percentile(bc, 2.5)), float(np.percentile(bc, 97.5))])
    print(f'\n=== {label} (n={n}) ===')
    print(f'  rho = {rho:.4f}  bootstrap 95% CI [{res["rho_ci"][0]:.4f}, {res["rho_ci"][1]:.4f}]')
    print(f'  concordance = {conc:.1f}%  CI [{res["conc_ci"][0]:.1f}, {res["conc_ci"][1]:.1f}]')
    print(f'  null rho {nr.mean():.4f} +/- {nr.std(ddof=1):.4f}  P(rho) = {p_rho:.2e}')
    print(f'  null conc {nc.mean():.1f} +/- {nc.std(ddof=1):.1f}  P(conc) = {p_conc:.2e}')
    return res


R = {}
R['oksa_vs_li95_all50'] = concordance_test(h['oksa_coef'].to_numpy(float),
                                           h['li95_mean_diff'].to_numpy(float),
                                           'Oksa EOI slow-vs-fast  vs  Li95 C1-vs-C2  (all 50)')
sub219 = h[h.in_primary_219]
R['oksa_vs_li95_in219'] = concordance_test(sub219['oksa_coef'].to_numpy(float),
                                           sub219['li95_mean_diff'].to_numpy(float),
                                           'Oksa vs Li95 restricted to the 24 genes in the primary 219')
R['li194_vs_li95'] = concordance_test(h['li194_mean_diff'].to_numpy(float),
                                      h['li95_mean_diff'].to_numpy(float),
                                      'Li194 discovery  vs  Li95 validation  (genes present in both)')

# three-way
b = h.dropna(subset=['oksa_coef', 'li194_mean_diff', 'li95_mean_diff'])
s3 = ((np.sign(b.oksa_coef) == np.sign(b.li194_mean_diff)) &
      (np.sign(b.li194_mean_diff) == np.sign(b.li95_mean_diff)))
R['three_way'] = dict(n=int(len(b)), concordant_pct=float(100 * s3.mean()))
print(f'\nthree-way sign concordance (Oksa = Li194 = Li95): {100*s3.mean():.1f}%  (n={len(b)})')

# Tier-1 genes present in the validation panel
tier1 = ['SHE', 'GATA2', 'MID1', 'MPV17L', 'ANPEP', 'SERPINI2', 'UCK2']
t = h[h.gene.isin(tier1)]
print('\nTier-1 genes measurable in Li95:')
print(t[['gene', 'oksa_coef', 'oksa_fdr', 'li194_mean_diff', 'li95_mean_diff',
         'li95_rank_biserial', 'li95_bh_q', 'oksa_li95_match']].to_string(index=False))
R['tier1_in_li95'] = t[['gene', 'oksa_coef', 'oksa_fdr', 'li194_mean_diff', 'li95_mean_diff',
                        'li95_rank_biserial', 'li95_bh_q', 'oksa_li95_match']].to_dict('records')

R['li95_n'] = 95; R['li95_C1'] = 61; R['li95_C2'] = 34
R['li95_n_genes'] = int(wide.shape[0])
json.dump(R, open(f'{OUT}/_li95_results.json', 'w'), indent=1, default=str)
print('\nWROTE', f'{OUT}/_li95_results.json')
