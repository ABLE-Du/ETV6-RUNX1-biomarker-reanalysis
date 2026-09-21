"""Module A: Li95 independent-validation analysis.

Questions
  1. Provenance check: recompute the frozen Oksa-Li194 219-gene correspondence
     from the OFFICIAL Nature Communications Source Data (not the project copy).
  2. Compute the Li95 (95-patient independent validation cohort) C1-vs-C2 effect
     vector for every gene with patient-level expression in the source data.
  3. For genes measurable in BOTH the primary 219 set and Li95, test whether the
     Oksa slow-response direction and the Li194 C1 direction persist in Li95.

Rules honoured
  - Li95 C1/C2 labels are the PUBLISHED PAM labels; not re-derived.
  - No gene re-selection, no threshold tuning.
  - Effect MAGNITUDES are never combined across platforms (discovery matrix is
    centred/scaled; validation matrix is log-scale) - only RANKS and SIGNS are
    compared across platforms.
  - Primary null: gene-identity permutation (isomorphic to the claim).
"""
import numpy as np, pandas as pd
from scipy import stats
import os, json

RNG = np.random.default_rng(20260919)
OUT = 'li95_out'
os.makedirs(OUT, exist_ok=True)

# ----------------------------------------------------------------- load Oksa x Li194 (frozen primary)
prim = pd.read_csv('../tables/08_OKSA_LI_gene_overlap.tsv', sep='\t', keep_default_na=False)
prim = prim.rename(columns={prim.columns[0]: 'ensembl'})
prim['oksa_coef'] = pd.to_numeric(prim['oksa_coefficient_slow_vs_fast_EOI'], errors='coerce')
prim['oksa_fdr'] = pd.to_numeric(prim['oksa_fdr_EOI'], errors='coerce')
print('primary 219 table rows:', len(prim))

# ----------------------------------------------------------------- Li194 from official source data
m = pd.read_csv(f'{OUT}/LI194_matrix.tsv', sep='\t', keep_default_na=False, index_col=0)
lab = pd.read_csv(f'{OUT}/LI194_labels.tsv', sep='\t', keep_default_na=False)
lab = lab.set_index('patient_id')
sub = lab.loc[m.columns, 'Subtype']
c1 = [c for c in m.columns if sub[c] == 'C1']
c2 = [c for c in m.columns if sub[c] == 'C2']
print('Li194 from source data: n_C1', len(c1), 'n_C2', len(c2), 'n_genes', m.shape[0])

M = m.apply(pd.to_numeric, errors='coerce')
li194_effect = M[c1].mean(axis=1) - M[c2].mean(axis=1)          # C1 minus C2, discovery scale
li194_median = M[c1].median(axis=1) - M[c2].median(axis=1)
u = np.array([stats.mannwhitneyu(M.loc[g, c1].dropna(), M.loc[g, c2].dropna(),
              alternative='two-sided').pvalue for g in M.index])
li194 = pd.DataFrame({'li194_mean_diff': li194_effect, 'li194_median_diff': li194_median,
                      'li194_mwu_p': u}, index=M.index)
li194.to_csv(f'{OUT}/LI194_gene_effects.tsv', sep='\t')

# ---------------- provenance check: recompute the frozen primary statistic
j = prim.set_index('oksa_symbol').join(li194, how='left')
j = j.dropna(subset=['oksa_coef', 'li194_mean_diff'])
rho_recomp = stats.spearmanr(j['oksa_coef'], j['li194_mean_diff']).statistic
rho_recomp_med = stats.spearmanr(j['oksa_coef'], j['li194_median_diff']).statistic
print(f'PROVENANCE CHECK  recomputed rho on n={len(j)} genes: '
      f'mean-diff {rho_recomp:.4f} | median-diff {rho_recomp_med:.4f}  (frozen = 0.7937)')

# ----------------------------------------------------------------- Li95 from official source data
e95 = pd.read_csv(f'{OUT}/LI95_expression_long.tsv', sep='\t', keep_default_na=False)
e95['expr'] = pd.to_numeric(e95['expression'], errors='coerce')
p95 = e95[e95.expression != '']
print('Li95: n_patients', p95.patient_id.nunique(), 'n_genes', p95.gene.nunique(),
      'n_cells', len(p95))
lab95 = p95[['patient_id', 'validation_subtype']].drop_duplicates()
print('Li95 labels:', lab95.validation_subtype.value_counts().to_dict())

wide = p95.pivot_table(index='gene', columns='patient_id', values='expr', aggfunc='first')
sub95 = lab95.set_index('patient_id')['validation_subtype']
c1_95 = [c for c in wide.columns if sub95[c] == 'C1']
c2_95 = [c for c in wide.columns if sub95[c] == 'C2']

li95 = pd.DataFrame(index=wide.index)
li95['li95_mean_diff'] = wide[c1_95].mean(axis=1) - wide[c2_95].mean(axis=1)
li95['li95_median_diff'] = wide[c1_95].median(axis=1) - wide[c2_95].median(axis=1)
li95['li95_mwu_p'] = [stats.mannwhitneyu(wide.loc[g, c1_95], wide.loc[g, c2_95],
                      alternative='two-sided').pvalue for g in wide.index]
# rank-biserial effect size (scale-free within platform)
li95['li95_rank_biserial'] = [2 * stats.mannwhitneyu(wide.loc[g, c1_95], wide.loc[g, c2_95],
                              alternative='two-sided').statistic / (len(c1_95) * len(c2_95)) - 1
                              for g in wide.index]
li95['published_direction'] = p95.groupby('gene')['published_direction_in_validation'].first()
li95.to_csv(f'{OUT}/LI95_gene_effects.tsv', sep='\t')

# ----------------------------------------------------------------- overlap + harmonisation
prim_sym = set(prim['oksa_symbol'])
li194_sym = set(M.index)
g95 = set(wide.index)
print('\nOverlap  validation50 & primary219 :', sorted(g95 & prim_sym))
print('Overlap  validation50 & Li194-252  :', len(g95 & li194_sym))
print('validation genes NOT in primary219 :', sorted(g95 - prim_sym))

harm = pd.DataFrame({'gene': sorted(g95)})
harm = harm.merge(prim[['oksa_symbol', 'ensembl', 'oksa_coef', 'oksa_fdr']],
                  left_on='gene', right_on='oksa_symbol', how='left').drop(columns=['oksa_symbol'])
harm = harm.merge(li95.reset_index(), on='gene', how='left')
harm = harm.merge(li194[['li194_mean_diff', 'li194_mwu_p']].reset_index().rename(
    columns={'index': 'gene'}), on='gene', how='left')
harm.to_csv(f'{OUT}/LI95_GENE_HARMONISATION.tsv', sep='\t', index=False)

# ----------------------------------------------------------------- concordance on the overlap
sub = harm.dropna(subset=['oksa_coef', 'li95_mean_diff']).copy()
n = len(sub)
print(f'\n=== Oksa vs Li95 on {n} overlapping genes ===')
rho = stats.spearmanr(sub['oksa_coef'], sub['li95_mean_diff']).statistic
conc = 100 * np.mean(np.sign(sub['oksa_coef']) == np.sign(sub['li95_mean_diff']))
sub['li194_95_both'] = sub['li194_mean_diff'].notna()
both = sub.dropna(subset=['li194_mean_diff'])
rho3 = stats.spearmanr(both['oksa_coef'], both['li95_mean_diff']).statistic
rho_li_li = stats.spearmanr(both['li194_mean_diff'], both['li95_mean_diff']).statistic
print(f'Spearman rho (Oksa vs Li95)              = {rho:.4f}')
print(f'direction concordance                    = {conc:.1f}%  ({int(conc*n/100)}/{n})')
print(f'rho on genes also in Li194 (n={len(both)})   = {rho3:.4f}')
print(f'rho Li194 vs Li95 (same genes)           = {rho_li_li:.4f}')

# three-way sign concordance
if len(both):
    s_ok = np.sign(both['oksa_coef']); s_l1 = np.sign(both['li194_mean_diff']); s_l9 = np.sign(both['li95_mean_diff'])
    three = 100 * np.mean((s_ok == s_l1) & (s_l1 == s_l9))
    print(f'three-way sign concordance (Oksa=Li194=Li95) = {three:.1f}%  (n={len(both)})')

# permutation null: gene identity
NPERM = 20000
ok = sub['oksa_coef'].to_numpy(float); l9 = sub['li95_mean_diff'].to_numpy(float)
null_rho = np.empty(NPERM); null_conc = np.empty(NPERM)
for i in range(NPERM):
    p = RNG.permutation(n)
    null_rho[i] = stats.spearmanr(ok, l9[p]).statistic
    null_conc[i] = 100 * np.mean(np.sign(ok) == np.sign(l9[p]))
p_rho = (np.sum(np.abs(null_rho) >= abs(rho)) + 1) / (NPERM + 1)
p_conc = (np.sum(null_conc >= conc) + 1) / (NPERM + 1)
print(f'\nnull rho  mean {null_rho.mean():.5f} sd {null_rho.std(ddof=1):.5f}  -> two-sided P = {p_rho:.2e}')
print(f'null conc mean {null_conc.mean():.2f} sd {null_conc.std(ddof=1):.2f}  -> upper P = {p_conc:.2e}')

# bootstrap CI over genes
NB = 10000
bs_r = np.empty(NB); bs_c = np.empty(NB)
for i in range(NB):
    idx = RNG.integers(0, n, n)
    bs_r[i] = stats.spearmanr(ok[idx], l9[idx]).statistic
    bs_c[i] = 100 * np.mean(np.sign(ok[idx]) == np.sign(l9[idx]))
print(f'bootstrap rho 95% CI [{np.percentile(bs_r,2.5):.4f}, {np.percentile(bs_r,97.5):.4f}]')
print(f'bootstrap conc 95% CI [{np.percentile(bs_c,2.5):.2f}, {np.percentile(bs_c,97.5):.2f}]')

# per-gene detail for the overlap
detail = sub[['gene', 'oksa_coef', 'oksa_fdr', 'li194_mean_diff', 'li95_mean_diff',
              'li95_median_diff', 'li95_rank_biserial', 'li95_mwu_p', 'published_direction']].copy()
detail['oksa_direction'] = np.where(detail.oksa_coef > 0, 'up_in_slow', 'down_in_slow')
detail['li95_direction'] = np.where(detail.li95_mean_diff > 0, 'up_in_C1', 'down_in_C1')
detail['direction_match'] = np.sign(detail.oksa_coef) == np.sign(detail.li95_mean_diff)
BH = detail['li95_mwu_p'].to_numpy(float)
order = np.argsort(BH); ranked = BH[order]; m_ = len(BH)
q = ranked * m_ / (np.arange(m_) + 1); q = np.minimum.accumulate(q[::-1])[::-1]
qq = np.empty(m_); qq[order] = np.clip(q, 0, 1)
detail['li95_bh_q'] = qq
detail = detail.sort_values('gene')
detail.to_csv(f'{OUT}/LI95_VALIDATION_RESULTS.tsv', sep='\t', index=False)
print('\nper-gene detail:')
print(detail[['gene', 'oksa_coef', 'oksa_fdr', 'li95_mean_diff', 'li95_bh_q',
              'direction_match']].to_string(index=False))

res = dict(
    n_genes_overlap=int(n),
    spearman_rho_oks_a_vs_li95=float(rho),
    direction_concordance_pct=float(conc),
    permutation_n=int(NPERM),
    null_rho_mean=float(null_rho.mean()), null_rho_sd=float(null_rho.std(ddof=1)),
    empirical_p_rho=float(p_rho), empirical_p_concordance=float(p_conc),
    bootstrap_rho_ci=[float(np.percentile(bs_r, 2.5)), float(np.percentile(bs_r, 97.5))],
    bootstrap_conc_ci=[float(np.percentile(bs_c, 2.5)), float(np.percentile(bs_c, 97.5))],
    provenance_rho_from_source_data=float(rho_recomp),
    provenance_rho_median_version=float(rho_recomp_med),
    frozen_primary_rho=0.7937,
    li95_n=95, li95_C1=61, li95_C2=34,
    n_genes_also_in_li194=int(len(both)),
    rho_on_li194_subset=float(rho3) if len(both) else None,
    rho_li194_vs_li95=float(rho_li_li) if len(both) else None,
    three_way_concordance_pct=float(three) if len(both) else None,
    overlap_genes=sorted(g95 & prim_sym),
)
json.dump(res, open(f'{OUT}/_li95_results.json', 'w'), indent=1)
print('\nWROTE _li95_results.json')
