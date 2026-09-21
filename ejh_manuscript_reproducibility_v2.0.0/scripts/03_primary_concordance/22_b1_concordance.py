#!/usr/bin/env python -u
# -*- coding: utf-8 -*-
"""22_b1_concordance.py -- B1: Oksa x Li transcriptional concordance.

Contrasts
    Oksa  : EOI slow vs fast coefficient (Supplementary Table 22, block
            'EOI slow VS fast').  Intermediate excluded by the source contrast.
    Li    : C1 - C2 effect, computed from the shipped row-standardised 252-gene
            matrix using the authors' ORIGINAL C1/C2 labels.

The gene set is the 219 genes that map onto BOTH cohorts
(`tables/08_OKSA_LI_gene_overlap.tsv`).  The 33 Li genes with no Oksa counterpart
are NOT_ASSESSABLE and are reported as such, never as negatives.

Statistics
    Spearman rho between the two effect vectors; direction concordance with an
    explicit positive-positive / negative-negative / discordant breakdown; a
    10,000-permutation test that permutes the **Li labels** (the Oksa contrast is
    fixed, because Oksa is the reference vector) and reports an empirical P; a
    rank-based concordance; and a robust (Huber) slope for the figure only.

    The permutation is the whole point: 219 genes selected because they are shared
    between a targeted panel and a transcriptome will not be an unbiased sample, so
    a bare rho is not interpretable on its own.  Permuting the cohort labels asks
    the question that matters -- is the agreement better than what the same two
    vectors would give if the Li patients were randomly relabelled?

Sensitivity (descriptions only; the primary result is all 219)
    A  all 219
    B  genes with adequate Li variability (|C1-C2| >= 0.5 row-SD units)
    C  Oksa FDR < 0.05 subset
    D  top |Oksa effect| quartile

Outputs
  results/B1_OKSA_LI_219_gene_concordance.tsv
  tables/B1_concordance_summary.tsv
  figures/B1_OKSA_vs_LI_effect_scatter.pdf / .png
  reports/B1_TRANSCRIPTIONAL_CONVERGENCE.md
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import _common as C  # noqa: E402

N_PERM = 10000
SEED = 20260918


def li_effects(mat: np.ndarray, labels: np.ndarray, rows: list[str]) -> np.ndarray:
    """C1 - C2 mean difference per gene."""
    c1 = labels == "C1"
    c2 = labels == "C2"
    return mat[:, c1].mean(axis=1) - mat[:, c2].mean(axis=1)


def main() -> int:
    C.ensure_dirs()
    rng = np.random.default_rng(SEED)

    ov = C.rd(C.TABLES / "08_OKSA_LI_gene_overlap.tsv")
    # the Oksa side is keyed by the Oksa symbol; the Li matrix by the Li original symbol
    genes = [g.strip() for g in ov["oksa_symbol"]]
    genes_li = [g.strip() for g in ov["li_original_symbol"]]
    ok_coef = pd.to_numeric(ov["oksa_coefficient_slow_vs_fast_EOI"], errors="coerce")
    ok_p = pd.to_numeric(ov["oksa_p_EOI"], errors="coerce")
    ok_fdr = pd.to_numeric(ov["oksa_fdr_EOI"], errors="coerce")

    li = C.rd(C.LI / "li_252gene_matrix.tsv")
    meta = C.rd(C.LI / "li_patient_meta.tsv")
    li_rows = [g.strip() for g in li["gene"]]
    patients = [c for c in li.columns if c != "gene"]
    labels = meta.set_index("patient_id").loc[patients, "subtype"].to_numpy()

    li_idx = {g.upper(): i for i, g in enumerate(li_rows)}
    ok_mask = ~ok_coef.isna().to_numpy()
    keep = [i for i, g in enumerate(genes_li) if ok_mask[i] and g.upper() in li_idx]
    n_mapped = len(genes_li)
    n_drop_no_li = int(sum(1 for i, g in enumerate(genes_li)
                           if ok_mask[i] and g.upper() not in li_idx))

    ok_vec = ok_coef.to_numpy()[keep]
    li_mat = li.drop(columns=["gene"]).to_numpy(dtype=float)
    li_rows_sel = np.array([li_idx[genes_li[i].upper()] for i in keep])
    sub = li_mat[li_rows_sel, :]
    genes = [genes[i] for i in keep]

    n_pat = len(patients)
    n_c1 = int((labels == "C1").sum())
    n_c2 = int((labels == "C2").sum())

    li_vec = li_effects(sub, labels, genes)

    # ---------------------------------------------------------- primary stats
    rho, rho_p = stats.spearmanr(ok_vec, li_vec)
    tau, _ = stats.kendalltau(ok_vec, li_vec)
    d = [C.direction(a, b) for a, b in zip(ok_vec, li_vec)]
    n_pp = d.count("positive_positive")
    n_nn = d.count("negative_negative")
    n_dis = d.count("discordant")
    n_zero = d.count("zero_effect")
    n_ass = n_pp + n_nn + n_dis
    conc_pct = 100 * (n_pp + n_nn) / n_ass if n_ass else float("nan")

    # ------------------------------------------------------ label permutation
    perm_rho = np.empty(N_PERM)
    perm_conc = np.empty(N_PERM)
    c1n = n_c1
    for i in range(N_PERM):
        pl = np.where(rng.permutation(n_pat) < c1n, "C1", "C2")
        v = li_effects(sub, pl, genes)
        perm_rho[i] = stats.spearmanr(ok_vec, v).statistic
        pd_ = [C.direction(a, b) for a, b in zip(ok_vec, v)]
        perm_conc[i] = 100 * (pd_.count("positive_positive") + pd_.count("negative_negative")) / n_ass
    emp_p_rho = float((np.abs(perm_rho) >= abs(rho)).mean())
    emp_p_conc = float((perm_conc >= conc_pct).mean())
    null_mean = float(np.nanmean(perm_rho))
    null_sd = float(np.nanstd(perm_rho))
    z_rho = float((rho - null_mean) / null_sd) if null_sd else float("nan")

    # robust slope for the figure only
    try:
        from statsmodels.api import RLM, add_constant
        fit = RLM(li_vec, add_constant(ok_vec)).fit()
        slope = float(fit.params[1])
        slope_ci = tuple(float(x) for x in fit.conf_int()[1])
    except Exception:  # noqa: BLE001
        slope, slope_ci = float("nan"), (float("nan"), float("nan"))

    # ------------------------------------------------------------ sensitivity
    ok_abs = np.abs(ok_vec)
    q75 = np.nanpercentile(ok_abs, 75)

    def stats_for(mask: np.ndarray) -> dict:
        if mask.sum() < 5:
            return dict(n=int(mask.sum()), rho=np.nan, concordant_pct=np.nan,
                        n_pp=int(np.nansum(mask & (ok_vec > 0) & (li_vec > 0))),
                        n_nn=int(np.nansum(mask & (ok_vec < 0) & (li_vec < 0))),
                        n_discordant=np.nan, empirical_p=np.nan)
        r = stats.spearmanr(ok_vec[mask], li_vec[mask]).statistic
        dd = [C.direction(a, b) for a, b in zip(ok_vec[mask], li_vec[mask])]
        pp, nn, di = dd.count("positive_positive"), dd.count("negative_negative"), dd.count("discordant")
        # same permutation design, restricted to the subset
        sub_m = sub[mask, :]
        cm = np.empty(N_PERM)
        for i in range(N_PERM):
            pl = np.where(rng.permutation(n_pat) < c1n, "C1", "C2")
            v = li_effects(sub_m, pl, genes)
            sp = stats.spearmanr(ok_vec[mask], v).statistic
            cm[i] = sp
        return dict(n=int(mask.sum()), rho=float(r), conc=100 * (pp + nn) / (pp + nn + di),
                    n_pp=pp, n_nn=nn, n_discordant=di,
                    empirical_p=float((np.abs(cm) >= abs(r)).mean()))

    sens = {}
    sA = stats_for(np.ones(len(ok_vec), dtype=bool))
    sens["A_all_219"] = sA
    sens["B_li_variability"] = stats_for(np.abs(li_vec) >= 0.5)
    sens["C_oksa_fdr_lt_0.05"] = stats_for((ok_fdr.to_numpy()[keep] < 0.05))
    sens["D_top_oksa_effect_quartile"] = stats_for(ok_abs >= q75)

    # RRHO feasibility: with 219 genes a rank-rank hypergeometric map is unstable,
    # so it is documented as not applicable rather than forced.
    rrho = dict(applicable=False,
                reason="n=219 is too small for a stable rank-rank hypergeometric overlap; "
                       "Spearman + permutation retained as the primary rank statistic")

    # ------------------------------------------------------------- outputs
    out = pd.DataFrame({
        "gene": [genes[i] for i in keep],
        "oksa_coefficient_slow_vs_fast_EOI": ok_vec,
        "oksa_p_EOI": ok_p.to_numpy()[keep],
        "oksa_fdr_EOI": ok_fdr.to_numpy()[keep],
        "li_C1_minus_C2": li_vec,
        "li_abs_effect": np.abs(li_vec),
        "direction_concordant": [C.direction(a, b) for a, b in zip(ok_vec, li_vec)],
        "rank_oksa": stats.rankdata(-ok_vec),
        "rank_li": stats.rankdata(-li_vec),
        "oksa_contrast": "EOI slow vs fast",
        "li_contrast": "C1 vs C2 (original labels)",
    }).sort_values("oksa_coefficient_slow_vs_fast_EOI")
    out.to_csv(C.RESULTS / "B1_OKSA_LI_219_gene_concordance.tsv", sep="\t", index=False)
    C.say(f"[write] results/B1_OKSA_LI_219_gene_concordance.tsv ({len(out)} rows)")

    summary = pd.DataFrame([
        dict(analysis="primary_all_219", n=sA["n"], spearman_rho=float(rho),
             kendall_tau=float(tau), concordant_pct=conc_pct, n_positive_positive=n_pp,
             n_negative_negative=n_nn, n_discordant=n_dis, n_zero_effect=n_zero,
             permutation_n=N_PERM, permutation_null_mean=null_mean,
             permutation_null_sd=null_sd, empirical_p_two_sided=emp_p_rho,
             empirical_p_concordance=emp_p_conc, z_vs_null=z_rho,
             robust_slope=slope, robust_slope_ci_low=slope_ci[0], robust_slope_ci_high=slope_ci[1],
             note="Oksa EOI slow vs fast; Li C1 vs C2; Li labels permuted"),
    ] + [
        dict(analysis=k, n=v["n"], spearman_rho=v.get("rho", np.nan), kendall_tau=np.nan,
             concordant_pct=v.get("conc", np.nan),
             n_positive_positive=v.get("n_pp", np.nan),
             n_negative_negative=v.get("n_nn", np.nan),
             n_discordant=v.get("n_discordant", np.nan), n_zero_effect=np.nan,
             permutation_n=N_PERM, permutation_null_mean=np.nan, permutation_null_sd=np.nan,
             empirical_p_two_sided=v.get("empirical_p", np.nan),
             empirical_p_concordance=np.nan, z_vs_null=np.nan,
             robust_slope=np.nan, robust_slope_ci_low=np.nan, robust_slope_ci_high=np.nan,
             note="sensitivity description; primary result is all 219")
        for k, v in sens.items() if k != "A_all_219"
    ] + [
        dict(analysis="RRHO_rank_rank_overlap", n=len(ok_vec), spearman_rho=np.nan,
             kendall_tau=np.nan, concordant_pct=np.nan, n_positive_positive=np.nan,
             n_negative_negative=np.nan, n_discordant=np.nan, n_zero_effect=np.nan,
             permutation_n=np.nan, permutation_null_mean=np.nan, permutation_null_sd=np.nan,
             empirical_p_two_sided=np.nan, empirical_p_concordance=np.nan, z_vs_null=np.nan,
             robust_slope=np.nan, robust_slope_ci_low=np.nan, robust_slope_ci_high=np.nan,
             note="NOT_APPLICABLE: " + rrho["reason"]),
    ])
    summary.to_csv(C.TABLES / "B1_concordance_summary.tsv", sep="\t", index=False)
    C.say(f"[write] tables/B1_concordance_summary.tsv")

    # ------------------------------------------------------------- figure
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(9.0, 3.9))
    ax = axes[0]
    col = {"positive_positive": "#B2182B", "negative_negative": "#2166AC",
           "discordant": "#A6A6A6", "zero_effect": "#D9D9D9"}
    for k, c in col.items():
        m = [x == k for x in out["direction_concordant"]]
        if any(m):
            ax.scatter(out.loc[m, "oksa_coefficient_slow_vs_fast_EOI"],
                       out.loc[m, "li_C1_minus_C2"], s=9, c=c, alpha=0.85,
                       linewidths=0.25, edgecolors="white", label=f"{k} (n={sum(m)})")
    xs = np.linspace(ok_vec.min(), ok_vec.max(), 50)
    ax.plot(xs, slope * xs + float(fit.params[0]), color="black", lw=1.0, ls="--",
            label=f"Huber slope {C.fmt(slope)}")
    ax.axhline(0, color="#CCCCCC", lw=0.6)
    ax.axvline(0, color="#CCCCCC", lw=0.6)
    ax.set_xlabel("Oksa: EOI slow vs fast coefficient")
    ax.set_ylabel("Li: C1 - C2 effect (row-SD units)")
    ax.set_title(f"219 shared genes\nSpearman rho = {C.fmt(rho)} (empirical P = {C.fmt(emp_p_rho)})",
                 fontsize=8)
    ax.legend(fontsize=5.4, loc="best", frameon=False)
    ax.tick_params(labelsize=6.5)

    ax2 = axes[1]
    ax2.hist(perm_rho, bins=40, color="#BFBFBF", edgecolor="white", linewidth=0.3)
    ax2.axvline(rho, color="#B2182B", lw=1.4,
                label=f"observed rho = {C.fmt(rho)}")
    ax2.axvline(null_mean, color="black", lw=0.9, ls=":", label=f"null mean = {C.fmt(null_mean)}")
    ax2.set_xlabel("Spearman rho under permuted Li labels")
    ax2.set_ylabel(f"count ({N_PERM} permutations)")
    ax2.set_title(f"Null distribution\nnull SD = {C.fmt(null_sd)}; z = {C.fmt(z_rho)}", fontsize=8)
    ax2.legend(fontsize=5.4, frameon=False)
    ax2.tick_params(labelsize=6.5)
    for a in axes:
        for s in ("top", "right"):
            a.spines[s].set_visible(False)
    fig.tight_layout()
    fig.savefig(C.FIGURES / "B1_OKSA_vs_LI_effect_scatter.pdf")
    fig.savefig(C.FIGURES / "B1_OKSA_vs_LI_effect_scatter.png", dpi=400)
    plt.close(fig)
    C.say("[write] figures/B1_OKSA_vs_LI_effect_scatter.pdf/.png")

    # ------------------------------------------------------------- report
    strongest = out.reindex(out["oksa_coefficient_slow_vs_fast_EOI"].abs()
                            .sort_values(ascending=False).index).head(15)
    verdict = ("the two independent ETV6::RUNX1 cohorts show cross-cohort convergence of an "
               "adverse-response transcriptional state within the 219 assessable genes"
               if (emp_p_rho < 0.05 and conc_pct > 50 and rho > 0) else
               "cross-cohort convergence is not statistically supported within the 219 assessable genes")

    md = f"""# B1 -- Oksa x Li transcriptional convergence

**Contrasts.** Oksa: **EOI slow vs fast** coefficient (Supplementary Table 22, block
`EOI slow VS fast`; the intermediate group is excluded by that contrast).
Li: **C1 - C2** effect, computed from the shipped row-standardised matrix using the
authors' **original** C1/C2 labels.

**Gene set.** The **{len(out)}** genes that map onto both cohorts
(`tables/08_OKSA_LI_gene_overlap.tsv`). Of Li's 252 genes, 33 have no Oksa counterpart;
they are **NOT_ASSESSABLE** and are listed as such in `08_GENE_HARMONISATION.tsv` — never
counted as negative results. {n_drop_no_li} mapped genes had no numeric Oksa coefficient
and were likewise excluded.

**Labels are used as published.** The Li cohort is never re-clustered and the C1/C2
assignment is never altered.

---

## 1. Primary result (all {len(out)} genes)

| quantity | value |
|---|---|
| genes analysed | {len(out)} |
| Li patients | {n_pat} (C1 = {n_c1}, C2 = {n_c2}) |
| Spearman rho | **{C.fmt(rho)}** (asymptotic P = {C.fmt(rho_p)}) |
| Kendall tau | {C.fmt(tau)} |
| direction concordant | **{C.fmt(conc_pct,1)}%** ({n_pp + n_nn} of {n_ass}) |
| positive-positive | {n_pp} |
| negative-negative | {n_nn} |
| discordant | {n_dis} |
| zero effect in either cohort | {n_zero} |
| Huber slope (figure only) | {C.fmt(slope)} ({C.fmt(slope_ci[0])} to {C.fmt(slope_ci[1])}) |

## 2. Permutation test — and why it is the number that matters

A bare rho is not interpretable here. The gene set was not sampled at random: it is the
intersection of a targeted panel's gene content and a curated 252-gene transcriptomic
panel, so the two effect vectors are not independent draws. The permutation therefore
holds the Oksa vector **fixed** and permutes the **Li labels** — randomly reassigning
which patients are C1 and which are C2, recomputing the Li effect vector each time — and
asks whether the observed agreement exceeds what the same two vectors produce under
relabelling.

| quantity | value |
|---|---|
| permutations | {N_PERM} |
| null mean rho | {C.fmt(null_mean)} |
| null SD | {C.fmt(null_sd)} |
| observed rho | {C.fmt(rho)} |
| z vs null | {C.fmt(z_rho)} |
| **empirical P (two-sided)** | **{C.fmt(emp_p_rho)}** |
| empirical P (concordance >= observed) | {C.fmt(emp_p_conc)} |

## 3. Sensitivity descriptions

The primary result is all {len(out)} genes. These subsets are descriptions only, and the
best-performing subset is **not** promoted to the headline.

| subset | n | rho | concordant % | pp | nn | discordant | empirical P |
|---|---|---|---|---|---|---|---|
"""
    for k, v in sens.items():
        md += (f"| {k} | {v['n']} | {C.fmt(v.get('rho'))} | {C.fmt(v.get('conc'),1)} | "
               f"{v.get('n_pp','')} | {v.get('n_nn','')} | {v.get('n_discordant','')} | "
               f"{C.fmt(v.get('empirical_p'))} |\n")

    md += f"""
## 4. RRHO / rank-rank overlap

**NOT_APPLICABLE.** {rrho['reason']}. Spearman rho plus the permutation test are retained as
the rank-based statistics.

## 5. Strongest genes by |Oksa effect|

Reported in full in `results/B1_OKSA_LI_219_gene_concordance.tsv`; the 15 largest |Oksa
effect| genes are shown here for orientation. This is a display, not a selection: the
per-gene table carries all {len(out)} genes with both effects, both P values and the
concordance label.

| gene | Oksa coef (slow-fast) | Oksa FDR | Li (C1-C2) | direction |
|---|---|---|---|---|
"""
    for _, r in strongest.iterrows():
        md += (f"| {r['gene']} | {C.fmt(r['oksa_coefficient_slow_vs_fast_EOI'])} | "
               f"{C.fmt(r['oksa_fdr_EOI'])} | {C.fmt(r['li_C1_minus_C2'])} | "
               f"{r['direction_concordant']} |\n")

    md += f"""
## 6. Interpretation rule applied

{"**Allowed and used:** " if emp_p_rho < 0.05 else "**Not allowed, not used:** "}
{verdict}.

The following statement is **not** made anywhere in this project, and must not be
extracted from this result:

> "a validated {len(out)}-gene signature was identified."

No score is constructed. `SIGNATURE_CONSTRUCTION_PERMITTED` remains FALSE, and this rung
produces no composite variable of any kind.

## 7. Limitations stated up front

1. **The Oksa side is a published effect vector, not our own fit.** Re-using the authors'
   coefficients is a re-analysis of their model, with their covariates and their
   normalisation. It is not an independent discovery.
2. **The Li side is a curated 252-gene panel**, row-standardised. Row standardisation
   removes absolute scale, so only direction and rank are comparable.
3. **The gene set is not a random sample** — see section 2.
4. **Only {len(out)} of 252 Li genes are assessable.** The rest are not assessable, which is
   different from negative.
5. **Different assays and different timepoints**: Oksa is day-29 EOI MRD in a panel/WGS
   cohort; Li is early-induction MRD in a St Jude cohort. The two phenotypes are not the
   same measurement, and this rung compares transcriptional states associated with them,
   not the phenotypes themselves.
"""
    (C.REPORTS / "B1_TRANSCRIPTIONAL_CONVERGENCE.md").write_text(md, encoding="utf-8")
    C.say(f"[write] reports/B1_TRANSCRIPTIONAL_CONVERGENCE.md")

    (C.TABLES / "_B1_CONCORDANCE.json").write_text(json.dumps({
        "n_genes": int(len(out)), "n_li_patients": n_pat, "n_C1": n_c1, "n_C2": n_c2,
        "spearman_rho": float(rho), "kendall_tau": float(tau),
        "concordant_pct": float(conc_pct), "n_pp": n_pp, "n_nn": n_nn,
        "n_discordant": n_dis, "n_zero": n_zero,
        "permutation_n": N_PERM, "null_mean": null_mean, "null_sd": null_sd,
        "empirical_p_two_sided": emp_p_rho, "z_vs_null": z_rho,
        "rrho_applicable": rrho["applicable"], "rrho_reason": rrho["reason"],
        "sensitivity": sens,
        "signature_constructed": False,
    }, indent=2, ensure_ascii=False), encoding="utf-8")

    C.say(f"[B1] n={len(out)} rho={rho:.3f} concordant={conc_pct:.1f}% "
          f"pp={n_pp} nn={n_nn} disc={n_dis} empP={emp_p_rho:.4f} z={z_rho:.2f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
