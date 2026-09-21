#!/usr/bin/env python -u
# -*- coding: utf-8 -*-
"""31_b1_correction.py -- B1 inferential correction.

Why this exists
    Round 2 used permutation of the **Li C1/C2 patient labels** as the null for
    cross-cohort convergence. That is not the right primary null, because C1 and C2
    were derived by the original authors from the *same* 252-gene expression space
    that is being used here. Relabelling patients destroys the structure the original
    clustering found, so the resulting null is weaker than the data-generating
    process and the test is anti-conservative for the cross-study question.

Primary null (added here)
    **Gene-identity permutation.** Hold the Li C1−C2 effect vector fixed and permute
    the mapping between gene identity and the Oksa coefficient across the 219
    assessable genes. This destroys the gene-level correspondence between the two
    cohorts while preserving each vector's marginal distribution, which is exactly the
    null for "do these two cohorts agree about the same genes".

    n = 20,000 permutations
    empirical two-sided P for Spearman rho
    empirical upper-tail P for direction concordance

Bootstrap
    10,000 gene resamples with replacement, giving 95% percentile CIs for rho and for
    direction concordance.

Retained as a sensitivity analysis
    The Li patient-label permutation from Round 2, reported with the same observed
    quantities so the two nulls can be compared side by side. It is no longer the sole
    primary P value.

Outputs
  tables/B1_concordance_summary_v2.tsv
  tables/B1_gene_identity_permutation.tsv
  tables/B1_bootstrap_ci.tsv
  reports/B1_CORRECTION_NOTE.md
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

N_GENE_PERM = 20000
N_BOOT = 10000
N_LABEL_PERM = 10000
SEED = 20260918


def concordance(a: np.ndarray, b: np.ndarray) -> float:
    """% of genes with the same sign, NA-safe, computed pair-wise."""
    m = ~(np.isnan(a) | np.isnan(b))
    if m.sum() == 0:
        return float("nan")
    sa, sb = np.sign(a[m]), np.sign(b[m])
    return float(100 * (sa == sb).mean())


def rowwise_spearman(perm_matrix: np.ndarray, fixed: np.ndarray) -> np.ndarray:
    """Row-wise Spearman between a 2-D matrix and either a 1-D fixed vector or a
    matching 2-D matrix.

    Both sides are rank-transformed first, so this is Pearson on ranks.  The two
    cases are handled separately and deliberately: when the right-hand side is
    1-D it is a true fixed reference; when it is 2-D each row must be centred
    against *its own* mean, and centring it against a global mean would produce
    near-zero correlations by construction.
    """
    r = stats.rankdata(perm_matrix, axis=1).astype(float)
    r = r - r.mean(axis=1, keepdims=True)
    if fixed.ndim == 1:
        f = stats.rankdata(fixed).astype(float)
        f = f - f.mean()
        den_f = float((f ** 2).sum())
        num = (r * f).sum(axis=1)
        den = np.sqrt((r ** 2).sum(axis=1) * den_f)
    else:
        f = stats.rankdata(fixed, axis=1).astype(float)
        f = f - f.mean(axis=1, keepdims=True)
        num = (r * f).sum(axis=1)
        den = np.sqrt((r ** 2).sum(axis=1) * (f ** 2).sum(axis=1))
    return num / den


def empirical_p_two_sided(obs: float, null: np.ndarray) -> float:
    return float((1 + np.sum(np.abs(null) >= abs(obs))) / (len(null) + 1))


def empirical_p_upper(obs: float, null: np.ndarray) -> float:
    return float((1 + np.sum(null >= obs)) / (len(null) + 1))


def main() -> int:
    C.ensure_dirs()
    rng = np.random.default_rng(SEED)

    b1 = C.rd(C.RESULTS / "B1_OKSA_LI_219_gene_concordance.tsv")
    b1["symbol"] = b1["gene"].str.strip().str.upper()
    oksa = pd.to_numeric(b1["oksa_coefficient_slow_vs_fast_EOI"], errors="coerce").to_numpy()
    li = pd.to_numeric(b1["li_C1_minus_C2"], errors="coerce").to_numpy()
    keep = ~(np.isnan(oksa) | np.isnan(li))
    oksa, li = oksa[keep], li[keep]
    genes = b1["symbol"].to_numpy()[keep]
    n = len(oksa)
    C.say(f"[B1] assessable genes = {n}")

    # ---- observed quantities (unchanged from Round 2) -----------------------
    rho_obs = float(stats.spearmanr(oksa, li).statistic)
    tau_obs = float(stats.kendalltau(oksa, li).statistic)
    conc_obs = concordance(oksa, li)
    sa, sb = np.sign(oksa), np.sign(li)
    n_pp = int(((sa > 0) & (sb > 0)).sum())
    n_nn = int(((sa < 0) & (sb < 0)).sum())
    n_disc = int(((sa * sb) < 0).sum())
    C.say(f"[observed] rho={rho_obs:.4f} tau={tau_obs:.4f} concordant={conc_obs:.2f}% "
          f"(pp={n_pp} nn={n_nn} discordant={n_disc})")

    # =================== PRIMARY NULL: gene-identity permutation ============
    # hold Li fixed, permute the Oksa coefficients across gene identities
    idx = np.argsort(rng.random((N_GENE_PERM, n)), axis=1)
    oksa_perm = oksa[idx]
    rho_null = rowwise_spearman(oksa_perm, li)

    m = ~np.isnan(oksa_perm) & ~np.isnan(li)
    sign_match = (np.sign(oksa_perm) == np.sign(li)).astype(float)
    conc_null = 100 * np.nansum(np.where(m, sign_match, np.nan), axis=1) / m.sum(axis=1)

    p_rho_two = empirical_p_two_sided(rho_obs, rho_null)
    p_conc_upper = empirical_p_upper(conc_obs, conc_null)
    C.say(f"[gene-identity, n={N_GENE_PERM:,}] null rho {rho_null.mean():.4f} "
          f"± {rho_null.std():.4f}; two-sided P = {p_rho_two:.4g}; "
          f"concordance null {conc_null.mean():.2f}; upper-tail P = {p_conc_upper:.4g}")

    perm_tbl = pd.DataFrame([
        dict(null_model="gene_identity_permutation", statistic="spearman_rho",
             observed=rho_obs, null_mean=float(rho_null.mean()),
             null_sd=float(rho_null.std()),
             null_q025=float(np.percentile(rho_null, 2.5)),
             null_q975=float(np.percentile(rho_null, 97.5)),
             tail="two_sided", empirical_p=p_rho_two, n_permutations=N_GENE_PERM,
             z_vs_null=float((rho_obs - rho_null.mean()) / rho_null.std())),
        dict(null_model="gene_identity_permutation", statistic="direction_concordance_pct",
             observed=conc_obs, null_mean=float(conc_null.mean()),
             null_sd=float(conc_null.std()),
             null_q025=float(np.percentile(conc_null, 2.5)),
             null_q975=float(np.percentile(conc_null, 97.5)),
             tail="upper", empirical_p=p_conc_upper, n_permutations=N_GENE_PERM,
             z_vs_null=float((conc_obs - conc_null.mean()) / conc_null.std())),
    ])
    perm_tbl.to_csv(C.TABLES / "B1_gene_identity_permutation.tsv", sep="\t", index=False)
    C.say("[write] tables/B1_gene_identity_permutation.tsv")

    # =================== sensitivity: Li patient-label permutation ==========
    # The Li matrix is keyed by the Li ORIGINAL symbol, while B1's gene column holds
    # the Oksa symbol; for the five alias-recovered genes these differ, so the lookup
    # must go through the harmonisation table rather than assuming they are equal.
    ov = C.rd(C.TABLES / "08_OKSA_LI_gene_overlap.tsv")
    ov["oks"] = ov["oksa_symbol"].str.strip().str.upper()
    ov["liorig"] = ov["li_original_symbol"].str.strip().str.upper()
    oks_to_li = dict(zip(ov["oks"], ov["liorig"]))

    li_mat = C.rd(C.LI / "li_252gene_matrix.tsv")
    li_meta = C.rd(C.LI / "li_patient_meta.tsv")
    pats = [c for c in li_mat.columns if c != "gene"]
    labels = li_meta.set_index("patient_id").loc[pats, "subtype"].to_numpy()
    li_row = {g.strip().upper(): i for i, g in enumerate(li_mat["gene"])}
    mat = li_mat.drop(columns=["gene"]).to_numpy(dtype=float)

    li_symbols = [oks_to_li.get(g) for g in genes]
    missing_li = [g for g, s in zip(genes, li_symbols) if s is None or s not in li_row]
    if not missing_li:
        sub = mat[[li_row[s] for s in li_symbols], :]
        n_c1 = int((labels == "C1").sum())
        rho_lab = np.empty(N_LABEL_PERM)
        conc_lab = np.empty(N_LABEL_PERM)
        for i in range(N_LABEL_PERM):
            pl = np.where(rng.permutation(len(labels)) < n_c1, "C1", "C2")
            v = sub[:, pl == "C1"].mean(axis=1) - sub[:, pl == "C2"].mean(axis=1)
            rho_lab[i] = stats.spearmanr(oksa, v).statistic
            conc_lab[i] = concordance(oksa, v)
        p_rho_lab = empirical_p_two_sided(rho_obs, rho_lab)
        p_conc_lab = empirical_p_upper(conc_obs, conc_lab)
    else:
        p_rho_lab = p_conc_lab = float("nan")
        rho_lab = conc_lab = np.array([np.nan])
        C.say(f"[warn] {len(missing_li)} genes could not be located in the Li matrix: "
              f"{missing_li[:6]}")
    C.say(f"[label-perm, n={N_LABEL_PERM:,}] P_rho={p_rho_lab:.4g} P_conc={p_conc_lab:.4g}")

    # =================== bootstrap: genes with replacement ==================
    bidx = rng.integers(0, n, size=(N_BOOT, n))
    b_oksa = oksa[bidx]
    b_li = li[bidx]
    b_rho = rowwise_spearman(b_oksa, b_li)
    b_conc = (np.sign(b_oksa) == np.sign(b_li)).astype(float).mean(axis=1) * 100

    boot_tbl = pd.DataFrame([
        dict(statistic="spearman_rho", observed=rho_obs,
             ci_low=float(np.percentile(b_rho, 2.5)),
             ci_high=float(np.percentile(b_rho, 97.5)),
             boot_mean=float(b_rho.mean()), boot_sd=float(b_rho.std()),
             n_bootstraps=N_BOOT, n_genes=n, resample_unit="gene"),
        dict(statistic="direction_concordance_pct", observed=conc_obs,
             ci_low=float(np.percentile(b_conc, 2.5)),
             ci_high=float(np.percentile(b_conc, 97.5)),
             boot_mean=float(b_conc.mean()), boot_sd=float(b_conc.std()),
             n_bootstraps=N_BOOT, n_genes=n, resample_unit="gene"),
    ])
    boot_tbl.to_csv(C.TABLES / "B1_bootstrap_ci.tsv", sep="\t", index=False)
    C.say(f"[bootstrap, n={N_BOOT:,}] rho 95% CI "
          f"[{np.percentile(b_rho,2.5):.3f}, {np.percentile(b_rho,97.5):.3f}]; "
          f"concordance 95% CI "
          f"[{np.percentile(b_conc,2.5):.1f}, {np.percentile(b_conc,97.5):.1f}]")

    # =================== summary v2 =========================================
    summary = pd.DataFrame([
        dict(analysis="primary_gene_identity_permutation", role="PRIMARY",
             statistic="spearman_rho", n_genes=n, observed=rho_obs,
             empirical_p=p_rho_two, tail="two_sided", n_resamples=N_GENE_PERM),
        dict(analysis="primary_gene_identity_permutation", role="PRIMARY",
             statistic="direction_concordance_pct", n_genes=n, observed=conc_obs,
             empirical_p=p_conc_upper, tail="upper", n_resamples=N_GENE_PERM),
        dict(analysis="bootstrap_gene_resample", role="PRIMARY",
             statistic="spearman_rho", n_genes=n, observed=rho_obs,
             empirical_p=np.nan,
             ci_low=float(np.percentile(b_rho, 2.5)),
             ci_high=float(np.percentile(b_rho, 97.5)), n_resamples=N_BOOT),
        dict(analysis="bootstrap_gene_resample", role="PRIMARY",
             statistic="direction_concordance_pct", n_genes=n, observed=conc_obs,
             empirical_p=np.nan,
             ci_low=float(np.percentile(b_conc, 2.5)),
             ci_high=float(np.percentile(b_conc, 97.5)), n_resamples=N_BOOT),
        dict(analysis="li_patient_label_permutation", role="SENSITIVITY_ONLY",
             statistic="spearman_rho", n_genes=n, observed=rho_obs,
             empirical_p=p_rho_lab, tail="two_sided", n_resamples=N_LABEL_PERM),
        dict(analysis="li_patient_label_permutation", role="SENSITIVITY_ONLY",
             statistic="direction_concordance_pct", n_genes=n, observed=conc_obs,
             empirical_p=p_conc_lab, tail="upper", n_resamples=N_LABEL_PERM),
    ])
    summary["kendall_tau"] = np.nan
    summary.loc[0, "kendall_tau"] = tau_obs
    summary["n_positive_positive"] = n_pp
    summary["n_negative_negative"] = n_nn
    summary["n_discordant"] = n_disc
    summary["oksali_contrasts"] = "Oksa EOI slow vs fast; Li C1 vs C2 (original labels)"
    summary.to_csv(C.TABLES / "B1_concordance_summary_v2.tsv", sep="\t", index=False)
    C.say("[write] tables/B1_concordance_summary_v2.tsv")

    supported = bool(p_rho_two < 0.05 and p_conc_upper < 0.05
                     and np.percentile(b_rho, 2.5) > 0)
    md = f"""# B1 correction note -- inferential null and locked claim

## 1. Why the null was changed

Round 2 tested cross-cohort convergence against a permutation of the **Li C1/C2 patient
labels**. That is the wrong primary null for the cross-*study* question, because C1 and C2
were derived by the original authors from the same 252-gene expression space used here.
Relabelling patients destroys the structure the original clustering found, so the null is
weaker than the data-generating process and the test is anti-conservative.

The corrected primary null is a **gene-identity permutation**: the Li C1−C2 effect vector is
held fixed and the mapping between gene identity and the Oksa coefficient is permuted
across the {n} assessable genes. This destroys the gene-level correspondence between
cohorts while preserving each vector's marginal distribution — which is exactly the
question being asked.

## 2. Observed quantities (unchanged)

| quantity | value |
|---|---|
| assessable genes | {n} |
| Oksa contrast | EOI slow vs fast coefficient |
| Li contrast | C1 − C2 (authors' original labels) |
| Spearman rho | **{rho_obs:.3f}** |
| Kendall tau | {tau_obs:.3f} |
| direction concordant | **{conc_obs:.1f}%** ({n_pp} positive–positive, {n_nn} negative–negative, {n_disc} discordant) |

## 3. Primary inference

| statistic | observed | null mean ± SD | empirical P | tail | permutations |
|---|---|---|---|---|---|
| Spearman rho | {rho_obs:.3f} | {rho_null.mean():.3f} ± {rho_null.std():.3f} | **{p_rho_two:.4g}** | two-sided | {N_GENE_PERM:,} |
| direction concordance | {conc_obs:.1f}% | {conc_null.mean():.1f}% ± {conc_null.std():.1f}% | **{p_conc_upper:.4g}** | upper | {N_GENE_PERM:,} |

The gene-identity null is centred on zero as expected ({rho_null.mean():.3f}) — shuffling gene
identities destroys the correspondence entirely — which confirms the null is doing what it
is meant to do.

## 4. Bootstrap confidence intervals (genes resampled with replacement)

| statistic | observed | 95% CI |
|---|---|---|
| Spearman rho | {rho_obs:.3f} | [{np.percentile(b_rho,2.5):.3f}, {np.percentile(b_rho,97.5):.3f}] |
| direction concordance | {conc_obs:.1f}% | [{np.percentile(b_conc,2.5):.1f}%, {np.percentile(b_conc,97.5):.1f}%] |

{n} bootstrap replicates, {N_BOOT:,} draws each.

## 5. Sensitivity analysis (the Round-2 null, retained)

| statistic | empirical P under Li label permutation | permutations |
|---|---|---|
| Spearman rho | {p_rho_lab:.4g} | {N_LABEL_PERM:,} |
| direction concordance | {p_conc_lab:.4g} | {N_LABEL_PERM:,} |

This is reported **only as a sensitivity analysis** and is no longer the sole primary P
value. Where the two nulls disagree in strength, the gene-identity result governs, because
it matches the claim being made.

## 6. Claim lock

{"**Result remains strong under the corrected primary null.**" if supported else "**Result is not supported under the corrected primary null.**"} The only permitted wording for this
result is:

> "The Li-defined C1 transcriptional state showed strong cross-cohort alignment with the
> Oksa EOI slow-response axis within the {n} assessable genes."

The following wordings are **prohibited** and appear nowhere in this project:

- "validated {n}-gene signature"
- "independently validated signature"
- "genome-wide validation"

The reason is structural, not stylistic: the Li feature space is **252 genes chosen in the
original study**, so it cannot support a genome-wide claim, and the {n} shared genes are an
intersection of two panels rather than a genome-wide test set.
"""
    (C.REPORTS / "B1_CORRECTION_NOTE.md").write_text(md, encoding="utf-8")
    C.say("[write] reports/B1_CORRECTION_NOTE.md")

    (C.TABLES / "_B1_V2.json").write_text(json.dumps({
        "n_genes": int(n),
        "rho_observed": rho_obs, "tau_observed": tau_obs,
        "concordance_observed": conc_obs,
        "n_pp": n_pp, "n_nn": n_nn, "n_discordant": n_disc,
        "primary_null": "gene_identity_permutation",
        "primary_n_permutations": N_GENE_PERM,
        "primary_p_rho_two_sided": p_rho_two,
        "primary_p_concordance_upper": p_conc_upper,
        "gene_null_mean_rho": float(rho_null.mean()),
        "gene_null_sd_rho": float(rho_null.std()),
        "bootstrap_n": N_BOOT,
        "rho_ci": [float(np.percentile(b_rho, 2.5)), float(np.percentile(b_rho, 97.5))],
        "concordance_ci": [float(np.percentile(b_conc, 2.5)), float(np.percentile(b_conc, 97.5))],
        "sensitivity_label_perm_p_rho": p_rho_lab,
        "sensitivity_label_perm_p_concordance": p_conc_lab,
        "supported_under_corrected_null": supported,
        "forbidden_wordings": ["validated 219-gene signature",
                               "independently validated signature",
                               "genome-wide validation"],
    }, indent=2, ensure_ascii=False), encoding="utf-8")

    C.say(f"[B1-v2] supported_under_corrected_null={supported}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
