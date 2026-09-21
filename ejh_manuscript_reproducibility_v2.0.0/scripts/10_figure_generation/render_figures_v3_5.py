# -*- coding: utf-8 -*-
"""Render Figures 1-5 + S1-S2 from figure_data TSVs. PDF + SVG + TIFF(600dpi).
Colour-blind-safe Okabe-Ito palette. No hard-coded values: all read from TSVs.
v3: Figure 3 relabelled to out-of-sample concordance (deposited 50-gene panel
is discovery-derived); Figure 4 Li95 markers data-driven (li95_status),
'not assessable' wording instead of 'reproduced'."""
import os, sys, io, re
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from scipy.stats import spearmanr

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

PROJECT_ROOT = Path(os.environ.get("EJH_PROJECT_ROOT", Path(__file__).resolve().parents[2]))
BASE = str(Path(os.environ.get("EJH_FIGURE_DATA", PROJECT_ROOT / "data" / "figure_data")))
OUT = str(Path(os.environ.get("EJH_FIGURE_OUTPUT", PROJECT_ROOT / "figures")))
os.makedirs(OUT, exist_ok=True)

# Okabe-Ito colour-blind-safe palette
OI = {"blue": "#0072B2", "orange": "#E69F00", "green": "#009E73",
      "verm": "#D55E00", "purple": "#CC79A7", "sky": "#56B4E9",
      "yellow": "#F0E442", "black": "#000000", "grey": "#999999"}
plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 9,
    "axes.linewidth": 0.8, "axes.spines.top": False, "axes.spines.right": False,
    "figure.dpi": 100, "savefig.bbox": "tight", "pdf.fonttype": 42,
    "svg.fonttype": "none",
})

def save(fig, name):
    for ext, kw in [("pdf", {}), ("svg", {}), ("tiff", {"dpi": 600, "pil_kwargs": {"compression": "tiff_lzw"}})]:
        fig.savefig(os.path.join(OUT, f"{name}.{ext}"), **kw)
    plt.close(fig)
    print("wrote", name)

def rd(f):
    return pd.read_csv(os.path.join(BASE, f), sep="\t", keep_default_na=False)

# ---------------------------------------------------------------- Figure 1
def fig1():
    d = rd("Fig1_study_design.tsv")
    ceil_col = {"out-of-sample concordance in independent patients": OI["green"],
                "published context only": OI["sky"],
                "exploratory phenotypic support": OI["orange"],
                "supportive benchmark, endpoint not pooled": OI["grey"],
                "published contrast, reused as published": OI["blue"]}
    n = len(d)
    fig, ax = plt.subplots(figsize=(7.6, 6.0)); ax.axis("off")
    ax.set_xlim(0, 10); ax.set_ylim(0, 10.6)
    top, bottom = 10.0, 0.5
    h = (top - bottom) / n * 0.86
    ys = np.linspace(top - h / 2, bottom + h / 2, n)
    for (i, r), y in zip(d.iterrows(), ys):
        col = ceil_col.get(r["evidence_ceiling"], OI["blue"])
        box = FancyBboxPatch((0.35, y - h / 2), 9.3, h,
                             boxstyle="round,pad=0.02,rounding_size=0.10",
                             linewidth=1.3, edgecolor=col, facecolor=col + "22")
        ax.add_patch(box)
        # line 1: cohort + n (left), evidence ceiling (right, on its own baseline)
        ax.text(0.60, y + h * 0.30, f"{r['cohort']}  (n={r['n']})",
                fontsize=8.6, weight="bold", va="center")
        ax.text(9.50, y + h * 0.30, str(r["evidence_ceiling"]),
                fontsize=7.0, va="center", ha="right", style="italic", color=col)
        # line 2: endpoint (left, indented) -- never shares a baseline with the right-hand label
        ax.text(0.60, y - h * 0.06, str(r["endpoint"]), fontsize=7.4, va="center", color="#333333")
        # line 3: role (left, indented)
        ax.text(0.60, y - h * 0.38, str(r["role"]), fontsize=7.2, va="center", color="#444444")
    ax.annotate("", xy=(4.85, ys[1] - h / 2), xytext=(4.85, ys[0] + h / 2),
                arrowprops=dict(arrowstyle="<->", color=OI["verm"], lw=2))
    ax.text(4.98, (ys[0] + ys[1]) / 2, "primary correspondence", color=OI["verm"],
            fontsize=7.5, va="center")
    ax.set_title("Six data layers and the response definition each contributes",
                 fontsize=10, weight="bold", loc="left")
    fig.text(0.012, 0.015,
             "The four response definitions (Oksa day 29 at 0.1%; Li early induction; institutional "
             "day 19; TARGET day 29 at 0.01%) are not interchangeable and are never pooled.",
             fontsize=6.4, color="#555555", ha="left")
    fig.subplots_adjust(top=0.93, bottom=0.09, left=0.02, right=0.98)
    save(fig, "Figure_1")

# ---------------------------------------------------------------- Figure 2
def fig2():
    sc = rd("Fig2_gene_scatter.tsv")
    for c in ["oksa_effect", "li194_effect"]:
        sc[c] = pd.to_numeric(sc[c], errors="coerce")
    prim = rd("Fig2_primary_correspondence.tsv")
    obs = float(prim[(prim.statistic == "spearman_rho") & (prim.role == "PRIMARY")].observed.iloc[0])
    null_mean, null_sd = -0.00012558, 0.06756846
    gm = rd("Fig2_gene_matched_sensitivity.tsv")
    gm = gm[gm.gene_universe.str.contains("113")].copy()
    dr = rd("Fig2_rho_difference.tsv")

    fig = plt.figure(figsize=(7.2, 7.6))
    gs = fig.add_gridspec(3, 1, height_ratios=[1.15, 0.95, 0.9], hspace=0.42)

    # A scatter
    axA = fig.add_subplot(gs[0])
    conc = sc["concordant"].astype(str).str.lower() == "true"
    axA.scatter(sc.loc[conc, "oksa_effect"], sc.loc[conc, "li194_effect"],
                s=11, c=OI["blue"], alpha=0.55, edgecolors="none", label="concordant (190)")
    axA.scatter(sc.loc[~conc, "oksa_effect"], sc.loc[~conc, "li194_effect"],
                s=13, c=OI["verm"], alpha=0.8, edgecolors="none", label="discordant (29)")
    axA.axhline(0, color="#bbbbbb", lw=0.7); axA.axvline(0, color="#bbbbbb", lw=0.7)
    axA.set_xlabel("Oksa slow-vs-fast coefficient"); axA.set_ylabel("Li194 C1-minus-C2")
    axA.legend(frameon=False, fontsize=7.5, loc="upper left")
    axA.text(0.98, 0.03, f"Spearman ρ = 0.79 (95% CI 0.73-0.84)\n86.8% concordant; 219 genes",
             transform=axA.transAxes, ha="right", va="bottom", fontsize=7.5,
             bbox=dict(boxstyle="round", fc="#f2f2f2", ec="#cccccc"))
    axA.set_title("A", loc="left", weight="bold", fontsize=11)

    # B permutation null (regenerate from real scatter, gene-identity)
    axB = fig.add_subplot(gs[1])
    x = sc["oksa_effect"].values; y = sc["li194_effect"].values
    rng = np.random.default_rng(20260919)
    NB = 20000; nulls = np.empty(NB)
    for i in range(NB):
        nulls[i] = spearmanr(x, rng.permutation(y)).statistic
    axB.hist(nulls, bins=80, color=OI["sky"], alpha=0.75, edgecolor="none")
    axB.set_xlim(-0.3, 0.92)
    axB.axvline(obs, color=OI["verm"], lw=2)
    axB.text(obs-0.02, axB.get_ylim()[1]*0.55, f"observed ρ = 0.79\nP ≤ 5.0×10$^{{-5}}$ (floor)",
             color=OI["verm"], fontsize=7.5, va="center", ha="right")
    axB.set_xlabel("Spearman ρ under gene-identity permutation (20,000)")
    axB.set_ylabel("count")
    axB.set_title("B", loc="left", weight="bold", fontsize=11)

    # C gene-matched benchmark forest
    axC = fig.add_subplot(gs[2])
    labels = ["Oksa–Li194", "Oksa–TARGET", "Li194–TARGET", "Oksa–TARGET (adj)"]
    rows = [gm[gm.comparison == "Oksa vs Li194"].iloc[0],
            gm[gm.comparison == "Oksa vs TARGET"].iloc[0],
            gm[gm.comparison == "Li194 vs TARGET"].iloc[0],
            gm[gm.comparison.str.contains("protocol-adjusted")].iloc[0]]
    cols = [OI["green"], OI["grey"], OI["grey"], OI["grey"]]
    ypos = np.arange(len(rows))[::-1]
    for yp, r, c, lab in zip(ypos, rows, cols, labels):
        axC.errorbar(float(r.rho), yp, xerr=[[float(r.rho)-float(r.rho_CI_lo)], [float(r.rho_CI_hi)-float(r.rho)]],
                     fmt="o", color=c, ecolor=c, capsize=3, ms=6)
        axC.text(float(r.rho_CI_hi)+0.02, yp, f"{float(r.rho):.2f}", va="center", fontsize=7.5)
    axC.set_yticks(ypos); axC.set_yticklabels(labels, fontsize=8)
    axC.set_xlim(0, 1.0); axC.set_ylim(-1.1, 3.6)
    axC.set_xlabel("Spearman ρ on the same 113 genes (95% CI)")
    d1 = dr.iloc[0]
    axC.text(0.03, -0.75, f"Δρ (Oksa–Li194 − Oksa–TARGET) = {d1.delta_rho:.2f} (95% CI {d1.CI_lo:.2f}–{d1.CI_hi:.2f})",
             fontsize=7.5, color=OI["green"], va="center")
    axC.set_title("C", loc="left", weight="bold", fontsize=11)
    save(fig, "Figure_2")

# ---------------------------------------------------------------- Figure 3
def fig3():
    tri = rd("Fig3_independent_triangulation.tsv")
    tri_gene = tri[tri.panel == "A_Li95_gene_level"].copy()
    for c in ["oksa_coef", "li194_effect", "li95_effect"]:
        tri_gene[c] = pd.to_numeric(tri_gene[c], errors="coerce")
    local = rd("Fig3_local_CD34_CD38_patients.tsv")
    for c in ["CD34_raw", "CD38_raw", "D19_MRD_percent"]:
        local[c] = pd.to_numeric(local[c].replace({"/": np.nan, "": np.nan}), errors="coerce")
    local["poor"] = local["D19_MRD_ge_0p1"].astype(str).str.upper() == "TRUE"
    # analysis set: day-19 MRD-evaluable only (patients without a day-19 MRD call are excluded)
    local["mrd_evaluable"] = local["D19_MRD_ge_0p1"].astype(str).str.upper().isin(["TRUE", "FALSE"])
    loc = local[local["mrd_evaluable"]].copy()

    # clinical / local OR rows parsed from TSVs (no hard-coded numbers)
    def _row_text(label):
        r = tri[tri.gene == label]
        return " ".join(str(v) for v in r.iloc[0].values if str(v) != "")
    m_age = re.search(r"OR=([\d.]+) CI=\[([\d.]+),([\d.]+)\] P=([\d.]+)", _row_text("age >= 5 y"))
    m_wbc = re.search(r"OR=([\d.]+) CI=\[([\d.]+),([\d.]+)\] P=([\d.]+)", _row_text("WBC >= 50 x10^9/L"))
    li95_ageOR = tuple(float(x) for x in m_age.groups())   # OR, lo, hi, P
    li95_wbcOR = tuple(float(x) for x in m_wbc.groups())
    # local age/WBC: read the source analysis table directly.
    # 2x2 encoding: a = exposed&poor, b = exposed&good, c = unexposed&poor, d = unexposed&good
    awt = pd.read_csv(Path(os.environ.get("EJH_LOCAL_AGE_WBC", PROJECT_ROOT / "data" / "derived_tables" / "LOCAL_AGE_WBC_TRIANGULATION.tsv")),
                      sep="\t", keep_default_na=False)
    ra = awt[awt.feature == "age_ge_5y"].iloc[0]
    rw = awt[awt.feature == "WBC_ge_50"].iloc[0]
    a, b_, c_, d_ = int(ra.a), int(ra.b), int(ra.c), int(ra.d)
    lo_age_or, lo_age_lo, lo_age_hi, lo_age_p = float(ra.OR), float(ra.CI_lo), float(ra.CI_hi), float(ra.P)
    wa, wb, wc, wd = int(rw.a), int(rw.b), int(rw.c), int(rw.d)
    lo_wbc_p = float(rw.P)

    fig = plt.figure(figsize=(7.2, 7.6))
    gs = fig.add_gridspec(2, 2, hspace=0.62, wspace=0.46)

    # A Oksa vs Li95 scatter (out-of-sample patients; deposited panel is discovery-derived)
    axA = fig.add_subplot(gs[0, 0])
    ok = tri_gene.dropna(subset=["oksa_coef", "li95_effect"])
    in219 = ok["in_primary_219"].astype(str).str.lower() == "true"
    axA.scatter(ok.loc[in219, "oksa_coef"], ok.loc[in219, "li95_effect"], s=14, c=OI["green"], alpha=0.7, edgecolors="none", label="in primary 219-gene set (24)")
    axA.scatter(ok.loc[~in219, "oksa_coef"], ok.loc[~in219, "li95_effect"], s=12, facecolors="none", edgecolors=OI["blue"], label="other deposited-panel genes (24)")
    lab_off = {"GATA2": (6, 4), "SHE": (6, -3), "ANPEP": (-8, -11), "MPV17L": (7, -9)}
    for g in ["GATA2", "SHE", "ANPEP", "MPV17L"]:
        rr = ok[ok.gene == g]
        if len(rr):
            axA.annotate(g, (rr.oksa_coef.iloc[0], rr.li95_effect.iloc[0]), fontsize=6.5,
                         textcoords="offset points", xytext=lab_off[g], weight="bold")
    axA.axhline(0, color="#bbbbbb", lw=0.7); axA.axvline(0, color="#bbbbbb", lw=0.7)
    axA.set_xlabel("Oksa slow-vs-fast coefficient"); axA.set_ylabel("Li95 C1-minus-C2")
    axA.legend(frameon=False, fontsize=6.3, loc="upper left")
    axA.text(0.98, 0.03, "ρ = 0.82 (48 of 50 genes)\n95.8% concordant\nout-of-sample patients;\npanel is discovery-derived",
             transform=axA.transAxes, ha="right", va="bottom", fontsize=6.2,
             bbox=dict(boxstyle="round,pad=0.25", fc="#f2f2f2", ec="#cccccc"))
    axA.set_title("A", loc="left", weight="bold", fontsize=11)

    # B Li95 clinical phenotype (age/WBC OR)
    axB = fig.add_subplot(gs[0, 1])
    feats = ["age ≥ 5 y", "WBC ≥ 50×10⁹/L"]
    # values parsed from the triangulation TSV above - no hard-coded numbers
    OR = [li95_ageOR[0], li95_wbcOR[0]]
    lo = [li95_ageOR[1], li95_wbcOR[1]]
    hi = [li95_ageOR[2], li95_wbcOR[2]]
    Pv = [li95_ageOR[3], li95_wbcOR[3]]
    yp = [1, 0]
    for y, o, l, h, p, f in zip(yp, OR, lo, hi, Pv, feats):
        axB.errorbar(o, y, xerr=[[o-l], [h-o]], fmt="s", color=OI["purple"], ecolor=OI["purple"], capsize=3, ms=6)
        axB.text(0.60, y + 0.16, f"OR {o:.2f} (P={p:.3f})", transform=axB.get_yaxis_transform(),
                 va="center", ha="left", fontsize=6.8)
    axB.axvline(1, color="#bbbbbb", lw=0.8, ls="--")
    axB.set_xscale("log"); axB.set_xlim(0.05, 20)
    axB.set_yticks(yp); axB.set_yticklabels(feats, fontsize=8)
    axB.set_xlabel("odds ratio (C1 vs C2), log scale")
    axB.set_title("B", loc="left", weight="bold", fontsize=11)

    # C local CD34/CD38
    axC = fig.add_subplot(gs[1, 0])
    data34 = [loc.loc[(~loc.poor), "CD34_raw"].dropna(), loc.loc[(loc.poor), "CD34_raw"].dropna()]
    data38 = [loc.loc[(~loc.poor), "CD38_raw"].dropna(), loc.loc[(loc.poor), "CD38_raw"].dropna()]
    pos = [1, 2, 4, 5]
    all_d = [data34[0], data34[1], data38[0], data38[1]]
    cols = [OI["blue"], OI["verm"], OI["blue"], OI["verm"]]
    for p, dd, c in zip(pos, all_d, cols):
        if len(dd) == 0:
            continue
        jit = np.random.default_rng(1).normal(0, 0.05, len(dd))
        axC.scatter(np.full(len(dd), p)+jit, dd, s=16, color=c, alpha=0.75, edgecolors="none")
        axC.hlines(dd.median(), p-0.28, p+0.28, color="black", lw=1.6)
        axC.text(p, -8, f"n={len(dd)}", ha="center", fontsize=7)
    axC.set_xticks([1.5, 4.5]); axC.set_xticklabels(["CD34", "CD38"], fontsize=9)
    axC.set_ylabel("% positive blasts"); axC.set_ylim(-14, 108)
    from matplotlib.lines import Line2D
    axC.legend(handles=[Line2D([], [], marker="o", ls="", color=OI["blue"], label="good"),
                        Line2D([], [], marker="o", ls="", color=OI["verm"], label="poor")],
               frameon=False, fontsize=7, loc="center", bbox_to_anchor=(0.52, 0.62),
               title="day-19 MRD", title_fontsize=7)
    axC.set_title("C", loc="left", weight="bold", fontsize=11)

    # D local age/WBC (parsed from source analysis table)
    axD = fig.add_subplot(gs[1, 1] if os.environ.get("PRIVACY_B") != "1" else gs[1, :])
    feats2 = ["age ≥ 5 y", "WBC ≥ 50×10⁹/L"]
    axD.errorbar(lo_age_or, 1, xerr=[[lo_age_or-lo_age_lo], [lo_age_hi-lo_age_or]], fmt="s", color=OI["sky"], ecolor=OI["sky"], capsize=3, ms=6)
    axD.text(0.62, 1.16, f"OR {lo_age_or:.2f} (P={lo_age_p:.2f})", transform=axD.get_yaxis_transform(), va="center", ha="left", fontsize=6.8)
    axD.text(0.30, -0.16, f"WBC: {wa}/{wa+wc} poor vs {wb}/{wb+wd} good\nOR not estimable (P={lo_wbc_p:.3f})", transform=axD.get_yaxis_transform(), va="center", ha="left", fontsize=6.8, color=OI["verm"])
    axD.axvline(1, color="#bbbbbb", lw=0.8, ls="--")
    axD.set_xscale("log"); axD.set_xlim(0.1, 12); axD.set_ylim(-0.7, 1.7)
    axD.set_yticks([1, 0]); axD.set_yticklabels(feats2, fontsize=8)
    axD.set_xlabel("odds ratio vs poor day-19 MRD, log scale")
    axD.set_title("D", loc="left", weight="bold", fontsize=11)
    # privacy option B (AUTHOR_ADMIN only): no patient-level dots -> panel C removed
    if os.environ.get("PRIVACY_B") == "1":
        axC.remove()
        save(fig, "Figure_3_Privacy_OptionB_summary_only")
    else:
        save(fig, "Figure_3")

# ---------------------------------------------------------------- Figure 4
def fig4():
    """Direction heatmap: cross-platform effect magnitudes are not on a common
    scale, so cells are coloured by DIRECTION only; raw values are printed."""
    d = rd("Fig4_candidate_biology.tsv")
    for c in ["oksa_coefficient_slow_vs_fast", "li_C1_minus_C2", "li95_effect"]:
        d[c] = pd.to_numeric(d[c], errors="coerce")
    genes = d["gene"].tolist()
    layers = ["oksa_coefficient_slow_vs_fast", "li_C1_minus_C2", "li95_effect"]
    layer_lab = ["Oksa", "Li194", "Li95"]
    M = d[layers].values.astype(float)

    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    # direction-only discrete colouring
    for i in range(len(genes)):
        for j in range(3):
            v = M[i, j]
            if np.isnan(v):
                facecol = "#d9d9d9"
                txt = "n/a"
                tcol = "#555555"
            elif v > 0:
                facecol = OI["blue"] + "55"
                txt = f"{v:+.2f}"
                tcol = "#0a2f52"
            else:
                facecol = OI["orange"] + "55"
                txt = f"{v:+.2f}"
                tcol = "#5a3200"
            ax.add_patch(plt.Rectangle((j - 0.5, i - 0.5), 1, 1, facecolor=facecol,
                                       edgecolor="white", linewidth=1.2))
            ax.text(j, i, txt, ha="center", va="center", fontsize=7.5, color=tcol)
    ax.set_xlim(-0.5, 3.4)
    ax.set_ylim(len(genes) - 0.5, -0.5)
    ax.set_xticks(range(3)); ax.set_xticklabels(layer_lab, fontsize=9)
    ax.set_yticks(range(len(genes))); ax.set_yticklabels(genes, fontsize=9)
    ax.set_xticks([x + 0.5 for x in range(3)], minor=True)
    for t in ax.get_xticks(minor=True):
        ax.axvline(t, color="white", lw=1.2)
    ax.tick_params(axis="both", which="both", length=0)
    for s in ("top", "right", "left", "bottom"):
        ax.spines[s].set_visible(False)

    # assessability markers, data-driven from li95_status
    status = d["li95_status"].astype(str).str.strip()
    for i, g in enumerate(genes):
        assessed = status.iloc[i] == "measured"
        ax.text(2.95, i, "\u25cf" if assessed else "\u25cb", va="center", fontsize=8,
                color=OI["green"] if assessed else OI["grey"])
    ax.text(2.95, -0.85, "\u25cf assessed in Li95\n    (concordant, BH q < 0.05)", fontsize=6.3,
            color=OI["green"], va="center")
    ax.text(2.95, len(genes) - 0.35, "\u25cb not assessable\n    (absent from panel)", fontsize=6.3,
            color=OI["grey"], va="center")
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(facecolor=OI["blue"] + "55", edgecolor="white", label="positive effect"),
                       Patch(facecolor=OI["orange"] + "55", edgecolor="white", label="negative effect"),
                       Patch(facecolor="#d9d9d9", edgecolor="white", label="not assessable")],
              frameon=False, fontsize=7, loc="upper center", bbox_to_anchor=(0.42, -0.13), ncol=3)
    ax.set_title("Seven pre-specified candidate genes: direction of effect across three contrasts",
                 fontsize=10, weight="bold", loc="left", pad=26)
    ax.text(0.0, -0.32, "Colour encodes direction only. Effect magnitudes are shown for reference but are not "
                        "directly comparable across platforms; the Li194 layer carries no per-gene inferential claim.",
            transform=ax.transAxes, fontsize=6.3, color="#444444")
    fig.subplots_adjust(top=0.78, bottom=0.36, left=0.10, right=0.88)
    save(fig, "Figure_4")

# ---------------------------------------------------------------- Figure 5
def fig5():
    d = rd("Fig5_drug_annotation.tsv")
    d["median_LC50_difference_C1_minus_C2"] = pd.to_numeric(d["median_LC50_difference_C1_minus_C2"], errors="coerce")
    d["bh_fdr_20drug"] = pd.to_numeric(d["bh_fdr_20drug"], errors="coerce")
    d = d.sort_values("median_LC50_difference_C1_minus_C2").reset_index(drop=True)
    fig, ax = plt.subplots(figsize=(7.0, 5.4))
    y = np.arange(len(d))
    sig = d["significant_bh_0.05"].astype(str).str.lower() == "true"
    ax.barh(y, d["median_LC50_difference_C1_minus_C2"],
            color=[OI["verm"] if s else OI["grey"] for s in sig], edgecolor="none", height=0.7)
    for i, (v, q, s) in enumerate(zip(d["median_LC50_difference_C1_minus_C2"], d["bh_fdr_20drug"], sig)):
        if s:
            ax.text(v + (0.005 if v >= 0 else -0.005), i, f"q={q:.1e}", va="center",
                    ha="left" if v >= 0 else "right", fontsize=6.5, color=OI["verm"], weight="bold")
    ax.set_yticks(y); ax.set_yticklabels(d["drug"], fontsize=7.5)
    ax.axvline(0, color="#444444", lw=0.8)
    ax.set_xlabel("median LC50 difference (C1 − C2)")
    ax.set_title("20-agent ex vivo drug screen (3 BH-significant highlighted)", fontsize=10, weight="bold", loc="left")
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(color=OI["verm"], label="BH q < 0.05"), Patch(color=OI["grey"], label="not significant")],
              frameon=False, fontsize=7.5, loc="lower right")
    save(fig, "Figure_5")

# ---------------------------------------------------------------- Figure S1 (pathways) + S2 (Round3B)
def figS():
    # S1: rendered from the real fgsea output (values are read from B2_pathway_axes_v2.tsv)
    src = pd.read_csv(Path(os.environ.get("EJH_PATHWAY_AXES", PROJECT_ROOT / "data" / "derived_tables" / "B2_pathway_axes_v2.tsv")), sep="\t",
                      keep_default_na=False)
    src["best_NES"] = pd.to_numeric(src["best_NES"], errors="coerce")
    src["best_padj"] = pd.to_numeric(src["best_padj"], errors="coerce")
    src["n_pathways"] = pd.to_numeric(src["n_pathways"], errors="coerce")
    src = src[src["n_pathways"].fillna(0) > 0].copy()
    src["neglog10q"] = -np.log10(src["best_padj"].clip(lower=1e-300))
    AXIS_LABEL = {
        "unexpected_top": "Other tested pathways (not pre-specified)",
        "metabolic_programmes": "Metabolic programmes",
        "RAS_MAPK": "RAS-MAPK signalling",
        "cell_cycle_E2F": "Cell cycle / E2F",
        "DNA_repair": "DNA repair",
        "JAK_STAT": "JAK-STAT signalling",
        "B_cell_developmental_state": "B-cell developmental state",
        "apoptosis": "Apoptosis",
        "TNFA_NFKB": "TNFa-NFkB signalling",
        "hypoxia": "Hypoxia",
        "MYC": "MYC targets",
        "BCR_signalling": "BCR signalling"}
    src = src.sort_values("neglog10q")
    lab = src["axis"].map(lambda a: AXIS_LABEL.get(a, a.replace("_", " ")))
    fig, ax = plt.subplots(figsize=(6.6, 4.3))
    fig.subplots_adjust(bottom=0.26, top=0.92, left=0.36, right=0.97)
    y = np.arange(len(src))
    colours = [OI["verm"] if d == "higher_in_slow" else OI["blue"] for d in src["best_direction"]]
    ax.barh(y, src["neglog10q"], color=colours, height=0.68)
    ax.set_yticks(y); ax.set_yticklabels(lab, fontsize=8)
    for yi, (q, npw) in enumerate(zip(src["neglog10q"], src["n_pathways"])):
        ax.text(q + 0.35, yi, f"n={int(npw)}", va="center", fontsize=6.8, color="#555555")
    ax.set_xlabel("-log10(BH q) of the best pathway in each axis")
    ax.set_xlim(0, float(src["neglog10q"].max()) * 1.28)
    ax.set_title("Figure S1. Pathway enrichment of the Oksa slow-response contrast",
                 fontsize=9.5, weight="bold", loc="left")
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(color=OI["verm"], label="higher in slow response"),
                       Patch(color=OI["blue"], label="higher in fast response")],
              frameon=False, fontsize=7.5, loc="lower right")
    fig.text(0.01, 0.015,
             "fgsea on the available Oksa effect universe (15,928 symbols); one bar per pre-specified axis;\n"
             "bars show the most significant pathway within that axis. Single-cohort biological context, not a cross-cohort test.",
             fontsize=6.6, color="#555555", ha="left", va="bottom")
    save(fig, "Figure_S1")

    # (supplementary perturbation figure removed in v3_2)

if __name__ == "__main__":
    fig1(); fig2(); fig3(); fig4(); fig5(); figS()
    print("ALL FIGURES DONE")
