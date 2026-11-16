"""Figures 1 and 2 of the manuscript from the aggregate result tables (no individual-level data).
Figure 1 (SPARK, SCQ): A within vs between family effects of the three exposures; B tetrachoric correlations of the
exposures; C sensitivity analyses of the within-family effects; D within-family effects of the variant-class tiers;
E within-family effects of the polygenic scores.
Figure 2 (SSC probands): A perinatal items x outcomes (M3); B genetic terms x outcomes (genetics-only model).
Usage: python 04_figures.py [results_root=results] [out_dir=figures]"""
import sys
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import numpy as np
import pandas as pd

ROOT = Path(sys.argv[1] if len(sys.argv) > 1 else "results"); OUT = Path(sys.argv[2] if len(sys.argv) > 2 else "figures"); OUT.mkdir(exist_ok=True, parents=True)
EXPO = ["nicu_oxygen", "premature", "weight_gain_diff"]
EL = {"nicu_oxygen": "NICU / oxygen at birth", "premature": "Premature birth", "weight_gain_diff": "Difficulty gaining weight"}
ES = {"nicu_oxygen": "NICU / oxygen", "premature": "Premature", "weight_gain_diff": "Weight gain"}
TIERS = ["dnm_lof_t1", "dnm_lof_t2", "dnm_mis_t1", "dnm_mis_t2", "rv_lof_t1", "rv_lof_t2", "rv_mis_t1", "rv_mis_t2", "dnm_lof_t3", "dnm_mis_t3", "dnm_syn", "rv_lof_t3", "rv_mis_t3"]
TL = {"dnm_lof_t1": "DNM LoF tier 1", "dnm_lof_t2": "DNM LoF tier 2", "dnm_mis_t1": "DNM missense tier 1", "dnm_mis_t2": "DNM missense tier 2",
      "rv_lof_t1": "Inherited LoF tier 1", "rv_lof_t2": "Inherited LoF tier 2", "rv_mis_t1": "Inherited missense tier 1", "rv_mis_t2": "Inherited missense tier 2",
      "dnm_lof_t3": "DNM LoF tier 3 (control)", "dnm_mis_t3": "DNM missense tier 3 (control)", "dnm_syn": "DNM synonymous (control)",
      "rv_lof_t3": "Inherited LoF tier 3 (control)", "rv_mis_t3": "Inherited missense tier 3 (control)"}
PGS = ["pgs", "pgs_gest_fetal", "pgs_pt_fetal"]; PL = {"pgs": "Autism PGS", "pgs_gest_fetal": "Gestational-duration PGS", "pgs_pt_fetal": "Preterm-birth PGS"}
SENS = [("main", "cont_eur", "Main analysis (EUR)"), ("main", "all_ancestry", "All ancestries"), ("male", "cont_eur", "Boys"), ("female", "cont_eur", "Girls"),
        ("no_twins", "cont_eur", "Multiple births excluded"), ("term_births", "cont_eur", "Term births only"), ("cases", "cont_eur", "Autistic children only"), ("controls", "cont_eur", "Non-autistic siblings only")]
C_W, C_B, C_DNM, C_RV, C_CTRL, C_PGS = "#1F5F8B", "#B45F06", "#1F5F8B", "#B45F06", "#8A8A8A", "#2E7D32"
pf = lambda p: "" if pd.isna(p) else (f"p = {p:.2f}" if p >= 0.01 else f"p = {p:.1e}")
def tcol(t): return C_CTRL if "t3" in t or t == "dnm_syn" else (C_DNM if t.startswith("dnm") else C_RV)
def style(ax):
    for s in ("top", "right"): ax.spines[s].set_visible(False)
    ax.axvline(0, color="grey", lw=0.8, zorder=0); ax.grid(axis="x", alpha=0.25)
def letter(ax, s): ax.text(-0.02, 1.04, s, transform=ax.transAxes, fontsize=14, fontweight="bold", ha="right", va="bottom")

meta = {sd: pd.read_csv(ROOT / sd / "fig1_meta.tsv", sep="\t") for sd in ["cont_eur", "all_ancestry"] if (ROOT / sd / "fig1_meta.tsv").exists()}
def m(sd, analysis, term, comp):
    d = meta[sd]; r = d[(d.analysis == analysis) & (d.term == term) & (d.component == comp)]; return r.iloc[0] if len(r) else None
tet = pd.read_csv(ROOT / "cont_eur" / "fig1_tetrachoric.tsv", sep="\t")

# ---------------- Figure 1 ----------------
fig = plt.figure(figsize=(18, 17)); gs = GridSpec(3, 6, figure=fig, height_ratios=[1, 1.05, 1.5], hspace=0.5, wspace=2.2)
def annot(ax, x, y, text, col, bold=False): ax.text(x, y, text, transform=ax.get_yaxis_transform(), va="center", ha="left", fontsize=8, color=col, fontweight="bold" if bold else "normal", clip_on=False)
# A: within vs between
ax = fig.add_subplot(gs[0, :3]); letter(ax, "A")
for i, X in enumerate(EXPO):
    for comp, off, col, mk in [("between", -0.2, C_B, "s"), ("within", 0.2, C_W, "o")]:
        r = m("cont_eur", "main", X, comp)
        ax.errorbar(r.est, i + off, xerr=1.96 * r.se, fmt=mk, color=col, mfc=col if r.pvalue < 0.05 else "white", ms=7, capsize=3, lw=1.4)
        annot(ax, 1.03, i + off, f"{r.est:+.3f} ({r['ci.lower']:+.3f}, {r['ci.upper']:+.3f})  {pf(r.pvalue)}", col)
    d = m("cont_eur", "main", X, "between_minus_within")
    annot(ax, 1.03, i + 0.5, f"between minus within: {d.est:+.3f}, {pf(d.pvalue)}", "#444444")
ax.set_yticks(range(3)); ax.set_yticklabels([EL[X] for X in EXPO], fontsize=9.5); ax.set_ylim(2.75, -0.6); style(ax); ax.set_xlim(-0.08, 0.3)
ax.set_xlabel("Effect on SCQ (SD), 95% CI", fontsize=9)
ax.set_title("Within- and between-family effects of perinatal exposures on SCQ\n(SPARK, hybrid mixed model with the three exposures entered jointly; three cohorts pooled)", fontsize=10, loc="left", fontweight="bold")
ax.legend(handles=[plt.Line2D([], [], marker="s", color=C_B, ls="", label="Between-family (family mean)"), plt.Line2D([], [], marker="o", color=C_W, ls="", label="Within-family (sibling contrast)")], loc="lower left", fontsize=8, frameon=False)
# B: tetrachoric heatmap
ax = fig.add_subplot(gs[0, 4:]); letter(ax, "B")
tp = tet[tet.cohort == "pooled"]; M = np.eye(3); S = [["" for _ in EXPO] for _ in EXPO]
for _, r in tp.iterrows():
    i, j = EXPO.index(r.x), EXPO.index(r.y); M[i, j] = M[j, i] = r.tetrachoric; S[i][j] = S[j][i] = f"{r.tetrachoric:.2f}\n(phi {r.phi:.2f})"
im = ax.imshow(M, cmap="Blues", vmin=0, vmax=1)
for i in range(3):
    for j in range(3): ax.text(j, i, "1" if i == j else S[i][j], ha="center", va="center", fontsize=9, color="white" if M[i, j] > 0.6 else "black")
ax.set_xticks(range(3)); ax.set_yticks(range(3)); ax.set_xticklabels([ES[X] for X in EXPO], fontsize=9, rotation=20, ha="right"); ax.set_yticklabels([ES[X] for X in EXPO], fontsize=9)
ax.set_title("Tetrachoric correlation of the exposures\n(children in the main analysis)", fontsize=10, loc="left", fontweight="bold")
plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04).set_label("tetrachoric r", fontsize=8)
# C: sensitivity (within-family), one panel per exposure
for k, X in enumerate(EXPO):
    ax = fig.add_subplot(gs[1, 2 * k:2 * k + 2]); labs = []; y = 0
    for a, sd, lab in SENS:
        if sd not in meta: continue
        r = m(sd, a, X, "within")
        if r is None: continue
        col = C_W if a == "main" and sd == "cont_eur" else "#5B8DB8"
        ax.errorbar(r.est, y, xerr=1.96 * r.se, fmt="o", color=col, mfc=col if r.pvalue < 0.05 else "white", ms=5.5, capsize=2.5, lw=1.2)
        annot(ax, 1.03, y, f"{r.est:+.3f} ({r.se:.3f})", col); labs.append(lab); y += 1
    ax.set_yticks(range(len(labs))); ax.set_yticklabels(labs if k == 0 else [], fontsize=8.5); ax.set_ylim(len(labs) - 0.5, -0.5); style(ax); ax.set_xlim(-0.45, 0.6)
    ax.set_xlabel("Within-family effect on SCQ (SD), 95% CI", fontsize=8.5); ax.set_title(EL[X], fontsize=10, loc="left", fontweight="bold")
    if k == 0: letter(ax, "C")
# D: tiers (within)
ax = fig.add_subplot(gs[2, :3]); letter(ax, "D"); rows = [t for t in TIERS if m("cont_eur", "main", t, "within") is not None]
for i, t in enumerate(rows):
    r = m("cont_eur", "main", t, "within"); col = tcol(t)
    ax.errorbar(r.est, i, xerr=1.96 * r.se, fmt="o", color=col, mfc=col if r.pvalue < 0.05 else "white", ms=6, capsize=3, lw=1.3)
    annot(ax, 1.03, i, f"{r.est:+.3f} ({r['ci.lower']:+.3f}, {r['ci.upper']:+.3f})  {pf(r.pvalue)}", col, r.pvalue < 0.05)
ax.set_yticks(range(len(rows))); ax.set_yticklabels([TL[t] for t in rows], fontsize=9); ax.set_ylim(len(rows) - 0.5, -0.5); ax.axhline(7.5, color="lightgrey", ls="--", lw=0.8); style(ax); ax.set_xlim(-0.2, 0.25)
ax.set_xlabel("Within-family effect on SCQ (SD) per variant, 95% CI", fontsize=9); ax.set_title("Within-family effects of rare-variant burden (same model as A)", fontsize=10, loc="left", fontweight="bold")
# E: PGS (within)
ax = fig.add_subplot(gs[2, 4:]); letter(ax, "E")
for i, t in enumerate(PGS):
    r = m("cont_eur", "main", t, "within")
    ax.errorbar(r.est, i, xerr=1.96 * r.se, fmt="o", color=C_PGS, mfc=C_PGS if r.pvalue < 0.05 else "white", ms=6, capsize=3, lw=1.3)
    annot(ax, 1.03, i, f"{r.est:+.3f} ({r['ci.lower']:+.3f}, {r['ci.upper']:+.3f})\n{pf(r.pvalue)}", C_PGS)
ax.set_yticks(range(len(PGS))); ax.set_yticklabels([PL[t] for t in PGS], fontsize=9); ax.set_ylim(len(PGS) - 0.5, -0.5); style(ax); ax.set_xlim(-0.04, 0.06)
ax.set_xlabel("Within-family effect on SCQ (SD)\nper SD of score, 95% CI", fontsize=9); ax.set_title("Polygenic scores\n(same model)", fontsize=10, loc="left", fontweight="bold")
fig.savefig(OUT / "figure1.png", dpi=200, bbox_inches="tight"); fig.savefig(OUT / "figure1.pdf", bbox_inches="tight"); plt.close(fig); print("wrote", OUT / "figure1.png")

# ---------------- Figure 2 ----------------
E2 = pd.read_csv(ROOT / "cont_eur" / "fig2_estimates.tsv", sep="\t"); G2 = pd.read_csv(ROOT / "cont_eur" / "fig2_genetics.tsv", sep="\t")
OUTS = ["scq", "srs", "vineland", "cbcl_int", "cbcl_ext"]; OL = {"scq": "SCQ", "srs": "SRS (parent)", "vineland": "Vineland-II composite", "cbcl_int": "CBCL 6-18 internalizing", "cbcl_ext": "CBCL 6-18 externalizing"}
IL = {"preg_preterm": "Preterm birth (< 37 wk)", "preg_gest_age_weeks": "Gestational age (per week)", "preg_birth_weight_sd": "Birth weight (per SD)", "preg_nicu": "NICU admission", "preg_oxygen": "Oxygen at birth",
      "preg_respiratory_distress": "Respiratory distress", "preg_hyperbilirubinemia": "Hyperbilirubinemia", "preg_apgar1_low": "1-min Apgar < 7", "preg_csection": "Caesarean section", "preg_labor_induced": "Induced labour"}
ITEMS = [i for i in IL if i in set(E2.item)]
fig, axes = plt.subplots(2, 5, figsize=(20, 11), gridspec_kw={"height_ratios": [10, 16], "hspace": 0.3, "wspace": 0.12})
for j, yv in enumerate(OUTS):
    ax = axes[0, j]; d = E2[(E2.outcome == yv) & (E2.model == "M3") & (E2.term == E2.item)].set_index("item")
    for i, it in enumerate(ITEMS):
        if it not in d.index: continue
        r = d.loc[it]; ax.errorbar(r.est, i, xerr=1.96 * r.se, fmt="o", color=C_B, mfc=C_B if r.p < 0.05 else "white", ms=5, capsize=2, lw=1.1)
    ax.set_yticks(range(len(ITEMS))); ax.set_yticklabels([IL[i] for i in ITEMS] if j == 0 else [], fontsize=8.5); ax.set_ylim(len(ITEMS) - 0.5, -0.5); style(ax); ax.set_xlim(-0.45, 0.45)
    ax.set_title(OL[yv], fontsize=10, fontweight="bold"); ax.set_xlabel("Effect (SD of outcome), 95% CI", fontsize=8)
    if j == 0: letter(ax, "A")
    n = d.n.max() if len(d) else 0; ax.text(0.98, 0.02, f"n up to {int(n):,}", transform=ax.transAxes, ha="right", fontsize=7.5, color="#444444")
for j, yv in enumerate(OUTS):
    ax = axes[1, j]; d = G2[G2.outcome == yv].set_index("term"); rows = [t for t in TIERS + PGS if t in d.index]
    for i, t in enumerate(rows):
        r = d.loc[t]; col = C_PGS if t in PGS else tcol(t)
        ax.errorbar(r.est, i, xerr=1.96 * r.se, fmt="o", color=col, mfc=col if r.p < 0.05 else "white", ms=5, capsize=2, lw=1.1)
    ax.set_yticks(range(len(rows))); ax.set_yticklabels([TL.get(t, PL.get(t)) for t in rows] if j == 0 else [], fontsize=8.5); ax.set_ylim(len(rows) - 0.5, -0.5)
    ax.axhline(7.5, color="lightgrey", ls="--", lw=0.8); ax.axhline(12.5, color="lightgrey", ls="--", lw=0.8); style(ax); ax.set_xlim(-0.45, 0.45)
    ax.set_xlabel("Effect (SD) per variant; PGS per SD", fontsize=8); ax.text(0.98, 0.02, f"n = {int(d.n.iloc[0]):,}", transform=ax.transAxes, ha="right", fontsize=7.5, color="#444444")
    if j == 0: letter(ax, "B")
fig.suptitle("SSC probands: perinatal items (A) and genetic burden (B) and five outcomes; population regressions with HC1 errors, adjusted for sex, birth year, PCs and coverage; A additionally adjusted for the genetic block (M3)", fontsize=10, y=0.93)
fig.savefig(OUT / "figure2.png", dpi=200, bbox_inches="tight"); fig.savefig(OUT / "figure2.pdf", bbox_inches="tight"); plt.close(fig); print("wrote", OUT / "figure2.png")
