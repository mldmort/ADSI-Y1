"""Manuscript tables from the aggregate result tables.
Table 1: analysed samples per cohort (children, families, multi-child families, autistic children, exposed counts).
Table S1: all Figure 1 coefficients (within, between, difference) for the main model and every sensitivity analysis.
Table S2: Figure 2 coefficients (M0 and M3) and the SSC genetics-only models.
Table S3: tetrachoric / phi correlations and co-occurrence counts.
Usage: python 05_tables.py [results_root=results] [out_dir=tables]"""
import sys
from pathlib import Path
import pandas as pd

ROOT = Path(sys.argv[1] if len(sys.argv) > 1 else "results"); OUT = Path(sys.argv[2] if len(sys.argv) > 2 else "tables"); OUT.mkdir(exist_ok=True, parents=True)
COH = {"spark_wes": "SPARK WES", "spark_wgs": "SPARK WGS (2023)", "spark_wgs_2026_08": "SPARK WGS (2026_08)"}
s = pd.read_csv(ROOT / "cont_eur" / "fig1_samples.tsv", sep="\t"); s = s[s.analysis == "main"]
t1 = pd.DataFrame({"Cohort": [COH[c] for c in s.cohort], "Children analysed": s.n, "Families": s.families, "Families with >= 2 analysed children": s.families_multi,
                   "Children in those families": s.children_in_multi, "Autistic children": s.cases, "NICU / oxygen": s.n_nicu, "Premature birth": s.n_premature, "Difficulty gaining weight": s.n_weight})
tot = t1.drop(columns="Cohort").sum(); tot["Cohort"] = "Total"; t1 = pd.concat([t1, tot.to_frame().T[t1.columns]], ignore_index=True)
if (ROOT / "all_ancestry" / "fig1_samples.tsv").exists():
    sa = pd.read_csv(ROOT / "all_ancestry" / "fig1_samples.tsv", sep="\t"); sa = sa[sa.analysis == "main"]
    t1["Children analysed, all ancestries"] = list(sa.n) + [sa.n.sum()]
t1.to_csv(OUT / "table1_samples.csv", index=False)
rows = []
for sd in ["cont_eur", "all_ancestry"]:
    p = ROOT / sd / "fig1_meta.tsv"
    if not p.exists(): continue
    m = pd.read_csv(p, sep="\t"); m["ancestry"] = "EUR" if sd == "cont_eur" else "all"
    rows.append(m)
m = pd.concat(rows); m["estimate (95% CI)"] = m.apply(lambda r: f"{r.est:+.3f} ({r['ci.lower']:+.3f}, {r['ci.upper']:+.3f})", axis=1)
m.pivot_table(index=["ancestry", "analysis", "term"], columns="component", values="estimate (95% CI)", aggfunc="first").reset_index().to_csv(OUT / "tableS1_fig1_all_models.csv", index=False)
m.to_csv(OUT / "tableS1_fig1_all_models_long.csv", index=False)
e2 = pd.read_csv(ROOT / "cont_eur" / "fig2_estimates.tsv", sep="\t"); e2[e2.term == e2.item].to_csv(OUT / "tableS2_fig2_items.csv", index=False)
pd.read_csv(ROOT / "cont_eur" / "fig2_genetics.tsv", sep="\t").to_csv(OUT / "tableS2_fig2_genetics.csv", index=False)
pd.read_csv(ROOT / "cont_eur" / "fig1_tetrachoric.tsv", sep="\t").to_csv(OUT / "tableS3_tetrachoric.csv", index=False)
print(t1.to_string(index=False)); print("wrote", OUT)
