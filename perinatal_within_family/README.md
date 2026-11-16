# Within-family effects of perinatal adversity on autism traits (SPARK, SSC)

Code for the manuscript on within- and between-family effects of NICU / oxygen at birth, premature birth and
difficulty gaining weight on the SCQ (SPARK), with the child's genetic burden (variant-class tiers of de novo and
inherited rare variants; autism, gestational-duration and preterm-birth polygenic scores) in the model, and the
proband-only perinatal items of SSC against SCQ, SRS, Vineland-II and CBCL 6-18.

## Pipeline

| Step | Script | What it does |
|---|---|---|
| 1 | `01_prepare_data.py` | Child-level tables per cohort (copy of `GxE_analysis/02_analysis/within_family/01_prepare_data.py`; run from that repo with `PYTHONPATH` set to it): outcomes, exposures, DNM / RV tier counts, residualised PGS, covariates. `--subdir cont_eur` (EUR, primary) and `--ancestry all --subdir all_ancestry`. |
| 2 | `02_fig1_models.R` | SPARK hybrid (Mundlak) mixed models: three exposures jointly, within / between decomposition, CR2 family-clustered SEs, IVW meta over cohorts; M0 (no genetic terms) and sensitivity analyses (sex, multiple births excluded, term births, cases / controls); tetrachoric correlations. |
| 3 | `03_fig2_ssc.R` | SSC probands: OLS + HC1 for ten perinatal items x five outcomes (M0, M3) and genetics-only models. |
| 4 | `04_figures.py` | Figures 1 and 2 from the aggregate tables. |
| 5 | `05_tables.py` | Table 1 and supplementary tables S1-S3. |
| | `run_models.sh` | SLURM launcher for steps 2-3 (EUR, all ancestries, SSC). |
| | `wf_common.R` | Shared helpers (`read_input`, `ivw_meta`). |

Inputs are the prepared tables in `GxE_analysis/results/within_family/{cont_eur,all_ancestry}/input_<cohort>.tsv`
(individual-level; not in this repository). `results/` holds only aggregate coefficient tables; `figures/` and
`tables/` the manuscript figures and tables.

## Model

SCQ_ij = sum_k [bW_k (x_kij - xbar_kj) + bB_k xbar_kj] + g' z_j + u_j + e_ij

x_k: NICU / oxygen, premature birth, difficulty gaining weight, genetic terms, sex, birth year, ASD;
z_j: 10 PCs, source / array, coverage. bW = within-family (sibling contrast), bB = between-family (family mean).

## Environment

R (micromamba env `sem`): lme4, clubSandwich, sandwich, data.table, pbivnorm. Python: pandas, matplotlib.
