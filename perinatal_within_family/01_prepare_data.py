#!/usr/bin/env python3
"""
Child-level input tables for the within-family analyses (one per cohort).

Rows: sequenced children (role == child) with case_status in {case, control}, restricted to
`--ancestry` (default: config analysis.sem.ancestry, EUR; "all" keeps every ancestry; cohorts whose table
carries the release `ancestry_call` are filtered on it, which leaves admixed / uncertain children out). Children
sequenced in two cohorts are kept once (config sem.overlap_exclusions). Every child is kept, not
only DNM-callset children: the DNM models restrict to `dnm_ok` (in the callset and passing the
zero / outlier filters of gxe.analysis.dnm_children), the other models use all children.

Columns (variable names used by 02_clogit_asd.R and 03_hybrid_scq.R):
  FID, IID         family and sample id
  asd              1 = case, 0 = control
  scq              SCQ total, z-scored within cohort (children)
  dnm_count, dnm   damaging de novo count (config sem.dnm / dnm_cohort_override) and its z-score
  dnm_ok           1 = usable in DNM models
  rv               inherited damaging count (config sem.rv_primary). Children outside the DNM callset
                   (no sequenced parents) have no parental filtering of their rare variants and, in the
                   2026_08 release, a heavy-tailed count (max 1,435 vs 14 in the callset), so the count is
                   standardized within cohort x dnm_ok group after dropping outliers above
                   median + 5 * 1.4826 * MAD of that group (rv_outlier = 1). Siblings share the group, so
                   the within-family contrasts are unaffected.
  pgs              ASD PGS (config sem.pgs_asd), z-scored, then residualized on the projected PCs
                   and the genotype-source indicators and re-standardized (as in the SEM's PGS
                   correction); pgs_raw = plain z-score
  pgs_gest, pgs_pt            child's own genotype scored with the EGG maternal gestational-duration /
                   preterm-delivery weights (gestdur_maternal_2023, preterm_maternal_2023); residualized like pgs
  pgs_gest_fetal, pgs_pt_fetal  child's genotype scored with the EGG fetal GWAS weights (Liu 2019; PRS-CS)
  pgs_gest_mother, pgs_pt_mother, pgs_gest_father   the parents' own scores (maternal weights), residualized on
                   the child's PCs; family constants, NA when the parent is not genotyped
  pgs_gest_nt      non-transmitted parental score = mother + father - child (Kong et al. 2018), maternal
                   gestational-duration weights; NA unless both parents are genotyped
  sex_male, father_age, mother_age (decades), birth_order, birth_year (decades, centred)
  mother_age_fam, father_age_fam, parent_age_diff_fam   family means of parental age (decades) and their difference
  premature        SPARK basic medical screening: 1 = "premature birth (< 37 weeks)" ticked,
                   0 = form completed and not ticked, NA = no form. gest_age: weeks if given.
                   nicu_oxygen (birth_oxygen) and weight_gain_diff (growth_low_wt) coded the same way.
                   SSC: NA (perinatal data exist for probands only, see preg_*).
  twin             SPARK multiple_birth: 1 = twin / triplet, 0 = singleton, NA = not recorded
  mother_alcohol   SPARK mother's ASR alcohol-user flag (family-level)
  mother_smoke, father_smoke   SSC parent ever-smoked flags (family-level)
  srs, srs_t, srs_teacher, vineland   SSC parent SRS raw total (z), SRS parent T score, teacher SRS raw
                   total (z), Vineland-II composite standard score (z); NA in SPARK
  preg_*           SSC proband-only maternal pregnancy / perinatal items from the Proband Data medhx
                   files (yes/no -> 1/0, not-sure -> NA; gestational age weeks, birth weight oz,
                   months breastfed numeric). Siblings have no values.
  pm25             prenatal PM2.5 mean, z-scored; pm25_enroll likewise; pm25_fam = family mean of
                   the prenatal value, z-scored (family-level moderator)
  adi              ADI national percentile (2019, else 2015), z-scored (SPARK only; NA in SSC)
  income, mother_edu, father_edu, ses   family-filled SES indicators, z-scored; ses = their mean
  coverage, mapped_pct                  child's sequencing covariates, z-scored
  src_*            genotype-source indicators (SPARK WES)

Output: results/within_family[/<subdir>]/input_<cohort>.tsv and input_summary.tsv
"""
import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

import numpy as np
import pandas as pd

from gxe import load_config, log
from gxe.analysis import dnm_children, load_cohort, results_path
import gxe.analysis as ga
import pyarrow.parquet as pq

COHORTS = ["spark_wes", "spark_wgs", "spark_wgs_2026_08", "ssc"]
# SSC proband medhx items: output name -> (file, column, kind). kind: yn (yes/no or true/false),
# num (first number in the string), csec (delivery_type == c-section), apgar (below-seven)
SSC_PREG = {
    "preg_smoke": ("medhx_mo_meds_substances.csv", "smoke", "yn"),
    "preg_alcohol": ("medhx_mo_meds_substances.csv", "drink", "yn"),
    "preg_antidepressants": ("medhx_mo_meds_substances.csv", "antidepressants_yes_no", "yn"),
    "preg_antibiotics": ("medhx_mo_meds_substances.csv", "antibiotics_yes_no", "yn"),
    "preg_painkillers": ("medhx_mo_meds_substances.csv", "pain_killers_yes_no", "yn"),
    "preg_folic_acid": ("medhx_preg_illneses_substances.csv", "folic_acid", "yn"),
    "preg_prenatal_vitamins": ("medhx_preg_illneses_substances.csv", "pre_prenatal_vits", "yn"),
    "preg_fever": ("medhx_preg_illneses_substances.csv", "mat_fever", "yn"),
    "preg_influenza": ("medhx_preg_illneses_substances.csv", "influenza", "yn"),
    "preg_uti": ("medhx_preg_illneses_substances.csv", "uti", "yn"),
    "preg_vaginal_infection": ("medhx_preg_illneses_substances.csv", "vaginal_infection", "yn"),
    "preg_gest_diabetes": ("medhx_mo_pregnancies.csv", "gestational_diabetes", "yn"),
    "preg_hypertension": ("medhx_mo_pregnancies.csv", "hypertension", "yn"),
    "preg_preeclampsia": ("medhx_mo_pregnancies.csv", "pre_eclampsia", "yn"),
    "preg_hyperemesis": ("medhx_mo_pregnancies.csv", "hyperemesis", "yn"),
    "preg_bedrest": ("medhx_mo_pregnancies.csv", "bedrest_2_months", "yn"),
    "preg_fertility_treatment": ("medhx_mo_pregnancies.csv", "fertility_treatment", "yn"),
    "preg_gest_age_weeks": ("medhx_perinatal.csv", "gestational_age_weeks", "num"),
    "preg_birth_weight_oz": ("medhx_perinatal.csv", "birth_weight_oz", "num"),
    "preg_csection": ("medhx_perinatal.csv", "delivery_type", "csec"),
    "preg_labor_induced": ("medhx_perinatal.csv", "labor_induced", "yn"),
    "preg_nicu": ("medhx_perinatal.csv", "nicu_admission", "yn"),
    "preg_oxygen": ("medhx_perinatal.csv", "o2_supplement", "yn"),
    "preg_hyperbilirubinemia": ("medhx_perinatal.csv", "hyperbilirubinemia", "yn"),
    "preg_respiratory_distress": ("medhx_perinatal.csv", "respiratory_distress", "yn"),
    "preg_apgar1_low": ("medhx_perinatal.csv", "apgar_score_at_1_minute", "apgar"),
    "preg_breastfed_months": ("medhx_perinatal.csv", "breast_total_months", "num"),
}
SES_COLS = {"income": "annual_household_income_ordinal_fam", "mother_edu": "mother_edu_code_fam",
            "father_edu": "father_edu_code_fam"}
# Gestational-duration / preterm PGS (EGG consortium; see config inputs.extra_pgs): output name -> assembled column.
GEST_COLS = {"pgs_gest": "gestdur_maternal_2023", "pgs_pt": "preterm_maternal_2023",
             "pgs_gest_fetal": "gestdur_fetal_2019", "pgs_pt_fetal": "preterm_fetal_2019",
             "pgs_gest_mother": "mother_gestdur_maternal_2023", "pgs_pt_mother": "mother_preterm_maternal_2023",
             "pgs_gest_father": "father_gestdur_maternal_2023"}


# Variant-class tiers for the split genetic block (raw counts; DNM prefix follows the cohort's dnm columns,
# e.g. dnm_snv_* for the 2026_08 release). LoF tiers by s_het (t1 >= 0.18, t2 0.03-0.18, t3 < 0.03);
# missense tiers by number of damaging predictors (t1 = 4, t2 = 3, t3 = 2); synonymous = negative control.
def tier_cols(dnm_cols):
    pre = dnm_cols[0].split("lof_t1")[0]
    return {"dnm_lof_t1": pre + "lof_t1", "dnm_lof_t2": pre + "lof_t2", "dnm_lof_t3": pre + "lof_t3",
            "dnm_mis_t1": pre + "miss_t1", "dnm_mis_t2": pre + "miss_t2", "dnm_mis_t3": pre + "miss_t3", "dnm_syn": pre + "synonymous"}
RV_TIERS = {"rv_lof_t1": "rv_lof_t1", "rv_lof_t2": "rv_lof_t2", "rv_lof_t3": "rv_lof_t3",
            "rv_mis_t1": "rv_miss_t1", "rv_mis_t2": "rv_miss_t2", "rv_mis_t3": "rv_miss_t3"}


def zscore(s):
    s = pd.to_numeric(s, errors="coerce").astype(float)
    return (s - s.mean()) / s.std()


def find_key(cfg, key):
    """First value of `key` anywhere in the nested config."""
    if isinstance(cfg, dict):
        if key in cfg:
            return cfg[key]
        for v in cfg.values():
            r = find_key(v, key)
            if r is not None:
                return r
    return None


def residualize(y, X):
    """OLS residuals of y on X (with intercept), NaN where y or X is missing; re-standardized."""
    X = X.copy()
    X.insert(0, "const", 1.0)
    ok = y.notna() & X.notna().all(axis=1)
    beta, *_ = np.linalg.lstsq(X[ok].to_numpy(float), y[ok].to_numpy(float), rcond=None)
    r = pd.Series(np.nan, index=y.index)
    r[ok] = y[ok] - X[ok].to_numpy(float) @ beta
    return (r - r.mean()) / r.std()


def spark_premature(cfg):
    """SP id -> perinatal flags (premature, nicu_oxygen, weight_gain_diff; 0 = form completed, not ticked)
    and gestational age from the basic medical screening form."""
    rel = find_key(cfg, "release_dir")
    tag = find_key(cfg, "release_tag")
    bms = pd.read_csv(f"{rel}/basic_medical_screening-{tag}.csv", dtype=str, low_memory=False,
                      usecols=["subject_sp_id", "birth_prem", "gest_age", "birth_oxygen", "growth_low_wt"]).drop_duplicates("subject_sp_id")
    return pd.DataFrame({"premature": (bms["birth_prem"] == "1").astype(float).values,
                         "nicu_oxygen": (bms["birth_oxygen"] == "1").astype(float).values,
                         "weight_gain_diff": (bms["growth_low_wt"] == "1").astype(float).values,
                         "gest_age": pd.to_numeric(bms["gest_age"], errors="coerce").values}, index=bms["subject_sp_id"].values)


def ssc_pregnancy(cfg):
    """SSC individual id (e.g. 11000.p1) -> proband-only pregnancy / perinatal items (SSC_PREG)."""
    d = Path(find_key(cfg, "pheno_dir")) / "Proband Data"
    files, out = {}, {}
    for name, (f, col, kind) in SSC_PREG.items():
        if f not in files:
            files[f] = pd.read_csv(d / f, dtype=str, low_memory=False).drop_duplicates("individual").set_index("individual")
        x = files[f][col].astype("string").str.strip().str.lower()
        if kind == "yn":
            v = x.map({"yes": 1.0, "no": 0.0, "true": 1.0, "false": 0.0})
        elif kind == "csec":
            v = x.map({"c-section": 1.0, "vaginally": 0.0})
        elif kind == "apgar":
            v = x.map({"below-seven": 1.0, "normal": 0.0})
        else:
            v = pd.to_numeric(x.str.extract(r"([0-9]+\.?[0-9]*)")[0], errors="coerce")
        out[name] = v.astype(float)
    return pd.DataFrame(out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ancestry", default=None, help="predicted_ancestry to keep (default: config sem.ancestry); 'all' keeps every ancestry")
    ap.add_argument("--cohorts", nargs="+", default=COHORTS)
    ap.add_argument("--subdir", default="", help="results subfolder under results/within_family (e.g. all_ancestry)")
    args = ap.parse_args()
    cfg = load_config()
    sem = cfg["analysis"]["sem"]
    ancestry = args.ancestry or sem["ancestry"]
    pcs = [f"proj_PC{i}" for i in range(1, sem["n_pcs"] + 1)]
    rv_cols = sem["rv_definitions"][sem["rv_primary"]]
    drop_ids = {}
    for rule in sem.get("overlap_exclusions", []):
        drop_ids.setdefault(rule["drop"], set()).update(load_cohort(rule["keep"], ["sample_id"])["sample_id"])
    bms = spark_premature(cfg)
    preg = ssc_pregnancy(cfg)
    log(f"  SSC pregnancy items: {len(preg):,} probands, {int(preg.notna().sum().min())}-{int(preg.notna().sum().max())} non-missing per item")
    out_dir = str(Path("within_family") / args.subdir) if args.subdir else "within_family"

    summary = []
    for cohort in args.cohorts:
        dnm_cols = sem.get("dnm_cohort_override", {}).get(cohort, sem["dnm"])
        spark = cohort.startswith("spark")
        need = ["sample_id", "FID", "source", "role", "case_status", "SEX", "predicted_ancestry", "birth_order",
                "birth_year_dec", "father_age_at_child_birth", "mother_age_at_child_birth", "scq_total_final_score",
                "prenatal_pm25_mean", "enroll_pm25_mean", "coverage_mean", "cram_mapped_pct", "in_dnm_callset",
                *SES_COLS.values(), *dnm_cols, *rv_cols, sem["pgs_asd"], *pcs, *GEST_COLS.values(),
                *tier_cols(dnm_cols).values(), *RV_TIERS.values()]
        if spark:
            need += ["adi_2019_national_pct", "adi_2015_national_pct", "multiple_birth", "mother_asr_alcohol_user"]
        else:
            need += ["ssc_individual", "srs_parent_raw_total", "srs_parent_t_score", "srs_teacher_raw_total",
                     "vineland_ii_composite_standard_score", "mother_ever_smoke", "father_ever_smoke",
                     "cbcl_6_18_internalizing_t_score", "cbcl_6_18_externalizing_t_score"]
        need = list(dict.fromkeys(need))   # tier columns overlap the summed dnm/rv columns
        df = load_cohort(cohort, need)
        df = df[(df["role"] == "child") & df["case_status"].isin(["case", "control"])]
        if ancestry != "all":
            anc_col = "predicted_ancestry"
            if "ancestry_call" in pq.read_schema(ga.assembled_path(f"gxe_{cohort}.parquet")).names:
                # release cohorts carry the release's own call, which labels admixed / uncertain children
                # separately instead of assigning them to the most likely group
                anc_col = "ancestry_call"
                df = load_cohort(cohort, list(dict.fromkeys(need + ["ancestry_call"])))
                df = df[(df["role"] == "child") & df["case_status"].isin(["case", "control"])]
            log(f"  {cohort}: keeping {anc_col} == {ancestry} ({int((df[anc_col] == ancestry).sum()):,} of {len(df):,} children)")
            df = df[df[anc_col] == ancestry]
        if cohort in drop_ids:
            overlap = df["sample_id"].isin(drop_ids[cohort])
            log(f"  {cohort}: dropping {int(overlap.sum()):,} children also in another cohort")
            df = df[~overlap]
        df = df.reset_index(drop=True)
        ok_ids = set(dnm_children(cohort, columns=["sample_id"], cfg=cfg)[0]["sample_id"])
        log(f"  {cohort}: DNM = {' + '.join(dnm_cols)}; RV = {' + '.join(rv_cols)}")

        out = pd.DataFrame({
            "FID": df["FID"].astype(str), "IID": df["sample_id"].astype(str),
            "asd": (df["case_status"] == "case").astype(int),
            "scq": zscore(df["scq_total_final_score"]),
            "dnm_count": df[dnm_cols].sum(axis=1, min_count=1),
        })
        out["dnm"] = zscore(out["dnm_count"])
        out["dnm_ok"] = df["sample_id"].isin(ok_ids).astype(int)
        out.loc[out["dnm_ok"] == 0, ["dnm_count", "dnm"]] = np.nan
        rv_count = df[rv_cols].sum(axis=1).astype(float)
        out["rv_count"] = rv_count
        out["rv"] = np.nan
        out["rv_outlier"] = 0
        for grp, idx in out.groupby("dnm_ok").groups.items():
            x = rv_count.loc[idx]
            cut = x.median() + 5 * 1.4826 * (x - x.median()).abs().median()
            keep = x <= cut
            out.loc[idx[~keep], "rv_outlier"] = 1
            out.loc[idx[keep], "rv"] = ((x[keep] - x[keep].mean()) / x[keep].std()).values
            log(f"  {cohort}: RV in dnm_ok={grp}: n {len(x):,}, mean {x.mean():.2f}, sd {x.std():.2f}, cutoff {cut:.0f}, outliers dropped {int((~keep).sum()):,}")
        for name, col in tier_cols(dnm_cols).items():
            out[name] = pd.to_numeric(df[col], errors="coerce").astype(float)
            out.loc[out["dnm_ok"] == 0, name] = np.nan
        for name, col in RV_TIERS.items():
            out[name] = pd.to_numeric(df[col], errors="coerce").astype(float)
            out.loc[out["rv_outlier"] == 1, name] = np.nan
        out["pgs_raw"] = zscore(df[sem["pgs_asd"]])
        out["sex_male"] = pd.to_numeric(df["SEX"], errors="coerce").map({1: 1, 2: 0})
        out["father_age"] = df["father_age_at_child_birth"] / 10
        out["mother_age"] = df["mother_age_at_child_birth"] / 10
        # family-level parental age (decades): family mean over the analysis children. Within a family
        # parental age is the same contrast as birth year, so only the family-level part is a new variable.
        out["mother_age_fam"] = out["mother_age"].groupby(out["FID"]).transform("mean")
        out["father_age_fam"] = out["father_age"].groupby(out["FID"]).transform("mean")
        out["parent_age_diff_fam"] = out["father_age_fam"] - out["mother_age_fam"]
        out["birth_order"] = pd.to_numeric(df["birth_order"], errors="coerce")
        by = pd.to_numeric(df["birth_year_dec"], errors="coerce")
        out["birth_year"] = (by - by.mean()) / 10
        for c in ["premature", "nicu_oxygen", "weight_gain_diff", "gest_age", "twin", "mother_alcohol", "adi",
                  "mother_smoke", "father_smoke", "srs", "srs_t", "srs_teacher", "vineland", "cbcl_int", "cbcl_ext", *SSC_PREG]:
            out[c] = np.nan
        if spark:
            m = bms.reindex(df["sample_id"].astype(str).values)
            for c in ["premature", "nicu_oxygen", "weight_gain_diff", "gest_age"]:
                out[c] = m[c].values
            out["twin"] = df["multiple_birth"].map({"No": 0.0, "Twin": 1.0, "Triplet": 1.0}).astype(float)
            out["mother_alcohol"] = pd.to_numeric(df["mother_asr_alcohol_user"], errors="coerce").astype(float)
            out["adi"] = zscore(df["adi_2019_national_pct"].fillna(df["adi_2015_national_pct"]))
        else:
            out["mother_smoke"] = pd.to_numeric(df["mother_ever_smoke"], errors="coerce").astype(float)
            out["father_smoke"] = pd.to_numeric(df["father_ever_smoke"], errors="coerce").astype(float)
            out["srs"] = zscore(df["srs_parent_raw_total"])
            out["srs_t"] = pd.to_numeric(df["srs_parent_t_score"], errors="coerce").astype(float)
            out["srs_teacher"] = zscore(df["srs_teacher_raw_total"])
            out["vineland"] = zscore(df["vineland_ii_composite_standard_score"])
            out["cbcl_int"] = zscore(df["cbcl_6_18_internalizing_t_score"])   # SSC probands and siblings, ages 6-18
            out["cbcl_ext"] = zscore(df["cbcl_6_18_externalizing_t_score"])
            m = preg.reindex(df["ssc_individual"].astype("string").values)
            for c in SSC_PREG:
                out[c] = m[c].values
        out["pm25"] = zscore(df["prenatal_pm25_mean"])
        out["pm25_enroll"] = zscore(df["enroll_pm25_mean"])
        fam_mean = pd.to_numeric(df["prenatal_pm25_mean"], errors="coerce").groupby(out["FID"]).transform("mean")
        out["pm25_fam"] = zscore(fam_mean)
        for c, col in SES_COLS.items():
            out[c] = zscore(df[col])
        out["ses"] = out[list(SES_COLS)].mean(axis=1)
        out["coverage"] = zscore(df["coverage_mean"])
        out["mapped_pct"] = zscore(df["cram_mapped_pct"])
        pcz = pd.DataFrame({f"PC{i}": zscore(df[c]) for i, c in enumerate(pcs, 1)})
        src_cols = []
        for src in df["source"].value_counts().index[1:]:
            out[f"src_{src}"] = (df["source"] == src).astype(int)
            src_cols.append(f"src_{src}")
        out["pgs"] = residualize(out["pgs_raw"], pd.concat([pcz, out[src_cols]], axis=1))
        covs = pd.concat([pcz, out[src_cols]], axis=1)
        for name, col in GEST_COLS.items():
            out[name] = residualize(zscore(df[col]), covs) if df[col].notna().sum() > 100 else np.nan
        # non-transmitted parental score: the assembled columns are z-scores of the same raw score within cohort,
        # so mother + father - child is on one scale; residualized on the child's PCs like the others
        nt = pd.to_numeric(df["mother_gestdur_maternal_2023"], errors="coerce") + pd.to_numeric(df["father_gestdur_maternal_2023"], errors="coerce") \
            - pd.to_numeric(df["gestdur_maternal_2023"], errors="coerce")
        out["pgs_gest_nt"] = residualize(zscore(nt), covs) if nt.notna().sum() > 100 else np.nan
        log(f"  {cohort}: gestational PGS present: child {int(out['pgs_gest'].notna().sum()):,}, mother {int(out['pgs_gest_mother'].notna().sum()):,}, "
            f"both parents (non-transmitted) {int(out['pgs_gest_nt'].notna().sum()):,}")
        out = pd.concat([out, pcz], axis=1)

        path = results_path(out_dir, f"input_{cohort}.tsv")
        out.to_csv(path, sep="\t", index=False, na_rep="NA")
        fam = out.groupby("FID")["asd"].agg(["min", "max", "size"])
        disc = fam[(fam["min"] == 0) & (fam["max"] == 1)]
        log(f"{cohort}: {len(out):,} children ({int(out['asd'].sum()):,} cases), {len(fam):,} families, "
            f"{len(disc):,} discordant, {int(out['dnm_ok'].sum()):,} DNM-usable -> {path}")
        summary.append({"cohort": cohort, "ancestry": ancestry, "n": len(out), "cases": int(out["asd"].sum()),
                        "families": len(fam), "discordant_families": len(disc), "dnm_ok": int(out["dnm_ok"].sum()),
                        **{f"{c}_n": int(out[c].notna().sum()) for c in ["scq", "srs", "vineland", "premature", "nicu_oxygen", "twin", "pm25", "pm25_enroll", "adi", "ses", "mother_alcohol", "mother_smoke", "preg_smoke", "mother_age_fam", "father_age_fam"]}})
    s = pd.DataFrame(summary)
    s.to_csv(results_path(out_dir, "input_summary.tsv"), sep="\t", index=False)
    log("\n" + s.to_string(index=False))


if __name__ == "__main__":
    main()
