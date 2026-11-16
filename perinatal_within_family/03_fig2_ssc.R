#!/usr/bin/env Rscript
# Figure 2 models: SSC probands. The SSC medical-history interview (perinatal items) exists for probands only,
# so exposure effects are population (one proband per family) OLS with HC1 standard errors:
#   M0  outcome ~ item + sex + birth year + 10 PCs + coverage
#   M3  M0 + the genetic block (variant-class tiers, raw counts; autism, fetal gestational-duration and fetal
#       preterm PGS per SD), on the same complete-case sample as M0
# plus one genetics-only model per outcome (M3 without an item) for the genetic-effects panel.
# Items: preterm (< 37 weeks), gestational age (weeks), birth weight (SD), NICU, oxygen, respiratory distress,
# hyperbilirubinemia, 1-min Apgar < 7, Caesarean section, induced labour (binary items need >= 40 exposed).
# Outcomes (z within cohort): SCQ, SRS parent raw total, Vineland-II composite, CBCL 6-18 internalizing / externalizing.
# Inputs : <gxe_repo>/results/within_family/<subdir>/input_ssc.tsv
# Outputs: results/<subdir>/fig2_{estimates,genetics,samples}.tsv
# Usage  : Rscript 03_fig2_ssc.R <gxe_repo> [subdir=cont_eur] [out_root=results]
suppressMessages({ library(sandwich); library(data.table) })
args <- commandArgs(trailingOnly = TRUE)
gxe <- args[1]; subdir <- if (length(args) >= 2) args[2] else "cont_eur"; out_root <- if (length(args) >= 3) args[3] else "results"
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))), "wf_common.R"))
in_dir <- file.path(gxe, "results", "within_family", subdir); out_dir <- file.path(out_root, subdir); dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
GEN <- c("dnm_lof_t1", "dnm_lof_t2", "dnm_mis_t1", "dnm_mis_t2", "rv_lof_t1", "rv_lof_t2", "rv_mis_t1", "rv_mis_t2",
         "dnm_lof_t3", "dnm_mis_t3", "dnm_syn", "rv_lof_t3", "rv_mis_t3", "pgs", "pgs_gest_fetal", "pgs_pt_fetal")
BASE <- c("sex_male", "birth_year"); OUT <- c("scq", "srs", "vineland", "cbcl_int", "cbcl_ext")
s <- read_input(in_dir, "ssc")[dnm_ok == 1 & asd == 1]
s[, preg_preterm := as.numeric(preg_gest_age_weeks < 37)]
s[, preg_birth_weight_sd := (preg_birth_weight_oz - mean(preg_birth_weight_oz, na.rm = TRUE)) / sd(preg_birth_weight_oz, na.rm = TRUE)]
ITEMS <- c("preg_preterm", "preg_gest_age_weeks", "preg_birth_weight_sd", "preg_nicu", "preg_oxygen", "preg_respiratory_distress",
           "preg_hyperbilirubinemia", "preg_apgar1_low", "preg_csection", "preg_labor_induced")
cs <- c(BASE, grep("^PC[0-9]+$", names(s), value = TRUE), "coverage"); cs <- cs[cs %in% names(s)]
gen <- GEN[sapply(GEN, function(v) v %in% names(s) && any(!is.na(s[[v]])))]
ols <- function(d, y, xs, cs) {
  keep <- Filter(function(v) sd(d[[v]]) > 0, c(xs, cs))
  fit <- lm(as.formula(paste(y, "~", paste(keep, collapse = " + "))), data = d); V <- vcovHC(fit, type = "HC1"); b <- coef(fit)
  rbindlist(lapply(xs, function(k) if (k %in% names(b) && !is.na(b[k])) data.table(term = k, est = unname(b[k]), se = sqrt(V[k, k]), p = 2 * pnorm(-abs(b[k] / sqrt(V[k, k])))) else NULL))
}
est <- list(); genr <- list(); samp <- list()
for (y in OUT) {
  for (it in ITEMS) {
    vars <- unique(c(y, it, gen, cs)); d <- s[complete.cases(s[, ..vars])]
    n_exp <- if (uniqueN(d[[it]]) == 2) sum(d[[it]] == 1) else NA_integer_
    if (nrow(d) < 200 || uniqueN(d[[it]]) < 2 || (is.finite(n_exp) && n_exp < 40)) next
    r0 <- ols(d, y, it, cs)[, model := "M0"]; r3 <- ols(d, y, c(it, gen), cs)[, model := "M3"]
    est[[length(est) + 1]] <- rbind(r0, r3)[, `:=`(outcome = y, item = it, n = nrow(d), n_exposed = n_exp)]
    samp[[length(samp) + 1]] <- data.table(outcome = y, item = it, n = nrow(d), n_exposed = n_exp)
  }
  vars <- unique(c(y, gen, cs)); d <- s[complete.cases(s[, ..vars])]
  genr[[length(genr) + 1]] <- ols(d, y, gen, cs)[, `:=`(outcome = y, n = nrow(d))]
  cat(sprintf("== %s: genetics-only model n %d\n", y, nrow(d)))
}
est <- rbindlist(est); genr <- rbindlist(genr); samp <- rbindlist(samp)
fwrite(est, file.path(out_dir, "fig2_estimates.tsv"), sep = "\t"); fwrite(genr, file.path(out_dir, "fig2_genetics.tsv"), sep = "\t"); fwrite(samp, file.path(out_dir, "fig2_samples.tsv"), sep = "\t")
print(est[model == "M3" & term == item, .(outcome, item, est = round(est, 3), se = round(se, 3), p = signif(p, 2), n, n_exposed)])
print(genr[term %in% c("dnm_lof_t1", "dnm_lof_t2", "pgs"), .(outcome, term, est = round(est, 3), se = round(se, 3), p = signif(p, 2))])
cat("Wrote", out_dir, "\n")
