#!/usr/bin/env Rscript
# Figure 1 models: within- and between-family effects of three perinatal exposures (NICU / oxygen at birth,
# premature birth, difficulty gaining weight; SPARK basic medical screening) on SCQ, with the genetic block
# split into variant-class tiers, in one hybrid (Mundlak) mixed model per cohort:
#
#   scq_ij = sum_k [bW_k (x_kij - xbar_kj) + bB_k xbar_kj] + g' z_j + u_j + e_ij,   u_j ~ N(0, tau^2)
#
#   x_k  the three exposures (entered jointly), the genetic terms, sex, birth year, ASD status
#   z_j  10 projected PCs, array / source indicators, sequencing coverage
#   bW   within-family effect (sibling contrast; = family fixed-effects estimate); bB between-family effect
#   SEs  CR2 family-clustered (clubSandwich) on the lmer fit
#
# Genetic block (raw counts, per variant): DNM LoF tier 1 (s_het >= 0.18) / tier 2 (0.03-0.18), DNM missense
# tier 1 (4 of 4 predictors) / tier 2 (3 of 4), the same four inherited tiers, negative controls (DNM LoF
# tier 3, DNM missense tier 3, synonymous DNMs, inherited LoF tier 3, inherited missense tier 3), and three
# polygenic scores per SD (autism; gestational duration and preterm birth, fetal EGG GWAS weights).
# The 2026_08 release has SNV-only DNM tiers and no synonymous count (dropped for that cohort).
#
# Analyses (all children with SCQ, DNM-usable, complete on every model variable; twin status needed only for no_twins):
#   main            pooled, full model                      M0: same without the genetic terms
#   male / female   stratified by child sex
#   no_twins        multiple births excluded
#   term_births     premature == 0 (premature term dropped)
#   cases / controls  affected / unaffected children only
# Ancestry sensitivity: run with subdir = all_ancestry.
# Also: tetrachoric and phi correlations among the three exposures (ML, pbivnorm) and co-occurrence counts.
#
# Inputs : <gxe_repo>/results/within_family/<subdir>/input_<cohort>.tsv   (01_prepare_data.py)
# Outputs: results/<subdir>/fig1_{estimates,meta,samples,tetrachoric}.tsv
# Usage  : Rscript 02_fig1_models.R <gxe_repo> [subdir=cont_eur] [out_root=results]
suppressMessages({ library(lme4); library(data.table); library(pbivnorm) })
args <- commandArgs(trailingOnly = TRUE)
gxe <- args[1]; subdir <- if (length(args) >= 2) args[2] else "cont_eur"; out_root <- if (length(args) >= 3) args[3] else "results"
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))), "wf_common.R"))
in_dir <- file.path(gxe, "results", "within_family", subdir); out_dir <- file.path(out_root, subdir); dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

EXPO <- c("nicu_oxygen", "premature", "weight_gain_diff")
GEN_MAIN <- c("dnm_lof_t1", "dnm_lof_t2", "dnm_mis_t1", "dnm_mis_t2", "rv_lof_t1", "rv_lof_t2", "rv_mis_t1", "rv_mis_t2")
GEN_CTRL <- c("dnm_lof_t3", "dnm_mis_t3", "dnm_syn", "rv_lof_t3", "rv_mis_t3")
PGS <- c("pgs", "pgs_gest_fetal", "pgs_pt_fetal")
BASE <- c("sex_male", "birth_year", "asd")
COHORTS <- c("spark_wes", "spark_wgs", "spark_wgs_2026_08")

fit_hybrid <- function(d, xs, zs) {
  for (x in xs) { d[, paste0(x, "_b") := mean(get(x)), by = FID]; d[, paste0(x, "_w") := get(x) - get(paste0(x, "_b"))] }
  rhs <- c(paste0(xs, "_w"), paste0(xs, "_b"), zs); rhs <- rhs[sapply(rhs, function(v) sd(d[[v]]) > 0)]
  fit <- lmer(as.formula(paste("scq ~", paste(rhs, collapse = " + "), "+ (1 | FID)")), data = d, REML = TRUE)
  b <- fixef(fit); V <- as.matrix(clubSandwich::vcovCR(fit, cluster = d$FID, type = "CR2")); V0 <- as.matrix(vcov(fit))
  ok <- function(k) k %in% names(b) && !is.na(b[k])
  rows <- lapply(xs, function(x) {
    w <- paste0(x, "_w"); bb <- paste0(x, "_b")
    bw <- if (ok(w)) unname(b[w]) else NA_real_; sw <- if (ok(w)) sqrt(V[w, w]) else NA_real_
    bB <- if (ok(bb)) unname(b[bb]) else NA_real_; sB <- if (ok(bb)) sqrt(V[bb, bb]) else NA_real_
    vd <- if (ok(w) && ok(bb)) V[bb, bb] + V[w, w] - 2 * V[w, bb] else NA_real_
    data.table(term = x, within_est = bw, within_se = sw, within_p = 2 * pnorm(-abs(bw / sw)),
               between_est = bB, between_se = sB, between_p = 2 * pnorm(-abs(bB / sB)),
               diff_est = bB - bw, diff_se = sqrt(vd), diff_p = 2 * pnorm(-abs((bB - bw) / sqrt(vd))))
  })
  list(est = rbindlist(rows), tau2 = as.numeric(VarCorr(fit)$FID[1]), sigma2 = sigma(fit)^2)
}

# Tetrachoric correlation by maximum likelihood (bivariate normal thresholds), with a numerical SE.
tetra <- function(x, y) {
  ok <- !is.na(x) & !is.na(y); x <- x[ok]; y <- y[ok]; n <- length(x)
  tx <- qnorm(1 - mean(x)); ty <- qnorm(1 - mean(y))
  n11 <- sum(x == 1 & y == 1); n10 <- sum(x == 1 & y == 0); n01 <- sum(x == 0 & y == 1); n00 <- n - n11 - n10 - n01
  ll <- function(r) { p11 <- pbivnorm(-tx, -ty, r); p1x <- 1 - pnorm(tx); px1 <- 1 - pnorm(ty)
    p10 <- p1x - p11; p01 <- px1 - p11; p00 <- 1 - p11 - p10 - p01
    n11 * log(max(p11, 1e-12)) + n10 * log(max(p10, 1e-12)) + n01 * log(max(p01, 1e-12)) + n00 * log(max(p00, 1e-12)) }
  o <- optimize(ll, c(-0.995, 0.995), maximum = TRUE); r <- o$maximum; h <- 1e-4
  info <- -(ll(r + h) - 2 * ll(r) + ll(r - h)) / h^2
  phi <- cor(x, y)
  list(rho = r, se = if (is.finite(info) && info > 0) 1 / sqrt(info) else NA_real_, phi = phi, n = n, n11 = n11, n10 = n10, n01 = n01, n00 = n00)
}

est_rows <- list(); samp_rows <- list(); tet_rows <- list()
for (cohort in COHORTS) {
  d0 <- read_input(in_dir, cohort)[dnm_ok == 1 & !is.na(scq)]
  gen <- c(GEN_MAIN, GEN_CTRL, PGS); gen <- gen[sapply(gen, function(v) v %in% names(d0) && any(!is.na(d0[[v]])))]
  zs <- c(grep("^PC[0-9]+$", names(d0), value = TRUE), grep("^src_", names(d0), value = TRUE), "coverage"); zs <- zs[zs %in% names(d0)]
  vars <- unique(c("FID", "scq", EXPO, gen, BASE, zs)); d <- d0[complete.cases(d0[, ..vars])]   # twin status is required only by the no_twins analysis (missing for ~11% of WES children)
  fam <- d[, .N, by = FID]
  cat(sprintf("== %s (%s): %d children, %d families, %d with >= 2 children; genetic terms: %s\n", cohort, subdir, nrow(d), nrow(fam), sum(fam$N >= 2), paste(gen, collapse = " ")))
  for (i in 1:2) for (j in (i + 1):3) if (j <= 3) { t <- tetra(d[[EXPO[i]]], d[[EXPO[j]]])
    tet_rows[[length(tet_rows) + 1]] <- data.table(cohort = cohort, x = EXPO[i], y = EXPO[j], tetrachoric = t$rho, se = t$se, phi = t$phi, n = t$n, both = t$n11, x_only = t$n10, y_only = t$n01, neither = t$n00) }
  analyses <- list(main = d, M0 = d, male = d[sex_male == 1], female = d[sex_male == 0], no_twins = d[!is.na(twin) & twin == 0], term_births = d[premature == 0], cases = d[asd == 1], controls = d[asd == 0])
  for (a in names(analyses)) {
    da <- copy(analyses[[a]]); fa <- da[, .N, by = FID]
    xs <- c(EXPO, if (a == "M0") character(0) else gen, BASE)
    if (a == "term_births") xs <- setdiff(xs, "premature")
    if (a %in% c("cases", "controls")) xs <- setdiff(xs, "asd")
    if (a %in% c("male", "female")) xs <- setdiff(xs, "sex_male")
    samp <- data.table(cohort = cohort, analysis = a, n = nrow(da), families = nrow(fa), families_multi = sum(fa$N >= 2), children_in_multi = sum(fa$N[fa$N >= 2]),
                       cases = sum(da$asd), n_nicu = sum(da$nicu_oxygen), n_premature = sum(da$premature), n_weight = sum(da$weight_gain_diff))
    if (sum(fa$N >= 2) < 20) { samp_rows[[length(samp_rows) + 1]] <- samp; next }
    r <- tryCatch(fit_hybrid(da, xs, zs), error = function(e) { cat("   ", a, "failed:", conditionMessage(e), "\n"); NULL })
    if (is.null(r)) { samp_rows[[length(samp_rows) + 1]] <- samp; next }
    samp[, `:=`(tau2 = r$tau2, sigma2 = r$sigma2)]; samp_rows[[length(samp_rows) + 1]] <- samp
    est_rows[[length(est_rows) + 1]] <- r$est[, `:=`(cohort = cohort, analysis = a)]
    cat(sprintf("   %-12s n %6d  multi-child families %5d  NICU within %+.3f (%.3f) between %+.3f (%.3f)\n", a, nrow(da), sum(fa$N >= 2),
                r$est[term == "nicu_oxygen", within_est], r$est[term == "nicu_oxygen", within_se], r$est[term == "nicu_oxygen", between_est], r$est[term == "nicu_oxygen", between_se]))
  }
}
est <- rbindlist(est_rows); samp <- rbindlist(samp_rows, fill = TRUE); tet <- rbindlist(tet_rows)
# pooled tetrachoric over the three cohorts (stacked children)
dd <- rbindlist(lapply(COHORTS, function(cohort) { d0 <- read_input(in_dir, cohort)[dnm_ok == 1 & !is.na(scq)]; d0[, ..EXPO] }))
for (i in 1:2) for (j in (i + 1):3) if (j <= 3) { t <- tetra(dd[[EXPO[i]]], dd[[EXPO[j]]])
  tet <- rbind(tet, data.table(cohort = "pooled", x = EXPO[i], y = EXPO[j], tetrachoric = t$rho, se = t$se, phi = t$phi, n = t$n, both = t$n11, x_only = t$n10, y_only = t$n01, neither = t$n00)) }
by <- c("analysis", "term")
meta <- rbind(ivw_meta(est, "within_est", "within_se", by = by)[, component := "within"],
              ivw_meta(est, "between_est", "between_se", by = by)[, component := "between"],
              ivw_meta(est, "diff_est", "diff_se", by = by)[, component := "between_minus_within"])
fwrite(est, file.path(out_dir, "fig1_estimates.tsv"), sep = "\t"); fwrite(meta, file.path(out_dir, "fig1_meta.tsv"), sep = "\t")
fwrite(samp, file.path(out_dir, "fig1_samples.tsv"), sep = "\t"); fwrite(tet, file.path(out_dir, "fig1_tetrachoric.tsv"), sep = "\t")
cat("\nMeta, main model, exposures:\n"); print(meta[analysis == "main" & term %in% EXPO & component != "between_minus_within", .(term, component, est = round(est, 3), se = round(se, 3), p = signif(pvalue, 2), k)])
cat("\nTetrachoric:\n"); print(tet[, .(cohort, x, y, tetrachoric = round(tetrachoric, 3), se = round(se, 3), phi = round(phi, 3), both)])
cat("Wrote", out_dir, "\n")
