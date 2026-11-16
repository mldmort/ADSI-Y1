# Shared pieces for the within-family analyses (02_clogit_asd.R, 03_hybrid_scq.R).
suppressMessages(library(data.table))

cohorts_all <- c("spark_wes", "spark_wgs", "spark_wgs_2026_08", "ssc")

read_input <- function(out_dir, cohort) fread(file.path(out_dir, sprintf("input_%s.tsv", cohort)), na.strings = "NA")

# Child-level predictors. Within a family, parental age, birth year and birth order are the same
# contrast (the parents' birth years are family constants), so birth_year (decades; available for
# nearly every child, unlike father_age) stands in for all three and is also the birth-year
# adjustment of prenatal PM2.5. Prematurity and prenatal PM2.5 are SPARK only (SSC has them for
# probands only / for 1,032 children) and are kept as a separate block because PM2.5 is missing for
# half of the 2026_08 children.
child_terms <- function(cohort) c("pgs", "rv", "sex_male", "birth_year")
env_terms <- function(cohort) if (startsWith(cohort, "spark")) c("premature", "pm25") else character(0)
# Family-level moderators (constant within family): main effects drop out of the conditional
# model and enter the hybrid model only through the family-mean part.
family_terms <- function(cohort) if (startsWith(cohort, "spark")) c("adi", "pm25_fam", "ses") else c("ses")
gxe_terms <- function(cohort) if (startsWith(cohort, "spark")) c("pgs:adi", "pgs:pm25_fam", "dnm:adi") else character(0)

# Fixed-effect inverse-variance meta of one estimate across cohorts.
ivw_meta <- function(dt, est = "est", se = "se", by) {   # by: character vector of grouping columns
  dt <- dt[is.finite(get(est)) & is.finite(get(se)) & get(se) > 0]
  dt[, {
    w <- 1 / get(se)^2; b <- sum(w * get(est)) / sum(w); s <- sqrt(1 / sum(w))
    q <- sum(w * (get(est) - b)^2); k <- .N
    .(est = b, se = s, pvalue = 2 * pnorm(-abs(b / s)), ci.lower = b - 1.96 * s, ci.upper = b + 1.96 * s,
      k = k, cohorts = paste(cohort, collapse = ","), het_p = if (k > 1) pchisq(q, k - 1, lower.tail = FALSE) else NA_real_)
  }, by = by]
}
