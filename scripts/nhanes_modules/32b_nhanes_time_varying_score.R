# ==============================================================================
# Time-varying effect sensitivity for the corrected NHANES continuous score
# ==============================================================================

required_packages <- c("dplyr", "readr", "tibble", "survival")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop("Missing required package(s): ", paste(missing_packages, collapse = ", "), call. = FALSE)
}

input_path <- paste0(
  "./output/nhanes_albumin_amendment17_2026-08-30/",
  "NHANES_albumin_benchmark_results.rds"
)
output_dir <- "./output/nhanes_albumin_amendment17_2026-08-30/cluster_bootstrap_PH"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

dat <- readRDS(input_path)$model_data
analysis_vars <- c(
  "PERMTH_INT", "MORTSTAT", "domain_balanced_score", "RIDAGEYR", "male",
  "race", "INDFMPIR", "Comorbidity_Score_Extended", "cycle", "smoking",
  "hypertension", "diabetes"
)
dat <- dat[stats::complete.cases(dat[, analysis_vars]), , drop = FALSE]
base_terms <- "RIDAGEYR + male + race + INDFMPIR + Comorbidity_Score_Extended + cycle + smoking + hypertension + diabetes"
ph_formula <- stats::as.formula(paste(
  "survival::Surv(PERMTH_INT, MORTSTAT) ~ domain_balanced_score +", base_terms
))
tv_formula <- stats::as.formula(paste(
  "survival::Surv(PERMTH_INT, MORTSTAT) ~ domain_balanced_score + tt(domain_balanced_score) +",
  base_terms
))

fit_ph <- survival::coxph(ph_formula, data = dat, ties = "efron")
fit_tv <- survival::coxph(
  tv_formula,
  data = dat,
  ties = "efron",
  tt = function(x, t, ...) x * log(pmax(t, 1))
)

lr_chisq <- 2 * (fit_tv$loglik[2] - fit_ph$loglik[2])
lr_df <- length(stats::coef(fit_tv)) - length(stats::coef(fit_ph))
lr_p <- stats::pchisq(lr_chisq, df = lr_df, lower.tail = FALSE)
co <- stats::coef(fit_tv)
vc <- stats::vcov(fit_tv)
main_term <- "domain_balanced_score"
time_term <- "tt(domain_balanced_score)"

evaluation_months <- c(12, 36, 60, 96)
time_effects <- dplyr::bind_rows(lapply(evaluation_months, function(month) {
  log_time <- log(month)
  log_hr <- co[main_term] + co[time_term] * log_time
  variance <- vc[main_term, main_term] + log_time^2 * vc[time_term, time_term] +
    2 * log_time * vc[main_term, time_term]
  se <- sqrt(variance)
  tibble::tibble(
    followup_month = month,
    HR_per_1_SD = exp(log_hr),
    lower_95 = exp(log_hr - 1.96 * se),
    upper_95 = exp(log_hr + 1.96 * se),
    hazard_ratio_95ci = sprintf(
      "%.2f (%.2f-%.2f)", exp(log_hr), exp(log_hr - 1.96 * se), exp(log_hr + 1.96 * se)
    )
  )
}))

model_comparison <- tibble::tibble(
  comparison = "Time-varying score effect vs constant score effect",
  n_constant_model = fit_ph$n,
  n_time_varying_model = fit_tv$n,
  likelihood_ratio_chisq = lr_chisq,
  df = lr_df,
  p_value = lr_p
)

readr::write_csv(time_effects, file.path(output_dir, "Table32E_time_varying_score_effects.csv"))
readr::write_csv(model_comparison, file.path(output_dir, "Table32F_time_varying_model_comparison.csv"))
saveRDS(
  list(fit_ph = fit_ph, fit_tv = fit_tv, time_effects = time_effects, comparison = model_comparison),
  file.path(output_dir, "NHANES_time_varying_score_results.rds")
)

summary_lines <- c(
  "Corrected NHANES time-varying continuous-score sensitivity",
  "",
  "Model comparison:",
  paste(capture.output(print(model_comparison)), collapse = "\n"),
  "",
  "Estimated score effect over follow-up:",
  paste(capture.output(print(time_effects)), collapse = "\n")
)
writeLines(summary_lines, file.path(output_dir, "NHANES_time_varying_score_summary.txt"))
message("Time-varying score analysis completed.")
cat(paste(summary_lines, collapse = "\n"), "\n")
