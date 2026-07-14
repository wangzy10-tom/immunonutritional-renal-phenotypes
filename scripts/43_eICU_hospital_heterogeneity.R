# ==============================================================================
# eICU hospital-level heterogeneity and leave-one-hospital-out audit
# ==============================================================================

required_packages <- c("dplyr", "readr", "tibble", "tidyr", "lme4", "metafor")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop("Missing required package(s): ", paste(missing_packages, collapse = ", "), call. = FALSE)
}

root <- normalizePath(
  Sys.getenv("PROJECT_ROOT", unset = getwd()), winslash = "/", mustWork = TRUE
)
source_path <- file.path(root, "output", "eicu_24h_robust", "eICU_24h_robust_results.rds")
output_dir <- file.path(root, "output", "eicu_hospital_heterogeneity")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

extract_glmer_p1 <- function(fit, model, data) {
  coefficients <- summary(fit)$coefficients
  term <- if ("p1" %in% rownames(coefficients)) "p1" else "phenotypeP1"
  beta <- coefficients[term, "Estimate"]
  se <- coefficients[term, "Std. Error"]
  tibble::tibble(
    model = model,
    n = nrow(data),
    hospitals = dplyr::n_distinct(data$hospitalid),
    events = sum(data$hospital_mortality == 1),
    comparison = "P1 vs P3",
    estimate = exp(beta),
    lower_95 = exp(beta - 1.96 * se),
    upper_95 = exp(beta + 1.96 * se),
    p_value = 2 * stats::pnorm(abs(beta / se), lower.tail = FALSE),
    effect_95ci = sprintf("%.2f (%.2f-%.2f)", exp(beta), exp(beta - 1.96 * se), exp(beta + 1.96 * se)),
    singular = lme4::isSingular(fit, tol = 1e-5)
  )
}

cat("Loading frozen strict first-24-hour eICU APACHE cohort...\n")
eicu_results <- readRDS(source_path)
data <- eicu_results$model_data |>
  dplyr::mutate(
    hospitalid = as.integer(hospitalid),
    phenotype = factor(as.character(phenotype), levels = c("P3", "P2", "P1")),
    p1 = as.integer(phenotype == "P1"),
    p2 = as.integer(phenotype == "P2"),
    male = as.integer(male)
  )
if (nrow(data) != 12548L || dplyr::n_distinct(data$hospitalid) != 166L) {
  stop("Unexpected eICU primary model cohort dimensions.", call. = FALSE)
}

hospital_info <- data |>
  dplyr::group_by(hospitalid) |>
  dplyr::summarise(
    total_n = dplyr::n(),
    p1_n = sum(phenotype == "P1"),
    p2_n = sum(phenotype == "P2"),
    p3_n = sum(phenotype == "P3"),
    p1_p3_n = sum(phenotype %in% c("P1", "P3")),
    p1_p3_deaths = sum(hospital_mortality == 1 & phenotype %in% c("P1", "P3")),
    p1_p3_survivors = sum(hospital_mortality == 0 & phenotype %in% c("P1", "P3")),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    meets_total_n = total_n >= 100,
    meets_p1_n = p1_n >= 20,
    meets_p3_n = p3_n >= 20,
    meets_deaths = p1_p3_deaths >= 20,
    meets_survivors = p1_p3_survivors >= 20,
    high_information = meets_total_n & meets_p1_n & meets_p3_n & meets_deaths & meets_survivors,
    exclusion_reason = dplyr::case_when(
      high_information ~ "Included",
      !meets_total_n ~ "Total n < 100",
      !meets_p1_n ~ "P1 n < 20",
      !meets_p3_n ~ "P3 n < 20",
      !meets_deaths ~ "P1/P3 deaths < 20",
      !meets_survivors ~ "P1/P3 survivors < 20",
      TRUE ~ "Other"
    )
  ) |>
  dplyr::arrange(dplyr::desc(high_information), dplyr::desc(total_n), hospitalid)

eligible_hospitals <- hospital_info |>
  dplyr::filter(high_information) |>
  dplyr::pull(hospitalid)
if (length(eligible_hospitals) < 10L) {
  stop("Fewer than 10 hospitals meet the pre-specified information criteria.", call. = FALSE)
}

cat("Fitting adjusted P1-versus-P3 models in ", length(eligible_hospitals), " high-information hospitals...\n", sep = "")
hospital_effects_list <- lapply(eligible_hospitals, function(hospital) {
  hospital_data <- data |>
    dplyr::filter(hospitalid == hospital, phenotype %in% c("P1", "P3")) |>
    dplyr::mutate(p1_local = as.integer(phenotype == "P1"))

  fit <- tryCatch(
    suppressWarnings(stats::glm(
      hospital_mortality ~ p1_local + age_num + male + apachescore,
      data = hospital_data,
      family = stats::binomial(),
      control = stats::glm.control(maxit = 100)
    )),
    error = function(error) error
  )

  if (inherits(fit, "error")) {
    return(tibble::tibble(
      hospitalid = hospital, model_n = nrow(hospital_data), deaths = sum(hospital_data$hospital_mortality),
      converged = FALSE, estimable = FALSE, log_or = NA_real_, standard_error = NA_real_,
      exclusion_reason_model = conditionMessage(fit)
    ))
  }

  coefficients <- summary(fit)$coefficients
  estimable <- "p1_local" %in% rownames(coefficients) &&
    all(is.finite(coefficients["p1_local", c("Estimate", "Std. Error")])) &&
    coefficients["p1_local", "Std. Error"] > 0
  tibble::tibble(
    hospitalid = hospital,
    model_n = nrow(hospital_data),
    deaths = sum(hospital_data$hospital_mortality == 1),
    converged = isTRUE(fit$converged),
    estimable = estimable,
    log_or = if (estimable) coefficients["p1_local", "Estimate"] else NA_real_,
    standard_error = if (estimable) coefficients["p1_local", "Std. Error"] else NA_real_,
    exclusion_reason_model = if (isTRUE(fit$converged) && estimable) "Included" else "Non-converged or non-estimable"
  )
})

hospital_effects <- dplyr::bind_rows(hospital_effects_list) |>
  dplyr::left_join(
    hospital_info |>
      dplyr::select(hospitalid, total_n, p1_n, p3_n, p1_p3_deaths, p1_p3_survivors),
    by = "hospitalid"
  ) |>
  dplyr::mutate(
    estimate = exp(log_or),
    lower_95 = exp(log_or - 1.96 * standard_error),
    upper_95 = exp(log_or + 1.96 * standard_error),
    effect_95ci = ifelse(
      estimable,
      sprintf("%.2f (%.2f-%.2f)", estimate, lower_95, upper_95),
      NA_character_
    )
  ) |>
  dplyr::arrange(hospitalid)

meta_data <- hospital_effects |>
  dplyr::filter(converged, estimable, is.finite(log_or), is.finite(standard_error))
if (nrow(meta_data) < 10L) {
  stop("Fewer than 10 hospital-specific effects are estimable.", call. = FALSE)
}

cat("Fitting REML Hartung-Knapp random-effects model...\n")
meta_fit <- metafor::rma.uni(
  yi = log_or,
  sei = standard_error,
  data = meta_data,
  method = "REML",
  test = "knha",
  slab = hospitalid
)
meta_prediction <- predict(meta_fit)
meta_summary <- tibble::tibble(
  model = "Hospital-specific APACHE-adjusted P1 vs P3 random-effects meta-analysis",
  hospitals = meta_fit$k,
  participants_p1_p3 = sum(meta_data$model_n),
  deaths_p1_p3 = sum(meta_data$deaths),
  pooled_or = exp(as.numeric(meta_fit$b)),
  lower_95 = exp(meta_fit$ci.lb),
  upper_95 = exp(meta_fit$ci.ub),
  prediction_lower_95 = exp(meta_prediction$pi.lb),
  prediction_upper_95 = exp(meta_prediction$pi.ub),
  p_value = meta_fit$pval,
  q_statistic = meta_fit$QE,
  q_df = meta_fit$k - 1,
  q_p_value = meta_fit$QEp,
  tau_squared = meta_fit$tau2,
  i_squared_percent = meta_fit$I2,
  effect_95ci = sprintf("%.2f (%.2f-%.2f)", exp(as.numeric(meta_fit$b)), exp(meta_fit$ci.lb), exp(meta_fit$ci.ub)),
  prediction_interval = sprintf("%.2f-%.2f", exp(meta_prediction$pi.lb), exp(meta_prediction$pi.ub))
)

leave_one_out <- dplyr::bind_rows(lapply(meta_data$hospitalid, function(omitted_hospital) {
  remaining <- meta_data |>
    dplyr::filter(hospitalid != omitted_hospital)
  fit <- metafor::rma.uni(
    yi = log_or,
    sei = standard_error,
    data = remaining,
    method = "REML",
    test = "knha"
  )
  prediction <- predict(fit)
  tibble::tibble(
    omitted_hospitalid = omitted_hospital,
    hospitals_remaining = fit$k,
    pooled_or = exp(as.numeric(fit$b)),
    lower_95 = exp(fit$ci.lb),
    upper_95 = exp(fit$ci.ub),
    prediction_lower_95 = exp(prediction$pi.lb),
    prediction_upper_95 = exp(prediction$pi.ub),
    p_value = fit$pval,
    i_squared_percent = fit$I2,
    tau_squared = fit$tau2,
    effect_95ci = sprintf("%.2f (%.2f-%.2f)", exp(as.numeric(fit$b)), exp(fit$ci.lb), exp(fit$ci.ub))
  )
}))

leave_one_out_summary <- tibble::tibble(
  omitted_models = nrow(leave_one_out),
  pooled_or_min = min(leave_one_out$pooled_or),
  pooled_or_max = max(leave_one_out$pooled_or),
  lower_95_min = min(leave_one_out$lower_95),
  lower_95_max = max(leave_one_out$lower_95),
  upper_95_min = min(leave_one_out$upper_95),
  upper_95_max = max(leave_one_out$upper_95),
  models_with_ci_above_one = sum(leave_one_out$lower_95 > 1),
  models_with_ci_crossing_one = sum(leave_one_out$lower_95 <= 1 & leave_one_out$upper_95 >= 1),
  i_squared_min = min(leave_one_out$i_squared_percent),
  i_squared_max = max(leave_one_out$i_squared_percent)
)

cat("Fitting full-cohort random-intercept and exploratory random-slope models...\n")
random_intercept_fit <- lme4::glmer(
  hospital_mortality ~ p1 + p2 + age_num + male + apachescore + (1 | hospitalid),
  family = stats::binomial(), data = data, nAGQ = 1,
  control = lme4::glmerControl(optimizer = "bobyqa")
)
random_slope_fit <- suppressWarnings(lme4::glmer(
  hospital_mortality ~ p1 + p2 + age_num + male + apachescore + (1 + p1 || hospitalid),
  family = stats::binomial(), data = data, nAGQ = 1,
  control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 200000))
))
model_comparison <- suppressWarnings(stats::anova(random_intercept_fit, random_slope_fit, test = "Chisq"))
random_slope_variance <- as.data.frame(lme4::VarCorr(random_slope_fit)) |>
  dplyr::filter(grp == "hospitalid", var1 == "p1", is.na(var2))

high_information_data <- data |>
  dplyr::filter(hospitalid %in% eligible_hospitals)
high_information_fit <- lme4::glmer(
  hospital_mortality ~ p1 + p2 + age_num + male + apachescore + (1 | hospitalid),
  family = stats::binomial(), data = high_information_data, nAGQ = 1,
  control = lme4::glmerControl(optimizer = "bobyqa")
)

mixed_model_effects <- dplyr::bind_rows(
  extract_glmer_p1(random_intercept_fit, "All hospitals, APACHE-adjusted random intercept", data),
  extract_glmer_p1(random_slope_fit, "All hospitals, APACHE-adjusted uncorrelated P1 random slope", data),
  extract_glmer_p1(high_information_fit, "High-information hospitals, APACHE-adjusted random intercept", high_information_data)
)

random_slope_diagnostics <- tibble::tibble(
  random_intercept_singular = lme4::isSingular(random_intercept_fit, tol = 1e-5),
  random_slope_singular = lme4::isSingular(random_slope_fit, tol = 1e-5),
  p1_random_slope_sd = if (nrow(random_slope_variance) == 1L) random_slope_variance$sdcor else NA_real_,
  likelihood_ratio_chisq = model_comparison$Chisq[2],
  likelihood_ratio_df = model_comparison$`Chi Df`[2],
  likelihood_ratio_p = model_comparison$`Pr(>Chisq)`[2]
)

direction_summary <- tibble::tibble(
  estimable_hospitals = nrow(meta_data),
  hospitals_or_above_one = sum(meta_data$estimate > 1),
  hospitals_or_below_one = sum(meta_data$estimate < 1),
  hospitals_ci_above_one = sum(meta_data$lower_95 > 1),
  hospitals_ci_crossing_one = sum(meta_data$lower_95 <= 1 & meta_data$upper_95 >= 1),
  hospitals_ci_below_one = sum(meta_data$upper_95 < 1)
)

readr::write_csv(hospital_info, file.path(output_dir, "Table43A_hospital_information_and_eligibility.csv"))
readr::write_csv(hospital_effects, file.path(output_dir, "Table43B_hospital_specific_P1_effects.csv"))
readr::write_csv(meta_summary, file.path(output_dir, "Table43C_random_effects_meta_analysis.csv"))
readr::write_csv(leave_one_out, file.path(output_dir, "Table43D_leave_one_hospital_out.csv"))
readr::write_csv(leave_one_out_summary, file.path(output_dir, "Table43E_leave_one_out_summary.csv"))
readr::write_csv(mixed_model_effects, file.path(output_dir, "Table43F_mixed_model_effects.csv"))
readr::write_csv(random_slope_diagnostics, file.path(output_dir, "Table43G_random_slope_diagnostics.csv"))
readr::write_csv(direction_summary, file.path(output_dir, "Table43H_hospital_effect_direction_summary.csv"))
saveRDS(
  list(
    hospital_info = hospital_info,
    hospital_effects = hospital_effects,
    meta_summary = meta_summary,
    leave_one_out = leave_one_out,
    leave_one_out_summary = leave_one_out_summary,
    mixed_model_effects = mixed_model_effects,
    random_slope_diagnostics = random_slope_diagnostics,
    direction_summary = direction_summary,
    meta_fit = meta_fit
  ),
  file.path(output_dir, "eICU_hospital_heterogeneity_results.rds")
)

summary_lines <- c(
  "eICU hospital-level heterogeneity audit",
  "",
  paste0("Primary model cohort: ", nrow(data), " patients across ", dplyr::n_distinct(data$hospitalid), " hospitals."),
  paste0("Pre-specified high-information hospitals: ", length(eligible_hospitals), "."),
  paste0("Estimable hospital-specific models: ", nrow(meta_data), "."),
  "",
  "Random-effects meta-analysis:",
  paste(capture.output(print(meta_summary)), collapse = "\n"),
  "",
  "Hospital effect directions:",
  paste(capture.output(print(direction_summary)), collapse = "\n"),
  "",
  "Leave-one-hospital-out summary:",
  paste(capture.output(print(leave_one_out_summary)), collapse = "\n"),
  "",
  "Mixed models:",
  paste(capture.output(print(mixed_model_effects)), collapse = "\n"),
  "",
  "Random-slope diagnostics:",
  paste(capture.output(print(random_slope_diagnostics)), collapse = "\n"),
  "",
  "Interpretation: this is a post-freeze transportability audit; non-significant individual hospitals are not evidence of effect absence."
)
writeLines(summary_lines, file.path(output_dir, "eICU_hospital_heterogeneity_summary.txt"))
cat(paste(summary_lines, collapse = "\n"), "\n")
