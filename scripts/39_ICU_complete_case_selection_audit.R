# ==============================================================================
# Post-freeze MIMIC-IV and eICU complete-case selection-bias audit
# ==============================================================================

required_packages <- c(
  "data.table", "dplyr", "readr", "tibble", "survival", "lme4"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop("Missing required package(s): ", paste(missing_packages, collapse = ", "), call. = FALSE)
}

root <- normalizePath(
  Sys.getenv("PROJECT_ROOT", unset = getwd()), winslash = "/", mustWork = TRUE
)
output_dir <- file.path(root, "output", "icu_selection_bias")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

mimic_denominator_path <- file.path(
  root, "output", "mimic_24h_validation",
  "MIMIC_first_ICU_feature_availability_dataset.csv"
)
mimic_result_path <- file.path(
  root, "output", "mimic_24h_robust_corrected_severity",
  "MIMIC_24h_robust_results.rds"
)
eicu_denominator_path <- file.path(
  root, "output", "eicu_24h_extraction",
  "eICU_first_ICU_feature_availability_dataset.csv"
)
eicu_result_path <- file.path(
  root, "output", "eicu_24h_robust",
  "eICU_24h_robust_results.rds"
)

format_effect <- function(estimate, lower, upper) {
  sprintf("%.2f (%.2f-%.2f)", estimate, lower, upper)
}

collapse_for_selection <- function(x, selected, min_total = 100L, min_selected = 10L) {
  value <- as.character(x)
  value[is.na(value) | trimws(value) == ""] <- "Missing"
  tab_total <- table(value)
  tab_selected <- table(value[selected == 1], useNA = "no")
  selected_count <- as.numeric(tab_selected[names(tab_total)])
  selected_count[is.na(selected_count)] <- 0
  excluded_count <- as.numeric(tab_total) - selected_count
  keep <- names(tab_total)[
    as.numeric(tab_total) >= min_total &
      selected_count >= min_selected &
      excluded_count >= min_selected
  ]
  factor(ifelse(value %in% keep, value, "Other"))
}

continuous_smd <- function(x, selected) {
  x <- suppressWarnings(as.numeric(x))
  included <- x[selected == 1]
  excluded <- x[selected == 0]
  included_mean <- mean(included, na.rm = TRUE)
  excluded_mean <- mean(excluded, na.rm = TRUE)
  included_sd <- stats::sd(included, na.rm = TRUE)
  excluded_sd <- stats::sd(excluded, na.rm = TRUE)
  pooled_sd <- sqrt(
    (included_sd^2 + excluded_sd^2) / 2
  )
  tibble::tibble(
    included = sprintf(
      "%.2f (%.2f)", included_mean, included_sd
    ),
    excluded = sprintf(
      "%.2f (%.2f)", excluded_mean, excluded_sd
    ),
    standardised_mean_difference = (included_mean - excluded_mean) / pooled_sd
  )
}

binary_smd <- function(x, selected) {
  p_included <- mean(x[selected == 1], na.rm = TRUE)
  p_excluded <- mean(x[selected == 0], na.rm = TRUE)
  pooled_p <- (p_included + p_excluded) / 2
  denominator <- sqrt(pooled_p * (1 - pooled_p))
  tibble::tibble(
    included = sprintf("%.1f%%", 100 * p_included),
    excluded = sprintf("%.1f%%", 100 * p_excluded),
    standardised_mean_difference = ifelse(
      denominator > 0, (p_included - p_excluded) / denominator, NA_real_
    )
  )
}

category_smd_rows <- function(data, dataset, variable, selected_variable) {
  selected <- data[[selected_variable]]
  value <- as.character(data[[variable]])
  levels <- sort(unique(value))
  dplyr::bind_rows(lapply(levels, function(level) {
    binary_smd(value == level, selected) |>
      dplyr::mutate(
        dataset = dataset,
        variable = variable,
        category = level,
        .before = 1
      )
  }))
}

derive_selection_weights <- function(data, model, selected_variable, dataset) {
  if (!isTRUE(model$converged)) {
    stop(dataset, " selection model did not converge.", call. = FALSE)
  }
  probability <- stats::predict(model, type = "response")
  probability <- pmin(pmax(probability, 0.001), 0.999)
  selected <- data[[selected_variable]]
  selection_rate <- mean(selected == 1)
  stabilised <- ifelse(selected == 1, selection_rate / probability, NA_real_)
  limits <- stats::quantile(stabilised[selected == 1], c(0.01, 0.99), na.rm = TRUE)
  trimmed <- ifelse(
    selected == 1,
    pmin(pmax(stabilised, limits[1]), limits[2]),
    NA_real_
  )
  effective_n <- sum(trimmed, na.rm = TRUE)^2 / sum(trimmed^2, na.rm = TRUE)

  list(
    probability = probability,
    weight = trimmed,
    diagnostics = tibble::tibble(
      dataset = dataset,
      denominator_n = nrow(data),
      feature_complete_n = sum(data$selected_feature == 1),
      primary_model_n = sum(selected == 1),
      primary_model_percent = 100 * selection_rate,
      selected_probability_min = min(probability[selected == 1]),
      selected_probability_median = stats::median(probability[selected == 1]),
      selected_probability_max = max(probability[selected == 1]),
      trimmed_weight_1st_percentile = limits[1],
      trimmed_weight_99th_percentile = limits[2],
      trimmed_weight_max = max(trimmed, na.rm = TRUE),
      effective_sample_size = effective_n
    )
  )
}

extract_cox_phenotype <- function(fit, dataset, model_label) {
  beta <- stats::coef(fit)
  se <- sqrt(diag(stats::vcov(fit)))
  terms <- names(beta)
  keep <- grepl("^phenotype", terms)
  tibble::tibble(
    dataset = dataset,
    model = model_label,
    effect_measure = "HR",
    comparison = dplyr::recode(
      terms[keep], phenotypeP2 = "P2 vs P3", phenotypeP1 = "P1 vs P3"
    ),
    estimate = exp(beta[keep]),
    lower_95 = exp(beta[keep] - 1.96 * se[keep]),
    upper_95 = exp(beta[keep] + 1.96 * se[keep]),
    p_value = 2 * stats::pnorm(abs(beta[keep] / se[keep]), lower.tail = FALSE)
  ) |>
    dplyr::mutate(effect_95ci = format_effect(estimate, lower_95, upper_95))
}

extract_glmer_phenotype <- function(fit, dataset, model_label) {
  coefficient_table <- summary(fit)$coefficients
  beta <- coefficient_table[, "Estimate"]
  se <- coefficient_table[, "Std. Error"]
  terms <- rownames(coefficient_table)
  keep <- grepl("^phenotype", terms)
  tibble::tibble(
    dataset = dataset,
    model = model_label,
    effect_measure = "OR",
    comparison = dplyr::recode(
      terms[keep], phenotypeP2 = "P2 vs P3", phenotypeP1 = "P1 vs P3"
    ),
    estimate = exp(beta[keep]),
    lower_95 = exp(beta[keep] - 1.96 * se[keep]),
    upper_95 = exp(beta[keep] + 1.96 * se[keep]),
    p_value = 2 * stats::pnorm(abs(beta[keep] / se[keep]), lower.tail = FALSE)
  ) |>
    dplyr::mutate(effect_95ci = format_effect(estimate, lower_95, upper_95))
}

# MIMIC-IV ---------------------------------------------------------------------
mimic_results <- readRDS(mimic_result_path)
mimic_denominator <- readr::read_csv(mimic_denominator_path, show_col_types = FALSE)
mimic_feature_ids <- mimic_results$analysis$stay_id
mimic_primary_ids <- mimic_results$analysis_oasis$stay_id

mimic_denominator <- mimic_denominator |>
  dplyr::mutate(
    selected_feature = as.integer(stay_id %in% mimic_feature_ids),
    selected_primary = as.integer(stay_id %in% mimic_primary_ids),
    gender_selection = collapse_for_selection(gender, selected_primary, 100, 10),
    race_selection = collapse_for_selection(race, selected_primary, 100, 10),
    careunit_selection = collapse_for_selection(first_careunit, selected_primary, 100, 10),
    year_selection = collapse_for_selection(anchor_year_group, selected_primary, 100, 10)
  )

if (sum(mimic_denominator$selected_feature) != nrow(mimic_results$analysis) ||
    sum(mimic_denominator$selected_primary) != nrow(mimic_results$analysis_oasis)) {
  stop("MIMIC selection indicators do not reproduce frozen cohort sizes.", call. = FALSE)
}

mimic_selection_model <- stats::glm(
  selected_primary ~ splines::ns(anchor_age, df = 3) + gender_selection +
    race_selection + careunit_selection + year_selection,
  data = mimic_denominator,
  family = stats::binomial(),
  control = stats::glm.control(maxit = 100)
)
mimic_weights <- derive_selection_weights(
  mimic_denominator, mimic_selection_model, "selected_primary", "MIMIC-IV"
)
mimic_denominator$selection_probability <- mimic_weights$probability
mimic_denominator$selection_ipw <- mimic_weights$weight

mimic_model_data <- mimic_results$analysis_oasis |>
  dplyr::left_join(
    mimic_denominator |>
      dplyr::select(stay_id, selection_probability, selection_ipw),
    by = "stay_id"
  ) |>
  dplyr::mutate(
    phenotype = factor(phenotype, levels = c("P3", "P2", "P1")),
    selection_ipw_scaled = selection_ipw / mean(selection_ipw)
  )

mimic_unweighted_fit <- survival::coxph(
  survival::Surv(survival_days_365, mortality_365d) ~
    phenotype + male + strata(oasis_quartile),
  data = mimic_model_data,
  ties = "efron"
)
mimic_ipw_fit <- survival::coxph(
  survival::Surv(survival_days_365, mortality_365d) ~
    phenotype + male + strata(oasis_quartile),
  data = mimic_model_data,
  weights = selection_ipw_scaled,
  robust = TRUE,
  cluster = stay_id,
  ties = "efron"
)

mimic_effects <- dplyr::bind_rows(
  extract_cox_phenotype(mimic_unweighted_fit, "MIMIC-IV", "Frozen OASIS-like stratified Cox"),
  extract_cox_phenotype(mimic_ipw_fit, "MIMIC-IV", "Selection-IPW OASIS-like stratified Cox")
) |>
  dplyr::mutate(
    n = nrow(mimic_model_data),
    events = sum(mimic_model_data$mortality_365d == 1)
  )

# eICU -------------------------------------------------------------------------
eicu_results <- readRDS(eicu_result_path)
eicu_denominator <- data.table::fread(eicu_denominator_path, data.table = FALSE) |>
  tibble::as_tibble()
eicu_feature_ids <- eicu_results$analysis$patientunitstayid
eicu_primary_ids <- eicu_results$model_data$patientunitstayid

eicu_denominator <- eicu_denominator |>
  dplyr::mutate(
    selected_feature = as.integer(patientunitstayid %in% eicu_feature_ids),
    selected_primary = as.integer(patientunitstayid %in% eicu_primary_ids),
    gender_selection = collapse_for_selection(gender, selected_primary, 100, 10),
    ethnicity_selection = collapse_for_selection(ethnicity, selected_primary, 100, 10),
    unittype_selection = collapse_for_selection(unittype, selected_primary, 100, 10),
    hospital_key = ifelse(is.na(hospitalid), "Missing", as.character(hospitalid))
  )

if (sum(eicu_denominator$selected_feature) != nrow(eicu_results$analysis) ||
    sum(eicu_denominator$selected_primary) != nrow(eicu_results$model_data)) {
  stop("eICU selection indicators do not reproduce frozen cohort sizes.", call. = FALSE)
}

eicu_overall_rate <- mean(eicu_denominator$selected_primary)
hospital_measurement <- eicu_denominator |>
  dplyr::group_by(hospital_key) |>
  dplyr::summarise(
    hospital_n = dplyr::n(),
    hospital_selected = sum(selected_primary),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    smoothed_selection_rate =
      (hospital_selected + 20 * eicu_overall_rate) / (hospital_n + 20),
    hospital_selection_logit = stats::qlogis(
      pmin(pmax(smoothed_selection_rate, 0.001), 0.999)
    )
  )
eicu_denominator <- eicu_denominator |>
  dplyr::left_join(hospital_measurement, by = "hospital_key")

eicu_selection_model <- stats::glm(
  selected_primary ~ splines::ns(age_num, df = 3) + gender_selection +
    ethnicity_selection + unittype_selection + hospital_selection_logit,
  data = eicu_denominator,
  family = stats::binomial(),
  control = stats::glm.control(maxit = 100)
)
eicu_weights <- derive_selection_weights(
  eicu_denominator, eicu_selection_model, "selected_primary", "eICU"
)
eicu_denominator$selection_probability <- eicu_weights$probability
eicu_denominator$selection_ipw <- eicu_weights$weight

eicu_model_data <- eicu_results$model_data |>
  dplyr::left_join(
    eicu_denominator |>
      dplyr::select(patientunitstayid, selection_probability, selection_ipw),
    by = "patientunitstayid"
  ) |>
  dplyr::mutate(
    phenotype = factor(phenotype, levels = c("P3", "P2", "P1")),
    selection_ipw_scaled = selection_ipw / mean(selection_ipw)
  )

eicu_unweighted_fit <- lme4::glmer(
  hospital_mortality ~ phenotype + age_num + male + apachescore + (1 | hospitalid),
  family = stats::binomial(),
  data = eicu_model_data,
  nAGQ = 1,
  control = lme4::glmerControl(optimizer = "bobyqa")
)
eicu_ipw_fit <- suppressWarnings(lme4::glmer(
  hospital_mortality ~ phenotype + age_num + male + apachescore + (1 | hospitalid),
  family = stats::binomial(),
  data = eicu_model_data,
  weights = selection_ipw_scaled,
  nAGQ = 1,
  control = lme4::glmerControl(optimizer = "bobyqa")
))

if (lme4::isSingular(eicu_ipw_fit, tol = 1e-5)) {
  warning("Selection-IPW eICU mixed model is singular; interpret cautiously.")
}

eicu_effects <- dplyr::bind_rows(
  extract_glmer_phenotype(eicu_unweighted_fit, "eICU", "Frozen APACHE-adjusted mixed model"),
  extract_glmer_phenotype(eicu_ipw_fit, "eICU", "Selection-IPW APACHE-adjusted mixed model")
) |>
  dplyr::mutate(
    n = nrow(eicu_model_data),
    events = sum(eicu_model_data$hospital_mortality == 1)
  )

# Descriptive selection comparisons -------------------------------------------
comparison_rows <- dplyr::bind_rows(
  continuous_smd(mimic_denominator$anchor_age, mimic_denominator$selected_primary) |>
    dplyr::mutate(dataset = "MIMIC-IV", variable = "Age, years", .before = 1),
  binary_smd(mimic_denominator$gender == "M", mimic_denominator$selected_primary) |>
    dplyr::mutate(dataset = "MIMIC-IV", variable = "Male", .before = 1),
  binary_smd(mimic_denominator$hospital_expire_flag == 1, mimic_denominator$selected_primary) |>
    dplyr::mutate(dataset = "MIMIC-IV", variable = "Hospital mortality", .before = 1),
  binary_smd(mimic_denominator$mortality_365d == 1, mimic_denominator$selected_primary) |>
    dplyr::mutate(dataset = "MIMIC-IV", variable = "365-day mortality", .before = 1),
  continuous_smd(eicu_denominator$age_num, eicu_denominator$selected_primary) |>
    dplyr::mutate(dataset = "eICU", variable = "Age, years", .before = 1),
  binary_smd(eicu_denominator$gender == "Male", eicu_denominator$selected_primary) |>
    dplyr::mutate(dataset = "eICU", variable = "Male", .before = 1),
  binary_smd(eicu_denominator$hospital_mortality == 1, eicu_denominator$selected_primary) |>
    dplyr::mutate(dataset = "eICU", variable = "Hospital mortality", .before = 1),
  binary_smd(eicu_denominator$icu_mortality == 1, eicu_denominator$selected_primary) |>
    dplyr::mutate(dataset = "eICU", variable = "ICU mortality", .before = 1)
) |>
  dplyr::mutate(
    absolute_smd = abs(standardised_mean_difference),
    potentially_important_imbalance = absolute_smd >= 0.10
  )

category_rows <- dplyr::bind_rows(
  category_smd_rows(mimic_denominator, "MIMIC-IV", "race_selection", "selected_primary"),
  category_smd_rows(mimic_denominator, "MIMIC-IV", "careunit_selection", "selected_primary"),
  category_smd_rows(mimic_denominator, "MIMIC-IV", "year_selection", "selected_primary"),
  category_smd_rows(eicu_denominator, "eICU", "ethnicity_selection", "selected_primary"),
  category_smd_rows(eicu_denominator, "eICU", "unittype_selection", "selected_primary")
) |>
  dplyr::mutate(
    absolute_smd = abs(standardised_mean_difference),
    potentially_important_imbalance = absolute_smd >= 0.10
  )

feature_availability <- dplyr::bind_rows(
  lapply(
    c("nlr", "sii", "haemoglobin", "protein_proxy", "bmi", "creatinine"),
    function(variable) {
      available <- is.finite(mimic_denominator[[variable]])
      tibble::tibble(
        dataset = "MIMIC-IV",
        feature = variable,
        denominator_n = nrow(mimic_denominator),
        available_n = sum(available),
        available_percent = 100 * mean(available)
      )
    }
  ),
  lapply(
    c("nlr", "sii_like", "haemoglobin", "total_protein", "bmi", "creatinine"),
    function(variable) {
      available <- is.finite(eicu_denominator[[variable]])
      tibble::tibble(
        dataset = "eICU",
        feature = variable,
        denominator_n = nrow(eicu_denominator),
        available_n = sum(available),
        available_percent = 100 * mean(available)
      )
    }
  )
)

flow <- dplyr::bind_rows(
  mimic_weights$diagnostics,
  eicu_weights$diagnostics
)
effects <- dplyr::bind_rows(mimic_effects, eicu_effects)

readr::write_csv(flow, file.path(output_dir, "Table39A_ICU_selection_flow_and_weights.csv"))
readr::write_csv(comparison_rows, file.path(output_dir, "Table39B_included_vs_excluded.csv"))
readr::write_csv(category_rows, file.path(output_dir, "Table39C_selection_category_distributions.csv"))
readr::write_csv(feature_availability, file.path(output_dir, "Table39D_feature_availability.csv"))
readr::write_csv(effects, file.path(output_dir, "Table39E_selection_IPW_outcome_models.csv"))
readr::write_csv(hospital_measurement, file.path(output_dir, "Table39F_eICU_hospital_measurement_rates.csv"))

saveRDS(
  list(
    flow = flow,
    comparison = comparison_rows,
    categories = category_rows,
    feature_availability = feature_availability,
    effects = effects,
    mimic_selection_model = mimic_selection_model,
    eicu_selection_model = eicu_selection_model
  ),
  file.path(output_dir, "ICU_selection_bias_results.rds")
)

summary_lines <- c(
  "MIMIC-IV and eICU complete-case selection-bias audit",
  "",
  "Selection flow and weight diagnostics:",
  paste(capture.output(print(flow)), collapse = "\n"),
  "",
  "Included versus excluded comparison:",
  paste(capture.output(print(comparison_rows)), collapse = "\n"),
  "",
  "Primary outcome models before and after selection IPW:",
  paste(capture.output(print(effects)), collapse = "\n"),
  "",
  paste0(
    "Interpretation: selection models exclude mortality outcomes, phenotype labels, and acute severity scores. ",
    "IPW addresses selection conditional on observed demographics and clinical setting only; unmeasured ",
    "measurement-related selection can remain."
  )
)
writeLines(summary_lines, file.path(output_dir, "ICU_selection_bias_summary.txt"))
cat(paste(summary_lines, collapse = "\n"), "\n")
