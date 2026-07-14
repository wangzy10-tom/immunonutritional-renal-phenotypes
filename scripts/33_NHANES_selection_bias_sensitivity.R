# ==============================================================================
# NHANES complete-feature selection-bias assessment
# ==============================================================================

required_packages <- c("dplyr", "readr", "tibble", "survey", "survival")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop("Missing required package(s): ", paste(missing_packages, collapse = ", "), call. = FALSE)
}

options(survey.lonely.psu = "adjust")

root <- normalizePath(
  Sys.getenv("PROJECT_ROOT", unset = getwd()), winslash = "/", mustWork = TRUE
)
eligible_path <- file.path(
  root, "output", "nhanes_2011_2018_rebuild",
  "NHANES_2011_2018_age65_mortality_eligible.rds"
)
benchmark_path <- file.path(
  root, "output", "nhanes_albumin_benchmarks",
  "NHANES_albumin_benchmark_results.rds"
)
output_dir <- file.path(root, "output", "nhanes_selection_bias")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

eligible <- readRDS(eligible_path) |>
  dplyr::mutate(
    selected = as.integer(complete_six_feature),
    male_selection = factor(
      dplyr::case_when(RIAGENDR == 1 ~ "Male", RIAGENDR == 2 ~ "Female", TRUE ~ "Missing")
    ),
    race_selection = factor(ifelse(is.na(RIDRETH3), "Missing", paste0("Race_", RIDRETH3))),
    education_selection = factor(ifelse(is.na(DMDEDUC2), "Missing", paste0("Education_", DMDEDUC2))),
    pir_missing = as.integer(is.na(INDFMPIR)),
    pir_imputed = ifelse(
      is.na(INDFMPIR), stats::median(INDFMPIR, na.rm = TRUE), INDFMPIR
    ),
    cycle_selection = factor(Cycle_ID),
    selection_sampling_weight = WTMEC8YR / mean(WTMEC8YR, na.rm = TRUE)
  ) |>
  droplevels()

if (!all(eligible$selected %in% c(0L, 1L))) {
  stop("Invalid complete-feature selection indicator.", call. = FALSE)
}

format_continuous <- function(x) {
  sprintf("%.2f (%.2f)", mean(x, na.rm = TRUE), stats::sd(x, na.rm = TRUE))
}

continuous_row <- function(variable, x, selected) {
  x1 <- x[selected == 1]
  x0 <- x[selected == 0]
  pooled_sd <- sqrt((stats::var(x1, na.rm = TRUE) + stats::var(x0, na.rm = TRUE)) / 2)
  tibble::tibble(
    variable = variable,
    included = format_continuous(x1),
    excluded = format_continuous(x0),
    standardised_mean_difference = (mean(x1, na.rm = TRUE) - mean(x0, na.rm = TRUE)) / pooled_sd
  )
}

binary_row <- function(variable, x, selected) {
  p1 <- mean(x[selected == 1], na.rm = TRUE)
  p0 <- mean(x[selected == 0], na.rm = TRUE)
  pooled_p <- (p1 + p0) / 2
  denominator <- sqrt(pooled_p * (1 - pooled_p))
  tibble::tibble(
    variable = variable,
    included = sprintf("%.1f%%", 100 * p1),
    excluded = sprintf("%.1f%%", 100 * p0),
    standardised_mean_difference = ifelse(denominator > 0, (p1 - p0) / denominator, NA_real_)
  )
}

comparison <- dplyr::bind_rows(
  continuous_row("Age, years", eligible$RIDAGEYR, eligible$selected),
  continuous_row("Family income-to-poverty ratio", eligible$INDFMPIR, eligible$selected),
  binary_row("Male", eligible$RIAGENDR == 1, eligible$selected),
  binary_row("Death during linked follow-up", eligible$MORTSTAT == 1, eligible$selected),
  lapply(sort(unique(eligible$Cycle_ID)), function(cycle) {
    binary_row(paste0("Survey cycle ", cycle), eligible$Cycle_ID == cycle, eligible$selected)
  }),
  lapply(sort(unique(stats::na.omit(eligible$RIDRETH3))), function(race_code) {
    binary_row(paste0("Race code ", race_code), eligible$RIDRETH3 == race_code, eligible$selected)
  })
) |>
  dplyr::mutate(
    absolute_smd = abs(standardised_mean_difference),
    potentially_important_imbalance = absolute_smd >= 0.10
  )

core_features <- c("NLR", "SII", "LBXHGB", "LBXSTP", "BMXBMI", "LBXSCR")
feature_missingness <- dplyr::bind_rows(lapply(core_features, function(variable) {
  missing <- is.na(eligible[[variable]]) | !is.finite(eligible[[variable]])
  tibble::tibble(
    feature = variable,
    missing_n = sum(missing),
    missing_percent = 100 * mean(missing)
  )
})) |>
  dplyr::arrange(dplyr::desc(missing_percent))

selection_model <- stats::glm(
  selected ~ RIDAGEYR + male_selection + race_selection + pir_imputed +
    pir_missing + cycle_selection,
  data = eligible,
  family = stats::quasibinomial(),
  weights = selection_sampling_weight,
  control = stats::glm.control(maxit = 100)
)

if (!isTRUE(selection_model$converged)) {
  stop("Selection model did not converge.", call. = FALSE)
}

eligible$selection_probability <- stats::predict(selection_model, type = "response")
eligible$selection_probability <- pmin(pmax(eligible$selection_probability, 0.01), 0.99)
selection_rate <- mean(eligible$selected == 1)
eligible$selection_ipw <- ifelse(
  eligible$selected == 1,
  selection_rate / eligible$selection_probability,
  NA_real_
)

selected_weights <- eligible$selection_ipw[eligible$selected == 1]
trim_limits <- stats::quantile(selected_weights, c(0.01, 0.99), na.rm = TRUE, names = FALSE)
eligible$selection_ipw_trimmed <- ifelse(
  eligible$selected == 1,
  pmin(pmax(eligible$selection_ipw, trim_limits[1]), trim_limits[2]),
  NA_real_
)

effective_sample_size <- function(w) sum(w, na.rm = TRUE)^2 / sum(w^2, na.rm = TRUE)
weight_diagnostics <- tibble::tibble(
  eligible_n = nrow(eligible),
  selected_n = sum(eligible$selected == 1),
  excluded_n = sum(eligible$selected == 0),
  selection_percent = 100 * mean(eligible$selected == 1),
  selected_probability_min = min(eligible$selection_probability[eligible$selected == 1]),
  selected_probability_median = stats::median(eligible$selection_probability[eligible$selected == 1]),
  selected_probability_max = max(eligible$selection_probability[eligible$selected == 1]),
  ipw_1st_percentile = trim_limits[1],
  ipw_99th_percentile = trim_limits[2],
  trimmed_ipw_max = max(eligible$selection_ipw_trimmed, na.rm = TRUE),
  trimmed_ipw_effective_sample_size = effective_sample_size(eligible$selection_ipw_trimmed)
)

model_data <- readRDS(benchmark_path)$model_data |>
  dplyr::left_join(
    eligible |>
      dplyr::select(SEQN, selection_probability, selection_ipw_trimmed),
    by = "SEQN"
  ) |>
  dplyr::mutate(
    phenotype_total_protein = factor(
      phenotype_total_protein, levels = c("P3", "P2", "P1")
    ),
    combined_selection_survey_weight = WTMEC8YR * selection_ipw_trimmed
  )

if (any(!is.finite(model_data$combined_selection_survey_weight))) {
  stop("Non-finite combined selection/survey weights detected.", call. = FALSE)
}

base_design <- survey::svydesign(
  ids = ~SDMVPSU,
  strata = ~SDMVSTRA,
  weights = ~WTMEC8YR,
  nest = TRUE,
  data = model_data
)
selection_design <- survey::svydesign(
  ids = ~SDMVPSU,
  strata = ~SDMVSTRA,
  weights = ~combined_selection_survey_weight,
  nest = TRUE,
  data = model_data
)

model_formula <- survival::Surv(PERMTH_INT, MORTSTAT) ~
  phenotype_total_protein + RIDAGEYR + male + race + INDFMPIR +
  Comorbidity_Score_Extended + cycle + smoking + hypertension + diabetes

fit_base <- survey::svycoxph(model_formula, design = base_design)
fit_selection <- survey::svycoxph(model_formula, design = selection_design)

extract_phenotype <- function(fit, model_label) {
  beta <- stats::coef(fit)
  standard_error <- sqrt(diag(stats::vcov(fit)))
  terms <- names(beta)
  keep <- grepl("^phenotype_total_protein", terms)
  tibble::tibble(
    model = model_label,
    comparison = dplyr::recode(
      terms[keep],
      phenotype_total_proteinP2 = "P2 vs P3",
      phenotype_total_proteinP1 = "P1 vs P3"
    ),
    HR = exp(beta[keep]),
    lower_95 = exp(beta[keep] - 1.96 * standard_error[keep]),
    upper_95 = exp(beta[keep] + 1.96 * standard_error[keep]),
    p_value = 2 * stats::pnorm(abs(beta[keep] / standard_error[keep]), lower.tail = FALSE)
  ) |>
    dplyr::mutate(
      hazard_ratio_95ci = sprintf("%.2f (%.2f-%.2f)", HR, lower_95, upper_95)
    )
}

cox_comparison <- dplyr::bind_rows(
  extract_phenotype(fit_base, "NHANES survey weights"),
  extract_phenotype(fit_selection, "NHANES survey weights × selection IPW")
)

readr::write_csv(comparison, file.path(output_dir, "Table33A_included_vs_excluded.csv"))
readr::write_csv(feature_missingness, file.path(output_dir, "Table33B_feature_missingness.csv"))
readr::write_csv(weight_diagnostics, file.path(output_dir, "Table33C_selection_weight_diagnostics.csv"))
readr::write_csv(cox_comparison, file.path(output_dir, "Table33D_selection_IPW_Cox.csv"))
saveRDS(
  list(
    comparison = comparison,
    feature_missingness = feature_missingness,
    weight_diagnostics = weight_diagnostics,
    cox_comparison = cox_comparison,
    selection_model = selection_model
  ),
  file.path(output_dir, "NHANES_selection_bias_results.rds")
)

summary_lines <- c(
  "NHANES complete-feature selection-bias sensitivity",
  "",
  "Selection-weight diagnostics:",
  paste(capture.output(print(weight_diagnostics)), collapse = "\n"),
  "",
  "Largest absolute standardised differences:",
  paste(
    capture.output(print(comparison |> dplyr::arrange(dplyr::desc(absolute_smd)) |> dplyr::slice_head(n = 8))),
    collapse = "\n"
  ),
  "",
  "Feature missingness:",
  paste(capture.output(print(feature_missingness)), collapse = "\n"),
  "",
  "Survey-weighted Cox comparison:",
  paste(capture.output(print(cox_comparison)), collapse = "\n"),
  "",
  paste0(
    "Interpretation: selection IPW addresses missingness in the six clustering features conditional on observed ",
    "demographics; it does not remove bias from unmeasured determinants of laboratory availability or from ",
    "complete-case restriction in later adjustment covariates."
  )
)
writeLines(summary_lines, file.path(output_dir, "NHANES_selection_bias_summary.txt"))
cat(paste(summary_lines, collapse = "\n"), "\n")
