# ==============================================================================
# Post-freeze clinical interpretation: adjusted absolute risk and RMST
# ==============================================================================

required_packages <- c("dplyr", "readr", "tibble", "tidyr", "survival", "survey")
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
input_path <- file.path(
  root, "output", "nhanes_albumin_benchmarks",
  "NHANES_albumin_benchmark_results.rds"
)
output_dir <- file.path(root, "output", "nhanes_adjusted_absolute_risk")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

n_boot <- as.integer(Sys.getenv("ABSOLUTE_RISK_BOOTSTRAP_REPS", unset = "1000"))
if (!is.finite(n_boot) || n_boot < 100L) {
  stop("ABSOLUTE_RISK_BOOTSTRAP_REPS must be at least 100.", call. = FALSE)
}

time_grid <- 0:60
risk_times <- c(36, 60)
phenotype_levels <- c("P3", "P2", "P1")

model_data <- readRDS(input_path)$model_data |>
  dplyr::mutate(
    phenotype_total_protein = factor(
      phenotype_total_protein, levels = phenotype_levels
    ),
    survey_weight_scaled = WTMEC8YR / mean(WTMEC8YR),
    survey_cluster = interaction(SDMVSTRA, SDMVPSU, drop = TRUE)
  )

full_formula <- survival::Surv(PERMTH_INT, MORTSTAT) ~
  phenotype_total_protein + RIDAGEYR + male + race + INDFMPIR +
  Comorbidity_Score_Extended + cycle + smoking + hypertension + diabetes

fit_weighted_cox <- function(data, analysis_weights, cluster = NULL) {
  data$.analysis_weight <- analysis_weights / mean(analysis_weights)
  if (is.null(cluster)) {
    survival::coxph(
      full_formula,
      data = data,
      weights = .analysis_weight,
      ties = "efron",
      robust = FALSE,
      model = TRUE
    )
  } else {
    data$.analysis_cluster <- cluster
    survival::coxph(
      full_formula,
      data = data,
      weights = .analysis_weight,
      ties = "efron",
      robust = TRUE,
      cluster = .analysis_cluster,
      model = TRUE
    )
  }
}

cumulative_hazard_at <- function(base_hazard, times) {
  vapply(times, function(current_time) {
    eligible <- which(base_hazard$time <= current_time)
    if (length(eligible) == 0L) return(0)
    base_hazard$hazard[max(eligible)]
  }, numeric(1))
}

standardise_survival <- function(fit, standardisation_data, standardisation_weights) {
  base_hazard <- survival::basehaz(fit, centered = TRUE)
  hazard_grid <- cumulative_hazard_at(base_hazard, time_grid)
  group_curves <- vector("list", length(phenotype_levels))
  group_metrics <- vector("list", length(phenotype_levels))

  for (i in seq_along(phenotype_levels)) {
    phenotype <- phenotype_levels[i]
    counterfactual <- standardisation_data
    counterfactual$phenotype_total_protein <- factor(
      phenotype, levels = phenotype_levels
    )
    linear_predictor <- stats::predict(
      fit, newdata = counterfactual, type = "lp", reference = "sample"
    )
    relative_hazard <- exp(linear_predictor)
    standardised_survival <- vapply(hazard_grid, function(hazard) {
      stats::weighted.mean(
        exp(-hazard * relative_hazard),
        w = standardisation_weights,
        na.rm = TRUE
      )
    }, numeric(1))

    group_curves[[i]] <- tibble::tibble(
      phenotype = phenotype,
      month = time_grid,
      adjusted_survival = standardised_survival,
      adjusted_risk = 1 - standardised_survival
    )
    rmst_60 <- sum(
      (standardised_survival[-length(standardised_survival)] +
         standardised_survival[-1]) / 2
    )
    group_metrics[[i]] <- tibble::tibble(
      phenotype = phenotype,
      adjusted_risk_36 = 1 - standardised_survival[time_grid == 36],
      adjusted_risk_60 = 1 - standardised_survival[time_grid == 60],
      adjusted_rmst_60_months = rmst_60
    )
  }

  metrics <- dplyr::bind_rows(group_metrics)
  contrasts <- dplyr::bind_rows(lapply(c("P2", "P1"), function(phenotype) {
    current <- metrics[metrics$phenotype == phenotype, , drop = FALSE]
    reference <- metrics[metrics$phenotype == "P3", , drop = FALSE]
    tibble::tibble(
      comparison = paste0(phenotype, " vs P3"),
      risk_difference_36_percentage_points =
        100 * (current$adjusted_risk_36 - reference$adjusted_risk_36),
      risk_ratio_36 = current$adjusted_risk_36 / reference$adjusted_risk_36,
      risk_difference_60_percentage_points =
        100 * (current$adjusted_risk_60 - reference$adjusted_risk_60),
      risk_ratio_60 = current$adjusted_risk_60 / reference$adjusted_risk_60,
      rmst_difference_60_months =
        current$adjusted_rmst_60_months - reference$adjusted_rmst_60_months
    )
  }))

  list(
    curves = dplyr::bind_rows(group_curves),
    metrics = metrics,
    contrasts = contrasts
  )
}

primary_fit <- fit_weighted_cox(
  model_data,
  model_data$WTMEC8YR,
  cluster = model_data$survey_cluster
)
primary_standardised <- standardise_survival(
  primary_fit,
  model_data,
  model_data$WTMEC8YR
)

primary_beta <- stats::coef(primary_fit)
weighted_cox_check <- tibble::tibble(
  comparison = c("P2 vs P3", "P1 vs P3"),
  HR = exp(primary_beta[c("phenotype_total_proteinP2", "phenotype_total_proteinP1")]),
  frozen_svycoxph_HR = c(0.968806754013275, 1.6830536847349)
) |>
  dplyr::mutate(absolute_difference = abs(HR - frozen_svycoxph_HR))

if (max(weighted_cox_check$absolute_difference) > 0.01) {
  stop("Weighted coxph coefficients do not reproduce the frozen survey Cox estimates.", call. = FALSE)
}

# Survey bootstrap replicate weights preserve the NHANES strata/PSU structure.
survey_design <- survey::svydesign(
  ids = ~SDMVPSU,
  strata = ~SDMVSTRA,
  weights = ~WTMEC8YR,
  nest = TRUE,
  data = model_data
)
set.seed(20260710)
replicate_design <- survey::as.svrepdesign(
  survey_design,
  type = "bootstrap",
  replicates = n_boot,
  mse = TRUE
)
replicate_weights <- stats::weights(replicate_design, type = "analysis")

replicate_results <- vector("list", n_boot)
replicate_group_results <- vector("list", n_boot)
failure_messages <- character()

for (b in seq_len(n_boot)) {
  current <- tryCatch(withCallingHandlers({
    weights_b <- replicate_weights[, b]
    keep <- is.finite(weights_b) & weights_b > 0
    data_b <- model_data[keep, , drop = FALSE]
    weights_b <- weights_b[keep]

    fit_b <- fit_weighted_cox(data_b, weights_b)
    if (any(!is.finite(stats::coef(fit_b)))) {
      stop("Non-finite Cox coefficient")
    }
    standardised_b <- standardise_survival(fit_b, data_b, weights_b)
    list(
      contrasts = standardised_b$contrasts |>
        dplyr::mutate(replicate = b, .before = 1),
      groups = standardised_b$metrics |>
        dplyr::mutate(replicate = b, .before = 1)
    )
  }, warning = function(w) {
    if (grepl("coefficient may be infinite", conditionMessage(w), fixed = TRUE)) {
      stop(conditionMessage(w), call. = FALSE)
    }
    invokeRestart("muffleWarning")
  }), error = function(e) {
    failure_messages <<- c(
      failure_messages, paste0("Replicate ", b, ": ", conditionMessage(e))
    )
    NULL
  })

  replicate_results[[b]] <- if (is.null(current)) NULL else current$contrasts
  replicate_group_results[[b]] <- if (is.null(current)) NULL else current$groups
  if (b %% 100L == 0L) {
    message("Completed ", b, "/", n_boot, " survey bootstrap replicates.")
  }
}

bootstrap_replicates <- dplyr::bind_rows(replicate_results)
bootstrap_group_replicates <- dplyr::bind_rows(replicate_group_results)
valid_by_comparison <- bootstrap_replicates |>
  dplyr::count(comparison, name = "valid_replicates")
if (
  nrow(valid_by_comparison) != 2L ||
    min(valid_by_comparison$valid_replicates) < 0.90 * n_boot
) {
  stop("Fewer than 90% valid survey bootstrap replicates.", call. = FALSE)
}

contrast_long <- primary_standardised$contrasts |>
  tidyr::pivot_longer(
    cols = -comparison,
    names_to = "metric",
    values_to = "estimate"
  )
bootstrap_long <- bootstrap_replicates |>
  tidyr::pivot_longer(
    cols = c(
      risk_difference_36_percentage_points,
      risk_ratio_36,
      risk_difference_60_percentage_points,
      risk_ratio_60,
      rmst_difference_60_months
    ),
    names_to = "metric",
    values_to = "bootstrap_estimate"
  )

contrast_summary <- bootstrap_long |>
  dplyr::group_by(comparison, metric) |>
  dplyr::summarise(
    valid_replicates = dplyr::n(),
    bootstrap_median = stats::median(bootstrap_estimate),
    lower_95 = stats::quantile(bootstrap_estimate, 0.025, names = FALSE),
    upper_95 = stats::quantile(bootstrap_estimate, 0.975, names = FALSE),
    .groups = "drop"
  ) |>
  dplyr::left_join(contrast_long, by = c("comparison", "metric")) |>
  dplyr::select(
    comparison, metric, estimate, lower_95, upper_95,
    bootstrap_median, valid_replicates
  )

group_metrics_point <- primary_standardised$metrics |>
  dplyr::mutate(
    adjusted_risk_36_percent = 100 * adjusted_risk_36,
    adjusted_risk_60_percent = 100 * adjusted_risk_60
  ) |>
  dplyr::select(
    phenotype,
    adjusted_risk_36_percent,
    adjusted_risk_60_percent,
    adjusted_rmst_60_months
  )

group_point_long <- group_metrics_point |>
  tidyr::pivot_longer(
    cols = -phenotype,
    names_to = "metric",
    values_to = "estimate"
  )
group_bootstrap_long <- bootstrap_group_replicates |>
  dplyr::mutate(
    adjusted_risk_36_percent = 100 * adjusted_risk_36,
    adjusted_risk_60_percent = 100 * adjusted_risk_60
  ) |>
  dplyr::select(
    replicate, phenotype, adjusted_risk_36_percent,
    adjusted_risk_60_percent, adjusted_rmst_60_months
  ) |>
  tidyr::pivot_longer(
    cols = -c(replicate, phenotype),
    names_to = "metric",
    values_to = "bootstrap_estimate"
  )
group_metrics_display <- group_bootstrap_long |>
  dplyr::group_by(phenotype, metric) |>
  dplyr::summarise(
    valid_replicates = dplyr::n(),
    lower_95 = stats::quantile(bootstrap_estimate, 0.025, names = FALSE),
    upper_95 = stats::quantile(bootstrap_estimate, 0.975, names = FALSE),
    .groups = "drop"
  ) |>
  dplyr::left_join(group_point_long, by = c("phenotype", "metric")) |>
  dplyr::select(phenotype, metric, estimate, lower_95, upper_95, valid_replicates) |>
  dplyr::mutate(
    estimate_95ci = sprintf("%.2f (%.2f-%.2f)", estimate, lower_95, upper_95)
  )

readr::write_csv(
  group_metrics_display,
  file.path(output_dir, "Table38A_adjusted_absolute_risk_and_RMST.csv")
)
readr::write_csv(
  contrast_summary,
  file.path(output_dir, "Table38B_adjusted_risk_and_RMST_contrasts.csv")
)
readr::write_csv(
  primary_standardised$curves,
  file.path(output_dir, "Table38C_monthly_adjusted_survival_curves.csv")
)
readr::write_csv(
  weighted_cox_check,
  file.path(output_dir, "Table38D_weighted_Cox_reproduction_check.csv")
)
readr::write_csv(
  bootstrap_replicates,
  file.path(output_dir, "Table38E_survey_bootstrap_replicates.csv")
)
readr::write_csv(
  bootstrap_group_replicates,
  file.path(output_dir, "Table38F_survey_bootstrap_group_estimates.csv")
)
if (length(failure_messages) > 0L) {
  writeLines(failure_messages, file.path(output_dir, "survey_bootstrap_failure_log.txt"))
}
saveRDS(
  list(
    fit = primary_fit,
    group_metrics = group_metrics_display,
    contrasts = contrast_summary,
    curves = primary_standardised$curves,
    bootstrap_replicates = bootstrap_replicates,
    bootstrap_group_replicates = bootstrap_group_replicates,
    weighted_cox_check = weighted_cox_check
  ),
  file.path(output_dir, "NHANES_adjusted_absolute_risk_RMST_results.rds")
)

summary_lines <- c(
  "NHANES post-freeze adjusted absolute risk and RMST analysis",
  "",
  paste0("Requested survey bootstrap replicates: ", n_boot),
  paste0("Failed replicates: ", length(failure_messages)),
  "",
  "Weighted Cox reproduction check:",
  paste(capture.output(print(weighted_cox_check)), collapse = "\n"),
  "",
  "Adjusted phenotype-specific estimates:",
  paste(capture.output(print(group_metrics_display)), collapse = "\n"),
  "",
  "Adjusted contrasts with survey bootstrap intervals:",
  paste(capture.output(print(contrast_summary)), collapse = "\n"),
  "",
  paste0(
    "Interpretation note: this post-freeze analysis improves clinical interpretation and does not replace ",
    "the frozen complex-survey Cox HR. RMST uses monthly trapezoidal integration through 60 months."
  )
)
writeLines(summary_lines, file.path(output_dir, "NHANES_adjusted_absolute_risk_RMST_summary.txt"))
cat(paste(summary_lines, collapse = "\n"), "\n")
