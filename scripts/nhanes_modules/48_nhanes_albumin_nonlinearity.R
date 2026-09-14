# ==============================================================================
# Script: 48_NHANES_nonlinearity_threshold_freeze_v2.R
# Purpose: Corrected-cohort nonlinearity and exploratory threshold analysis for
#          the NHANES domain-balanced vulnerability score.
#
# Safeguards:
#   - Corrected common model sample only (n = 3979; 720 deaths).
#   - Complex-survey Cox models for average shape and hinge scans.
#   - Fixed weighted percentile search range with BH correction.
#   - Leave-one-cycle-out threshold stability audit.
#   - Weighted robust time-varying spline curves because PH is violated.
#   - No clinical cutoff is claimed.
# ==============================================================================

required_packages <- c("survey", "survival", "splines", "dplyr", "readr", "tibble")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop("Missing required package(s): ", paste(missing_packages, collapse = ", "), call. = FALSE)
}

options(survey.lonely.psu = "adjust")

input_path <- paste0(
  "./output/nhanes_albumin_amendment17_2026-08-30/",
  "NHANES_albumin_benchmark_results.rds"
)
time_effect_path <- paste0(
  "./output/nhanes_albumin_amendment17_2026-08-30/cluster_bootstrap_PH/",
  "Table32E_time_varying_score_effects.csv"
)
time_test_path <- paste0(
  "./output/nhanes_albumin_amendment17_2026-08-30/cluster_bootstrap_PH/",
  "Table32F_time_varying_model_comparison.csv"
)
output_dir <- "./output/nhanes_albumin_amendment17_2026-08-30/nonlinearity"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

for (path in c(input_path, time_effect_path, time_test_path)) {
  if (!file.exists(path)) stop("Missing input: ", path, call. = FALSE)
}

analysis_vars <- c(
  "SEQN", "PERMTH_INT", "MORTSTAT", "domain_balanced_score",
  "RIDAGEYR", "male", "race", "INDFMPIR", "Comorbidity_Score_Extended",
  "cycle", "smoking", "hypertension", "diabetes",
  "WTMEC8YR", "SDMVPSU", "SDMVSTRA", "Cycle_ID"
)

data <- readRDS(input_path)[["model_data"]]
missing_vars <- setdiff(analysis_vars, names(data))
if (length(missing_vars) > 0L) {
  stop("Missing variable(s): ", paste(missing_vars, collapse = ", "), call. = FALSE)
}
data <- data[stats::complete.cases(data[, analysis_vars]), , drop = FALSE]
data <- data |>
  dplyr::mutate(
    psu_unique = interaction(SDMVSTRA, SDMVPSU, drop = TRUE),
    survey_weight_scaled = WTMEC8YR / mean(WTMEC8YR)
  )

stopifnot(nrow(data) == 3979L, sum(data$MORTSTAT == 1) == 720L)

base_terms <- paste(
  c(
    "RIDAGEYR", "male", "race", "INDFMPIR", "Comorbidity_Score_Extended",
    "cycle", "smoking", "hypertension", "diabetes"
  ),
  collapse = " + "
)

# SDMVSTRA absorbs survey-cycle baseline differences in the weighted robust
# time-varying model, so adding cycle again would create exact collinearity.
time_base_terms <- paste(
  c(
    "RIDAGEYR", "male", "race", "INDFMPIR", "Comorbidity_Score_Extended",
    "smoking", "hypertension", "diabetes"
  ),
  collapse = " + "
)

make_design <- function(df) {
  survey::svydesign(
    ids = ~SDMVPSU,
    strata = ~SDMVSTRA,
    weights = ~WTMEC8YR,
    nest = TRUE,
    data = df
  )
}

design_initial <- make_design(data)

weighted_quantiles <- function(design, probabilities) {
  as.numeric(stats::coef(survey::svyquantile(
    ~domain_balanced_score,
    design,
    quantiles = probabilities,
    ci = FALSE
  )))
}

knot_probabilities <- c(0.05, 0.275, 0.50, 0.725, 0.95)
knots_all <- weighted_quantiles(design_initial, knot_probabilities)
boundary_knots <- knots_all[c(1, 5)]
inner_knots <- knots_all[2:4]

spline_basis <- splines::ns(
  data$domain_balanced_score,
  knots = inner_knots,
  Boundary.knots = boundary_knots
)
colnames(spline_basis) <- paste0("score_ns", seq_len(ncol(spline_basis)))
data <- cbind(data, spline_basis)
design <- make_design(data)

linear_formula <- stats::as.formula(paste(
  "survival::Surv(PERMTH_INT, MORTSTAT) ~ domain_balanced_score +", base_terms
))
spline_formula <- stats::as.formula(paste(
  "survival::Surv(PERMTH_INT, MORTSTAT) ~",
  paste(colnames(spline_basis), collapse = " + "), "+", base_terms
))

fit_linear <- survey::svycoxph(linear_formula, design = design)
fit_spline <- survey::svycoxph(spline_formula, design = design)

nonlinearity_test <- stats::anova(
  fit_linear,
  fit_spline,
  test = "F",
  method = "Wald",
  force = TRUE
)

nonlinearity_summary <- tibble::tibble(
  analysis = "Five-knot natural spline versus linear score",
  n = nrow(data),
  events = sum(data$MORTSTAT == 1),
  weighted_knot_probabilities = paste(knot_probabilities, collapse = ";"),
  weighted_knots = paste(sprintf("%.4f", knots_all), collapse = ";"),
  numerator_df = as.numeric(nonlinearity_test$df),
  denominator_df = as.numeric(nonlinearity_test$ddf),
  F_statistic = as.numeric(nonlinearity_test$Ftest),
  p_nonlinearity = as.numeric(nonlinearity_test$p),
  p_nonlinearity_formatted = ifelse(
    as.numeric(nonlinearity_test$p) < 0.001,
    "<0.001",
    sprintf("%.3f", as.numeric(nonlinearity_test$p))
  ),
  interpretation = "Design-based average hazard shape; PH-violating score requires time-varying sensitivity"
)

relative_spline_curve <- function(fit, grid, reference, inner, boundary) {
  grid_basis <- splines::ns(grid, knots = inner, Boundary.knots = boundary)
  ref_basis <- splines::ns(reference, knots = inner, Boundary.knots = boundary)
  colnames(grid_basis) <- paste0("score_ns", seq_len(ncol(grid_basis)))
  colnames(ref_basis) <- colnames(grid_basis)
  difference <- sweep(grid_basis, 2, as.numeric(ref_basis[1, ]), "-")

  coefficients <- stats::coef(fit)[colnames(grid_basis)]
  covariance <- as.matrix(stats::vcov(fit))[colnames(grid_basis), colnames(grid_basis), drop = FALSE]
  log_effect <- as.numeric(difference %*% coefficients)
  se <- sqrt(rowSums((difference %*% covariance) * difference))

  tibble::tibble(
    score = grid,
    reference_score = reference,
    HR = exp(log_effect),
    lower_95 = exp(log_effect - 1.96 * se),
    upper_95 = exp(log_effect + 1.96 * se),
    curve_type = "Complex-survey average-hazard natural spline"
  )
}

curve_limits <- weighted_quantiles(design, c(0.01, 0.99))
curve_grid <- seq(curve_limits[1], curve_limits[2], length.out = 161L)
average_curve <- relative_spline_curve(
  fit_spline, curve_grid, reference = 0,
  inner = inner_knots, boundary = boundary_knots
)

extract_regtest <- function(test) {
  tibble::tibble(
    F_statistic = as.numeric(test$Ftest),
    numerator_df = as.numeric(test$df),
    denominator_df = as.numeric(test$ddf),
    p_raw = as.numeric(test$p)
  )
}

combined_slope <- function(fit) {
  coefficients <- stats::coef(fit)
  covariance <- as.matrix(stats::vcov(fit))
  below_beta <- coefficients[["domain_balanced_score"]]
  hinge_beta <- coefficients[["score_hinge"]]
  below_variance <- covariance["domain_balanced_score", "domain_balanced_score"]
  above_variance <-
    covariance["domain_balanced_score", "domain_balanced_score"] +
    covariance["score_hinge", "score_hinge"] +
    2 * covariance["domain_balanced_score", "score_hinge"]

  tibble::tibble(
    segment = c("Below candidate", "Above candidate", "Slope change"),
    estimate = exp(c(below_beta, below_beta + hinge_beta, hinge_beta)),
    lower_95 = exp(c(
      below_beta - 1.96 * sqrt(below_variance),
      below_beta + hinge_beta - 1.96 * sqrt(above_variance),
      hinge_beta - 1.96 * sqrt(covariance["score_hinge", "score_hinge"])
    )),
    upper_95 = exp(c(
      below_beta + 1.96 * sqrt(below_variance),
      below_beta + hinge_beta + 1.96 * sqrt(above_variance),
      hinge_beta + 1.96 * sqrt(covariance["score_hinge", "score_hinge"])
    ))
  ) |>
    dplyr::mutate(effect_95ci = sprintf("%.2f (%.2f-%.2f)", estimate, lower_95, upper_95))
}

scan_thresholds <- function(df, candidates, analysis_label) {
  rows <- vector("list", length(candidates))
  fits <- vector("list", length(candidates))

  for (i in seq_along(candidates)) {
    candidate <- candidates[[i]]
    message(
      "Threshold scan: ", analysis_label,
      ", candidate ", i, "/", length(candidates),
      ", score=", sprintf("%.4f", candidate)
    )
    temporary <- df |>
      dplyr::mutate(
        score_hinge = pmax(domain_balanced_score - candidate, 0),
        above_candidate = as.integer(domain_balanced_score > candidate)
      )
    current_design <- make_design(temporary)
    formula <- stats::as.formula(paste(
      "survival::Surv(PERMTH_INT, MORTSTAT) ~ domain_balanced_score + score_hinge +",
      base_terms
    ))
    fit <- survey::svycoxph(formula, design = current_design)
    test <- survey::regTermTest(fit, ~score_hinge, method = "Wald")
    test_values <- extract_regtest(test)
    weighted_above <- as.numeric(stats::coef(survey::svymean(~above_candidate, current_design)))

    rows[[i]] <- tibble::tibble(
      analysis = analysis_label,
      candidate_percentile = seq(20, 80, by = 5)[[i]],
      candidate_score = candidate,
      n_below_or_equal = sum(temporary$domain_balanced_score <= candidate),
      n_above = sum(temporary$domain_balanced_score > candidate),
      events_below_or_equal = sum(
        temporary$MORTSTAT[temporary$domain_balanced_score <= candidate] == 1
      ),
      events_above = sum(temporary$MORTSTAT[temporary$domain_balanced_score > candidate] == 1),
      weighted_percent_above = 100 * weighted_above,
      F_statistic = test_values$F_statistic,
      numerator_df = test_values$numerator_df,
      denominator_df = test_values$denominator_df,
      p_raw = test_values$p_raw
    )
    fits[[i]] <- fit
  }

  table <- dplyr::bind_rows(rows) |>
    dplyr::mutate(p_BH = stats::p.adjust(p_raw, method = "BH"))
  best_index <- which.min(table$p_raw)
  list(table = table, best_fit = fits[[best_index]], best_index = best_index)
}

candidate_probabilities <- seq(0.20, 0.80, by = 0.05)
candidate_scores <- weighted_quantiles(design, candidate_probabilities)
full_scan <- scan_thresholds(data, candidate_scores, "Full corrected cohort")
threshold_scan <- full_scan$table
best_threshold <- threshold_scan[full_scan$best_index, , drop = FALSE]
best_slopes <- combined_slope(full_scan$best_fit) |>
  dplyr::mutate(
    candidate_score = best_threshold$candidate_score,
    candidate_percentile = best_threshold$candidate_percentile,
    hinge_p_raw = best_threshold$p_raw,
    hinge_p_BH = best_threshold$p_BH
  )

cycle_values <- sort(unique(as.character(data$Cycle_ID)))
loco_rows <- vector("list", length(cycle_values))
for (i in seq_along(cycle_values)) {
  held_out <- cycle_values[[i]]
  training <- data |>
    dplyr::filter(as.character(Cycle_ID) != held_out) |>
    droplevels()
  current <- scan_thresholds(
    training,
    candidate_scores,
    paste0("Leave out cycle ", held_out)
  )$table
  best <- current[which.min(current$p_raw), , drop = FALSE]
  loco_rows[[i]] <- best |>
    dplyr::mutate(
      held_out_cycle = held_out,
      n_training = nrow(training),
      events_training = sum(training$MORTSTAT == 1)
    ) |>
    dplyr::select(
      held_out_cycle, n_training, events_training,
      candidate_percentile, candidate_score, p_raw, p_BH
    )
}
loco_thresholds <- dplyr::bind_rows(loco_rows)

loco_range <- diff(range(loco_thresholds$candidate_score))
loco_significant_count <- sum(loco_thresholds$p_BH < 0.05)
fixed_threshold_supported <-
  nonlinearity_summary$p_nonlinearity < 0.05 &&
  best_threshold$p_BH < 0.05 &&
  loco_range <= 0.50 &&
  loco_significant_count >= 3L

threshold_conclusion <- if (fixed_threshold_supported) {
  paste0(
    "A statistically supported but exploratory average-hazard hinge candidate was identified at score ",
    sprintf("%.2f", best_threshold$candidate_score),
    "; it is not a validated clinical cutoff because the score violates PH."
  )
} else {
  paste0(
    "No stable fixed threshold was supported after design-based nonlinearity testing, BH correction, ",
    "leave-one-cycle-out stability, and recognition of the PH violation."
  )
}

threshold_decision <- tibble::tibble(
  best_candidate_score = best_threshold$candidate_score,
  best_candidate_percentile = best_threshold$candidate_percentile,
  best_hinge_p_raw = best_threshold$p_raw,
  best_hinge_p_BH = best_threshold$p_BH,
  spline_p_nonlinearity = nonlinearity_summary$p_nonlinearity,
  loco_selected_score_min = min(loco_thresholds$candidate_score),
  loco_selected_score_max = max(loco_thresholds$candidate_score),
  loco_selected_score_range = loco_range,
  loco_BH_significant_count = loco_significant_count,
  loco_total = nrow(loco_thresholds),
  fixed_threshold_supported = fixed_threshold_supported,
  clinical_cutoff_validated = FALSE,
  conclusion = threshold_conclusion
)

# Weighted robust time-varying spline sensitivity. This is not a replacement
# for the complex-survey model; it shows how the spline shape changes over time.
time_varying_formula <- stats::as.formula(paste(
  "survival::Surv(PERMTH_INT, MORTSTAT) ~",
  paste(colnames(spline_basis), collapse = " + "), "+",
  paste(paste0("tt(", colnames(spline_basis), ")"), collapse = " + "), "+",
  time_base_terms, "+ strata(SDMVSTRA) + cluster(psu_unique)"
))
fit_time_varying_spline <- survival::coxph(
  time_varying_formula,
  data = data,
  weights = survey_weight_scaled,
  ties = "efron",
  tt = function(x, t, ...) x * log(pmax(t, 1))
)

time_specific_curve <- function(fit, grid, reference, times, inner, boundary) {
  grid_basis <- splines::ns(grid, knots = inner, Boundary.knots = boundary)
  ref_basis <- splines::ns(reference, knots = inner, Boundary.knots = boundary)
  column_names <- paste0("score_ns", seq_len(ncol(grid_basis)))
  colnames(grid_basis) <- column_names
  colnames(ref_basis) <- column_names
  difference <- sweep(grid_basis, 2, as.numeric(ref_basis[1, ]), "-")

  coefficients <- stats::coef(fit)
  covariance <- as.matrix(stats::vcov(fit))
  rows <- vector("list", length(times))
  for (i in seq_along(times)) {
    month <- times[[i]]
    time_names <- paste0("tt(", column_names, ")")
    score_terms <- c(column_names, time_names)
    contrast <- matrix(
      0,
      nrow = nrow(difference),
      ncol = length(score_terms),
      dimnames = list(NULL, score_terms)
    )
    contrast[, column_names] <- difference
    contrast[, time_names] <- difference * log(month)
    score_coefficients <- coefficients[score_terms]
    score_covariance <- covariance[score_terms, score_terms, drop = FALSE]
    log_effect <- as.numeric(contrast %*% score_coefficients)
    se <- sqrt(rowSums((contrast %*% score_covariance) * contrast))
    rows[[i]] <- tibble::tibble(
      followup_month = month,
      score = grid,
      reference_score = reference,
      HR = exp(log_effect),
      lower_95 = exp(log_effect - 1.96 * se),
      upper_95 = exp(log_effect + 1.96 * se),
      curve_type = "Weighted robust log-time-varying natural spline"
    )
  }
  dplyr::bind_rows(rows)
}

evaluation_months <- c(12, 36, 60, 96)
time_specific_curves <- time_specific_curve(
  fit_time_varying_spline,
  curve_grid,
  reference = 0,
  times = evaluation_months,
  inner = inner_knots,
  boundary = boundary_knots
)

linear_time_effects <- readr::read_csv(time_effect_path, show_col_types = FALSE)
linear_time_test <- readr::read_csv(time_test_path, show_col_types = FALSE)

qa <- tibble::tribble(
  ~check, ~passed, ~detail,
  "Corrected common model n", nrow(data) == 3979L, paste0("n=", nrow(data)),
  "Corrected common model events", sum(data$MORTSTAT == 1) == 720L,
    paste0("events=", sum(data$MORTSTAT == 1)),
  "Five weighted spline knots ordered", length(knots_all) == 5L && all(diff(knots_all) > 0),
    paste(sprintf("%.3f", knots_all), collapse = ","),
  "Design-based nonlinearity test finite", is.finite(nonlinearity_summary$p_nonlinearity),
    paste0("p=", signif(nonlinearity_summary$p_nonlinearity, 4)),
  "Thirteen fixed threshold candidates", nrow(threshold_scan) == 13L,
    paste0("rows=", nrow(threshold_scan)),
  "Threshold multiplicity correction valid", all(
    is.finite(threshold_scan$p_BH) & threshold_scan$p_BH >= 0 & threshold_scan$p_BH <= 1
  ), "BH-adjusted p values in [0,1]",
  "Threshold sides retain events", all(
    threshold_scan$events_below_or_equal >= 50 & threshold_scan$events_above >= 50
  ), "At least 50 events on both sides for every candidate",
  "Four leave-one-cycle-out audits", nrow(loco_thresholds) == 4L,
    paste0("rows=", nrow(loco_thresholds)),
  "Four time-specific spline curves", setequal(
    unique(time_specific_curves$followup_month), evaluation_months
  ), paste(evaluation_months, collapse = ","),
  "Time-varying spline coefficients finite",
    all(is.finite(stats::coef(fit_time_varying_spline))),
    paste0("coefficients=", length(stats::coef(fit_time_varying_spline))),
  "Time-specific spline curve finite",
    all(is.finite(as.matrix(
      time_specific_curves[, c("HR", "lower_95", "upper_95")]
    ))),
    paste0("curve rows=", nrow(time_specific_curves)),
  "No validated clinical cutoff claim", !threshold_decision$clinical_cutoff_validated,
    threshold_decision$conclusion
)

readr::write_csv(
  nonlinearity_summary,
  file.path(output_dir, "Table48A_survey_nonlinearity_test.csv")
)
readr::write_csv(
  threshold_scan,
  file.path(output_dir, "Table48B_hinge_threshold_scan.csv")
)
readr::write_csv(
  best_slopes,
  file.path(output_dir, "Table48C_best_candidate_piecewise_slopes.csv")
)
readr::write_csv(
  loco_thresholds,
  file.path(output_dir, "Table48D_leave_one_cycle_threshold_stability.csv")
)
readr::write_csv(
  threshold_decision,
  file.path(output_dir, "Table48E_threshold_decision.csv")
)
readr::write_csv(
  linear_time_effects,
  file.path(output_dir, "Table48F_time_varying_linear_effects.csv")
)
readr::write_csv(
  linear_time_test,
  file.path(output_dir, "Table48G_time_varying_linear_test.csv")
)
readr::write_csv(
  average_curve,
  file.path(output_dir, "FigureS2A_average_spline_curve_data.csv")
)
readr::write_csv(
  time_specific_curves,
  file.path(output_dir, "FigureS2B_time_specific_spline_curve_data.csv")
)
readr::write_csv(qa, file.path(output_dir, "Table48H_QA.csv"))

saveRDS(
  list(
    fit_linear = fit_linear,
    fit_spline = fit_spline,
    fit_time_varying_spline = fit_time_varying_spline,
    nonlinearity_summary = nonlinearity_summary,
    threshold_scan = threshold_scan,
    best_slopes = best_slopes,
    loco_thresholds = loco_thresholds,
    threshold_decision = threshold_decision,
    average_curve = average_curve,
    time_specific_curves = time_specific_curves,
    qa = qa,
    knots = list(probabilities = knot_probabilities, values = knots_all)
  ),
  file.path(output_dir, "NHANES_nonlinearity_threshold_freeze_v2_results.rds")
)

summary_lines <- c(
  "Freeze V2 NHANES nonlinearity and threshold analysis",
  "",
  paste0("QA checks passed: ", sum(qa$passed), "/", nrow(qa)),
  "",
  paste0(
    "Design-based spline nonlinearity: F=",
    sprintf("%.3f", nonlinearity_summary$F_statistic),
    ", p=", signif(nonlinearity_summary$p_nonlinearity, 4)
  ),
  paste0(
    "Best exploratory hinge: score=", sprintf("%.3f", best_threshold$candidate_score),
    " (weighted percentile ", best_threshold$candidate_percentile, "), raw p=",
    signif(best_threshold$p_raw, 4), ", BH p=", signif(best_threshold$p_BH, 4)
  ),
  paste0(
    "Leave-one-cycle selected score range: ",
    sprintf("%.3f", min(loco_thresholds$candidate_score)), " to ",
    sprintf("%.3f", max(loco_thresholds$candidate_score)),
    "; BH-significant cycles=", loco_significant_count, "/4"
  ),
  paste0("Fixed threshold supported: ", fixed_threshold_supported),
  threshold_conclusion,
  "",
  "Known time-varying linear score effects:",
  paste(capture.output(print(linear_time_effects, n = Inf)), collapse = "\n"),
  "",
  "Interpretation boundary:",
  "Spline and hinge results describe exploratory average-hazard shape.",
  "Time-specific curves must be considered because the score violates PH.",
  "No result validates a clinical cutoff."
)
writeLines(
  summary_lines,
  file.path(output_dir, "NHANES_nonlinearity_threshold_freeze_v2_summary.txt")
)

if (!all(qa$passed)) {
  stop("Nonlinearity/threshold QA failed. Review Table48H_QA.csv.", call. = FALSE)
}

cat(paste(summary_lines, collapse = "\n"), "\n")
