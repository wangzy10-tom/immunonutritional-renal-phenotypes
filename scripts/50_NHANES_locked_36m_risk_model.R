# ==============================================================================
# Locked NHANES 36-month risk model used in exploratory Table S19.
#
# This model uses NHANES data only. Its seven-category comorbidity definition is
# chosen because the same categories are available in the historical NHANES
# transport period. It is analytically separate from the phenotype analysis.
# ==============================================================================

required_packages <- c("glmnet", "survey", "dplyr", "readr", "tibble")
missing_packages <- required_packages[!vapply(
  required_packages, requireNamespace, logical(1), quietly = TRUE
)]
if (length(missing_packages) > 0L) {
  stop("Missing required package(s): ", paste(missing_packages, collapse = ", "), call. = FALSE)
}

set.seed(20260710)
options(survey.lonely.psu = "adjust")

root <- normalizePath(
  Sys.getenv("PROJECT_ROOT", unset = getwd()), winslash = "/", mustWork = TRUE
)
input_path <- file.path(
  root, "output", "nhanes_albumin_benchmarks",
  "NHANES_albumin_benchmark_results.rds"
)
output_dir <- file.path(root, "output", "nhanes_locked_36m_risk_model")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

horizon_months <- 36L
bootstrap_repetitions <- 1000L

raw_data <- readRDS(input_path)[["model_data"]]
required_vars <- c(
  "SEQN", "Cycle_ID", "MORTSTAT", "PERMTH_INT", "WTMEC8YR",
  "SDMVSTRA", "SDMVPSU", "RIDAGEYR", "male",
  "MCQ220", "MCQ160A", "MCQ160B", "MCQ160C", "MCQ160D", "MCQ160E",
  "MCQ160F", "MCQ160G", "MCQ160K", "BPQ020", "DIQ010",
  "NLR", "SII", "LBXHGB", "LBXSTP", "BMXBMI", "LBXSCR"
)
missing_vars <- setdiff(required_vars, names(raw_data))
if (length(missing_vars) > 0L) {
  stop("Missing variable(s): ", paste(missing_vars, collapse = ", "), call. = FALSE)
}

yes_no_flag <- function(x) {
  dplyr::case_when(x == 1 ~ 1, x == 2 ~ 0, TRUE ~ NA_real_)
}

collapse_positive <- function(...) {
  values <- cbind(...)
  apply(values, 1L, function(row) {
    if (any(row == 1, na.rm = TRUE)) return(1)
    if (all(row == 0)) return(0)
    NA_real_
  })
}

data <- raw_data |>
  dplyr::select(dplyr::all_of(required_vars)) |>
  dplyr::filter(Cycle_ID %in% c("G", "H", "I")) |>
  dplyr::filter(stats::complete.cases(dplyr::across(dplyr::all_of(required_vars)))) |>
  dplyr::mutate(
    cancer_flag = yes_no_flag(MCQ220),
    arthritis_flag = yes_no_flag(MCQ160A),
    heart_disease_flag = collapse_positive(
      yes_no_flag(MCQ160B), yes_no_flag(MCQ160C),
      yes_no_flag(MCQ160D), yes_no_flag(MCQ160E)
    ),
    stroke_flag = yes_no_flag(MCQ160F),
    chronic_lung_flag = collapse_positive(
      yes_no_flag(MCQ160G), yes_no_flag(MCQ160K)
    ),
    hypertension_flag = yes_no_flag(BPQ020),
    diabetes_flag = dplyr::case_when(
      DIQ010 == 1 ~ 1,
      DIQ010 %in% c(2, 3) ~ 0,
      TRUE ~ NA_real_
    ),
    common_comorbidity_count = cancer_flag + arthritis_flag +
      heart_disease_flag + stroke_flag + chronic_lung_flag +
      hypertension_flag + diabetes_flag,
    event_36m = as.integer(MORTSTAT == 1 & PERMTH_INT <= horizon_months),
    horizon_observed = event_36m == 1 | PERMTH_INT >= horizon_months,
    strata_unique = interaction(Cycle_ID, SDMVSTRA, drop = TRUE),
    psu_unique = interaction(Cycle_ID, SDMVSTRA, SDMVPSU, drop = TRUE)
  )

development <- data |>
  dplyr::filter(Cycle_ID %in% c("G", "H"), horizon_observed) |>
  droplevels()
validation <- data |>
  dplyr::filter(Cycle_ID == "I", horizon_observed) |>
  droplevels()

stopifnot(
  nrow(development) == 1954L,
  sum(development$event_36m) == 177L,
  nrow(validation) == 988L,
  sum(validation$event_36m) == 93L,
  all(development$horizon_observed),
  all(validation$horizon_observed)
)

weighted_quantile <- function(x, w, probabilities) {
  ordering <- order(x)
  x <- x[ordering]
  w <- w[ordering]
  cumulative <- cumsum(w) / sum(w)
  vapply(probabilities, function(probability) {
    x[which(cumulative >= probability)[1L]]
  }, numeric(1))
}

biomarker_vars <- c("NLR", "SII", "LBXHGB", "LBXSTP", "BMXBMI", "LBXSCR")
winsor_bounds <- lapply(biomarker_vars, function(variable) {
  weighted_quantile(
    development[[variable]], development$WTMEC8YR, c(0.01, 0.99)
  )
})
names(winsor_bounds) <- biomarker_vars

winsorise <- function(x, bounds) {
  pmin(pmax(x, bounds[[1L]]), bounds[[2L]])
}

prepare_features <- function(df) {
  nlr <- winsorise(df$NLR, winsor_bounds$NLR)
  sii <- winsorise(df$SII, winsor_bounds$SII)
  haemoglobin <- winsorise(df$LBXHGB, winsor_bounds$LBXHGB)
  total_protein <- winsorise(df$LBXSTP, winsor_bounds$LBXSTP)
  bmi <- winsorise(df$BMXBMI, winsor_bounds$BMXBMI)
  creatinine <- winsorise(df$LBXSCR, winsor_bounds$LBXSCR)

  dplyr::mutate(
    df,
    age_5y = (RIDAGEYR - 75) / 5,
    comorbidity_count = common_comorbidity_count,
    log_nlr = log(nlr),
    log_sii = log(sii),
    haemoglobin = haemoglobin,
    total_protein = total_protein,
    bmi_5 = (bmi - 27) / 5,
    bmi_5_sq = bmi_5^2,
    log_creatinine = log(creatinine)
  )
}

development <- prepare_features(development)
validation <- prepare_features(validation)

base_features <- c("age_5y", "male", "comorbidity_count")
full_features <- c(
  base_features, "log_nlr", "log_sii", "haemoglobin", "total_protein",
  "bmi_5", "bmi_5_sq", "log_creatinine"
)

matrix_from <- function(df, features) {
  as.matrix(df[, features, drop = FALSE])
}

make_cluster_folds <- function(df, folds = 5L) {
  cluster_summary <- df |>
    dplyr::group_by(psu_unique) |>
    dplyr::summarise(
      weighted_events = sum(WTMEC8YR * event_36m),
      .groups = "drop"
    ) |>
    dplyr::arrange(dplyr::desc(weighted_events)) |>
    dplyr::mutate(fold = rep(seq_len(folds), length.out = dplyr::n()))
  fold_map <- stats::setNames(cluster_summary$fold, as.character(cluster_summary$psu_unique))
  unname(fold_map[as.character(df$psu_unique)])
}

fit_ridge_model <- function(df, features, fold_id = NULL, fixed_lambda = NULL) {
  x <- matrix_from(df, features)
  y <- df$event_36m
  weights <- df$WTMEC8YR / mean(df$WTMEC8YR)

  if (is.null(fixed_lambda)) {
    cv_fit <- glmnet::cv.glmnet(
      x = x,
      y = y,
      weights = weights,
      family = "binomial",
      alpha = 0,
      foldid = fold_id,
      type.measure = "deviance",
      standardize = TRUE,
      intercept = TRUE
    )
    # Prediction, rather than variable selection, is the objective. Use the
    # development-only lambda that minimises cross-validated deviance.
    lambda <- cv_fit$lambda.min
  } else {
    cv_fit <- NULL
    lambda <- fixed_lambda
  }

  fit <- glmnet::glmnet(
    x = x,
    y = y,
    weights = weights,
    family = "binomial",
    alpha = 0,
    lambda = lambda,
    standardize = TRUE,
    intercept = TRUE
  )
  list(fit = fit, lambda = lambda, cv_fit = cv_fit, features = features)
}

predict_ridge <- function(model, df) {
  as.numeric(stats::predict(
    model$fit,
    newx = matrix_from(df, model$features),
    s = model$lambda,
    type = "response"
  ))
}

fold_id <- make_cluster_folds(development, folds = 5L)
base_model <- fit_ridge_model(development, base_features, fold_id = fold_id)
full_model <- fit_ridge_model(development, full_features, fold_id = fold_id)

development$pred_base <- predict_ridge(base_model, development)
development$pred_full <- predict_ridge(full_model, development)
validation$pred_base <- predict_ridge(base_model, validation)
validation$pred_full <- predict_ridge(full_model, validation)

weighted_auc <- function(outcome, prediction, weight) {
  case_index <- outcome == 1L
  control_index <- outcome == 0L
  case_prediction <- prediction[case_index]
  control_prediction <- prediction[control_index]
  case_weight <- weight[case_index]
  control_weight <- weight[control_index]
  comparisons <- outer(case_prediction, control_prediction, "-")
  pair_weights <- outer(case_weight, control_weight, "*")
  sum(pair_weights * ((comparisons > 0) + 0.5 * (comparisons == 0))) /
    sum(pair_weights)
}

weighted_brier <- function(outcome, prediction, weight) {
  stats::weighted.mean((outcome - prediction)^2, weight)
}

threshold_metrics <- function(outcome, prediction, weight, threshold) {
  classified_high <- prediction >= threshold
  case_weight <- sum(weight[outcome == 1L])
  control_weight <- sum(weight[outcome == 0L])
  true_positive <- sum(weight[outcome == 1L & classified_high])
  false_negative <- sum(weight[outcome == 1L & !classified_high])
  true_negative <- sum(weight[outcome == 0L & !classified_high])
  false_positive <- sum(weight[outcome == 0L & classified_high])

  tibble::tibble(
    threshold = threshold,
    sensitivity = true_positive / case_weight,
    specificity = true_negative / control_weight,
    ppv = true_positive / (true_positive + false_positive),
    npv = true_negative / (true_negative + false_negative),
    classified_high_percent = 100 * sum(weight[classified_high]) / sum(weight),
    youden = sensitivity + specificity - 1
  )
}

select_threshold <- function(outcome, prediction, weight) {
  candidates <- unique(as.numeric(stats::quantile(
    prediction, probs = seq(0.02, 0.98, length.out = 301L), na.rm = TRUE
  )))
  scan <- dplyr::bind_rows(lapply(candidates, function(candidate) {
    threshold_metrics(outcome, prediction, weight, candidate)
  }))
  scan[which.max(scan$youden), , drop = FALSE]
}

locked_threshold <- select_threshold(
  development$event_36m, development$pred_full, development$WTMEC8YR
)
development_threshold_metrics <- threshold_metrics(
  development$event_36m, development$pred_full, development$WTMEC8YR,
  locked_threshold$threshold
)
validation_threshold_metrics <- threshold_metrics(
  validation$event_36m, validation$pred_full, validation$WTMEC8YR,
  locked_threshold$threshold
)

make_design <- function(df) {
  survey::svydesign(
    ids = ~psu_unique,
    strata = ~strata_unique,
    weights = ~WTMEC8YR,
    nest = TRUE,
    data = df
  )
}

validation$lp_full <- stats::qlogis(pmin(pmax(validation$pred_full, 1e-6), 1 - 1e-6))
validation_design <- make_design(validation)
calibration_intercept_fit <- survey::svyglm(
  event_36m ~ 1,
  offset = lp_full,
  design = validation_design,
  family = stats::quasibinomial()
)
calibration_slope_fit <- survey::svyglm(
  event_36m ~ lp_full,
  design = validation_design,
  family = stats::quasibinomial()
)

performance_row <- function(df, prediction_name, cohort, model_name) {
  prediction <- df[[prediction_name]]
  prevalence <- stats::weighted.mean(df$event_36m, df$WTMEC8YR)
  tibble::tibble(
    cohort = cohort,
    model = model_name,
    n = nrow(df),
    events = sum(df$event_36m),
    weighted_event_rate = prevalence,
    AUC = weighted_auc(df$event_36m, prediction, df$WTMEC8YR),
    Brier = weighted_brier(df$event_36m, prediction, df$WTMEC8YR),
    null_Brier = prevalence * (1 - prevalence)
  )
}

performance <- dplyr::bind_rows(
  performance_row(development, "pred_base", "Development G/H", "Clinical base"),
  performance_row(development, "pred_full", "Development G/H", "Full bedside"),
  performance_row(validation, "pred_base", "Temporal validation I", "Clinical base"),
  performance_row(validation, "pred_full", "Temporal validation I", "Full bedside")
)

decision_curve <- function(df, prediction_name, thresholds) {
  prediction <- df[[prediction_name]]
  outcome <- df$event_36m
  weight <- df$WTMEC8YR
  prevalence <- stats::weighted.mean(outcome, weight)
  dplyr::bind_rows(lapply(thresholds, function(threshold) {
    high <- prediction >= threshold
    true_positive_rate_population <- sum(weight[outcome == 1L & high]) / sum(weight)
    false_positive_rate_population <- sum(weight[outcome == 0L & high]) / sum(weight)
    tibble::tibble(
      threshold_probability = threshold,
      net_benefit = true_positive_rate_population -
        false_positive_rate_population * threshold / (1 - threshold),
      net_benefit_all = prevalence - (1 - prevalence) * threshold / (1 - threshold),
      net_benefit_none = 0
    )
  }))
}

dca_thresholds <- seq(0.02, 0.25, by = 0.01)
dca <- dplyr::bind_rows(
  dplyr::mutate(
    decision_curve(validation, "pred_base", dca_thresholds),
    model = "Clinical base"
  ),
  dplyr::mutate(
    decision_curve(validation, "pred_full", dca_thresholds),
    model = "Full bedside"
  )
)

cluster_bootstrap <- function(df) {
  strata_values <- unique(as.character(df$strata_unique))
  pieces <- list()
  piece_index <- 1L
  for (stratum in strata_values) {
    stratum_data <- df[as.character(df$strata_unique) == stratum, , drop = FALSE]
    units <- unique(as.character(stratum_data$psu_unique))
    sampled_units <- sample(units, length(units), replace = TRUE)
    for (draw in seq_along(sampled_units)) {
      piece <- stratum_data[
        as.character(stratum_data$psu_unique) == sampled_units[[draw]],
        , drop = FALSE
      ]
      piece$bootstrap_psu <- paste0(stratum, "_", draw)
      pieces[[piece_index]] <- piece
      piece_index <- piece_index + 1L
    }
  }
  dplyr::bind_rows(pieces)
}

message("Running ", bootstrap_repetitions, " development bootstrap repetitions...")
development_bootstrap <- vector("list", bootstrap_repetitions)
for (iteration in seq_len(bootstrap_repetitions)) {
  if (iteration %% 100L == 0L) message("Development bootstrap: ", iteration)
  development_bootstrap[[iteration]] <- tryCatch({
    boot_data <- cluster_bootstrap(development)
    boot_model <- fit_ridge_model(
      boot_data,
      full_features,
      fixed_lambda = full_model$lambda
    )
    boot_data$prediction <- predict_ridge(boot_model, boot_data)
    boot_threshold <- select_threshold(
      boot_data$event_36m, boot_data$prediction, boot_data$WTMEC8YR
    )$threshold
    validation_prediction <- predict_ridge(boot_model, validation)
    validation_metrics <- threshold_metrics(
      validation$event_36m, validation_prediction, validation$WTMEC8YR,
      boot_threshold
    )
    tibble::tibble(
      iteration = iteration,
      success = TRUE,
      threshold = boot_threshold,
      validation_AUC = weighted_auc(
        validation$event_36m, validation_prediction, validation$WTMEC8YR
      ),
      validation_sensitivity = validation_metrics$sensitivity,
      validation_specificity = validation_metrics$specificity
    )
  }, error = function(error) {
    tibble::tibble(
      iteration = iteration,
      success = FALSE,
      threshold = NA_real_,
      validation_AUC = NA_real_,
      validation_sensitivity = NA_real_,
      validation_specificity = NA_real_
    )
  })
}
development_bootstrap <- dplyr::bind_rows(development_bootstrap)

message("Running ", bootstrap_repetitions, " validation bootstrap repetitions...")
validation_bootstrap <- vector("list", bootstrap_repetitions)
for (iteration in seq_len(bootstrap_repetitions)) {
  if (iteration %% 100L == 0L) message("Validation bootstrap: ", iteration)
  validation_bootstrap[[iteration]] <- tryCatch({
    boot_data <- cluster_bootstrap(validation)
    base_auc <- weighted_auc(
      boot_data$event_36m, boot_data$pred_base, boot_data$WTMEC8YR
    )
    full_auc <- weighted_auc(
      boot_data$event_36m, boot_data$pred_full, boot_data$WTMEC8YR
    )
    metrics <- threshold_metrics(
      boot_data$event_36m, boot_data$pred_full, boot_data$WTMEC8YR,
      locked_threshold$threshold
    )
    tibble::tibble(
      iteration = iteration,
      success = TRUE,
      base_AUC = base_auc,
      full_AUC = full_auc,
      delta_AUC = full_auc - base_auc,
      full_Brier = weighted_brier(
        boot_data$event_36m, boot_data$pred_full, boot_data$WTMEC8YR
      ),
      sensitivity = metrics$sensitivity,
      specificity = metrics$specificity,
      ppv = metrics$ppv,
      npv = metrics$npv
    )
  }, error = function(error) {
    tibble::tibble(
      iteration = iteration,
      success = FALSE,
      base_AUC = NA_real_,
      full_AUC = NA_real_,
      delta_AUC = NA_real_,
      full_Brier = NA_real_,
      sensitivity = NA_real_,
      specificity = NA_real_,
      ppv = NA_real_,
      npv = NA_real_
    )
  })
}
validation_bootstrap <- dplyr::bind_rows(validation_bootstrap)

interval_row <- function(data, variable, label) {
  values <- data[[variable]][data$success & is.finite(data[[variable]])]
  intervals <- stats::quantile(values, c(0.025, 0.5, 0.975), na.rm = TRUE)
  tibble::tibble(
    metric = label,
    estimate = intervals[[2L]],
    lower_95 = intervals[[1L]],
    upper_95 = intervals[[3L]],
    successful_repetitions = length(values)
  )
}

validation_intervals <- dplyr::bind_rows(
  interval_row(validation_bootstrap, "base_AUC", "Validation AUC: clinical base"),
  interval_row(validation_bootstrap, "full_AUC", "Validation AUC: full bedside"),
  interval_row(validation_bootstrap, "delta_AUC", "Validation delta AUC"),
  interval_row(validation_bootstrap, "full_Brier", "Validation Brier: full bedside"),
  interval_row(validation_bootstrap, "sensitivity", "Locked threshold sensitivity"),
  interval_row(validation_bootstrap, "specificity", "Locked threshold specificity"),
  interval_row(validation_bootstrap, "ppv", "Locked threshold PPV"),
  interval_row(validation_bootstrap, "npv", "Locked threshold NPV")
)

successful_thresholds <- development_bootstrap$threshold[
  development_bootstrap$success & is.finite(development_bootstrap$threshold)
]
threshold_quantiles <- stats::quantile(
  successful_thresholds, c(0.025, 0.10, 0.25, 0.50, 0.75, 0.90, 0.975),
  na.rm = TRUE
)
threshold_stability <- tibble::tibble(
  locked_development_threshold = locked_threshold$threshold,
  bootstrap_median = threshold_quantiles[[4L]],
  bootstrap_q25 = threshold_quantiles[[3L]],
  bootstrap_q75 = threshold_quantiles[[5L]],
  bootstrap_p10 = threshold_quantiles[[2L]],
  bootstrap_p90 = threshold_quantiles[[6L]],
  bootstrap_lower_95 = threshold_quantiles[[1L]],
  bootstrap_upper_95 = threshold_quantiles[[7L]],
  bootstrap_iqr_width = threshold_quantiles[[5L]] - threshold_quantiles[[3L]],
  bootstrap_p10_p90_width = threshold_quantiles[[6L]] - threshold_quantiles[[2L]],
  successful_repetitions = length(successful_thresholds)
)

calibration <- tibble::tibble(
  cohort = "Temporal validation I",
  calibration_intercept = as.numeric(stats::coef(calibration_intercept_fit)[[1L]]),
  calibration_intercept_se = sqrt(stats::vcov(calibration_intercept_fit)[1L, 1L]),
  calibration_slope = as.numeric(stats::coef(calibration_slope_fit)[[2L]]),
  calibration_slope_se = sqrt(stats::vcov(calibration_slope_fit)[2L, 2L])
) |>
  dplyr::mutate(
    calibration_intercept_lower_95 = calibration_intercept - 1.96 * calibration_intercept_se,
    calibration_intercept_upper_95 = calibration_intercept + 1.96 * calibration_intercept_se,
    calibration_slope_lower_95 = calibration_slope - 1.96 * calibration_slope_se,
    calibration_slope_upper_95 = calibration_slope + 1.96 * calibration_slope_se
  )

validation_full_row <- performance |>
  dplyr::filter(cohort == "Temporal validation I", model == "Full bedside")
validation_base_row <- performance |>
  dplyr::filter(cohort == "Temporal validation I", model == "Clinical base")
locked_dca <- dca |>
  dplyr::filter(model == "Full bedside") |>
  dplyr::slice(which.min(abs(threshold_probability - locked_threshold$threshold)))

gate_components <- tibble::tribble(
  ~criterion, ~passed, ~value,
  "Temporal validation AUC at least 0.70",
    validation_full_row$AUC >= 0.70,
    validation_full_row$AUC,
  "Full model validation AUC exceeds base model",
    validation_full_row$AUC > validation_base_row$AUC,
    validation_full_row$AUC - validation_base_row$AUC,
  "Validation Brier improves on null model",
    validation_full_row$Brier < validation_full_row$null_Brier,
    validation_full_row$Brier - validation_full_row$null_Brier,
  "Calibration slope between 0.70 and 1.30",
    calibration$calibration_slope >= 0.70 & calibration$calibration_slope <= 1.30,
    calibration$calibration_slope,
  "Absolute calibration intercept no greater than 0.50",
    abs(calibration$calibration_intercept) <= 0.50,
    calibration$calibration_intercept,
  "Locked threshold sensitivity at least 0.60",
    validation_threshold_metrics$sensitivity >= 0.60,
    validation_threshold_metrics$sensitivity,
  "Locked threshold specificity at least 0.60",
    validation_threshold_metrics$specificity >= 0.60,
    validation_threshold_metrics$specificity,
  "Bootstrap threshold 10th-90th width no greater than 0.08",
    threshold_stability$bootstrap_p10_p90_width <= 0.08,
    threshold_stability$bootstrap_p10_p90_width,
  "Full model net benefit exceeds treat-all at locked threshold",
    locked_dca$net_benefit > locked_dca$net_benefit_all,
    locked_dca$net_benefit - locked_dca$net_benefit_all,
  "At least 950 successful development bootstrap repetitions",
    threshold_stability$successful_repetitions >= 950L,
    threshold_stability$successful_repetitions
)

candidate_bedside_supported <- all(gate_components$passed)
final_decision <- tibble::tibble(
  prediction_target = "36-month all-cause mortality in community-dwelling older adults",
  development_cycles = "NHANES G/H (2011-2014)",
  temporal_validation_cycle = "NHANES I (2015-2016)",
  excluded_from_36m_validation = "NHANES J because administrative follow-up is incomplete at 36 months",
  comorbidity_definition = paste0(
    "Seven harmonised categories: cancer, arthritis, heart disease, stroke, ",
    "chronic lung disease, hypertension, diabetes"
  ),
  intended_transport_evaluation = "Locked-coefficient transport to historical NHANES 2005-2010",
  locked_risk_threshold = locked_threshold$threshold,
  candidate_bedside_supported = candidate_bedside_supported,
  clinical_deployment_validated = FALSE,
  conclusion = if (candidate_bedside_supported) {
    paste0(
      "The model passes the defined feasibility gate as a temporally validated ",
      "candidate bedside risk tool; independent external validation and impact evaluation remain required."
    )
  } else {
    paste0(
      "The model does not pass all prespecified feasibility criteria and must remain a ",
      "research-only risk model rather than a bedside decision tool."
    )
  }
)

coefficient_matrix <- as.matrix(stats::coef(full_model$fit, s = full_model$lambda))
coefficients <- tibble::tibble(
  term = rownames(coefficient_matrix),
  coefficient = as.numeric(coefficient_matrix[, 1L])
)

manual_prediction <- function(model, df) {
  coefficient_matrix <- as.matrix(stats::coef(model$fit, s = model$lambda))
  intercept <- coefficient_matrix["(Intercept)", 1L]
  beta <- coefficient_matrix[model$features, 1L]
  stats::plogis(intercept + as.numeric(matrix_from(df, model$features) %*% beta))
}
manual_development_prediction <- manual_prediction(full_model, development)
manual_validation_prediction <- manual_prediction(full_model, validation)
manual_prediction_max_difference <- max(abs(c(
  manual_development_prediction - development$pred_full,
  manual_validation_prediction - validation$pred_full
)))

tuning_details <- tibble::tibble(
  model = c("Clinical base", "Full bedside"),
  alpha = 0,
  cross_validation_folds = 5L,
  lambda_rule = "Minimum development cross-validated binomial deviance",
  lambda_selected = c(base_model$lambda, full_model$lambda),
  lambda_min = c(base_model$cv_fit$lambda.min, full_model$cv_fit$lambda.min),
  lambda_1se = c(base_model$cv_fit$lambda.1se, full_model$cv_fit$lambda.1se)
)

fold_details <- development |>
  dplyr::mutate(cv_fold = fold_id) |>
  dplyr::group_by(cv_fold) |>
  dplyr::summarise(
    n = dplyr::n(),
    events = sum(event_36m),
    weighted_events = sum(WTMEC8YR * event_36m) / sum(WTMEC8YR),
    psu_count = dplyr::n_distinct(psu_unique),
    .groups = "drop"
  )

feature_dictionary <- tibble::tribble(
  ~feature, ~bedside_input, ~transformation,
  "age_5y", "Age (years)", "(age - 75) / 5",
  "male", "Sex", "male=1; female=0",
  "comorbidity_count", "Seven-category harmonised comorbidity count",
    "sum of cancer, arthritis, heart disease, stroke, chronic lung disease, hypertension and diabetes",
  "log_nlr", "Neutrophil-to-lymphocyte ratio", "development-winsorised natural log",
  "log_sii", "Systemic immune-inflammation index", "development-winsorised natural log",
  "haemoglobin", "Haemoglobin (g/dL)", "development-winsorised",
  "total_protein", "Serum total protein (g/dL)", "development-winsorised",
  "bmi_5", "Body mass index (kg/m2)", "(development-winsorised BMI - 27) / 5",
  "bmi_5_sq", "Body mass index (kg/m2)", "square of bmi_5",
  "log_creatinine", "Serum creatinine (mg/dL)", "development-winsorised natural log"
)

winsor_table <- dplyr::bind_rows(lapply(names(winsor_bounds), function(variable) {
  tibble::tibble(
    variable = variable,
    development_weighted_p01 = winsor_bounds[[variable]][[1L]],
    development_weighted_p99 = winsor_bounds[[variable]][[2L]]
  )
}))

qa <- tibble::tribble(
  ~check, ~passed, ~detail,
  "Development cohort fixed", nrow(development) == 1954L,
    paste0("n=", nrow(development)),
  "Development events fixed", sum(development$event_36m) == 177L,
    paste0("events=", sum(development$event_36m)),
  "Temporal validation cohort fixed", nrow(validation) == 988L,
    paste0("n=", nrow(validation)),
  "Temporal validation events fixed", sum(validation$event_36m) == 93L,
    paste0("events=", sum(validation$event_36m)),
  "Complete 36-month ascertainment", all(development$horizon_observed) & all(validation$horizon_observed),
    "G/H development and I validation only",
  "Cycle J excluded", !any(as.character(validation$Cycle_ID) == "J"),
    "J has incomplete administrative follow-up at 36 months",
  "Harmonised comorbidity count valid",
    all(c(development$comorbidity_count, validation$comorbidity_count) %in% 0:7),
    "seven disease categories available in the development and historical NHANES periods",
  "Predictions finite", all(is.finite(c(development$pred_full, validation$pred_full))),
    "development and validation",
  "Predictions are probabilities", all(c(development$pred_full, validation$pred_full) > 0 &
    c(development$pred_full, validation$pred_full) < 1),
    "all predictions in (0,1)",
  "Development bootstrap success", sum(development_bootstrap$success) >= 950L,
    paste0(sum(development_bootstrap$success), "/", bootstrap_repetitions),
  "Validation bootstrap success", sum(validation_bootstrap$success) >= 950L,
    paste0(sum(validation_bootstrap$success), "/", bootstrap_repetitions),
  "Published formula reproduces predictions",
    manual_prediction_max_difference < 1e-12,
    paste0("maximum absolute difference=", signif(manual_prediction_max_difference, 4)),
  "No deployment claim", !final_decision$clinical_deployment_validated,
    final_decision$conclusion
)

readr::write_csv(performance, file.path(output_dir, "Table50A_model_performance.csv"))
readr::write_csv(calibration, file.path(output_dir, "Table50B_temporal_calibration.csv"))
readr::write_csv(
  dplyr::bind_rows(
    dplyr::mutate(development_threshold_metrics, cohort = "Development G/H"),
    dplyr::mutate(validation_threshold_metrics, cohort = "Temporal validation I")
  ),
  file.path(output_dir, "Table50C_locked_threshold_performance.csv")
)
readr::write_csv(threshold_stability, file.path(output_dir, "Table50D_threshold_stability.csv"))
readr::write_csv(validation_intervals, file.path(output_dir, "Table50E_validation_bootstrap_intervals.csv"))
readr::write_csv(dca, file.path(output_dir, "Table50F_decision_curve_data.csv"))
readr::write_csv(gate_components, file.path(output_dir, "Table50G_feasibility_gate.csv"))
readr::write_csv(final_decision, file.path(output_dir, "Table50H_final_decision.csv"))
readr::write_csv(coefficients, file.path(output_dir, "Table50I_model_coefficients.csv"))
readr::write_csv(feature_dictionary, file.path(output_dir, "Table50J_feature_dictionary.csv"))
readr::write_csv(winsor_table, file.path(output_dir, "Table50K_winsor_bounds.csv"))
readr::write_csv(qa, file.path(output_dir, "Table50L_QA.csv"))
readr::write_csv(
  development_bootstrap,
  file.path(output_dir, "Table50M_development_bootstrap_repetitions.csv")
)
readr::write_csv(
  validation_bootstrap,
  file.path(output_dir, "Table50N_validation_bootstrap_repetitions.csv")
)
readr::write_csv(tuning_details, file.path(output_dir, "Table50O_tuning_details.csv"))
readr::write_csv(fold_details, file.path(output_dir, "Table50P_cross_validation_folds.csv"))

saveRDS(
  list(
    base_model = base_model,
    full_model = full_model,
    development = development,
    validation = validation,
    performance = performance,
    calibration = calibration,
    locked_threshold = locked_threshold,
    threshold_stability = threshold_stability,
    validation_intervals = validation_intervals,
    dca = dca,
    gate_components = gate_components,
    final_decision = final_decision,
    qa = qa
  ),
  file.path(output_dir, "NHANES_locked_36m_risk_model_results.rds")
)

summary_lines <- c(
  "Locked NHANES 36-month risk model",
  "",
  paste0("Development: n=", nrow(development), ", events=", sum(development$event_36m)),
  paste0("Temporal validation: n=", nrow(validation), ", events=", sum(validation$event_36m)),
  paste0("Locked 36-month risk threshold: ", sprintf("%.4f", locked_threshold$threshold)),
  paste0("Validation full-model AUC: ", sprintf("%.3f", validation_full_row$AUC)),
  paste0("Validation calibration slope: ", sprintf("%.3f", calibration$calibration_slope)),
  paste0("Validation sensitivity: ", sprintf("%.3f", validation_threshold_metrics$sensitivity)),
  paste0("Validation specificity: ", sprintf("%.3f", validation_threshold_metrics$specificity)),
  paste0("Threshold bootstrap P10-P90 width: ", sprintf("%.4f", threshold_stability$bootstrap_p10_p90_width)),
  paste0("Feasibility criteria passed: ", sum(gate_components$passed), "/", nrow(gate_components)),
  paste0("Candidate bedside supported: ", candidate_bedside_supported),
  final_decision$conclusion,
  "",
  paste0("QA checks passed: ", sum(qa$passed), "/", nrow(qa))
)
writeLines(summary_lines, file.path(output_dir, "NHANES_locked_36m_risk_model_summary.txt"))

if (!all(qa$passed)) {
  stop("Bedside prediction QA failed. Review Table50L_QA.csv.", call. = FALSE)
}

cat(paste(summary_lines, collapse = "\n"), "\n")

