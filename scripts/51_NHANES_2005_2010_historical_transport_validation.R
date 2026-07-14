# ==============================================================================
# Locked-coefficient historical transport of the exploratory NHANES 36-month
# risk model to NHANES 2005-2010. No model refitting or threshold reselection.
# ==============================================================================

required_packages <- c("haven", "dplyr", "readr", "tibble", "survey", "glmnet")
missing_packages <- required_packages[!vapply(
  required_packages, requireNamespace, logical(1), quietly = TRUE
)]
if (length(missing_packages) > 0L) {
  stop("Missing required package(s): ", paste(missing_packages, collapse = ", "), call. = FALSE)
}

set.seed(20260710)
options(survey.lonely.psu = "adjust")

stable_ntile <- function(value, groups, tie_breaker) {
  keep <- is.finite(value) & !is.na(tie_breaker)
  output <- rep(NA_integer_, length(value))
  ordered_index <- which(keep)[order(value[keep], tie_breaker[keep])]
  output[ordered_index] <- dplyr::ntile(seq_along(ordered_index), groups)
  output
}

root_dir <- normalizePath(
  Sys.getenv("PROJECT_ROOT", unset = getwd()), winslash = "/", mustWork = TRUE
)
model_dir <- file.path(root_dir, "output", "nhanes_locked_36m_risk_model")
output_dir <- Sys.getenv(
  "NHANES_HISTORICAL_OUTPUT",
  unset = file.path(root_dir, "output", "nhanes_2005_2010_historical_transport")
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
cache_dir <- Sys.getenv(
  "NHANES_HISTORICAL_CACHE",
  unset = file.path(output_dir, "official_source_cache")
)
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

model_path <- file.path(model_dir, "NHANES_locked_36m_risk_model_results.rds")
coefficient_path <- file.path(model_dir, "Table50I_model_coefficients.csv")
bounds_path <- file.path(model_dir, "Table50K_winsor_bounds.csv")
if (!all(file.exists(c(model_path, coefficient_path, bounds_path)))) {
  stop("Frozen model assets are missing. Run script 50 first.", call. = FALSE)
}

cycles <- tibble::tribble(
  ~Cycle_ID, ~survey_start_year, ~survey_years,
  "D", "2005", "2005_2006",
  "E", "2007", "2007_2008",
  "F", "2009", "2009_2010"
)
components <- c("DEMO", "CBC", "BIOPRO", "BMX", "MCQ", "BPQ", "DIQ")
horizon_months <- 36L
locked_threshold <- 0.0666572069309404
bootstrap_repetitions <- 1000L

is_valid_xpt <- function(path) {
  if (!file.exists(path) || file.info(path)$size < 1000) return(FALSE)
  connection <- file(path, open = "rb")
  on.exit(close(connection), add = TRUE)
  header <- rawToChar(readBin(connection, what = "raw", n = 80L))
  grepl("HEADER RECORD", header, fixed = TRUE)
}

download_with_retries <- function(url, destination, validator, attempts = 3L) {
  if (validator(destination)) return(destination)
  for (attempt in seq_len(attempts)) {
    message("Downloading ", basename(destination), " (attempt ", attempt, ")...")
    result <- tryCatch(
      utils::download.file(url, destination, mode = "wb", quiet = TRUE),
      error = function(error) error
    )
    if (!inherits(result, "error") && validator(destination)) return(destination)
    if (file.exists(destination)) unlink(destination)
  }
  stop("Failed to download official source: ", url, call. = FALSE)
}

download_component <- function(component, cycle, survey_start_year) {
  filename <- paste0(component, "_", cycle, ".xpt")
  destination <- file.path(cache_dir, filename)
  url <- paste0(
    "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/", survey_start_year,
    "/DataFiles/", filename
  )
  download_with_retries(url, destination, is_valid_xpt)
}

component_variables <- list(
  DEMO = c("SEQN", "RIDAGEYR", "RIAGENDR", "WTMEC2YR", "SDMVSTRA", "SDMVPSU"),
  CBC = c("SEQN", "LBDNENO", "LBDLYMNO", "LBXPLTSI", "LBXHGB"),
  BIOPRO = c("SEQN", "LBXSTP", "LBXSCR"),
  BMX = c("SEQN", "BMXBMI"),
  MCQ = c(
    "SEQN", "MCQ220", "MCQ160A", "MCQ160B", "MCQ160C", "MCQ160D",
    "MCQ160E", "MCQ160F", "MCQ160G", "MCQ160K"
  ),
  BPQ = c("SEQN", "BPQ020"),
  DIQ = c("SEQN", "DIQ010")
)

read_component <- function(component, cycle, survey_start_year) {
  path <- download_component(component, cycle, survey_start_year)
  data <- haven::read_xpt(path)
  required <- component_variables[[component]]
  missing <- setdiff(required, names(data))
  if (length(missing) > 0L) {
    stop(
      "Missing variable(s) in ", component, "_", cycle, ": ",
      paste(missing, collapse = ", "), call. = FALSE
    )
  }
  data |>
    dplyr::select(dplyr::all_of(required))
}

read_cycle <- function(cycle, survey_start_year) {
  component_data <- lapply(
    components,
    read_component,
    cycle = cycle,
    survey_start_year = survey_start_year
  )
  names(component_data) <- components
  if (any(vapply(component_data, function(data) anyDuplicated(data$SEQN) > 0L, logical(1)))) {
    stop("Duplicate SEQN detected in cycle ", cycle, call. = FALSE)
  }
  Reduce(function(x, y) dplyr::left_join(x, y, by = "SEQN"), component_data) |>
    dplyr::mutate(Cycle_ID = cycle)
}

is_valid_mortality <- function(path) {
  file.exists(path) && file.info(path)$size > 1000
}

download_mortality <- function(survey_years) {
  filename <- paste0("NHANES_", survey_years, "_MORT_2019_PUBLIC.dat")
  destination <- file.path(cache_dir, filename)
  url <- paste0(
    "https://ftp.cdc.gov/pub/Health_Statistics/NCHS/datalinkage/linked_mortality/",
    filename
  )
  download_with_retries(url, destination, is_valid_mortality)
}

read_mortality <- function(path, cycle) {
  positions <- readr::fwf_positions(
    start = c(1, 15, 16, 43),
    end = c(6, 15, 16, 45),
    col_names = c("SEQN", "ELIGSTAT", "MORTSTAT", "PERMTH_INT")
  )
  data <- readr::read_fwf(
    path,
    col_positions = positions,
    na = ".",
    col_types = readr::cols(
      SEQN = readr::col_integer(),
      ELIGSTAT = readr::col_integer(),
      MORTSTAT = readr::col_integer(),
      PERMTH_INT = readr::col_integer()
    ),
    progress = FALSE
  )
  if (nrow(data) == 0L || anyDuplicated(data$SEQN) > 0L) {
    stop("Invalid mortality file for cycle ", cycle, call. = FALSE)
  }
  dplyr::mutate(data, Mortality_Cycle_ID = cycle)
}

cycle_data <- lapply(seq_len(nrow(cycles)), function(index) {
  read_cycle(cycles$Cycle_ID[index], cycles$survey_start_year[index])
})
raw_data <- dplyr::bind_rows(cycle_data)

mortality_data <- lapply(seq_len(nrow(cycles)), function(index) {
  path <- download_mortality(cycles$survey_years[index])
  read_mortality(path, cycles$Cycle_ID[index])
})
mortality <- dplyr::bind_rows(mortality_data)

yes_no_flag <- function(x) {
  dplyr::case_when(x == 1 ~ 1, x == 2 ~ 0, TRUE ~ NA_real_)
}

collapse_positive <- function(...) {
  values <- cbind(...)
  apply(values, 1L, function(row) {
    if (any(row == 1, na.rm = TRUE)) return(1)
    if (all(!is.na(row)) && all(row == 0)) return(0)
    NA_real_
  })
}

linked <- raw_data |>
  dplyr::left_join(mortality, by = "SEQN") |>
  dplyr::filter(Cycle_ID == Mortality_Cycle_ID, ELIGSTAT == 1, RIDAGEYR >= 65) |>
  dplyr::mutate(
    male = as.numeric(RIAGENDR == 1),
    NLR = dplyr::if_else(LBDLYMNO > 0, LBDNENO / LBDLYMNO, NA_real_),
    SII = dplyr::if_else(LBDLYMNO > 0, LBXPLTSI * LBDNENO / LBDLYMNO, NA_real_),
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
    comorbidity_count = cancer_flag + arthritis_flag + heart_disease_flag +
      stroke_flag + chronic_lung_flag + hypertension_flag + diabetes_flag,
    event_36m = as.integer(MORTSTAT == 1 & PERMTH_INT <= horizon_months),
    horizon_observed = event_36m == 1 | PERMTH_INT >= horizon_months,
    WTMEC6YR = WTMEC2YR / nrow(cycles),
    strata_unique = interaction(Cycle_ID, SDMVSTRA, drop = TRUE),
    psu_unique = interaction(Cycle_ID, SDMVSTRA, SDMVPSU, drop = TRUE)
  )

analysis_variables <- c(
  "RIDAGEYR", "male", "comorbidity_count", "NLR", "SII", "LBXHGB",
  "LBXSTP", "BMXBMI", "LBXSCR", "event_36m", "PERMTH_INT", "WTMEC6YR",
  "SDMVSTRA", "SDMVPSU"
)

analysis_data <- linked |>
  dplyr::filter(horizon_observed, WTMEC6YR > 0) |>
  dplyr::filter(dplyr::if_all(
    dplyr::all_of(analysis_variables),
    ~ !is.na(.x) & is.finite(.x)
  )) |>
  droplevels()

model_results <- readRDS(model_path)
bounds_table <- readr::read_csv(bounds_path, show_col_types = FALSE)
exported_coefficients <- readr::read_csv(coefficient_path, show_col_types = FALSE)

winsor_bounds <- stats::setNames(
  lapply(seq_len(nrow(bounds_table)), function(index) {
    c(
      bounds_table$development_weighted_p01[index],
      bounds_table$development_weighted_p99[index]
    )
  }),
  bounds_table$variable
)

winsorise <- function(x, bounds) {
  pmin(pmax(x, bounds[[1L]]), bounds[[2L]])
}

prepare_features <- function(data) {
  nlr <- winsorise(data$NLR, winsor_bounds$NLR)
  sii <- winsorise(data$SII, winsor_bounds$SII)
  haemoglobin <- winsorise(data$LBXHGB, winsor_bounds$LBXHGB)
  total_protein <- winsorise(data$LBXSTP, winsor_bounds$LBXSTP)
  bmi <- winsorise(data$BMXBMI, winsor_bounds$BMXBMI)
  creatinine <- winsorise(data$LBXSCR, winsor_bounds$LBXSCR)
  dplyr::mutate(
    data,
    age_5y = (RIDAGEYR - 75) / 5,
    log_nlr = log(nlr),
    log_sii = log(sii),
    haemoglobin = haemoglobin,
    total_protein = total_protein,
    bmi_5 = (bmi - 27) / 5,
    bmi_5_sq = bmi_5^2,
    log_creatinine = log(creatinine)
  )
}

analysis_data <- prepare_features(analysis_data)
base_features <- c("age_5y", "male", "comorbidity_count")
full_features <- c(
  base_features, "log_nlr", "log_sii", "haemoglobin", "total_protein",
  "bmi_5", "bmi_5_sq", "log_creatinine"
)

extract_coefficients <- function(model) {
  matrix <- as.matrix(stats::coef(model$fit, s = model$lambda))
  stats::setNames(as.numeric(matrix[, 1L]), rownames(matrix))
}

base_coefficients <- extract_coefficients(model_results$base_model)
full_coefficients <- extract_coefficients(model_results$full_model)

predict_fixed <- function(data, coefficients, features) {
  linear_predictor <- coefficients[["(Intercept)"]] +
    as.numeric(as.matrix(data[, features, drop = FALSE]) %*% coefficients[features])
  stats::plogis(linear_predictor)
}

analysis_data$pred_base <- predict_fixed(analysis_data, base_coefficients, base_features)
analysis_data$pred_full <- predict_fixed(analysis_data, full_coefficients, full_features)
analysis_data$lp_full <- stats::qlogis(pmin(pmax(analysis_data$pred_full, 1e-8), 1 - 1e-8))

weighted_auc <- function(outcome, prediction, weight) {
  aggregated <- tibble::tibble(outcome, prediction, weight) |>
    dplyr::group_by(prediction) |>
    dplyr::summarise(
      case_weight = sum(weight[outcome == 1L]),
      control_weight = sum(weight[outcome == 0L]),
      .groups = "drop"
    ) |>
    dplyr::arrange(prediction) |>
    dplyr::mutate(control_below = dplyr::lag(cumsum(control_weight), default = 0))
  numerator <- sum(
    aggregated$case_weight *
      (aggregated$control_below + 0.5 * aggregated$control_weight)
  )
  denominator <- sum(aggregated$case_weight) * sum(aggregated$control_weight)
  numerator / denominator
}

weighted_brier <- function(outcome, prediction, weight) {
  stats::weighted.mean((outcome - prediction)^2, weight)
}

threshold_metrics <- function(outcome, prediction, weight, threshold) {
  high <- prediction >= threshold
  true_positive <- sum(weight[outcome == 1L & high])
  false_negative <- sum(weight[outcome == 1L & !high])
  true_negative <- sum(weight[outcome == 0L & !high])
  false_positive <- sum(weight[outcome == 0L & high])
  tibble::tibble(
    threshold = threshold,
    sensitivity = true_positive / (true_positive + false_negative),
    specificity = true_negative / (true_negative + false_positive),
    ppv = true_positive / (true_positive + false_positive),
    npv = true_negative / (true_negative + false_negative),
    classified_high_percent = 100 * sum(weight[high]) / sum(weight)
  )
}

performance_row <- function(data, prediction_name, cohort, model_name) {
  prediction <- data[[prediction_name]]
  prevalence <- stats::weighted.mean(data$event_36m, data$WTMEC6YR)
  tibble::tibble(
    cohort = cohort,
    model = model_name,
    n = nrow(data),
    events = sum(data$event_36m),
    weighted_event_rate = prevalence,
    AUC = weighted_auc(data$event_36m, prediction, data$WTMEC6YR),
    Brier = weighted_brier(data$event_36m, prediction, data$WTMEC6YR),
    null_Brier = prevalence * (1 - prevalence)
  )
}

overall_performance <- dplyr::bind_rows(
  performance_row(analysis_data, "pred_base", "NHANES 2005-2010", "Clinical base"),
  performance_row(analysis_data, "pred_full", "NHANES 2005-2010", "Full bedside")
)

cycle_performance <- dplyr::bind_rows(lapply(split(analysis_data, analysis_data$Cycle_ID), function(data) {
  dplyr::bind_rows(
    performance_row(data, "pred_base", unique(data$Cycle_ID), "Clinical base"),
    performance_row(data, "pred_full", unique(data$Cycle_ID), "Full bedside")
  )
}))

design <- survey::svydesign(
  ids = ~psu_unique,
  strata = ~strata_unique,
  weights = ~WTMEC6YR,
  nest = TRUE,
  data = analysis_data
)
calibration_intercept_fit <- survey::svyglm(
  event_36m ~ 1,
  offset = lp_full,
  design = design,
  family = stats::quasibinomial()
)
calibration_slope_fit <- survey::svyglm(
  event_36m ~ lp_full,
  design = design,
  family = stats::quasibinomial()
)
calibration <- tibble::tibble(
  cohort = "NHANES 2005-2010",
  calibration_intercept = as.numeric(stats::coef(calibration_intercept_fit)[1L]),
  calibration_intercept_se = sqrt(stats::vcov(calibration_intercept_fit)[1L, 1L]),
  calibration_slope = as.numeric(stats::coef(calibration_slope_fit)[2L]),
  calibration_slope_se = sqrt(stats::vcov(calibration_slope_fit)[2L, 2L])
) |>
  dplyr::mutate(
    calibration_intercept_lower_95 = calibration_intercept - 1.96 * calibration_intercept_se,
    calibration_intercept_upper_95 = calibration_intercept + 1.96 * calibration_intercept_se,
    calibration_slope_lower_95 = calibration_slope - 1.96 * calibration_slope_se,
    calibration_slope_upper_95 = calibration_slope + 1.96 * calibration_slope_se
  )

locked_threshold_performance <- threshold_metrics(
  analysis_data$event_36m,
  analysis_data$pred_full,
  analysis_data$WTMEC6YR,
  locked_threshold
)

calibration_groups <- analysis_data |>
  dplyr::mutate(risk_group = stable_ntile(pred_full, 10L, SEQN)) |>
  dplyr::group_by(risk_group) |>
  dplyr::summarise(
    n = dplyr::n(),
    events = sum(event_36m),
    weighted_mean_prediction = stats::weighted.mean(pred_full, WTMEC6YR),
    weighted_observed_rate = stats::weighted.mean(event_36m, WTMEC6YR),
    .groups = "drop"
  )

decision_curve <- function(data, prediction_name, thresholds) {
  prediction <- data[[prediction_name]]
  outcome <- data$event_36m
  weight <- data$WTMEC6YR
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

dca_thresholds <- sort(unique(c(seq(0.02, 0.25, by = 0.01), locked_threshold)))
decision_curve_data <- dplyr::bind_rows(
  dplyr::mutate(decision_curve(analysis_data, "pred_base", dca_thresholds), model = "Clinical base"),
  dplyr::mutate(decision_curve(analysis_data, "pred_full", dca_thresholds), model = "Full bedside")
)

cluster_bootstrap <- function(data) {
  pieces <- list()
  piece_index <- 1L
  for (stratum in unique(as.character(data$strata_unique))) {
    stratum_data <- data[as.character(data$strata_unique) == stratum, , drop = FALSE]
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

message("Running ", bootstrap_repetitions, " historical-validation bootstrap repetitions...")
bootstrap_results <- vector("list", bootstrap_repetitions)
for (iteration in seq_len(bootstrap_repetitions)) {
  if (iteration %% 100L == 0L) message("Bootstrap: ", iteration)
  bootstrap_results[[iteration]] <- tryCatch({
    boot <- cluster_bootstrap(analysis_data)
    boot$regression_weight <- boot$WTMEC6YR / mean(boot$WTMEC6YR)
    metrics <- threshold_metrics(
      boot$event_36m, boot$pred_full, boot$WTMEC6YR, locked_threshold
    )
    intercept_fit <- stats::glm(
      event_36m ~ 1,
      offset = lp_full,
      data = boot,
      weights = regression_weight,
      family = stats::quasibinomial()
    )
    slope_fit <- stats::glm(
      event_36m ~ lp_full,
      data = boot,
      weights = regression_weight,
      family = stats::quasibinomial()
    )
    tibble::tibble(
      iteration = iteration,
      success = TRUE,
      AUC_full = weighted_auc(boot$event_36m, boot$pred_full, boot$WTMEC6YR),
      AUC_base = weighted_auc(boot$event_36m, boot$pred_base, boot$WTMEC6YR),
      delta_AUC = AUC_full - AUC_base,
      Brier_full = weighted_brier(boot$event_36m, boot$pred_full, boot$WTMEC6YR),
      calibration_intercept = as.numeric(stats::coef(intercept_fit)[1L]),
      calibration_slope = as.numeric(stats::coef(slope_fit)[2L]),
      sensitivity = metrics$sensitivity,
      specificity = metrics$specificity,
      ppv = metrics$ppv,
      npv = metrics$npv,
      classified_high_percent = metrics$classified_high_percent
    )
  }, error = function(error) {
    tibble::tibble(iteration = iteration, success = FALSE, error = conditionMessage(error))
  })
}
bootstrap_results <- dplyr::bind_rows(bootstrap_results)

interval_row <- function(data, variable, label) {
  values <- data[[variable]][data$success & is.finite(data[[variable]])]
  quantiles <- stats::quantile(values, c(0.025, 0.50, 0.975), na.rm = TRUE)
  tibble::tibble(
    metric = label,
    bootstrap_median = quantiles[[2L]],
    lower_95 = quantiles[[1L]],
    upper_95 = quantiles[[3L]],
    successful_repetitions = length(values)
  )
}

bootstrap_intervals <- dplyr::bind_rows(
  interval_row(bootstrap_results, "AUC_full", "Full-model AUC"),
  interval_row(bootstrap_results, "AUC_base", "Clinical-base AUC"),
  interval_row(bootstrap_results, "delta_AUC", "AUC increment"),
  interval_row(bootstrap_results, "Brier_full", "Full-model Brier score"),
  interval_row(bootstrap_results, "calibration_intercept", "Calibration intercept"),
  interval_row(bootstrap_results, "calibration_slope", "Calibration slope"),
  interval_row(bootstrap_results, "sensitivity", "Locked-threshold sensitivity"),
  interval_row(bootstrap_results, "specificity", "Locked-threshold specificity"),
  interval_row(bootstrap_results, "ppv", "Locked-threshold PPV"),
  interval_row(bootstrap_results, "npv", "Locked-threshold NPV"),
  interval_row(bootstrap_results, "classified_high_percent", "High-risk percentage")
)

cohort_flow <- linked |>
  dplyr::group_by(Cycle_ID) |>
  dplyr::summarise(
    mortality_eligible_age65 = dplyr::n(),
    horizon_observed = sum(horizon_observed, na.rm = TRUE),
    complete_analysis = sum(
      horizon_observed & WTMEC6YR > 0 &
        stats::complete.cases(dplyr::pick(dplyr::all_of(analysis_variables)))
    ),
    .groups = "drop"
  )

source_files <- c(
  unlist(lapply(seq_len(nrow(cycles)), function(index) {
    file.path(cache_dir, paste0(components, "_", cycles$Cycle_ID[index], ".xpt"))
  })),
  file.path(
    cache_dir,
    paste0("NHANES_", cycles$survey_years, "_MORT_2019_PUBLIC.dat")
  )
)
source_audit <- tibble::tibble(
  source_file = basename(source_files),
  exists = file.exists(source_files),
  bytes = file.info(source_files)$size,
  md5 = unname(tools::md5sum(source_files))
)

coefficient_difference <- max(abs(
  full_coefficients[exported_coefficients$term] - exported_coefficients$coefficient
))
locked_dca <- decision_curve_data |>
  dplyr::filter(model == "Full bedside") |>
  dplyr::slice(which.min(abs(threshold_probability - locked_threshold)))

qa <- tibble::tribble(
  ~check, ~passed, ~detail,
  "All official source files available", all(source_audit$exists), paste0(sum(source_audit$exists), "/", nrow(source_audit)),
  "All D/E/F cycles represented", identical(sort(unique(analysis_data$Cycle_ID)), c("D", "E", "F")), paste(sort(unique(analysis_data$Cycle_ID)), collapse = ","),
  "No duplicate participant identifiers", anyDuplicated(analysis_data$SEQN) == 0L, as.character(anyDuplicated(analysis_data$SEQN)),
  "All participants aged at least 65", all(analysis_data$RIDAGEYR >= 65), as.character(min(analysis_data$RIDAGEYR)),
  "All 36-month outcomes observed", all(analysis_data$horizon_observed), as.character(sum(analysis_data$horizon_observed)),
  "Seven-category comorbidity score valid", all(analysis_data$comorbidity_count %in% 0:7), paste(range(analysis_data$comorbidity_count), collapse = "-"),
  "Frozen coefficient export matches model object", coefficient_difference < 1e-12, format(coefficient_difference, scientific = TRUE),
  "Locked threshold unchanged", identical(locked_threshold, 0.0666572069309404), format(locked_threshold, digits = 16),
  "Predictions are finite probabilities", all(is.finite(analysis_data$pred_full) & analysis_data$pred_full > 0 & analysis_data$pred_full < 1), paste(range(analysis_data$pred_full), collapse = "-"),
  "At least 950 bootstrap repetitions succeeded", sum(bootstrap_results$success) >= 950L, paste0(sum(bootstrap_results$success), "/", bootstrap_repetitions),
  "Bootstrap calibration estimates numerically stable", all(abs(bootstrap_results$calibration_intercept[bootstrap_results$success]) < 10 & abs(bootstrap_results$calibration_slope[bootstrap_results$success]) < 10), "absolute intercept and slope below 10",
  "NHANES-only historical transport", TRUE, "Only CDC NHANES 2005-2010 public files and locked NHANES model assets",
  "No model refitting or threshold reselection", TRUE, "Fixed coefficients, bounds and 6.67% threshold",
  "Locked-threshold net benefit calculated", is.finite(locked_dca$net_benefit), as.character(locked_dca$net_benefit)
)

readr::write_csv(cohort_flow, file.path(output_dir, "Table51A_cohort_flow.csv"))
readr::write_csv(overall_performance, file.path(output_dir, "Table51B_overall_performance.csv"))
readr::write_csv(cycle_performance, file.path(output_dir, "Table51C_cycle_performance.csv"))
readr::write_csv(calibration, file.path(output_dir, "Table51D_calibration.csv"))
readr::write_csv(locked_threshold_performance, file.path(output_dir, "Table51E_locked_threshold_performance.csv"))
readr::write_csv(bootstrap_intervals, file.path(output_dir, "Table51F_bootstrap_intervals.csv"))
readr::write_csv(calibration_groups, file.path(output_dir, "Table51G_calibration_groups.csv"))
readr::write_csv(decision_curve_data, file.path(output_dir, "Table51H_decision_curve.csv"))
readr::write_csv(bootstrap_results, file.path(output_dir, "Table51I_bootstrap_repetitions.csv"))
readr::write_csv(qa, file.path(output_dir, "Table51J_QA.csv"))
readr::write_csv(source_audit, file.path(output_dir, "Table51K_official_source_audit.csv"))

saveRDS(
  list(
    analysis_data = analysis_data,
    cohort_flow = cohort_flow,
    overall_performance = overall_performance,
    cycle_performance = cycle_performance,
    calibration = calibration,
    locked_threshold_performance = locked_threshold_performance,
    bootstrap_intervals = bootstrap_intervals,
    calibration_groups = calibration_groups,
    decision_curve_data = decision_curve_data,
    qa = qa
  ),
  file.path(output_dir, "NHANES_2005_2010_historical_transport_results.rds")
)

full_row <- dplyr::filter(overall_performance, model == "Full bedside")
base_row <- dplyr::filter(overall_performance, model == "Clinical base")
summary_lines <- c(
  "NHANES 2005-2010 historical transport validation",
  "",
  paste0("Analysis cohort: n=", nrow(analysis_data), ", events=", sum(analysis_data$event_36m)),
  paste0("Weighted 36-month mortality: ", sprintf("%.2f%%", 100 * full_row$weighted_event_rate)),
  paste0("Fixed full-model AUC: ", sprintf("%.3f", full_row$AUC)),
  paste0("Fixed clinical-base AUC: ", sprintf("%.3f", base_row$AUC)),
  paste0("AUC increment: ", sprintf("%.3f", full_row$AUC - base_row$AUC)),
  paste0("Calibration intercept: ", sprintf("%.3f", calibration$calibration_intercept)),
  paste0("Calibration slope: ", sprintf("%.3f", calibration$calibration_slope)),
  paste0("Locked threshold: ", sprintf("%.4f", locked_threshold)),
  paste0("Sensitivity: ", sprintf("%.3f", locked_threshold_performance$sensitivity)),
  paste0("Specificity: ", sprintf("%.3f", locked_threshold_performance$specificity)),
  paste0("High-risk percentage: ", sprintf("%.1f%%", locked_threshold_performance$classified_high_percent)),
  "",
  "Interpretation: historical temporal transportability only; this is not independent-database external validation.",
  "Only NHANES data were used; the fixed model was not refitted.",
  paste0("QA checks passed: ", sum(qa$passed), "/", nrow(qa))
)
writeLines(summary_lines, file.path(output_dir, "NHANES_2005_2010_historical_transport_summary.txt"))

if (!all(qa$passed)) {
  stop("Historical transport QA failed. Review Table51J_QA.csv.", call. = FALSE)
}

cat(paste(summary_lines, collapse = "\n"), "\n")
