# ==============================================================================
# MIT-LCP MIMIC Code v3.0.1 official OASIS reconstruction and outcome analysis
# ==============================================================================

required_packages <- c("DBI", "duckdb", "dplyr", "readr", "tibble", "survival")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop("Missing required package(s): ", paste(missing_packages, collapse = ", "), call. = FALSE)
}

root <- normalizePath(
  Sys.getenv("PROJECT_ROOT", unset = getwd()), winslash = "/", mustWork = TRUE
)
mimic_root <- Sys.getenv(
  "MIMIC_ROOT", unset = file.path(root, "data", "mimic-iv-3.1")
)
official_root <- Sys.getenv(
  "MIMIC_CODE_ROOT",
  unset = file.path(root, "vendor", "mimic-code-3.0.1", "mimic-iv", "concepts_duckdb")
)
output_dir <- file.path(root, "output", "mimic_official_oasis_v301")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

source_rds <- file.path(
  root, "output", "mimic_24h_robust_corrected_severity", "MIMIC_24h_robust_results.rds"
)
old_oasis_path <- file.path(
  root, "output", "mimic_oasislike_severity_corrected", "MIMIC_denovo_with_OASISlike_score.csv"
)
denominator_path <- file.path(
  root, "output", "mimic_24h_validation",
  "MIMIC_first_ICU_feature_availability_dataset.csv"
)

raw_paths <- list(
  admissions = file.path(mimic_root, "hosp", "admissions.csv.gz"),
  patients = file.path(mimic_root, "hosp", "patients.csv.gz"),
  services = file.path(mimic_root, "hosp", "services.csv.gz"),
  icustays = file.path(mimic_root, "icu", "icustays.csv.gz"),
  chartevents = file.path(mimic_root, "icu", "chartevents.csv.gz"),
  outputevents = file.path(mimic_root, "icu", "outputevents.csv.gz")
)

official_scripts <- c(
  "demographics/age.sql",
  "measurement/gcs.sql",
  "measurement/vitalsign.sql",
  "measurement/urine_output.sql",
  "measurement/ventilator_setting.sql",
  "measurement/oxygen_delivery.sql",
  "treatment/ventilation.sql",
  "firstday/first_day_gcs.sql",
  "firstday/first_day_vitalsign.sql",
  "firstday/first_day_urine_output.sql",
  "score/oasis.sql"
)
official_paths <- file.path(official_root, official_scripts)

missing_inputs <- c(
  unlist(raw_paths)[!file.exists(unlist(raw_paths))],
  official_paths[!file.exists(official_paths)],
  c(source_rds, old_oasis_path, denominator_path)[
    !file.exists(c(source_rds, old_oasis_path, denominator_path))
  ]
)

stable_ntile <- function(value, groups, tie_breaker) {
  keep <- is.finite(value) & !is.na(tie_breaker)
  output <- rep(NA_integer_, length(value))
  ordered_index <- which(keep)[order(value[keep], tie_breaker[keep])]
  output[ordered_index] <- dplyr::ntile(seq_along(ordered_index), groups)
  output
}
if (length(missing_inputs) > 0L) {
  stop("Missing required input(s): ", paste(missing_inputs, collapse = "; "), call. = FALSE)
}

sql_path <- function(path) {
  value <- normalizePath(path, winslash = "/", mustWork = TRUE)
  gsub("'", "''", value, fixed = TRUE)
}

execute_official_sql <- function(connection, path) {
  cat("Executing official concept: ", basename(path), "\n", sep = "")
  sql <- paste(readLines(path, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  DBI::dbExecute(connection, sql)
}

collapse_for_selection <- function(x, selected, min_total = 100L, min_selected = 10L) {
  value <- as.character(x)
  value[is.na(value) | trimws(value) == ""] <- "Missing"
  tab_total <- table(value)
  tab_selected <- table(value[selected == 1])
  selected_count <- as.numeric(tab_selected[names(tab_total)])
  selected_count[is.na(selected_count)] <- 0
  excluded_count <- as.numeric(tab_total) - selected_count
  keep <- names(tab_total)[
    as.numeric(tab_total) >= min_total & selected_count >= min_selected & excluded_count >= min_selected
  ]
  factor(ifelse(value %in% keep, value, "Other"))
}

extract_cox <- function(fit, model, n, events) {
  beta <- stats::coef(fit)
  se <- sqrt(diag(stats::vcov(fit)))
  terms <- names(beta)
  keep <- grepl("^phenotype", terms)
  tibble::tibble(
    model = model,
    effect_measure = "HR",
    comparison = dplyr::recode(
      terms[keep], phenotypeP2 = "P2 vs P3", phenotypeP1 = "P1 vs P3"
    ),
    estimate = exp(beta[keep]),
    lower_95 = exp(beta[keep] - 1.96 * se[keep]),
    upper_95 = exp(beta[keep] + 1.96 * se[keep]),
    p_value = 2 * stats::pnorm(abs(beta[keep] / se[keep]), lower.tail = FALSE),
    n = n,
    events = events
  ) |>
    dplyr::mutate(effect_95ci = sprintf("%.2f (%.2f-%.2f)", estimate, lower_95, upper_95))
}

extract_logistic <- function(fit, model, n, events) {
  coefficients <- summary(fit)$coefficients
  beta <- coefficients[, "Estimate"]
  se <- coefficients[, "Std. Error"]
  terms <- rownames(coefficients)
  keep <- grepl("^phenotype", terms)
  tibble::tibble(
    model = model,
    effect_measure = "OR",
    comparison = dplyr::recode(
      terms[keep], phenotypeP2 = "P2 vs P3", phenotypeP1 = "P1 vs P3"
    ),
    estimate = exp(beta[keep]),
    lower_95 = exp(beta[keep] - 1.96 * se[keep]),
    upper_95 = exp(beta[keep] + 1.96 * se[keep]),
    p_value = 2 * stats::pnorm(abs(beta[keep] / se[keep]), lower.tail = FALSE),
    n = n,
    events = events
  ) |>
    dplyr::mutate(effect_95ci = sprintf("%.2f (%.2f-%.2f)", estimate, lower_95, upper_95))
}

cat("Loading strict first-24-hour MIMIC-IV phenotype cohort...\n")
mimic_results <- readRDS(source_rds)
analysis <- mimic_results$analysis |>
  dplyr::mutate(
    stay_id = as.integer(stay_id),
    hadm_id = as.integer(hadm_id),
    subject_id = as.integer(subject_id),
    phenotype = factor(as.character(phenotype), levels = c("P3", "P2", "P1")),
    male = as.integer(gender == "M")
  )
if (nrow(analysis) != 1145L || dplyr::n_distinct(analysis$stay_id) != nrow(analysis)) {
  stop("Unexpected strict first-24-hour cohort size or duplicate stay_id.", call. = FALSE)
}

cohort_ids <- analysis |>
  dplyr::distinct(stay_id, hadm_id, subject_id)

cat("Creating cohort-restricted read-only MIMIC-IV views...\n")
con <- DBI::dbConnect(duckdb::duckdb(), dbdir = ":memory:")
on.exit(try(DBI::dbDisconnect(con, shutdown = TRUE), silent = TRUE), add = TRUE)
DBI::dbExecute(con, "SELECT 1")
DBI::dbExecute(con, "CREATE SCHEMA mimiciv_hosp")
DBI::dbExecute(con, "CREATE SCHEMA mimiciv_icu")
DBI::dbExecute(con, "CREATE SCHEMA mimiciv_derived")
DBI::dbWriteTable(con, "cohort_ids", cohort_ids, overwrite = TRUE)

view_sql <- c(
  sprintf(
    "CREATE VIEW mimiciv_icu.icustays AS SELECT x.* FROM read_csv_auto('%s', union_by_name=true) x INNER JOIN cohort_ids c ON x.stay_id=c.stay_id",
    sql_path(raw_paths$icustays)
  ),
  sprintf(
    "CREATE VIEW mimiciv_icu.chartevents AS SELECT x.* FROM read_csv_auto('%s', union_by_name=true) x INNER JOIN cohort_ids c ON x.stay_id=c.stay_id",
    sql_path(raw_paths$chartevents)
  ),
  sprintf(
    "CREATE VIEW mimiciv_icu.outputevents AS SELECT x.* FROM read_csv_auto('%s', union_by_name=true) x INNER JOIN cohort_ids c ON x.stay_id=c.stay_id",
    sql_path(raw_paths$outputevents)
  ),
  sprintf(
    "CREATE VIEW mimiciv_hosp.admissions AS SELECT x.* FROM read_csv_auto('%s', union_by_name=true) x WHERE x.hadm_id IN (SELECT hadm_id FROM cohort_ids)",
    sql_path(raw_paths$admissions)
  ),
  sprintf(
    "CREATE VIEW mimiciv_hosp.patients AS SELECT x.* FROM read_csv_auto('%s', union_by_name=true) x WHERE x.subject_id IN (SELECT subject_id FROM cohort_ids)",
    sql_path(raw_paths$patients)
  ),
  sprintf(
    "CREATE VIEW mimiciv_hosp.services AS SELECT x.* FROM read_csv_auto('%s', union_by_name=true) x WHERE x.hadm_id IN (SELECT hadm_id FROM cohort_ids)",
    sql_path(raw_paths$services)
  )
)
invisible(lapply(view_sql, function(sql) DBI::dbExecute(con, sql)))

cat("Running MIT-LCP MIMIC Code v3.0.1 OASIS dependency chain...\n")
invisible(lapply(official_paths, function(path) execute_official_sql(con, path)))
official_oasis <- DBI::dbGetQuery(
  con,
  "SELECT * FROM mimiciv_derived.oasis ORDER BY stay_id"
) |>
  tibble::as_tibble()

if (nrow(official_oasis) != nrow(analysis) || anyDuplicated(official_oasis$stay_id)) {
  stop(
    "Official OASIS output does not contain exactly one row per analysis stay: ",
    nrow(official_oasis), " rows.", call. = FALSE
  )
}

provenance <- tibble::tibble(
  mimic_code_release = "v3.0.1",
  relative_script = official_scripts,
  md5 = unname(tools::md5sum(official_paths))
)

analysis_official <- analysis |>
  dplyr::left_join(official_oasis, by = c("subject_id", "hadm_id", "stay_id"))
if (anyNA(analysis_official$oasis)) {
  stop("Official OASIS score is missing after the cohort join.", call. = FALSE)
}
analysis_official <- analysis_official |>
  dplyr::mutate(
    oasis_quartile = factor(
      stable_ntile(oasis, 4, stay_id),
      levels = 1:4, labels = paste0("Q", 1:4)
    )
  )
reordered_oasis <- analysis_official[order(analysis_official$stay_id, decreasing = TRUE), ]
reordered_oasis$quartile_check <- factor(
  stable_ntile(reordered_oasis$oasis, 4, reordered_oasis$stay_id),
  levels = 1:4, labels = paste0("Q", 1:4)
)
oasis_order_invariant <- identical(
  as.character(analysis_official$oasis_quartile),
  as.character(reordered_oasis$quartile_check[
    match(analysis_official$stay_id, reordered_oasis$stay_id)
  ])
)
if (!oasis_order_invariant) {
  stop("Deterministic official OASIS quartile assignment failed the row-order audit.", call. = FALSE)
}
oasis_quartile_qa <- tibble::tibble(
  rule = "Ascending official OASIS; stay_id breaks score ties",
  row_order_invariant = oasis_order_invariant,
  q1_n = sum(analysis_official$oasis_quartile == "Q1"),
  q2_n = sum(analysis_official$oasis_quartile == "Q2"),
  q3_n = sum(analysis_official$oasis_quartile == "Q3"),
  q4_n = sum(analysis_official$oasis_quartile == "Q4")
)

component_availability <- tibble::tibble(
  component = c(
    "age", "pre-ICU length of stay", "GCS", "heart rate", "mean blood pressure",
    "respiratory rate", "temperature", "urine output", "mechanical ventilation",
    "elective surgery", "official OASIS"
  ),
  available_n = c(
    sum(!is.na(analysis_official$age)),
    sum(!is.na(analysis_official$preiculos)),
    sum(!is.na(analysis_official$gcs)),
    sum(!is.na(analysis_official$heartrate)),
    sum(!is.na(analysis_official$meanbp)),
    sum(!is.na(analysis_official$resprate)),
    sum(!is.na(analysis_official$temp)),
    sum(!is.na(analysis_official$urineoutput)),
    sum(!is.na(analysis_official$mechvent)),
    sum(!is.na(analysis_official$electivesurgery)),
    sum(!is.na(analysis_official$oasis))
  ),
  denominator_n = nrow(analysis_official)
) |>
  dplyr::mutate(available_percent = 100 * available_n / denominator_n)

old_oasis <- readr::read_csv(old_oasis_path, show_col_types = FALSE) |>
  dplyr::select(stay_id, oasis_like_score, oasis_like_complete) |>
  dplyr::distinct(stay_id, .keep_all = TRUE)
score_comparison_data <- analysis_official |>
  dplyr::select(stay_id, official_oasis = oasis) |>
  dplyr::left_join(old_oasis, by = "stay_id") |>
  dplyr::filter(oasis_like_complete == 1, is.finite(oasis_like_score)) |>
  dplyr::mutate(
    official_quartile = stable_ntile(official_oasis, 4, stay_id),
    oasis_like_quartile = stable_ntile(oasis_like_score, 4, stay_id)
  )
score_comparison <- tibble::tibble(
  official_n = nrow(analysis_official),
  oasis_like_complete_n = nrow(score_comparison_data),
  official_median = stats::median(analysis_official$oasis),
  official_q1 = stats::quantile(analysis_official$oasis, 0.25),
  official_q3 = stats::quantile(analysis_official$oasis, 0.75),
  oasis_like_median = stats::median(score_comparison_data$oasis_like_score),
  spearman_correlation = stats::cor(
    score_comparison_data$official_oasis, score_comparison_data$oasis_like_score,
    method = "spearman"
  ),
  pearson_correlation = stats::cor(
    score_comparison_data$official_oasis, score_comparison_data$oasis_like_score
  ),
  exact_quartile_agreement = mean(
    score_comparison_data$official_quartile == score_comparison_data$oasis_like_quartile
  )
)

cat("Fitting official-OASIS adjusted outcome models...\n")
fit_continuous <- survival::coxph(
  survival::Surv(survival_days_365, mortality_365d) ~ phenotype + male + oasis,
  data = analysis_official, ties = "efron"
)
fit_stratified <- survival::coxph(
  survival::Surv(survival_days_365, mortality_365d) ~ phenotype + male + strata(oasis_quartile),
  data = analysis_official, ties = "efron"
)
landmark <- analysis_official |>
  dplyr::filter(survival_days_365 > 1) |>
  dplyr::mutate(
    landmark_time = survival_days_365 - 1,
    landmark_death = mortality_365d
  )
fit_landmark <- survival::coxph(
  survival::Surv(landmark_time, landmark_death) ~ phenotype + male + strata(oasis_quartile),
  data = landmark, ties = "efron"
)
fit_hospital <- stats::glm(
  hospital_expire_flag ~ phenotype + male + oasis,
  data = analysis_official, family = stats::binomial()
)

denominator <- readr::read_csv(denominator_path, show_col_types = FALSE) |>
  dplyr::mutate(
    selected_primary = as.integer(stay_id %in% analysis_official$stay_id),
    gender_selection = collapse_for_selection(gender, selected_primary),
    race_selection = collapse_for_selection(race, selected_primary),
    careunit_selection = collapse_for_selection(first_careunit, selected_primary),
    year_selection = collapse_for_selection(anchor_year_group, selected_primary)
  )
selection_fit <- stats::glm(
  selected_primary ~ splines::ns(anchor_age, df = 3) + gender_selection +
    race_selection + careunit_selection + year_selection,
  data = denominator,
  family = stats::binomial(),
  control = stats::glm.control(maxit = 100)
)
if (!isTRUE(selection_fit$converged)) {
  stop("Official-OASIS selection model did not converge.", call. = FALSE)
}
probability <- pmin(pmax(stats::predict(selection_fit, type = "response"), 0.001), 0.999)
selection_rate <- mean(denominator$selected_primary == 1)
selection_weight <- ifelse(
  denominator$selected_primary == 1, selection_rate / probability, NA_real_
)
trim_limits <- stats::quantile(
  selection_weight[denominator$selected_primary == 1], c(0.01, 0.99), na.rm = TRUE
)
denominator$selection_ipw <- ifelse(
  denominator$selected_primary == 1,
  pmin(pmax(selection_weight, trim_limits[1]), trim_limits[2]),
  NA_real_
)
analysis_official <- analysis_official |>
  dplyr::left_join(
    denominator |>
      dplyr::select(stay_id, selection_ipw),
    by = "stay_id"
  ) |>
  dplyr::mutate(selection_ipw_scaled = selection_ipw / mean(selection_ipw))
fit_selection_ipw <- survival::coxph(
  survival::Surv(survival_days_365, mortality_365d) ~ phenotype + male + strata(oasis_quartile),
  data = analysis_official,
  weights = selection_ipw_scaled,
  robust = TRUE,
  cluster = stay_id,
  ties = "efron"
)

outcome_models <- dplyr::bind_rows(
  extract_cox(
    fit_continuous, "Official OASIS continuous adjusted Cox",
    nrow(analysis_official), sum(analysis_official$mortality_365d == 1)
  ),
  extract_cox(
    fit_stratified, "Official OASIS quartile-stratified Cox",
    nrow(analysis_official), sum(analysis_official$mortality_365d == 1)
  ),
  extract_cox(
    fit_selection_ipw, "Selection-IPW official OASIS quartile-stratified Cox",
    nrow(analysis_official), sum(analysis_official$mortality_365d == 1)
  ),
  extract_cox(
    fit_landmark, "24-hour landmark official OASIS quartile-stratified Cox",
    nrow(landmark), sum(landmark$landmark_death == 1)
  ),
  extract_logistic(
    fit_hospital, "Hospital mortality official OASIS adjusted logistic model",
    nrow(analysis_official), sum(analysis_official$hospital_expire_flag == 1)
  )
)

ph_continuous <- survival::cox.zph(fit_continuous)$table |>
  as.data.frame() |>
  tibble::rownames_to_column("term") |>
  dplyr::mutate(model = "Official OASIS continuous adjusted Cox", .before = 1)
ph_stratified <- survival::cox.zph(fit_stratified)$table |>
  as.data.frame() |>
  tibble::rownames_to_column("term") |>
  dplyr::mutate(model = "Official OASIS quartile-stratified Cox", .before = 1)
ph_checks <- dplyr::bind_rows(ph_continuous, ph_stratified)

selection_diagnostics <- tibble::tibble(
  denominator_n = nrow(denominator),
  feature_complete_and_official_oasis_n = nrow(analysis_official),
  selected_percent = 100 * selection_rate,
  selected_probability_min = min(probability[denominator$selected_primary == 1]),
  selected_probability_median = stats::median(probability[denominator$selected_primary == 1]),
  selected_probability_max = max(probability[denominator$selected_primary == 1]),
  weight_1st_percentile = trim_limits[1],
  weight_99th_percentile = trim_limits[2],
  weight_max = max(denominator$selection_ipw, na.rm = TRUE),
  effective_sample_size = sum(analysis_official$selection_ipw)^2 /
    sum(analysis_official$selection_ipw^2)
)

readr::write_csv(provenance, file.path(output_dir, "Table41A_official_SQL_provenance.csv"))
readr::write_csv(component_availability, file.path(output_dir, "Table41B_official_OASIS_component_availability.csv"))
readr::write_csv(score_comparison, file.path(output_dir, "Table41C_official_vs_OASISlike_comparison.csv"))
readr::write_csv(outcome_models, file.path(output_dir, "Table41D_official_OASIS_outcome_models.csv"))
readr::write_csv(ph_checks, file.path(output_dir, "Table41E_official_OASIS_PH_checks.csv"))
readr::write_csv(selection_diagnostics, file.path(output_dir, "Table41F_official_OASIS_selection_weights.csv"))
readr::write_csv(oasis_quartile_qa, file.path(output_dir, "Table41G_official_OASIS_quartile_QA.csv"))
readr::write_csv(
  official_oasis,
  file.path(output_dir, "MIMIC_official_OASIS_v301_scores.csv")
)
saveRDS(
  list(
    analysis = analysis_official,
    official_oasis = official_oasis,
    provenance = provenance,
    component_availability = component_availability,
    score_comparison = score_comparison,
    outcome_models = outcome_models,
    ph_checks = ph_checks,
    selection_diagnostics = selection_diagnostics
    , oasis_quartile_qa = oasis_quartile_qa
  ),
  file.path(output_dir, "MIMIC_official_OASIS_v301_results.rds")
)

summary_lines <- c(
  "MIMIC-IV official OASIS reconstruction using MIT-LCP MIMIC Code v3.0.1",
  "",
  paste0("Strict first-24-hour phenotype cohort: n = ", nrow(analysis_official)),
  paste0("Official OASIS median (IQR): ",
         stats::median(analysis_official$oasis), " (",
         stats::quantile(analysis_official$oasis, 0.25), "-",
         stats::quantile(analysis_official$oasis, 0.75), ")"),
  "",
  "Official versus previous OASIS-like score:",
  paste(capture.output(print(score_comparison)), collapse = "\n"),
  "",
  "Outcome models:",
  paste(capture.output(print(outcome_models)), collapse = "\n"),
  "",
  "Proportional-hazards checks:",
  paste(capture.output(print(ph_checks)), collapse = "\n"),
  "",
  "Selection diagnostics:",
  paste(capture.output(print(selection_diagnostics)), collapse = "\n"),
  "",
  "Deterministic OASIS quartile audit:",
  paste(capture.output(print(oasis_quartile_qa)), collapse = "\n"),
  "",
  "Interpretation: the official MIT-LCP SQL score supersedes the locally mapped OASIS-like score for MIMIC severity adjustment, subject to successful QA."
)
writeLines(summary_lines, file.path(output_dir, "MIMIC_official_OASIS_v301_summary.txt"))
cat(paste(summary_lines, collapse = "\n"), "\n")
