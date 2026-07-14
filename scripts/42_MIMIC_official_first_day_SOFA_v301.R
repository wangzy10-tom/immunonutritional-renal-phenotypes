# ==============================================================================
# MIT-LCP MIMIC Code v3.0.1 official first-day SOFA sensitivity analysis
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
output_dir <- file.path(root, "output", "mimic_official_sofa_v301")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

oasis_results_path <- file.path(
  root, "output", "mimic_official_oasis_v301", "MIMIC_official_OASIS_v301_results.rds"
)
raw_paths <- list(
  labevents = file.path(mimic_root, "hosp", "labevents.csv.gz"),
  icustays = file.path(mimic_root, "icu", "icustays.csv.gz"),
  chartevents = file.path(mimic_root, "icu", "chartevents.csv.gz"),
  outputevents = file.path(mimic_root, "icu", "outputevents.csv.gz"),
  inputevents = file.path(mimic_root, "icu", "inputevents.csv.gz")
)

official_scripts <- c(
  "measurement/gcs.sql",
  "measurement/vitalsign.sql",
  "measurement/urine_output.sql",
  "measurement/ventilator_setting.sql",
  "measurement/oxygen_delivery.sql",
  "measurement/bg.sql",
  "measurement/blood_differential.sql",
  "measurement/chemistry.sql",
  "measurement/coagulation.sql",
  "measurement/complete_blood_count.sql",
  "measurement/enzyme.sql",
  "medication/dobutamine.sql",
  "medication/dopamine.sql",
  "medication/epinephrine.sql",
  "medication/norepinephrine.sql",
  "treatment/ventilation.sql",
  "firstday/first_day_gcs.sql",
  "firstday/first_day_vitalsign.sql",
  "firstday/first_day_urine_output.sql",
  "firstday/first_day_lab.sql",
  "firstday/first_day_sofa.sql"
)
official_paths <- file.path(official_root, official_scripts)
required_inputs <- c(oasis_results_path, unlist(raw_paths), official_paths)
if (any(!file.exists(required_inputs))) {
  stop(
    "Missing required input(s): ",
    paste(required_inputs[!file.exists(required_inputs)], collapse = "; "),
    call. = FALSE
  )
}

stable_ntile <- function(value, groups, tie_breaker) {
  keep <- is.finite(value) & !is.na(tie_breaker)
  output <- rep(NA_integer_, length(value))
  ordered_index <- which(keep)[order(value[keep], tie_breaker[keep])]
  output[ordered_index] <- dplyr::ntile(seq_along(ordered_index), groups)
  output
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

cat("Loading official-OASIS MIMIC cohort and fixed selection weights...\n")
oasis_results <- readRDS(oasis_results_path)
analysis <- oasis_results$analysis |>
  dplyr::mutate(
    stay_id = as.integer(stay_id),
    hadm_id = as.integer(hadm_id),
    subject_id = as.integer(subject_id),
    phenotype = factor(as.character(phenotype), levels = c("P3", "P2", "P1")),
    male = as.integer(gender == "M")
  )
if (nrow(analysis) != 1145L || anyDuplicated(analysis$stay_id)) {
  stop("Unexpected MIMIC strict first-24-hour analysis cohort.", call. = FALSE)
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
    "CREATE VIEW mimiciv_icu.inputevents AS SELECT x.* FROM read_csv_auto('%s', union_by_name=true) x INNER JOIN cohort_ids c ON x.stay_id=c.stay_id",
    sql_path(raw_paths$inputevents)
  ),
  sprintf(
    "CREATE VIEW mimiciv_hosp.labevents AS SELECT x.* FROM read_csv_auto('%s', union_by_name=true) x WHERE x.subject_id IN (SELECT subject_id FROM cohort_ids)",
    sql_path(raw_paths$labevents)
  )
)
invisible(lapply(view_sql, function(sql) DBI::dbExecute(con, sql)))

cat("Running MIT-LCP MIMIC Code v3.0.1 first-day SOFA dependency chain...\n")
invisible(lapply(official_paths, function(path) execute_official_sql(con, path)))
official_sofa <- DBI::dbGetQuery(
  con,
  "SELECT * FROM mimiciv_derived.first_day_sofa ORDER BY stay_id"
) |>
  tibble::as_tibble()
if (nrow(official_sofa) != nrow(analysis) || anyDuplicated(official_sofa$stay_id)) {
  stop(
    "Official first-day SOFA output does not contain one row per stay: ",
    nrow(official_sofa), " rows.", call. = FALSE
  )
}

provenance <- tibble::tibble(
  mimic_code_release = "v3.0.1",
  relative_script = official_scripts,
  md5 = unname(tools::md5sum(official_paths))
)

analysis_sofa <- analysis |>
  dplyr::left_join(official_sofa, by = c("subject_id", "hadm_id", "stay_id"))
if (anyNA(analysis_sofa$sofa)) {
  stop("Official first-day SOFA score is missing after the cohort join.", call. = FALSE)
}
analysis_sofa <- analysis_sofa |>
  dplyr::mutate(
    sofa_quartile = factor(
      stable_ntile(sofa, 4, stay_id),
      levels = 1:4, labels = paste0("Q", 1:4)
    ),
    selection_ipw_scaled = selection_ipw / mean(selection_ipw)
  )
reordered_sofa <- analysis_sofa[order(analysis_sofa$stay_id, decreasing = TRUE), ]
reordered_sofa$quartile_check <- factor(
  stable_ntile(reordered_sofa$sofa, 4, reordered_sofa$stay_id),
  levels = 1:4, labels = paste0("Q", 1:4)
)
sofa_order_invariant <- identical(
  as.character(analysis_sofa$sofa_quartile),
  as.character(reordered_sofa$quartile_check[
    match(analysis_sofa$stay_id, reordered_sofa$stay_id)
  ])
)
if (!sofa_order_invariant) {
  stop("Deterministic first-day SOFA quartile assignment failed the row-order audit.", call. = FALSE)
}
sofa_quartile_qa <- tibble::tibble(
  rule = "Ascending official first-day SOFA; stay_id breaks score ties",
  row_order_invariant = sofa_order_invariant,
  q1_n = sum(analysis_sofa$sofa_quartile == "Q1"),
  q2_n = sum(analysis_sofa$sofa_quartile == "Q2"),
  q3_n = sum(analysis_sofa$sofa_quartile == "Q3"),
  q4_n = sum(analysis_sofa$sofa_quartile == "Q4")
)

component_availability <- tibble::tibble(
  component = c("respiration", "coagulation", "liver", "cardiovascular", "CNS", "renal", "official first-day SOFA"),
  available_n = c(
    sum(!is.na(analysis_sofa$respiration)),
    sum(!is.na(analysis_sofa$coagulation)),
    sum(!is.na(analysis_sofa$liver)),
    sum(!is.na(analysis_sofa$cardiovascular)),
    sum(!is.na(analysis_sofa$cns)),
    sum(!is.na(analysis_sofa$renal)),
    sum(!is.na(analysis_sofa$sofa))
  ),
  denominator_n = nrow(analysis_sofa)
) |>
  dplyr::mutate(available_percent = 100 * available_n / denominator_n)

score_summary <- analysis_sofa |>
  dplyr::group_by(phenotype) |>
  dplyr::summarise(
    n = dplyr::n(),
    sofa_median = stats::median(sofa),
    sofa_q1 = stats::quantile(sofa, 0.25),
    sofa_q3 = stats::quantile(sofa, 0.75),
    oasis_median = stats::median(oasis),
    mortality_365_percent = 100 * mean(mortality_365d == 1),
    .groups = "drop"
  )
score_comparison <- tibble::tibble(
  n = nrow(analysis_sofa),
  sofa_median = stats::median(analysis_sofa$sofa),
  sofa_q1 = stats::quantile(analysis_sofa$sofa, 0.25),
  sofa_q3 = stats::quantile(analysis_sofa$sofa, 0.75),
  oasis_median = stats::median(analysis_sofa$oasis),
  spearman_sofa_oasis = stats::cor(analysis_sofa$sofa, analysis_sofa$oasis, method = "spearman"),
  pearson_sofa_oasis = stats::cor(analysis_sofa$sofa, analysis_sofa$oasis)
)

cat("Fitting official first-day SOFA sensitivity models...\n")
fit_continuous <- survival::coxph(
  survival::Surv(survival_days_365, mortality_365d) ~
    phenotype + anchor_age + male + sofa,
  data = analysis_sofa, ties = "efron"
)
fit_stratified <- survival::coxph(
  survival::Surv(survival_days_365, mortality_365d) ~
    phenotype + anchor_age + male + strata(sofa_quartile),
  data = analysis_sofa, ties = "efron"
)
fit_selection_ipw <- survival::coxph(
  survival::Surv(survival_days_365, mortality_365d) ~
    phenotype + anchor_age + male + strata(sofa_quartile),
  data = analysis_sofa,
  weights = selection_ipw_scaled,
  robust = TRUE,
  cluster = stay_id,
  ties = "efron"
)
landmark <- analysis_sofa |>
  dplyr::filter(survival_days_365 > 1) |>
  dplyr::mutate(
    landmark_time = survival_days_365 - 1,
    landmark_death = mortality_365d
  )
fit_landmark <- survival::coxph(
  survival::Surv(landmark_time, landmark_death) ~
    phenotype + anchor_age + male + strata(sofa_quartile),
  data = landmark, ties = "efron"
)
fit_hospital <- stats::glm(
  hospital_expire_flag ~ phenotype + anchor_age + male + sofa,
  data = analysis_sofa, family = stats::binomial()
)

outcome_models <- dplyr::bind_rows(
  extract_cox(
    fit_continuous, "Official first-day SOFA continuous adjusted Cox",
    nrow(analysis_sofa), sum(analysis_sofa$mortality_365d == 1)
  ),
  extract_cox(
    fit_stratified, "Official first-day SOFA quartile-stratified Cox",
    nrow(analysis_sofa), sum(analysis_sofa$mortality_365d == 1)
  ),
  extract_cox(
    fit_selection_ipw, "Selection-IPW official first-day SOFA quartile-stratified Cox",
    nrow(analysis_sofa), sum(analysis_sofa$mortality_365d == 1)
  ),
  extract_cox(
    fit_landmark, "24-hour landmark official first-day SOFA quartile-stratified Cox",
    nrow(landmark), sum(landmark$landmark_death == 1)
  ),
  extract_logistic(
    fit_hospital, "Hospital mortality official first-day SOFA adjusted logistic model",
    nrow(analysis_sofa), sum(analysis_sofa$hospital_expire_flag == 1)
  )
)

ph_continuous <- survival::cox.zph(fit_continuous)$table |>
  as.data.frame() |>
  tibble::rownames_to_column("term") |>
  dplyr::mutate(model = "Official first-day SOFA continuous adjusted Cox", .before = 1)
ph_stratified <- survival::cox.zph(fit_stratified)$table |>
  as.data.frame() |>
  tibble::rownames_to_column("term") |>
  dplyr::mutate(model = "Official first-day SOFA quartile-stratified Cox", .before = 1)
ph_checks <- dplyr::bind_rows(ph_continuous, ph_stratified)

readr::write_csv(provenance, file.path(output_dir, "Table42A_official_SOFA_SQL_provenance.csv"))
readr::write_csv(component_availability, file.path(output_dir, "Table42B_official_SOFA_component_availability.csv"))
readr::write_csv(score_summary, file.path(output_dir, "Table42C_official_SOFA_by_phenotype.csv"))
readr::write_csv(score_comparison, file.path(output_dir, "Table42D_SOFA_vs_OASIS_comparison.csv"))
readr::write_csv(outcome_models, file.path(output_dir, "Table42E_official_SOFA_outcome_models.csv"))
readr::write_csv(ph_checks, file.path(output_dir, "Table42F_official_SOFA_PH_checks.csv"))
readr::write_csv(sofa_quartile_qa, file.path(output_dir, "Table42H_official_SOFA_quartile_QA.csv"))
readr::write_csv(official_sofa, file.path(output_dir, "MIMIC_official_first_day_SOFA_v301_scores.csv"))
saveRDS(
  list(
    analysis = analysis_sofa,
    official_sofa = official_sofa,
    provenance = provenance,
    component_availability = component_availability,
    score_summary = score_summary,
    score_comparison = score_comparison,
    outcome_models = outcome_models,
    sofa_quartile_qa = sofa_quartile_qa,
    ph_checks = ph_checks
  ),
  file.path(output_dir, "MIMIC_official_first_day_SOFA_v301_results.rds")
)

summary_lines <- c(
  "MIMIC-IV official first-day SOFA sensitivity using MIT-LCP MIMIC Code v3.0.1",
  "",
  paste0("Strict first-24-hour phenotype cohort: n = ", nrow(analysis_sofa)),
  paste0(
    "Official first-day SOFA median (IQR): ", stats::median(analysis_sofa$sofa),
    " (", stats::quantile(analysis_sofa$sofa, 0.25), "-",
    stats::quantile(analysis_sofa$sofa, 0.75), ")"
  ),
  "",
  "SOFA versus OASIS:",
  paste(capture.output(print(score_comparison)), collapse = "\n"),
  "",
  "SOFA by phenotype:",
  paste(capture.output(print(score_summary)), collapse = "\n"),
  "",
  "Outcome models:",
  paste(capture.output(print(outcome_models)), collapse = "\n"),
  "",
  "Deterministic SOFA quartile audit:",
  paste(capture.output(print(sofa_quartile_qa)), collapse = "\n"),
  "",
  "Proportional-hazards checks:",
  paste(capture.output(print(ph_checks)), collapse = "\n"),
  "",
  "Interpretation: official first-day SOFA is a post-freeze overadjustment-sensitive robustness analysis and does not replace official OASIS."
)
writeLines(summary_lines, file.path(output_dir, "MIMIC_official_first_day_SOFA_v301_summary.txt"))
cat(paste(summary_lines, collapse = "\n"), "\n")
