# ==============================================================================
# Build corrected, analysis-freeze-compliant publication tables
# ==============================================================================

required_packages <- c("dplyr", "readr", "tibble", "survival", "survey")
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
output_dir <- file.path(root, "output", "corrected_publication_tables")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cohort_path <- file.path(
  root, "output", "nhanes_covariate_upgrade",
  "NHANES_2011_2018_covariate_augmented.rds"
)
assignment_path <- file.path(
  root, "output", "nhanes_robust_reanalysis",
  "NHANES_2011_2018_variant_assignments.csv"
)
benchmark_path <- file.path(
  root, "output", "nhanes_albumin_benchmarks",
  "NHANES_albumin_benchmark_results.rds"
)

format_p <- function(p) {
  ifelse(is.na(p), NA_character_, ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
}

cohort <- readRDS(cohort_path)
assignments <- readr::read_csv(assignment_path, show_col_types = FALSE) |>
  dplyr::select(SEQN, phenotype_log6_winsor)

table1_data <- cohort |>
  dplyr::left_join(assignments, by = "SEQN") |>
  dplyr::mutate(
    phenotype = factor(phenotype_log6_winsor, levels = c("P3", "P2", "P1")),
    p1_indicator = as.numeric(phenotype == "P1"),
    p2_indicator = as.numeric(phenotype == "P2"),
    p3_indicator = as.numeric(phenotype == "P3"),
    male_indicator = as.numeric(RIAGENDR == 1),
    race_mexican = as.numeric(RIDRETH3 == 1),
    race_other_hispanic = as.numeric(RIDRETH3 == 2),
    race_nhw = as.numeric(RIDRETH3 == 3),
    race_nhb = as.numeric(RIDRETH3 == 4),
    race_nha = as.numeric(RIDRETH3 == 6),
    race_other = as.numeric(RIDRETH3 == 7),
    smoking_former = as.numeric(smoking == "Former"),
    smoking_current = as.numeric(smoking == "Current"),
    hypertension_yes = as.numeric(hypertension == "Yes"),
    diabetes_borderline = as.numeric(diabetes == "Borderline"),
    diabetes_yes = as.numeric(diabetes == "Yes"),
    death_indicator = as.numeric(MORTSTAT == 1)
  )

table1_design <- survey::svydesign(
  ids = ~SDMVPSU,
  strata = ~SDMVSTRA,
  weights = ~WTMEC8YR,
  nest = TRUE,
  data = table1_data
)

group_order <- c("Overall", "P1", "P2", "P3")
designs <- list(
  Overall = table1_design,
  P1 = subset(table1_design, phenotype == "P1"),
  P2 = subset(table1_design, phenotype == "P2"),
  P3 = subset(table1_design, phenotype == "P3")
)
raw_groups <- list(
  Overall = table1_data,
  P1 = table1_data |> dplyr::filter(phenotype == "P1"),
  P2 = table1_data |> dplyr::filter(phenotype == "P2"),
  P3 = table1_data |> dplyr::filter(phenotype == "P3")
)

weighted_mean_sd <- function(design, variable, digits = 2) {
  form <- stats::as.formula(paste0("~", variable))
  mean_value <- unname(stats::coef(survey::svymean(form, design, na.rm = TRUE))[1])
  variance <- unname(stats::coef(survey::svyvar(form, design, na.rm = TRUE))[1])
  sprintf(paste0("%.", digits, "f (%.", digits, "f)"), mean_value, sqrt(variance))
}

weighted_median_iqr <- function(design, variable, digits = 2) {
  form <- stats::as.formula(paste0("~", variable))
  values <- as.numeric(
    survey::svyquantile(form, design, c(0.25, 0.50, 0.75), ci = FALSE, na.rm = TRUE)[[1]]
  )
  sprintf(
    paste0("%.", digits, "f (%.", digits, "f-%.", digits, "f)"),
    values[2], values[1], values[3]
  )
}

weighted_percent <- function(design, variable) {
  form <- stats::as.formula(paste0("~", variable))
  value <- unname(stats::coef(survey::svymean(form, design, na.rm = TRUE))[1])
  sprintf("%.1f%%", 100 * value)
}

count_weighted_percent <- function(raw_data, design, variable) {
  count <- sum(raw_data[[variable]] == 1, na.rm = TRUE)
  percent <- unname(stats::coef(
    survey::svymean(stats::as.formula(paste0("~", variable)), design, na.rm = TRUE)
  )[1])
  sprintf("%d (%.1f%%)", count, 100 * percent)
}

population_percent <- c(
  Overall = "100.0%",
  P1 = weighted_percent(table1_design, "p1_indicator"),
  P2 = weighted_percent(table1_design, "p2_indicator"),
  P3 = weighted_percent(table1_design, "p3_indicator")
)

table1_rows <- list(
  list("Participants, unweighted n", function(g) as.character(nrow(raw_groups[[g]]))),
  list("Weighted population proportion", function(g) population_percent[[g]]),
  list("Age, years, mean (SD)", function(g) weighted_mean_sd(designs[[g]], "RIDAGEYR", 1)),
  list("Male, weighted %", function(g) weighted_percent(designs[[g]], "male_indicator")),
  list("Mexican American, weighted %", function(g) weighted_percent(designs[[g]], "race_mexican")),
  list("Other Hispanic, weighted %", function(g) weighted_percent(designs[[g]], "race_other_hispanic")),
  list("Non-Hispanic White, weighted %", function(g) weighted_percent(designs[[g]], "race_nhw")),
  list("Non-Hispanic Black, weighted %", function(g) weighted_percent(designs[[g]], "race_nhb")),
  list("Non-Hispanic Asian, weighted %", function(g) weighted_percent(designs[[g]], "race_nha")),
  list("Other/multiracial, weighted %", function(g) weighted_percent(designs[[g]], "race_other")),
  list("Former smoker, weighted %", function(g) weighted_percent(designs[[g]], "smoking_former")),
  list("Current smoker, weighted %", function(g) weighted_percent(designs[[g]], "smoking_current")),
  list("Hypertension, weighted %", function(g) weighted_percent(designs[[g]], "hypertension_yes")),
  list("Borderline diabetes, weighted %", function(g) weighted_percent(designs[[g]], "diabetes_borderline")),
  list("Diabetes, weighted %", function(g) weighted_percent(designs[[g]], "diabetes_yes")),
  list("Extended comorbidity score, mean (SD)", function(g) weighted_mean_sd(designs[[g]], "Comorbidity_Score_Extended", 2)),
  list("NLR, median (IQR)", function(g) weighted_median_iqr(designs[[g]], "NLR", 2)),
  list("SII, median (IQR)", function(g) weighted_median_iqr(designs[[g]], "SII", 1)),
  list("Haemoglobin, g/dL, mean (SD)", function(g) weighted_mean_sd(designs[[g]], "LBXHGB", 2)),
  list("Serum total protein, g/dL, mean (SD)", function(g) weighted_mean_sd(designs[[g]], "LBXSTP", 2)),
  list("BMI, kg/m2, mean (SD)", function(g) weighted_mean_sd(designs[[g]], "BMXBMI", 2)),
  list("Serum creatinine, mg/dL, median (IQR)", function(g) weighted_median_iqr(designs[[g]], "LBXSCR", 2)),
  list("Deaths during follow-up, n (weighted %)", function(g) count_weighted_percent(raw_groups[[g]], designs[[g]], "death_indicator")),
  list("Follow-up, months, median (IQR)", function(g) weighted_median_iqr(designs[[g]], "PERMTH_INT", 1))
)

table1 <- dplyr::bind_rows(lapply(table1_rows, function(specification) {
  values <- vapply(group_order, specification[[2]], character(1))
  tibble::tibble(
    Characteristic = specification[[1]],
    Overall = values[["Overall"]],
    P1 = values[["P1"]],
    P2 = values[["P2"]],
    P3 = values[["P3"]]
  )
}))

# Main NHANES model table -------------------------------------------------------
model_data <- readRDS(benchmark_path)$model_data |>
  dplyr::mutate(
    phenotype_total_protein = factor(phenotype_total_protein, levels = c("P3", "P2", "P1"))
  )

minimum_data <- table1_data |>
  dplyr::filter(!is.na(phenotype), !is.na(PERMTH_INT), !is.na(MORTSTAT), !is.na(RIDAGEYR), !is.na(RIAGENDR)) |>
  dplyr::mutate(male = as.integer(RIAGENDR == 1), phenotype_total_protein = phenotype)

minimum_fit <- survival::coxph(
  survival::Surv(PERMTH_INT, MORTSTAT) ~ phenotype_total_protein + RIDAGEYR + male,
  data = minimum_data,
  ties = "efron"
)
full_formula <- survival::Surv(PERMTH_INT, MORTSTAT) ~
  phenotype_total_protein + RIDAGEYR + male + race + INDFMPIR +
  Comorbidity_Score_Extended + cycle + smoking + hypertension + diabetes
full_fit <- survival::coxph(full_formula, data = model_data, ties = "efron")
primary_design <- survey::svydesign(
  ids = ~SDMVPSU,
  strata = ~SDMVSTRA,
  weights = ~WTMEC8YR,
  nest = TRUE,
  data = model_data
)
primary_fit <- survey::svycoxph(full_formula, design = primary_design)

extract_phenotype <- function(fit, model_label, n_model, events, model_role) {
  beta <- stats::coef(fit)
  standard_error <- sqrt(diag(stats::vcov(fit)))
  terms <- names(beta)
  keep <- grepl("^phenotype_total_protein", terms)
  tibble::tibble(
    model_role = model_role,
    model = model_label,
    n = n_model,
    events = events,
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
      effect_95ci = sprintf("%.2f (%.2f-%.2f)", HR, lower_95, upper_95),
      p_value_formatted = format_p(p_value)
    )
}

table2_core <- dplyr::bind_rows(
  extract_phenotype(
    minimum_fit, "Age- and sex-adjusted Cox", nrow(minimum_data), sum(minimum_data$MORTSTAT == 1),
    "Secondary"
  ),
  extract_phenotype(
    full_fit, "Fully adjusted Cox", nrow(model_data), sum(model_data$MORTSTAT == 1),
    "Key consistency"
  ),
  extract_phenotype(
    primary_fit, "Complex-survey fully adjusted Cox", nrow(model_data), sum(model_data$MORTSTAT == 1),
    "Primary"
  )
)

selection_rows <- readr::read_csv(
  file.path(root, "output", "nhanes_selection_bias", "Table33D_selection_IPW_Cox.csv"),
  show_col_types = FALSE
) |>
  dplyr::filter(model == "NHANES survey weights × selection IPW") |>
  dplyr::transmute(
    model_role = "Sensitivity",
    model = "Complex survey x selection IPW",
    n = nrow(model_data),
    events = sum(model_data$MORTSTAT == 1),
    comparison,
    HR,
    lower_95,
    upper_95,
    p_value,
    effect_95ci = sprintf("%.2f (%.2f-%.2f)", HR, lower_95, upper_95),
    p_value_formatted = format_p(p_value)
  )

loco_rows <- readr::read_csv(
  file.path(root, "output", "nhanes_leave_one_cycle_out", "Table35D_pooled_LOCO_Cox.csv"),
  show_col_types = FALSE
) |>
  dplyr::filter(model == "Survey-weighted pooled leave-one-cycle-out model") |>
  dplyr::transmute(
    model_role = "Sensitivity",
    model = "Leave-one-cycle-out complex-survey Cox",
    n = n_model,
    events,
    comparison,
    HR,
    lower_95,
    upper_95,
    p_value,
    effect_95ci = sprintf("%.2f (%.2f-%.2f)", HR, lower_95, upper_95),
    p_value_formatted = format_p(p_value)
  )

table2 <- dplyr::bind_rows(table2_core, selection_rows, loco_rows)

bootstrap_hr <- readr::read_csv(
  file.path(root, "output", "nhanes_corrected_bootstrap", "Table32C_cluster_aware_bootstrap_summary.csv"),
  show_col_types = FALSE
)
ph_global <- readr::read_csv(
  file.path(root, "output", "nhanes_corrected_bootstrap", "Table32A_corrected_PH_checks.csv"),
  show_col_types = FALSE
) |>
  dplyr::filter(model == "Total-protein phenotype", variable == "GLOBAL")
table2_validation <- tibble::tibble(
  validation = c("Cluster-aware bootstrap P1 HR", "Fully adjusted phenotype model global PH test"),
  repetitions_or_df = c(
    paste0(bootstrap_hr$valid_replicates, "/", bootstrap_hr$requested_replicates, " valid"),
    paste0("df = ", ph_global$df)
  ),
  result = c(
    sprintf(
      "Median HR %.2f (percentile 95%% interval %.2f-%.2f)",
      bootstrap_hr$bootstrap_median_HR, bootstrap_hr$bootstrap_lower_95, bootstrap_hr$bootstrap_upper_95
    ),
    sprintf("Global P = %.3f", ph_global$p_value)
  )
)

# Cross-database primary and sensitivity tables --------------------------------
mimic_rds <- readRDS(file.path(
  root, "output", "mimic_24h_robust_corrected_severity", "MIMIC_24h_robust_results.rds"
))
mimic_official_rds <- readRDS(file.path(
  root, "output", "mimic_official_oasis_v301", "MIMIC_official_OASIS_v301_results.rds"
))
eicu_rds <- readRDS(file.path(
  root, "output", "eicu_24h_robust", "eICU_24h_robust_results.rds"
))

primary_nhanes <- table2 |>
  dplyr::filter(model_role == "Primary")
primary_mimic <- mimic_official_rds$outcome_models |>
  dplyr::filter(model == "Official OASIS quartile-stratified Cox")
primary_eicu <- eicu_rds$models |>
  dplyr::filter(model == "Hospital mortality, APACHE-adjusted mixed effects")

effect_for <- function(data, comparison, effect_column) {
  row <- data[data$comparison == comparison, , drop = FALSE]
  effect <- row[[effect_column]][1]
  sprintf("%.2f (%.2f-%.2f)", effect, row$lower_95[1], row$upper_95[1])
}

table3_primary <- tibble::tibble(
  database = c("NHANES 2011-2018", "MIMIC-IV v3.1", "eICU v2.0"),
  setting = c("Community-dwelling older adults", "Single-centre critical care", "Multicentre critical care"),
  biomarker_window = c("Baseline examination", "ICU admission 0-24 h", "ICU admission 0-24 h"),
  n = c(nrow(model_data), nrow(mimic_official_rds$analysis), nrow(eicu_rds$model_data)),
  events = c(
    sum(model_data$MORTSTAT == 1),
    sum(mimic_official_rds$analysis$mortality_365d == 1),
    sum(eicu_rds$model_data$hospital_mortality == 1)
  ),
  outcome = c("All-cause mortality", "365-day all-cause mortality", "In-hospital mortality"),
  model = c(
    "Complex-survey fully adjusted Cox",
    "MIT-LCP v3.0.1 official OASIS quartile-stratified Cox",
    "APACHE-adjusted logistic mixed model with hospital random intercept"
  ),
  effect_measure = c("HR", "HR", "OR"),
  P1_vs_P3 = c(
    effect_for(primary_nhanes, "P1 vs P3", "HR"),
    effect_for(primary_mimic, "P1 vs P3", "estimate"),
    effect_for(primary_eicu, "P1 vs P3", "OR")
  ),
  P2_vs_P3 = c(
    effect_for(primary_nhanes, "P2 vs P3", "HR"),
    effect_for(primary_mimic, "P2 vs P3", "estimate"),
    effect_for(primary_eicu, "P2 vs P3", "OR")
  )
)

mimic_landmark <- mimic_official_rds$outcome_models |>
  dplyr::filter(model == "24-hour landmark official OASIS quartile-stratified Cox")
mimic_hospital <- mimic_official_rds$outcome_models |>
  dplyr::filter(model == "Hospital mortality official OASIS adjusted logistic model")
eicu_icu <- eicu_rds$models |>
  dplyr::filter(model == "ICU mortality, APACHE-adjusted mixed effects")
eicu_hospital_landmark <- eicu_rds$models |>
  dplyr::filter(model == "Hospital mortality, 24-hour landmark and APACHE-adjusted")
eicu_icu_landmark <- eicu_rds$models |>
  dplyr::filter(model == "ICU mortality, 24-hour landmark and APACHE-adjusted")

mimic_landmark_data <- mimic_official_rds$analysis |>
  dplyr::filter(survival_days_365 > 1)
eicu_hospital_landmark_data <- eicu_rds$model_data |>
  dplyr::filter(survival_days_hosp > 1)
eicu_icu_landmark_data <- eicu_rds$model_data |>
  dplyr::filter(icu_los_days > 1)

table3_sensitivity <- tibble::tibble(
  database = c("MIMIC-IV", "MIMIC-IV", "eICU", "eICU", "eICU"),
  analysis = c(
    "365-day mortality, 24-hour landmark",
    "In-hospital mortality",
    "ICU mortality",
    "In-hospital mortality, 24-hour landmark",
    "ICU mortality, 24-hour landmark"
  ),
  n = c(
    nrow(mimic_landmark_data),
    nrow(mimic_official_rds$analysis),
    nrow(eicu_rds$model_data),
    nrow(eicu_hospital_landmark_data),
    nrow(eicu_icu_landmark_data)
  ),
  events = c(
    sum(mimic_landmark_data$mortality_365d == 1),
    sum(mimic_official_rds$analysis$hospital_expire_flag == 1),
    sum(eicu_rds$model_data$icu_mortality == 1),
    sum(eicu_hospital_landmark_data$hospital_mortality == 1),
    sum(eicu_icu_landmark_data$icu_mortality == 1)
  ),
  effect_measure = c("HR", "OR", "OR", "OR", "OR"),
  P1_vs_P3 = c(
    effect_for(mimic_landmark, "P1 vs P3", "estimate"),
    effect_for(mimic_hospital, "P1 vs P3", "estimate"),
    effect_for(eicu_icu, "P1 vs P3", "OR"),
    effect_for(eicu_hospital_landmark, "P1 vs P3", "OR"),
    effect_for(eicu_icu_landmark, "P1 vs P3", "OR")
  ),
  P2_vs_P3 = c(
    effect_for(mimic_landmark, "P2 vs P3", "estimate"),
    effect_for(mimic_hospital, "P2 vs P3", "estimate"),
    effect_for(eicu_icu, "P2 vs P3", "OR"),
    effect_for(eicu_hospital_landmark, "P2 vs P3", "OR"),
    effect_for(eicu_icu_landmark, "P2 vs P3", "OR")
  )
)

absolute_risk_dir <- file.path(root, "output", "nhanes_adjusted_absolute_risk")
table4_by_phenotype <- readr::read_csv(
  file.path(absolute_risk_dir, "Table38A_adjusted_absolute_risk_and_RMST.csv"),
  show_col_types = FALSE
) |>
  dplyr::mutate(
    metric = dplyr::recode(
      metric,
      adjusted_risk_36_percent = "Adjusted 36-month mortality risk, %",
      adjusted_risk_60_percent = "Adjusted 60-month mortality risk, %",
      adjusted_rmst_60_months = "Adjusted RMST through 60 months, months"
    )
  ) |>
  dplyr::select(phenotype, metric, estimate, lower_95, upper_95, estimate_95ci)

table4_contrasts <- readr::read_csv(
  file.path(absolute_risk_dir, "Table38B_adjusted_risk_and_RMST_contrasts.csv"),
  show_col_types = FALSE
) |>
  dplyr::mutate(
    metric = dplyr::recode(
      metric,
      risk_difference_36_percentage_points = "36-month risk difference, percentage points",
      risk_ratio_36 = "36-month risk ratio",
      risk_difference_60_percentage_points = "60-month risk difference, percentage points",
      risk_ratio_60 = "60-month risk ratio",
      rmst_difference_60_months = "RMST difference through 60 months, months"
    ),
    estimate_95ci = sprintf("%.2f (%.2f-%.2f)", estimate, lower_95, upper_95)
  ) |>
  dplyr::select(comparison, metric, estimate, lower_95, upper_95, estimate_95ci, valid_replicates)

icu_selection_dir <- file.path(root, "output", "icu_selection_bias")
table5_selection_flow_legacy <- readr::read_csv(
  file.path(icu_selection_dir, "Table39A_ICU_selection_flow_and_weights.csv"),
  show_col_types = FALSE
) 
table5_selection_flow_mimic <- mimic_official_rds$selection_diagnostics |>
  dplyr::transmute(
    dataset = "MIMIC-IV",
    denominator_n,
    feature_complete_n = feature_complete_and_official_oasis_n,
    primary_model_n = feature_complete_and_official_oasis_n,
    primary_model_percent = selected_percent,
    trimmed_weight_99th_percentile = weight_99th_percentile,
    effective_sample_size
  )
table5_selection_flow <- dplyr::bind_rows(
  table5_selection_flow_mimic,
  table5_selection_flow_legacy |>
    dplyr::filter(dataset == "eICU") |>
    dplyr::select(
      dataset, denominator_n, feature_complete_n, primary_model_n,
      primary_model_percent, trimmed_weight_99th_percentile,
      effective_sample_size
    )
)
table5_selection_effects_eicu <- readr::read_csv(
  file.path(icu_selection_dir, "Table39E_selection_IPW_outcome_models.csv"),
  show_col_types = FALSE
) |>
  dplyr::filter(dataset == "eICU") |>
  dplyr::select(
    dataset, model, effect_measure, comparison, estimate, lower_95, upper_95,
    effect_95ci, p_value, n, events
  )
table5_selection_effects_mimic <- mimic_official_rds$outcome_models |>
  dplyr::filter(model %in% c(
    "Official OASIS quartile-stratified Cox",
    "Selection-IPW official OASIS quartile-stratified Cox"
  )) |>
  dplyr::mutate(dataset = "MIMIC-IV", .before = 1) |>
  dplyr::select(
    dataset, model, effect_measure, comparison, estimate, lower_95, upper_95,
    effect_95ci, p_value, n, events
  )
table5_selection_effects <- dplyr::bind_rows(
  table5_selection_effects_mimic,
  table5_selection_effects_eicu
)

albumin_dir <- file.path(root, "output", "harmonised_albumin_sensitivity")
table6_albumin_effects <- readr::read_csv(
  file.path(albumin_dir, "Table40C_harmonised_albumin_outcome_models.csv"),
  show_col_types = FALSE
) |>
  dplyr::select(
    dataset, model, effect_measure, comparison, estimate, lower_95, upper_95,
    effect_95ci, p_value, n, events
  )
table6_albumin_agreement <- readr::read_csv(
  file.path(albumin_dir, "Table40D_albumin_vs_current_agreement.csv"),
  show_col_types = FALSE
) |>
  dplyr::select(dataset, common_n, exact_label_agreement, adjusted_rand_index)
table6_albumin_selection <- readr::read_csv(
  file.path(albumin_dir, "Table40E_albumin_selection_weights.csv"),
  show_col_types = FALSE
) |>
  dplyr::select(
    dataset, denominator_n, albumin_feature_complete_n, primary_model_n,
    primary_model_percent, weight_99th_percentile, effective_sample_size
  )

sofa_dir <- file.path(root, "output", "mimic_official_sofa_v301")
table7_sofa_components <- readr::read_csv(
  file.path(sofa_dir, "Table42B_official_SOFA_component_availability.csv"),
  show_col_types = FALSE
)
table7_sofa_summary <- readr::read_csv(
  file.path(sofa_dir, "Table42C_official_SOFA_by_phenotype.csv"),
  show_col_types = FALSE
)
table7_sofa_comparison <- readr::read_csv(
  file.path(sofa_dir, "Table42D_SOFA_vs_OASIS_comparison.csv"),
  show_col_types = FALSE
)
table7_sofa_models <- readr::read_csv(
  file.path(sofa_dir, "Table42E_official_SOFA_outcome_models.csv"),
  show_col_types = FALSE
)

eicu_heterogeneity_dir <- file.path(root, "output", "eicu_hospital_heterogeneity")
table8_meta <- readr::read_csv(
  file.path(eicu_heterogeneity_dir, "Table43C_random_effects_meta_analysis.csv"),
  show_col_types = FALSE
)
table8_leave_one_out <- readr::read_csv(
  file.path(eicu_heterogeneity_dir, "Table43E_leave_one_out_summary.csv"),
  show_col_types = FALSE
)
table8_mixed_models <- readr::read_csv(
  file.path(eicu_heterogeneity_dir, "Table43F_mixed_model_effects.csv"),
  show_col_types = FALSE
)
table8_random_slope <- readr::read_csv(
  file.path(eicu_heterogeneity_dir, "Table43G_random_slope_diagnostics.csv"),
  show_col_types = FALSE
)
table8_direction <- readr::read_csv(
  file.path(eicu_heterogeneity_dir, "Table43H_hospital_effect_direction_summary.csv"),
  show_col_types = FALSE
)

effect_modification_dir <- file.path(root, "output", "age_sex_effect_modification")
table9_global_interactions <- readr::read_csv(
  file.path(effect_modification_dir, "Table44A_global_interaction_tests.csv"),
  show_col_types = FALSE
)
table9_interaction_ratios <- readr::read_csv(
  file.path(effect_modification_dir, "Table44B_P1_interaction_ratios.csv"),
  show_col_types = FALSE
)
table9_conditional_effects <- readr::read_csv(
  file.path(effect_modification_dir, "Table44C_P1_prespecified_conditional_effects.csv"),
  show_col_types = FALSE
)
table9_diagnostics <- readr::read_csv(
  file.path(effect_modification_dir, "Table44D_model_diagnostics.csv"),
  show_col_types = FALSE
)

supplement_index <- tibble::tribble(
  ~item, ~title, ~primary_source, ~role,
  "Table S1", "NHANES cohort flow and cycle-specific inclusion", "output/nhanes_2011_2018_rebuild/Table27C_cohort_flow.csv; Table27D_cycle_specific_flow.csv", "Required",
  "Table S2", "Robust preprocessing and feature-set sensitivity", "output/nhanes_robust_reanalysis/Table28E_Cox_models_all_variants.csv", "Required",
  "Table S3", "K-number diagnostics and bootstrap stability", "output/nhanes_robust_reanalysis/Table28A_cluster_number_diagnostics.csv", "Required",
  "Table S4", "Albumin substitution sensitivity", "output/nhanes_albumin_benchmarks/Table29A-Table29C", "Required",
  "Table S5", "Conventional immunonutritional benchmark models", "output/nhanes_albumin_benchmarks/Table29D-Table29F", "Required",
  "Table S6", "Proportional-hazards and time-varying score analyses", "output/nhanes_corrected_bootstrap/Table32A; Table32E-Table32F", "Required",
  "Table S7", "Complete-case selection and selection-IPW analysis", "output/nhanes_selection_bias/Table33A-Table33D", "Required",
  "Table S8", "Cluster-aware out-of-bag C-index bootstrap", "output/nhanes_bootstrap_cindex/Table34A-Table34C", "Required",
  "Table S9", "Leave-one-cycle-out transportability", "output/nhanes_leave_one_cycle_out/Table35A-Table35F", "Required",
  "Table S10", "MIMIC-IV first-24-hour profiles and official OASIS-adjusted models", "output/mimic_24h_robust_corrected_severity/Table30A-Table30F; output/mimic_official_oasis_v301/Table41B-Table41F", "Required",
  "Table S11", "eICU first-24-hour multicentre profiles and sensitivity models", "output/eicu_24h_robust/Table31A-Table31D", "Required",
  "Table S12", "Exploratory direct centroid projection and distribution shift", "output/validation_results", "Exploratory only",
  "Table S13", "Post-freeze adjusted absolute risk, survival curves, and RMST", "output/nhanes_adjusted_absolute_risk/Table38A-Table38F", "Clinical interpretation",
  "Table S14", "MIMIC-IV and eICU complete-case selection audit and selection-IPW models", "output/icu_selection_bias/Table39A-Table39F", "Required",
  "Table S15", "Harmonised albumin phenotype sensitivity across NHANES, MIMIC-IV, and eICU", "output/harmonised_albumin_sensitivity/Table40A-Table40E", "Post-freeze sensitivity",
  "Table S16", "MIT-LCP v3.0.1 official OASIS SQL provenance and comparison with OASIS-like score", "output/mimic_official_oasis_v301/Table41A-Table41F", "Required",
  "Table S17", "MIT-LCP v3.0.1 official first-day SOFA sensitivity and component availability", "output/mimic_official_sofa_v301/Table42A-Table42F", "Post-freeze sensitivity",
  "Table S18", "eICU hospital-specific effects, random-effects meta-analysis, and leave-one-hospital-out audit", "output/eicu_hospital_heterogeneity/Table43A-Table43H", "Post-freeze transportability audit",
  "Table S19", "Limited age and sex effect-modification analysis across three databases", "output/age_sex_effect_modification/Table44A-Table44D", "Post-freeze exploratory interaction audit"
)

write_markdown_table <- function(data, path) {
  display <- data
  display[] <- lapply(display, function(x) {
    x <- as.character(x)
    x[is.na(x)] <- ""
    gsub("\\|", "\\\\|", x)
  })
  header <- paste0("| ", paste(names(display), collapse = " | "), " |")
  separator <- paste0("| ", paste(rep("---", ncol(display)), collapse = " | "), " |")
  rows <- apply(display, 1, function(row) paste0("| ", paste(row, collapse = " | "), " |"))
  writeLines(c(header, separator, rows), path)
}

readr::write_csv(table1, file.path(output_dir, "Table1_corrected_weighted_baseline.csv"))
readr::write_csv(table2, file.path(output_dir, "Table2_corrected_NHANES_models.csv"))
readr::write_csv(table2_validation, file.path(output_dir, "Table2B_internal_validation.csv"))
readr::write_csv(table3_primary, file.path(output_dir, "Table3_cross_database_primary.csv"))
readr::write_csv(table3_sensitivity, file.path(output_dir, "Table3B_cross_database_sensitivity.csv"))
readr::write_csv(table4_by_phenotype, file.path(output_dir, "Table4A_adjusted_absolute_risk_RMST.csv"))
readr::write_csv(table4_contrasts, file.path(output_dir, "Table4B_adjusted_risk_RMST_contrasts.csv"))
readr::write_csv(table5_selection_flow, file.path(output_dir, "Table5A_ICU_selection_flow.csv"))
readr::write_csv(table5_selection_effects, file.path(output_dir, "Table5B_ICU_selection_IPW_models.csv"))
readr::write_csv(table6_albumin_effects, file.path(output_dir, "Table6A_harmonised_albumin_models.csv"))
readr::write_csv(table6_albumin_agreement, file.path(output_dir, "Table6B_harmonised_albumin_agreement.csv"))
readr::write_csv(table6_albumin_selection, file.path(output_dir, "Table6C_harmonised_albumin_selection.csv"))
readr::write_csv(table7_sofa_components, file.path(output_dir, "Table7A_official_SOFA_components.csv"))
readr::write_csv(table7_sofa_summary, file.path(output_dir, "Table7B_official_SOFA_by_phenotype.csv"))
readr::write_csv(table7_sofa_comparison, file.path(output_dir, "Table7C_SOFA_vs_OASIS.csv"))
readr::write_csv(table7_sofa_models, file.path(output_dir, "Table7D_official_SOFA_models.csv"))
readr::write_csv(table8_meta, file.path(output_dir, "Table8A_eICU_hospital_meta_analysis.csv"))
readr::write_csv(table8_leave_one_out, file.path(output_dir, "Table8B_eICU_leave_one_hospital_out.csv"))
readr::write_csv(table8_mixed_models, file.path(output_dir, "Table8C_eICU_hospital_mixed_models.csv"))
readr::write_csv(table8_random_slope, file.path(output_dir, "Table8D_eICU_random_slope_diagnostics.csv"))
readr::write_csv(table8_direction, file.path(output_dir, "Table8E_eICU_hospital_effect_directions.csv"))
readr::write_csv(table9_global_interactions, file.path(output_dir, "Table9A_age_sex_global_interactions.csv"))
readr::write_csv(table9_interaction_ratios, file.path(output_dir, "Table9B_P1_interaction_ratios.csv"))
readr::write_csv(table9_conditional_effects, file.path(output_dir, "Table9C_P1_conditional_effects.csv"))
readr::write_csv(table9_diagnostics, file.path(output_dir, "Table9D_effect_modification_diagnostics.csv"))
readr::write_csv(supplement_index, file.path(output_dir, "Supplementary_materials_index.csv"))

write_markdown_table(table1, file.path(output_dir, "Table1_corrected_weighted_baseline.md"))
write_markdown_table(table2, file.path(output_dir, "Table2_corrected_NHANES_models.md"))
write_markdown_table(table2_validation, file.path(output_dir, "Table2B_internal_validation.md"))
write_markdown_table(table3_primary, file.path(output_dir, "Table3_cross_database_primary.md"))
write_markdown_table(table3_sensitivity, file.path(output_dir, "Table3B_cross_database_sensitivity.md"))
write_markdown_table(table4_by_phenotype, file.path(output_dir, "Table4A_adjusted_absolute_risk_RMST.md"))
write_markdown_table(table4_contrasts, file.path(output_dir, "Table4B_adjusted_risk_RMST_contrasts.md"))
write_markdown_table(table5_selection_flow, file.path(output_dir, "Table5A_ICU_selection_flow.md"))
write_markdown_table(table5_selection_effects, file.path(output_dir, "Table5B_ICU_selection_IPW_models.md"))
write_markdown_table(table6_albumin_effects, file.path(output_dir, "Table6A_harmonised_albumin_models.md"))
write_markdown_table(table6_albumin_agreement, file.path(output_dir, "Table6B_harmonised_albumin_agreement.md"))
write_markdown_table(table6_albumin_selection, file.path(output_dir, "Table6C_harmonised_albumin_selection.md"))
write_markdown_table(table7_sofa_components, file.path(output_dir, "Table7A_official_SOFA_components.md"))
write_markdown_table(table7_sofa_summary, file.path(output_dir, "Table7B_official_SOFA_by_phenotype.md"))
write_markdown_table(table7_sofa_comparison, file.path(output_dir, "Table7C_SOFA_vs_OASIS.md"))
write_markdown_table(table7_sofa_models, file.path(output_dir, "Table7D_official_SOFA_models.md"))
write_markdown_table(table8_meta, file.path(output_dir, "Table8A_eICU_hospital_meta_analysis.md"))
write_markdown_table(table8_leave_one_out, file.path(output_dir, "Table8B_eICU_leave_one_hospital_out.md"))
write_markdown_table(table8_mixed_models, file.path(output_dir, "Table8C_eICU_hospital_mixed_models.md"))
write_markdown_table(table8_random_slope, file.path(output_dir, "Table8D_eICU_random_slope_diagnostics.md"))
write_markdown_table(table8_direction, file.path(output_dir, "Table8E_eICU_hospital_effect_directions.md"))
write_markdown_table(table9_global_interactions, file.path(output_dir, "Table9A_age_sex_global_interactions.md"))
write_markdown_table(table9_interaction_ratios, file.path(output_dir, "Table9B_P1_interaction_ratios.md"))
write_markdown_table(table9_conditional_effects, file.path(output_dir, "Table9C_P1_conditional_effects.md"))
write_markdown_table(table9_diagnostics, file.path(output_dir, "Table9D_effect_modification_diagnostics.md"))
write_markdown_table(supplement_index, file.path(output_dir, "Supplementary_materials_index.md"))

notes <- c(
  "Corrected publication table notes",
  "",
  "Table 1 reports unweighted participant counts and NHANES complex-survey weighted summaries.",
  "No P values are shown for clustering variables because the groups were defined using those variables.",
  "Table 2 designates the complex-survey fully adjusted Cox model as the primary NHANES inference.",
  "Table 3 uses the MIT-LCP MIMIC Code v3.0.1 official OASIS implementation; it does not pool HRs and ORs or imply identical phenotype membership.",
  "Table 4 is a post-freeze clinical interpretation analysis and does not replace the primary survey-weighted HR.",
  "Table 5 uses official OASIS for MIMIC selection adjustment; eICU P1 is directionally consistent but attenuated after selection IPW.",
  "Table 6 is a post-freeze harmonised albumin sensitivity; its de novo labels and effects do not replace the total-protein primary phenotype.",
  "Table 7 is a post-freeze official first-day SOFA sensitivity; respiratory SOFA availability is reported because missing components are scored as zero by the official concept.",
  "Table 8 is a post-freeze eICU hospital transportability audit; the hospital-specific random-effects interval crosses one and must not be described as uniformly significant replication.",
  "Table 9 is a post-freeze limited effect-modification audit; none of six global age/sex interactions remains significant after BH correction, and conditional estimates must not be described as subgroup differences.",
  "P1 is the primary comparison; P2 is retained as a secondary comparison.",
  "All values supersede tables based on the old 3,386-person NHANES cohort or 72-hour ICU windows."
)
writeLines(notes, file.path(output_dir, "PUBLICATION_TABLE_NOTES.txt"))

cat("Corrected publication tables generated in: ", output_dir, "\n", sep = "")
cat(paste(notes, collapse = "\n"), "\n")
