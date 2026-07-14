# ==============================================================================
# Automated QA for corrected publication tables
# ==============================================================================

required_packages <- c("dplyr", "readr", "tibble")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop("Missing required package(s): ", paste(missing_packages, collapse = ", "), call. = FALSE)
}

root <- normalizePath(
  Sys.getenv("PROJECT_ROOT", unset = getwd()), winslash = "/", mustWork = TRUE
)
table_dir <- file.path(root, "output", "corrected_publication_tables")

table1 <- readr::read_csv(file.path(table_dir, "Table1_corrected_weighted_baseline.csv"), show_col_types = FALSE)
table2 <- readr::read_csv(file.path(table_dir, "Table2_corrected_NHANES_models.csv"), show_col_types = FALSE)
table3 <- readr::read_csv(file.path(table_dir, "Table3_cross_database_primary.csv"), show_col_types = FALSE)
table4 <- readr::read_csv(file.path(table_dir, "Table4B_adjusted_risk_RMST_contrasts.csv"), show_col_types = FALSE)
table5_flow <- readr::read_csv(file.path(table_dir, "Table5A_ICU_selection_flow.csv"), show_col_types = FALSE)
table5_effects <- readr::read_csv(file.path(table_dir, "Table5B_ICU_selection_IPW_models.csv"), show_col_types = FALSE)
table6_effects <- readr::read_csv(file.path(table_dir, "Table6A_harmonised_albumin_models.csv"), show_col_types = FALSE)
table6_agreement <- readr::read_csv(file.path(table_dir, "Table6B_harmonised_albumin_agreement.csv"), show_col_types = FALSE)
table6_selection <- readr::read_csv(file.path(table_dir, "Table6C_harmonised_albumin_selection.csv"), show_col_types = FALSE)
table7_sofa_components <- readr::read_csv(file.path(table_dir, "Table7A_official_SOFA_components.csv"), show_col_types = FALSE)
table7_sofa_comparison <- readr::read_csv(file.path(table_dir, "Table7C_SOFA_vs_OASIS.csv"), show_col_types = FALSE)
table7_sofa_models <- readr::read_csv(file.path(table_dir, "Table7D_official_SOFA_models.csv"), show_col_types = FALSE)
table8_meta <- readr::read_csv(file.path(table_dir, "Table8A_eICU_hospital_meta_analysis.csv"), show_col_types = FALSE)
table8_leave_one_out <- readr::read_csv(file.path(table_dir, "Table8B_eICU_leave_one_hospital_out.csv"), show_col_types = FALSE)
table8_mixed_models <- readr::read_csv(file.path(table_dir, "Table8C_eICU_hospital_mixed_models.csv"), show_col_types = FALSE)
table8_random_slope <- readr::read_csv(file.path(table_dir, "Table8D_eICU_random_slope_diagnostics.csv"), show_col_types = FALSE)
table8_direction <- readr::read_csv(file.path(table_dir, "Table8E_eICU_hospital_effect_directions.csv"), show_col_types = FALSE)
table9_global <- readr::read_csv(file.path(table_dir, "Table9A_age_sex_global_interactions.csv"), show_col_types = FALSE)
table9_ratios <- readr::read_csv(file.path(table_dir, "Table9B_P1_interaction_ratios.csv"), show_col_types = FALSE)
table9_conditional <- readr::read_csv(file.path(table_dir, "Table9C_P1_conditional_effects.csv"), show_col_types = FALSE)
table9_diagnostics <- readr::read_csv(file.path(table_dir, "Table9D_effect_modification_diagnostics.csv"), show_col_types = FALSE)
official_oasis_dir <- file.path(root, "output", "mimic_official_oasis_v301")
official_oasis_provenance <- readr::read_csv(
  file.path(official_oasis_dir, "Table41A_official_SQL_provenance.csv"), show_col_types = FALSE
)
official_oasis_comparison <- readr::read_csv(
  file.path(official_oasis_dir, "Table41C_official_vs_OASISlike_comparison.csv"), show_col_types = FALSE
)
official_oasis_ph <- readr::read_csv(
  file.path(official_oasis_dir, "Table41E_official_OASIS_PH_checks.csv"), show_col_types = FALSE
)
official_sofa_dir <- file.path(root, "output", "mimic_official_sofa_v301")
official_sofa_provenance <- readr::read_csv(
  file.path(official_sofa_dir, "Table42A_official_SOFA_SQL_provenance.csv"), show_col_types = FALSE
)
official_sofa_ph <- readr::read_csv(
  file.path(official_sofa_dir, "Table42F_official_SOFA_PH_checks.csv"), show_col_types = FALSE
)

checks <- list()
add_check <- function(name, passed, detail) {
  checks[[length(checks) + 1L]] <<- tibble::tibble(
    check = name,
    passed = isTRUE(passed),
    detail = as.character(detail)
  )
}

participant_row <- table1 |>
  dplyr::filter(Characteristic == "Participants, unweighted n")
participant_counts <- as.numeric(participant_row[c("P1", "P2", "P3")])
add_check(
  "Table 1 phenotype counts sum to overall n",
  sum(participant_counts) == as.numeric(participant_row$Overall),
  paste(participant_counts, collapse = " + ")
)
add_check(
  "Table 1 overall n equals corrected cohort",
  as.numeric(participant_row$Overall) == 4636,
  participant_row$Overall
)

population_row <- table1 |>
  dplyr::filter(Characteristic == "Weighted population proportion")
population_percent <- as.numeric(sub("%", "", unlist(population_row[c("P1", "P2", "P3")])))
add_check(
  "Weighted phenotype proportions sum to 100%",
  abs(sum(population_percent) - 100) <= 0.2,
  paste0(sum(population_percent), "%")
)
add_check(
  "Weighted phenotype proportions are valid",
  all(population_percent > 0 & population_percent < 100),
  paste(population_percent, collapse = ", ")
)
add_check(
  "Table 1 contains no P-value column",
  !any(grepl("p_value|P value", names(table1), ignore.case = TRUE)),
  paste(names(table1), collapse = ", ")
)

primary_rows <- table2 |>
  dplyr::filter(model_role == "Primary")
primary_p1 <- primary_rows |>
  dplyr::filter(comparison == "P1 vs P3")
primary_p2 <- primary_rows |>
  dplyr::filter(comparison == "P2 vs P3")
add_check(
  "Table 2 has one primary row per phenotype comparison",
  nrow(primary_p1) == 1L && nrow(primary_p2) == 1L,
  paste0("P1 rows=", nrow(primary_p1), "; P2 rows=", nrow(primary_p2))
)
add_check(
  "Primary NHANES P1 estimate matches frozen result",
  abs(primary_p1$HR - 1.6830536847349) < 1e-10,
  primary_p1$effect_95ci
)
add_check(
  "Primary NHANES model uses corrected common cohort",
  primary_p1$n == 3979 && primary_p1$events == 720,
  paste0("n=", primary_p1$n, "; events=", primary_p1$events)
)

full_p1 <- table2 |>
  dplyr::filter(model == "Fully adjusted Cox", comparison == "P1 vs P3")
add_check(
  "Fully adjusted Cox matches cluster-aware bootstrap source model",
  abs(full_p1$HR - 1.6996655682007529) < 1e-10,
  full_p1$effect_95ci
)

add_check(
  "Table 3 contains exactly the three frozen databases",
  identical(table3$database, c("NHANES 2011-2018", "MIMIC-IV v3.1", "eICU v2.0")),
  paste(table3$database, collapse = "; ")
)
add_check(
  "Table 3 event counts are positive and below cohort sizes",
  all(table3$events > 0 & table3$events < table3$n),
  paste0(table3$database, ": ", table3$events, "/", table3$n, collapse = "; ")
)
add_check(
  "Table 3 preserves HR/OR distinction",
  identical(table3$effect_measure, c("HR", "HR", "OR")),
  paste(table3$effect_measure, collapse = ", ")
)
add_check(
  "MIMIC primary effect matches official OASIS first-24-hour source",
  table3$P1_vs_P3[table3$database == "MIMIC-IV v3.1"] == "1.67 (1.36-2.05)",
  table3$P1_vs_P3[table3$database == "MIMIC-IV v3.1"]
)
add_check(
  "MIMIC primary row uses official OASIS and complete strict cohort",
  table3$n[table3$database == "MIMIC-IV v3.1"] == 1145 &&
    table3$events[table3$database == "MIMIC-IV v3.1"] == 574 &&
    grepl("official OASIS", table3$model[table3$database == "MIMIC-IV v3.1"]),
  paste0(
    "n=", table3$n[table3$database == "MIMIC-IV v3.1"],
    "; events=", table3$events[table3$database == "MIMIC-IV v3.1"]
  )
)
add_check(
  "eICU primary effect matches APACHE mixed model",
  table3$P1_vs_P3[table3$database == "eICU v2.0"] == "1.28 (1.11-1.48)",
  table3$P1_vs_P3[table3$database == "eICU v2.0"]
)

p1_risk_60 <- table4 |>
  dplyr::filter(
    comparison == "P1 vs P3",
    metric == "60-month risk difference, percentage points"
  )
p1_rmst_60 <- table4 |>
  dplyr::filter(
    comparison == "P1 vs P3",
    metric == "RMST difference through 60 months, months"
  )
p2_risk_60 <- table4 |>
  dplyr::filter(
    comparison == "P2 vs P3",
    metric == "60-month risk difference, percentage points"
  )
add_check(
  "P1 adjusted 60-month absolute risk difference is positive with CI above zero",
  nrow(p1_risk_60) == 1L && p1_risk_60$estimate > 0 && p1_risk_60$lower_95 > 0,
  p1_risk_60$estimate_95ci
)
add_check(
  "P1 adjusted 60-month RMST difference is negative with CI below zero",
  nrow(p1_rmst_60) == 1L && p1_rmst_60$estimate < 0 && p1_rmst_60$upper_95 < 0,
  p1_rmst_60$estimate_95ci
)
add_check(
  "P2 adjusted 60-month risk difference interval includes zero",
  nrow(p2_risk_60) == 1L && p2_risk_60$lower_95 <= 0 && p2_risk_60$upper_95 >= 0,
  p2_risk_60$estimate_95ci
)

mimic_flow <- table5_flow |>
  dplyr::filter(dataset == "MIMIC-IV")
eicu_flow <- table5_flow |>
  dplyr::filter(dataset == "eICU")
add_check(
  "MIMIC selection audit reproduces denominator and primary model size",
  nrow(mimic_flow) == 1L && mimic_flow$denominator_n == 33671 && mimic_flow$primary_model_n == 1145,
  paste0(mimic_flow$primary_model_n, "/", mimic_flow$denominator_n)
)
add_check(
  "eICU selection audit reproduces denominator and primary model size",
  nrow(eicu_flow) == 1L && eicu_flow$denominator_n == 68798 && eicu_flow$primary_model_n == 12548,
  paste0(eicu_flow$primary_model_n, "/", eicu_flow$denominator_n)
)

mimic_ipw_p1 <- table5_effects |>
  dplyr::filter(
    dataset == "MIMIC-IV",
    model == "Selection-IPW official OASIS quartile-stratified Cox",
    comparison == "P1 vs P3"
  )
eicu_ipw_p1 <- table5_effects |>
  dplyr::filter(
    dataset == "eICU",
    model == "Selection-IPW APACHE-adjusted mixed model",
    comparison == "P1 vs P3"
  )
add_check(
  "MIMIC selection-IPW P1 interval remains above one",
  nrow(mimic_ipw_p1) == 1L && mimic_ipw_p1$lower_95 > 1 &&
    mimic_ipw_p1$effect_95ci == "1.71 (1.29-2.25)",
  mimic_ipw_p1$effect_95ci
)
add_check(
  "eICU selection-IPW P1 interval includes one",
  nrow(eicu_ipw_p1) == 1L && eicu_ipw_p1$lower_95 <= 1 && eicu_ipw_p1$upper_95 >= 1,
  eicu_ipw_p1$effect_95ci
)

albumin_p1 <- table6_effects |>
  dplyr::filter(comparison == "P1 vs P3")
albumin_nhanes_p1 <- albumin_p1 |>
  dplyr::filter(dataset == "NHANES")
albumin_mimic_ipw_p1 <- albumin_p1 |>
  dplyr::filter(
    dataset == "MIMIC-IV",
    model == "Selection-IPW official OASIS albumin phenotype Cox"
  )
albumin_eicu_ipw_p1 <- albumin_p1 |>
  dplyr::filter(
    dataset == "eICU",
    model == "Selection-IPW APACHE-adjusted albumin phenotype mixed model"
  )
add_check(
  "Harmonised albumin NHANES P1 interval is above one",
  nrow(albumin_nhanes_p1) == 1L && albumin_nhanes_p1$lower_95 > 1 && albumin_nhanes_p1$n == 3979,
  albumin_nhanes_p1$effect_95ci
)
add_check(
  "Harmonised albumin MIMIC selection-IPW P1 interval is above one",
  nrow(albumin_mimic_ipw_p1) == 1L && albumin_mimic_ipw_p1$lower_95 > 1 && albumin_mimic_ipw_p1$n == 1100,
  albumin_mimic_ipw_p1$effect_95ci
)
add_check(
  "Harmonised albumin eICU selection-IPW P1 interval is above one",
  nrow(albumin_eicu_ipw_p1) == 1L && albumin_eicu_ipw_p1$lower_95 > 1 && albumin_eicu_ipw_p1$n == 13234,
  albumin_eicu_ipw_p1$effect_95ci
)
add_check(
  "Harmonised albumin ICU labels retain moderate agreement with current variants",
  all(table6_agreement$adjusted_rand_index[table6_agreement$dataset %in% c("MIMIC-IV", "eICU")] > 0.6),
  paste0(table6_agreement$dataset, " ARI=", round(table6_agreement$adjusted_rand_index, 3), collapse = "; ")
)
add_check(
  "Harmonised albumin selection denominators match audited ICU cohorts",
  identical(table6_selection$denominator_n, c(33671, 68798)),
  paste0(table6_selection$dataset, ": ", table6_selection$denominator_n, collapse = "; ")
)

official_stratified_global <- official_oasis_ph |>
  dplyr::filter(
    model == "Official OASIS quartile-stratified Cox",
    term == "GLOBAL"
  )
add_check(
  "Official OASIS provenance contains the frozen MIT-LCP v3.0.1 dependency chain",
  nrow(official_oasis_provenance) == 11L &&
    all(official_oasis_provenance$mimic_code_release == "v3.0.1") &&
    all(nchar(official_oasis_provenance$md5) == 32L),
  paste0(nrow(official_oasis_provenance), " scripts")
)
add_check(
  "Official OASIS covers the complete strict MIMIC cohort",
  official_oasis_comparison$official_n == 1145 &&
    official_oasis_comparison$oasis_like_complete_n == 1085,
  paste0(
    "official=", official_oasis_comparison$official_n,
    "; OASIS-like=", official_oasis_comparison$oasis_like_complete_n
  )
)
add_check(
  "Official and OASIS-like scores retain strong rank correlation",
  official_oasis_comparison$spearman_correlation > 0.8,
  round(official_oasis_comparison$spearman_correlation, 3)
)
add_check(
  "Official OASIS quartile-stratified model passes the global PH test",
  nrow(official_stratified_global) == 1L && official_stratified_global$p > 0.05,
  official_stratified_global$p
)

sofa_stratified_p1 <- table7_sofa_models |>
  dplyr::filter(
    model == "Official first-day SOFA quartile-stratified Cox",
    comparison == "P1 vs P3"
  )
sofa_ipw_p1 <- table7_sofa_models |>
  dplyr::filter(
    model == "Selection-IPW official first-day SOFA quartile-stratified Cox",
    comparison == "P1 vs P3"
  )
sofa_stratified_global <- official_sofa_ph |>
  dplyr::filter(
    model == "Official first-day SOFA quartile-stratified Cox",
    term == "GLOBAL"
  )
sofa_stratified_phenotype <- official_sofa_ph |>
  dplyr::filter(
    model == "Official first-day SOFA quartile-stratified Cox",
    term == "phenotype"
  )
sofa_respiration <- table7_sofa_components |>
  dplyr::filter(component == "respiration")
sofa_total <- table7_sofa_components |>
  dplyr::filter(component == "official first-day SOFA")
add_check(
  "Official first-day SOFA P1 sensitivity remains above one",
  nrow(sofa_stratified_p1) == 1L &&
    sofa_stratified_p1$effect_95ci == "1.56 (1.27-1.92)" &&
    sofa_stratified_p1$n == 1145,
  sofa_stratified_p1$effect_95ci
)
add_check(
  "Official first-day SOFA selection-IPW P1 sensitivity remains above one",
  nrow(sofa_ipw_p1) == 1L &&
    sofa_ipw_p1$effect_95ci == "1.65 (1.25-2.17)",
  sofa_ipw_p1$effect_95ci
)
add_check(
  "Official first-day SOFA provenance contains 21 MIT-LCP v3.0.1 scripts",
  nrow(official_sofa_provenance) == 21L &&
    all(official_sofa_provenance$mimic_code_release == "v3.0.1") &&
    all(nchar(official_sofa_provenance$md5) == 32L),
  paste0(nrow(official_sofa_provenance), " scripts")
)
add_check(
  "Official SOFA and OASIS provide related but nonredundant severity information",
  table7_sofa_comparison$spearman_sofa_oasis > 0.3 &&
    table7_sofa_comparison$spearman_sofa_oasis < 0.8,
  round(table7_sofa_comparison$spearman_sofa_oasis, 3)
)
add_check(
  "Official SOFA quartile-stratified phenotype effect passes PH diagnostic",
  nrow(sofa_stratified_global) == 1L && sofa_stratified_global$p > 0.05 &&
    nrow(sofa_stratified_phenotype) == 1L && sofa_stratified_phenotype$p > 0.05,
  paste0(
    "global P=", round(sofa_stratified_global$p, 3),
    "; phenotype P=", round(sofa_stratified_phenotype$p, 3)
  )
)
add_check(
  "SOFA respiratory missingness remains explicitly visible while total score covers all stays",
  nrow(sofa_respiration) == 1L && sofa_respiration$available_percent < 50 &&
    nrow(sofa_total) == 1L && sofa_total$available_percent == 100,
  paste0(
    "respiration=", round(sofa_respiration$available_percent, 1),
    "%; total=", round(sofa_total$available_percent, 1), "%"
  )
)

eicu_all_hospitals_mixed <- table8_mixed_models |>
  dplyr::filter(model == "All hospitals, APACHE-adjusted random intercept")
eicu_high_information_mixed <- table8_mixed_models |>
  dplyr::filter(model == "High-information hospitals, APACHE-adjusted random intercept")
add_check(
  "eICU high-information hospital meta-analysis retains direction but interval crosses one",
  table8_meta$hospitals == 17 && table8_meta$pooled_or > 1 &&
    table8_meta$lower_95 < 1 && table8_meta$upper_95 > 1,
  table8_meta$effect_95ci
)
add_check(
  "eICU hospital heterogeneity is low-to-moderate with a wide prediction interval",
  table8_meta$i_squared_percent >= 0 && table8_meta$i_squared_percent < 50 &&
    table8_meta$q_p_value > 0.05 && table8_meta$prediction_lower_95 < 1 &&
    table8_meta$prediction_upper_95 > 1,
  paste0(
    "I2=", round(table8_meta$i_squared_percent, 1),
    "%; PI=", table8_meta$prediction_interval
  )
)
add_check(
  "Most high-information hospitals have P1 odds ratios above one",
  table8_direction$estimable_hospitals == 17 &&
    table8_direction$hospitals_or_above_one == 13 &&
    table8_direction$hospitals_or_below_one == 4,
  paste0(
    table8_direction$hospitals_or_above_one, "/",
    table8_direction$estimable_hospitals, " above one"
  )
)
add_check(
  "Leave-one-hospital-out estimates do not reverse direction but usually cross one",
  table8_leave_one_out$omitted_models == 17 &&
    table8_leave_one_out$pooled_or_min > 1 &&
    table8_leave_one_out$models_with_ci_crossing_one == 16 &&
    table8_leave_one_out$models_with_ci_above_one == 1,
  paste0(
    "OR range ", round(table8_leave_one_out$pooled_or_min, 2), "-",
    round(table8_leave_one_out$pooled_or_max, 2)
  )
)
add_check(
  "All-hospital mixed model reproduces the frozen eICU P1 estimate",
  nrow(eicu_all_hospitals_mixed) == 1L &&
    eicu_all_hospitals_mixed$effect_95ci == "1.28 (1.11-1.48)",
  eicu_all_hospitals_mixed$effect_95ci
)
add_check(
  "High-information hospital mixed-model interval includes one",
  nrow(eicu_high_information_mixed) == 1L &&
    eicu_high_information_mixed$lower_95 <= 1 && eicu_high_information_mixed$upper_95 >= 1,
  eicu_high_information_mixed$effect_95ci
)
add_check(
  "Exploratory eICU P1 random-slope model is singular and not overinterpreted",
  isTRUE(table8_random_slope$random_slope_singular) &&
    table8_random_slope$likelihood_ratio_p > 0.9,
  paste0("LRT P=", round(table8_random_slope$likelihood_ratio_p, 3))
)

add_check(
  "Age and sex audit contains exactly six global interaction tests",
  nrow(table9_global) == 6L &&
    setequal(table9_global$dataset, c("NHANES 2011-2018", "MIMIC-IV v3.1", "eICU v2.0")) &&
    setequal(table9_global$modifier, c("Age", "Sex")),
  paste0(nrow(table9_global), " tests")
)
add_check(
  "No global age or sex interaction survives BH correction",
  all(!table9_global$significant_bh_0_05) && all(table9_global$p_value_bh > 0.05),
  paste0("minimum BH P=", round(min(table9_global$p_value_bh), 3))
)
mimic_age_interaction <- table9_global |>
  dplyr::filter(dataset == "MIMIC-IV v3.1", modifier == "Age")
add_check(
  "MIMIC age interaction is transparently attenuated by multiplicity correction",
  nrow(mimic_age_interaction) == 1L && mimic_age_interaction$p_value < 0.05 &&
    mimic_age_interaction$p_value_bh > 0.05,
  paste0(
    "raw P=", round(mimic_age_interaction$p_value, 3),
    "; BH P=", round(mimic_age_interaction$p_value_bh, 3)
  )
)
add_check(
  "All P1 interaction-ratio intervals include one",
  nrow(table9_ratios) == 6L && all(table9_ratios$lower_95 <= 1 & table9_ratios$upper_95 >= 1),
  paste(table9_ratios$effect_95ci, collapse = "; ")
)
add_check(
  "Pre-specified conditional effects have valid ordered confidence intervals",
  nrow(table9_conditional) == 12L &&
    all(table9_conditional$lower_95 <= table9_conditional$estimate) &&
    all(table9_conditional$estimate <= table9_conditional$upper_95),
  paste0(nrow(table9_conditional), " conditional estimates")
)
add_check(
  "Pre-specified P1 conditional effects retain direction across age and sex",
  all(table9_conditional$estimate > 1) && all(table9_conditional$lower_95 > 1),
  paste(table9_conditional$effect_95ci, collapse = "; ")
)
add_check(
  "eICU interaction mixed models converge without singularity",
  nrow(table9_diagnostics) == 2L && all(table9_diagnostics$converged) &&
    all(!table9_diagnostics$singular),
  paste0(table9_diagnostics$modifier, ": converged=", table9_diagnostics$converged,
         ", singular=", table9_diagnostics$singular, collapse = "; ")
)

required_files <- c(
  "Table1_corrected_weighted_baseline.csv",
  "Table2_corrected_NHANES_models.csv",
  "Table2B_internal_validation.csv",
  "Table3_cross_database_primary.csv",
  "Table3B_cross_database_sensitivity.csv",
  "Table4A_adjusted_absolute_risk_RMST.csv",
  "Table4B_adjusted_risk_RMST_contrasts.csv",
  "Table5A_ICU_selection_flow.csv",
  "Table5B_ICU_selection_IPW_models.csv",
  "Table6A_harmonised_albumin_models.csv",
  "Table6B_harmonised_albumin_agreement.csv",
  "Table6C_harmonised_albumin_selection.csv",
  "Table7A_official_SOFA_components.csv",
  "Table7B_official_SOFA_by_phenotype.csv",
  "Table7C_SOFA_vs_OASIS.csv",
  "Table7D_official_SOFA_models.csv",
  "Table8A_eICU_hospital_meta_analysis.csv",
  "Table8B_eICU_leave_one_hospital_out.csv",
  "Table8C_eICU_hospital_mixed_models.csv",
  "Table8D_eICU_random_slope_diagnostics.csv",
  "Table8E_eICU_hospital_effect_directions.csv",
  "Table9A_age_sex_global_interactions.csv",
  "Table9B_P1_interaction_ratios.csv",
  "Table9C_P1_conditional_effects.csv",
  "Table9D_effect_modification_diagnostics.csv",
  "Supplementary_materials_index.csv",
  "Table1_corrected_weighted_baseline.md",
  "Table2_corrected_NHANES_models.md",
  "Table3_cross_database_primary.md",
  "Table4A_adjusted_absolute_risk_RMST.md",
  "Table4B_adjusted_risk_RMST_contrasts.md"
  ,"Table5A_ICU_selection_flow.md"
  ,"Table5B_ICU_selection_IPW_models.md"
  ,"Table6A_harmonised_albumin_models.md"
  ,"Table6B_harmonised_albumin_agreement.md"
  ,"Table6C_harmonised_albumin_selection.md"
  ,"Table7A_official_SOFA_components.md"
  ,"Table7B_official_SOFA_by_phenotype.md"
  ,"Table7C_SOFA_vs_OASIS.md"
  ,"Table7D_official_SOFA_models.md"
  ,"Table8A_eICU_hospital_meta_analysis.md"
  ,"Table8B_eICU_leave_one_hospital_out.md"
  ,"Table8C_eICU_hospital_mixed_models.md"
  ,"Table8D_eICU_random_slope_diagnostics.md"
  ,"Table8E_eICU_hospital_effect_directions.md"
  ,"Table9A_age_sex_global_interactions.md"
  ,"Table9B_P1_interaction_ratios.md"
  ,"Table9C_P1_conditional_effects.md"
  ,"Table9D_effect_modification_diagnostics.md"
)
file_status <- file.exists(file.path(table_dir, required_files))
add_check(
  "All frozen publication table files exist",
  all(file_status),
  paste(required_files[!file_status], collapse = "; ")
)

markdown_files <- list.files(table_dir, pattern = "\\.md$", full.names = TRUE)
markdown_text <- unlist(lapply(markdown_files, readLines, warn = FALSE, encoding = "UTF-8"))
add_check(
  "Markdown outputs contain no replacement-character encoding errors",
  !any(grepl("\uFFFD|脳", markdown_text)),
  paste(basename(markdown_files), collapse = "; ")
)

qa <- dplyr::bind_rows(checks)
readr::write_csv(qa, file.path(table_dir, "TABLE_QA_REPORT.csv"))

summary_lines <- c(
  "Corrected publication table QA",
  "",
  paste0("Checks passed: ", sum(qa$passed), "/", nrow(qa)),
  "",
  paste(capture.output(print(qa)), collapse = "\n")
)
writeLines(summary_lines, file.path(table_dir, "TABLE_QA_REPORT.txt"))
cat(paste(summary_lines, collapse = "\n"), "\n")

if (!all(qa$passed)) {
  stop("Publication table QA failed. Review TABLE_QA_REPORT.csv.", call. = FALSE)
}
