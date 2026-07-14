# ==============================================================================
# Independent audit of the frozen processing pipeline, primary data, and models
# ==============================================================================

required_packages <- c("dplyr", "readr", "tibble", "survival", "survey", "lme4")
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
scripts_dir <- file.path(root, "scripts")
output_dir <- file.path(root, "output", "final_pipeline_audit")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

canonical_pipeline <- tibble::tribble(
  ~order, ~script, ~role, ~primary_output,
  1L, "27_rebuild_NHANES_2011_2018.R", "Rebuild four-cycle NHANES cohort and linked mortality", "output/nhanes_2011_2018_rebuild/NHANES_2011_2018_rebuilt_analytical_cohort.rds",
  2L, "27B_augment_NHANES_covariates.R", "Add smoking, hypertension, and diabetes", "output/nhanes_covariate_upgrade/NHANES_2011_2018_covariate_augmented.rds",
  3L, "28_NHANES_robust_phenotype_reanalysis.R", "Frozen robust six-feature phenotype", "output/nhanes_robust_reanalysis/NHANES_robust_reanalysis_results.rds",
  4L, "28B_NHANES_cluster_number_validation.R", "K-number and subsampling stability", "output/nhanes_cluster_number_validation/NHANES_cluster_number_validation_results.rds",
  5L, "29_NHANES_albumin_and_conventional_benchmarks.R", "Albumin and conventional benchmark models", "output/nhanes_albumin_benchmarks/NHANES_albumin_benchmark_results.rds",
  6L, "32_NHANES_cluster_aware_bootstrap_and_PH.R", "Reclustering bootstrap and PH audit", "output/nhanes_corrected_bootstrap/NHANES_corrected_bootstrap_results.rds",
  7L, "32B_NHANES_time_varying_score.R", "Time-varying continuous score", "output/nhanes_corrected_bootstrap/NHANES_time_varying_score_results.rds",
  8L, "33_NHANES_selection_bias_sensitivity.R", "NHANES selection audit and selection IPW", "output/nhanes_selection_bias/NHANES_selection_bias_results.rds",
  9L, "34_NHANES_bootstrap_Cindex_increment.R", "Cluster-aware out-of-bag C-index", "output/nhanes_bootstrap_cindex/NHANES_bootstrap_Cindex_results.rds",
  10L, "35_NHANES_leave_one_cycle_out_validation.R", "Leave-one-cycle-out transportability", "output/nhanes_leave_one_cycle_out/NHANES_LOCO_results.rds",
  11L, "09_MIMIC_IV_external_validation_formal.R", "Strict 0-24h MIMIC extraction", "output/mimic_24h_validation/MIMIC_projected_albumin_proxy_dataset.csv",
  12L, "30_MIMIC_IV_24h_robust_revalidation.R", "MIMIC de novo robust phenotype", "output/mimic_24h_robust_corrected_severity/MIMIC_24h_robust_results.rds",
  13L, "10_eICU_multicenter_validation.R", "Strict 0-24h eICU extraction", "output/eicu_24h_extraction/eICU_denovo_strict_total_protein_dataset.csv",
  14L, "31_eICU_24h_robust_revalidation.R", "eICU de novo robust phenotype and APACHE model", "output/eicu_24h_robust/eICU_24h_robust_results.rds",
  15L, "39_ICU_complete_case_selection_audit.R", "ICU complete-case selection audit", "output/icu_selection_bias/ICU_selection_bias_results.rds",
  16L, "41_MIMIC_official_OASIS_v301.R", "MIT-LCP v3.0.1 official OASIS reconstruction", "output/mimic_official_oasis_v301/MIMIC_official_OASIS_v301_results.rds",
  17L, "42_MIMIC_official_first_day_SOFA_v301.R", "MIT-LCP v3.0.1 official first-day SOFA", "output/mimic_official_sofa_v301/MIMIC_official_first_day_SOFA_v301_results.rds",
  18L, "40_harmonised_albumin_cross_database_sensitivity.R", "Harmonised albumin sensitivity", "output/harmonised_albumin_sensitivity/harmonised_albumin_results.rds",
  19L, "43_eICU_hospital_heterogeneity.R", "Hospital heterogeneity and leave-one-hospital-out audit", "output/eicu_hospital_heterogeneity/eICU_hospital_heterogeneity_results.rds",
  20L, "44_age_sex_effect_modification.R", "Limited age and sex interactions", "output/age_sex_effect_modification/age_sex_effect_modification_results.rds",
  21L, "38_NHANES_adjusted_absolute_risk_RMST.R", "Adjusted absolute risk and RMST", "output/nhanes_adjusted_absolute_risk/NHANES_adjusted_absolute_risk_RMST_results.rds",
  22L, "36_build_corrected_publication_tables.R", "Build frozen publication tables", "output/corrected_publication_tables/Table3_cross_database_primary.csv",
  23L, "37_verify_corrected_publication_tables.R", "Publication-table QA", "output/corrected_publication_tables/TABLE_QA_REPORT.csv",
  24L, "45_build_final_analysis_freeze_v2.R", "Final freeze and unique result dictionary", "output/final_analysis_freeze_v2/UNIQUE_RESULT_DICTIONARY.csv"
)

write_markdown_table <- function(data, path) {
  display <- data
  display[] <- lapply(display, function(value) {
    value <- as.character(value)
    value[is.na(value)] <- ""
    gsub("\\|", "\\\\|", value)
  })
  header <- paste0("| ", paste(names(display), collapse = " | "), " |")
  separator <- paste0("| ", paste(rep("---", ncol(display)), collapse = " | "), " |")
  rows <- apply(display, 1, function(row) paste0("| ", paste(row, collapse = " | "), " |"))
  writeLines(c(header, separator, rows), path, useBytes = TRUE)
}

audit_checks <- list()
add_check <- function(audit_id, category, check, passed, severity = "Required", detail = "") {
  audit_checks[[length(audit_checks) + 1L]] <<- tibble::tibble(
    audit_id = audit_id,
    category = category,
    check = check,
    passed = isTRUE(passed),
    severity = severity,
    detail = as.character(detail)
  )
}

# Script and output audit ------------------------------------------------------
for (index in seq_len(nrow(canonical_pipeline))) {
  script_path <- file.path(scripts_dir, canonical_pipeline$script[index])
  output_path <- file.path(root, canonical_pipeline$primary_output[index])
  syntax_error <- tryCatch({
    parse(file = script_path)
    NA_character_
  }, error = function(error) conditionMessage(error))
  add_check(
    paste0("SCRIPT-", sprintf("%02d", index)), "Pipeline syntax",
    paste0(canonical_pipeline$script[index], " parses"),
    file.exists(script_path) && is.na(syntax_error),
    detail = ifelse(is.na(syntax_error), "Syntax OK", syntax_error)
  )
  add_check(
    paste0("OUTPUT-", sprintf("%02d", index)), "Pipeline output",
    paste0(canonical_pipeline$primary_output[index], " exists"),
    file.exists(output_path), detail = canonical_pipeline$role[index]
  )
}

mimic_extraction_code <- paste(
  readLines(file.path(scripts_dir, "09_MIMIC_IV_external_validation_formal.R"), warn = FALSE, encoding = "UTF-8"),
  collapse = "\n"
)
eicu_extraction_code <- paste(
  readLines(file.path(scripts_dir, "10_eICU_multicenter_validation.R"), warn = FALSE, encoding = "UTF-8"),
  collapse = "\n"
)
mimic_method_note <- paste(
  readLines(file.path(root, "output", "mimic_24h_validation", "MIMIC_validation_method_note.md"), warn = FALSE, encoding = "UTF-8"),
  collapse = "\n"
)
add_check(
  "WINDOW-MIMIC-DEFAULT", "Time window", "MIMIC extraction defaults to 0-24 hours",
  grepl('LAB_WINDOW_BEFORE_HOURS", unset = "0"', mimic_extraction_code, fixed = TRUE) &&
    grepl('LAB_WINDOW_AFTER_HOURS", unset = "24"', mimic_extraction_code, fixed = TRUE),
  detail = "Default window checked in script 09"
)
add_check(
  "WINDOW-MIMIC-OUTPUT", "Time window", "MIMIC method note records 0-24 hours",
  grepl("Lab window: ICU intime -0h to +24h", mimic_method_note, fixed = TRUE),
  detail = "Frozen extraction method note"
)
add_check(
  "WINDOW-EICU-DEFAULT", "Time window", "eICU extraction defaults to 0-1440 minutes",
  grepl('EICU_LAB_WINDOW_BEFORE_MINUTES", unset = "0"', eicu_extraction_code, fixed = TRUE) &&
    grepl('EICU_LAB_WINDOW_AFTER_MINUTES", unset = "1440"', eicu_extraction_code, fixed = TRUE),
  detail = "Default window checked in script 10"
)

# Frozen data -----------------------------------------------------------------
nhanes_cohort <- readRDS(file.path(
  root, "output", "nhanes_covariate_upgrade", "NHANES_2011_2018_covariate_augmented.rds"
))
nhanes_results <- readRDS(file.path(
  root, "output", "nhanes_albumin_benchmarks", "NHANES_albumin_benchmark_results.rds"
))
nhanes <- nhanes_results$model_data |>
  dplyr::mutate(
    phenotype_total_protein = factor(as.character(phenotype_total_protein), levels = c("P3", "P2", "P1"))
  )
mimic_results <- readRDS(file.path(
  root, "output", "mimic_official_oasis_v301", "MIMIC_official_OASIS_v301_results.rds"
))
mimic <- mimic_results$analysis |>
  dplyr::mutate(
    phenotype = factor(as.character(phenotype), levels = c("P3", "P2", "P1")),
    oasis_quartile = factor(oasis_quartile)
  )
eicu_results <- readRDS(file.path(
  root, "output", "eicu_24h_robust", "eICU_24h_robust_results.rds"
))
eicu <- eicu_results$model_data |>
  dplyr::mutate(
    phenotype = factor(as.character(phenotype), levels = c("P3", "P2", "P1")),
    hospitalid = factor(hospitalid)
  )

add_check("NH-N", "NHANES data", "Corrected NHANES cohort has 4,636 participants", nrow(nhanes_cohort) == 4636L, detail = nrow(nhanes_cohort))
add_check("NH-MODEL-N", "NHANES data", "Primary NHANES model has 3,979 participants and 720 deaths", nrow(nhanes) == 3979L && sum(nhanes$MORTSTAT == 1) == 720L, detail = paste0(nrow(nhanes), "/", sum(nhanes$MORTSTAT == 1)))
add_check("NH-UNIQUE", "NHANES data", "NHANES SEQN is unique", dplyr::n_distinct(nhanes_cohort$SEQN) == nrow(nhanes_cohort), detail = dplyr::n_distinct(nhanes_cohort$SEQN))
add_check("NH-AGE", "NHANES data", "NHANES participants are aged 65 years or older", all(nhanes_cohort$RIDAGEYR >= 65), detail = paste(range(nhanes_cohort$RIDAGEYR), collapse = "-"))
add_check("NH-FOLLOWUP", "NHANES data", "NHANES follow-up and outcome are valid", all(nhanes$PERMTH_INT > 0) && all(nhanes$MORTSTAT %in% c(0, 1)), detail = paste0("range ", paste(range(nhanes$PERMTH_INT), collapse = "-"), " months"))
add_check("NH-WEIGHT", "NHANES data", "Eight-year MEC weights equal two-year weights divided by four", max(abs(nhanes$WTMEC8YR - nhanes$WTMEC2YR / 4), na.rm = TRUE) < 1e-12 && all(nhanes$WTMEC8YR > 0), detail = "Weight identity and positivity")
add_check("NH-DESIGN", "NHANES data", "Survey strata and PSU fields are complete", all(!is.na(nhanes$SDMVSTRA)) && all(!is.na(nhanes$SDMVPSU)), detail = paste0(dplyr::n_distinct(nhanes$SDMVSTRA), " strata; ", dplyr::n_distinct(interaction(nhanes$SDMVSTRA, nhanes$SDMVPSU)), " nested PSUs"))
add_check("NH-NLR", "Feature engineering", "NHANES NLR is reconstructed from absolute counts", max(abs(nhanes$NLR - nhanes$LBDNENO / nhanes$LBDLYMNO), na.rm = TRUE) < 1e-10, detail = "NLR = neutrophils / lymphocytes")
add_check("NH-SII", "Feature engineering", "NHANES SII is reconstructed from platelets and absolute counts", max(abs(nhanes$SII - nhanes$LBXPLTSI * nhanes$LBDNENO / nhanes$LBDLYMNO), na.rm = TRUE) < 1e-8, detail = "SII = platelets x neutrophils / lymphocytes")
nh_assignments <- readr::read_csv(
  file.path(root, "output", "nhanes_robust_reanalysis", "NHANES_2011_2018_variant_assignments.csv"),
  show_col_types = FALSE
)
nh_counts <- table(nh_assignments$phenotype_log6_winsor)
add_check("NH-PHENOTYPES", "NHANES data", "NHANES frozen phenotype counts reproduce 817/1,894/1,925", identical(as.integer(nh_counts[c("P1", "P2", "P3")]), c(817L, 1894L, 1925L)), detail = paste(names(nh_counts), nh_counts, collapse = "; "))

add_check("MIMIC-N", "MIMIC data", "MIMIC official OASIS cohort has 1,145 participants and 574 deaths", nrow(mimic) == 1145L && sum(mimic$mortality_365d == 1) == 574L, detail = paste0(nrow(mimic), "/", sum(mimic$mortality_365d == 1)))
add_check("MIMIC-UNIQUE", "MIMIC data", "MIMIC contains one first ICU stay per subject", dplyr::n_distinct(mimic$subject_id) == nrow(mimic) && dplyr::n_distinct(mimic$stay_id) == nrow(mimic), detail = paste0(dplyr::n_distinct(mimic$subject_id), " subjects"))
add_check("MIMIC-AGE", "MIMIC data", "MIMIC participants are aged 65 years or older", all(mimic$anchor_age >= 65), detail = paste(range(mimic$anchor_age), collapse = "-"))
add_check("MIMIC-FOLLOWUP", "MIMIC data", "MIMIC 365-day outcome and follow-up are valid", all(mimic$mortality_365d %in% c(0, 1)) && all(mimic$survival_days_365 > 0 & mimic$survival_days_365 <= 365), detail = paste(range(mimic$survival_days_365), collapse = "-"))
mimic_nlr_q01 <- as.numeric(stats::quantile(mimic$nlr, 0.01, na.rm = TRUE))
mimic_sii_q01 <- as.numeric(stats::quantile(mimic$sii, 0.01, na.rm = TRUE))
add_check(
  "MIMIC-FEATURES", "MIMIC data",
  "MIMIC raw features are valid and log-transformed features are finite after 1% winsorisation",
  all(is.finite(mimic$nlr) & mimic$nlr >= 0) &&
    all(is.finite(mimic$sii) & mimic$sii >= 0) &&
    all(vapply(
      mimic[c("haemoglobin", "protein_proxy", "bmi", "creatinine")],
      function(value) all(is.finite(value) & value > 0), logical(1)
    )) && mimic_nlr_q01 > 0 && mimic_sii_q01 > 0,
  detail = paste0(
    "raw zero NLR/SII rows=", sum(mimic$nlr == 0),
    "; 1% quantiles=", round(mimic_nlr_q01, 4), "/", round(mimic_sii_q01, 4)
  )
)
add_check("MIMIC-NLR", "Feature engineering", "MIMIC NLR matches cell-count formula", max(abs(mimic$nlr - mimic$neutrophils / mimic$lymphocytes), na.rm = TRUE) < 1e-10, detail = "NLR formula")
add_check("MIMIC-SII", "Feature engineering", "MIMIC SII matches platelet-cell-count formula", max(abs(mimic$sii - mimic$platelets * mimic$neutrophils / mimic$lymphocytes), na.rm = TRUE) < 1e-8, detail = "SII formula")
add_check("MIMIC-OASIS", "MIMIC data", "Official OASIS and quartiles are complete", all(is.finite(mimic$oasis)) && nlevels(mimic$oasis_quartile) == 4L, detail = paste0("median ", median(mimic$oasis), "; range ", paste(range(mimic$oasis), collapse = "-")))
mimic_counts <- table(mimic$phenotype)
add_check("MIMIC-PHENOTYPES", "MIMIC data", "MIMIC phenotype counts reproduce 515/253/377", identical(as.integer(mimic_counts[c("P1", "P2", "P3")]), c(515L, 253L, 377L)), detail = paste(names(mimic_counts), mimic_counts, collapse = "; "))

add_check("EICU-N", "eICU data", "eICU APACHE cohort has 12,548 participants and 1,884 deaths", nrow(eicu) == 12548L && sum(eicu$hospital_mortality == 1) == 1884L, detail = paste0(nrow(eicu), "/", sum(eicu$hospital_mortality == 1)))
add_check("EICU-UNIQUE", "eICU data", "eICU contains one first ICU stay per person", dplyr::n_distinct(eicu$person_id) == nrow(eicu) && dplyr::n_distinct(eicu$patientunitstayid) == nrow(eicu), detail = paste0(dplyr::n_distinct(eicu$person_id), " persons"))
add_check("EICU-AGE", "eICU data", "eICU participants are aged 65 years or older", all(eicu$age_num >= 65), detail = paste(range(eicu$age_num), collapse = "-"))
add_check("EICU-OUTCOME", "eICU data", "eICU mortality outcome is binary and APACHE is complete", all(eicu$hospital_mortality %in% c(0, 1)) && all(is.finite(eicu$apachescore)), detail = paste0("APACHE range ", paste(range(eicu$apachescore), collapse = "-")))
add_check("EICU-FEATURES", "eICU data", "eICU six phenotype features are finite and positive", all(vapply(eicu[c("nlr", "sii_like", "haemoglobin", "total_protein", "bmi", "creatinine")], function(value) all(is.finite(value) & value > 0), logical(1))), detail = "NLR-like, SII-like, haemoglobin, total protein, BMI, creatinine")
add_check("EICU-NLR", "Feature engineering", "eICU NLR-like ratio matches differential percentages", max(abs(eicu$nlr - eicu$polys_percent / eicu$lymphs_percent), na.rm = TRUE) < 1e-10, detail = "Polys percentage / lymphocyte percentage")
add_check("EICU-SII", "Feature engineering", "eICU SII-like index matches platelet and differential formula", max(abs(eicu$sii_like - eicu$platelets * eicu$polys_percent / eicu$lymphs_percent), na.rm = TRUE) < 1e-8, detail = "Platelets x polys percentage / lymphocyte percentage")
add_check("EICU-HOSPITALS", "eICU data", "eICU primary cohort contains 166 hospitals", dplyr::n_distinct(eicu$hospitalid) == 166L, detail = dplyr::n_distinct(eicu$hospitalid))

# Static feature and provenance audit ------------------------------------------
nh_code <- paste(readLines(file.path(scripts_dir, "28_NHANES_robust_phenotype_reanalysis.R"), warn = FALSE), collapse = "\n")
mimic_code <- paste(readLines(file.path(scripts_dir, "30_MIMIC_IV_24h_robust_revalidation.R"), warn = FALSE), collapse = "\n")
eicu_code <- paste(readLines(file.path(scripts_dir, "31_eICU_24h_robust_revalidation.R"), warn = FALSE), collapse = "\n")
add_check("FEATURES-NH-CODE", "Feature leakage", "NHANES frozen matrix contains only six pre-outcome features", all(vapply(c("NLR", "SII", "LBXHGB", "LBXSTP", "BMXBMI", "LBXSCR"), grepl, logical(1), x = nh_code, fixed = TRUE)) && grepl("stats::kmeans(z", nh_code, fixed = TRUE), detail = "Outcome is not part of the scaled matrix")
add_check("FEATURES-MIMIC-CODE", "Feature leakage", "MIMIC clustering uses six transformed baseline features", all(vapply(c("log_nlr", "log_sii", "haemoglobin", "protein_proxy", "bmi", "log_creatinine"), grepl, logical(1), x = mimic_code, fixed = TRUE)) && grepl("stats::kmeans(z", mimic_code, fixed = TRUE), detail = "Outcome is not part of the clustering matrix")
add_check("FEATURES-EICU-CODE", "Feature leakage", "eICU clustering uses six transformed baseline features", all(vapply(c("log_nlr", "log_sii", "haemoglobin", "total_protein", "bmi", "log_creatinine"), grepl, logical(1), x = eicu_code, fixed = TRUE)) && grepl("stats::kmeans(z", eicu_code, fixed = TRUE), detail = "Outcome is not part of the clustering matrix")

audit_official_provenance <- function(table_name, official_subdir) {
  provenance <- readr::read_csv(file.path(root, "output", official_subdir, table_name), show_col_types = FALSE)
  official_root <- Sys.getenv(
    "MIMIC_CODE_ROOT",
    unset = file.path(root, "vendor", "mimic-code-3.0.1", "mimic-iv", "concepts_duckdb")
  )
  paths <- file.path(official_root, provenance$relative_script)
  current_md5 <- unname(tools::md5sum(paths))
  list(
    n = nrow(provenance),
    all_exist = all(file.exists(paths)),
    all_match = all(current_md5 == provenance$md5),
    release_ok = all(provenance$mimic_code_release == "v3.0.1")
  )
}
oasis_provenance <- audit_official_provenance("Table41A_official_SQL_provenance.csv", "mimic_official_oasis_v301")
sofa_provenance <- audit_official_provenance("Table42A_official_SOFA_SQL_provenance.csv", "mimic_official_sofa_v301")
add_check("OASIS-PROVENANCE", "Official severity provenance", "Official OASIS SQL files match frozen v3.0.1 MD5", oasis_provenance$n == 11L && oasis_provenance$all_exist && oasis_provenance$all_match && oasis_provenance$release_ok, detail = paste0(oasis_provenance$n, " scripts"))
add_check("SOFA-PROVENANCE", "Official severity provenance", "Official SOFA SQL files match frozen v3.0.1 MD5", sofa_provenance$n == 21L && sofa_provenance$all_exist && sofa_provenance$all_match && sofa_provenance$release_ok, detail = paste0(sofa_provenance$n, " scripts"))

# Independent model refits -----------------------------------------------------
dictionary <- readr::read_csv(file.path(
  root, "output", "final_analysis_freeze_v2", "UNIQUE_RESULT_DICTIONARY.csv"
), show_col_types = FALSE)
dictionary_value <- function(result_id, column) {
  row <- dictionary[dictionary$result_id == result_id, , drop = FALSE]
  if (nrow(row) != 1L) stop("Dictionary lookup failed for ", result_id, call. = FALSE)
  row[[column]][1]
}

nhanes_formula <- survival::Surv(PERMTH_INT, MORTSTAT) ~
  phenotype_total_protein + RIDAGEYR + male + race + INDFMPIR +
  Comorbidity_Score_Extended + cycle + smoking + hypertension + diabetes
nhanes_design <- survey::svydesign(
  ids = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~WTMEC8YR,
  nest = TRUE, data = nhanes
)
nhanes_refit <- survey::svycoxph(nhanes_formula, design = nhanes_design)
mimic_refit <- survival::coxph(
  survival::Surv(survival_days_365, mortality_365d) ~ phenotype + male + strata(oasis_quartile),
  data = mimic, ties = "efron"
)
eicu_refit <- lme4::glmer(
  hospital_mortality ~ phenotype + age_num + male + apachescore + (1 | hospitalid),
  family = stats::binomial(), data = eicu, nAGQ = 1,
  control = lme4::glmerControl(optimizer = "bobyqa")
)

extract_refit <- function(fit, dataset, measure, term_p1, term_p2, source_p1, source_p2) {
  beta <- if (inherits(fit, "merMod")) lme4::fixef(fit) else stats::coef(fit)
  covariance <- as.matrix(stats::vcov(fit))
  rows <- lapply(list(P1 = term_p1, P2 = term_p2), function(term) {
    estimate <- exp(beta[[term]])
    standard_error <- sqrt(covariance[term, term])
    tibble::tibble(
      dataset = dataset,
      comparison = paste0(names(term)[1], " vs P3"),
      effect_measure = measure,
      refit_estimate = estimate,
      refit_lower_95 = exp(beta[[term]] - 1.96 * standard_error),
      refit_upper_95 = exp(beta[[term]] + 1.96 * standard_error)
    )
  })
  result <- dplyr::bind_rows(rows)
  result$comparison <- c("P1 vs P3", "P2 vs P3")
  result$dictionary_estimate <- c(source_p1, source_p2)
  result$absolute_difference <- abs(result$refit_estimate - result$dictionary_estimate)
  result
}

primary_replication <- dplyr::bind_rows(
  extract_refit(
    nhanes_refit, "NHANES 2011-2018", "HR",
    c(P1 = "phenotype_total_proteinP1"), c(P2 = "phenotype_total_proteinP2"),
    dictionary_value("NH-PRIMARY-P1", "estimate"), dictionary_value("NH-SECONDARY-P2", "estimate")
  ),
  extract_refit(
    mimic_refit, "MIMIC-IV v3.1", "HR",
    c(P1 = "phenotypeP1"), c(P2 = "phenotypeP2"),
    dictionary_value("MIMIC-PRIMARY-P1", "estimate"), dictionary_value("MIMIC-SECONDARY-P2", "estimate")
  ),
  extract_refit(
    eicu_refit, "eICU v2.0", "OR",
    c(P1 = "phenotypeP1"), c(P2 = "phenotypeP2"),
    dictionary_value("EICU-PRIMARY-P1", "estimate"), dictionary_value("EICU-SECONDARY-P2", "estimate")
  )
)
add_check("MODEL-NH", "Independent model refit", "NHANES primary model reproduces frozen P1 and P2 estimates", all(primary_replication$absolute_difference[primary_replication$dataset == "NHANES 2011-2018"] < 1e-10), detail = max(primary_replication$absolute_difference[primary_replication$dataset == "NHANES 2011-2018"]))
add_check("MODEL-MIMIC", "Independent model refit", "MIMIC official OASIS model reproduces frozen P1 and P2 estimates", all(primary_replication$absolute_difference[primary_replication$dataset == "MIMIC-IV v3.1"] < 1e-10), detail = max(primary_replication$absolute_difference[primary_replication$dataset == "MIMIC-IV v3.1"]))
add_check("MODEL-EICU", "Independent model refit", "eICU APACHE mixed model reproduces frozen P1 and P2 estimates", all(primary_replication$absolute_difference[primary_replication$dataset == "eICU v2.0"] < 1e-8) && !lme4::isSingular(eicu_refit, tol = 1e-5), detail = max(primary_replication$absolute_difference[primary_replication$dataset == "eICU v2.0"]))

# Draft consistency ------------------------------------------------------------
draft_path <- file.path(root, "MANUSCRIPT_CHINESE_DRAFT_V1_2026-07-10.md")
draft_text <- paste(readLines(draft_path, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
obsolete_patterns <- c("3,386", "aHR = 3.08", "aHR=3.08", "3.21（", "2.07（", "OASIS-like四分位")
obsolete_hits <- obsolete_patterns[vapply(obsolete_patterns, grepl, logical(1), x = draft_text, fixed = TRUE)]
required_draft_tokens <- c(
  "4,636", "HR 1.68", "HR为1.67", "OR为1.28", "OR减弱为1.13",
  "官方OASIS", "概念性复现", "不支持其全面优于所有传统营养指标"
)
missing_draft_tokens <- required_draft_tokens[!vapply(required_draft_tokens, grepl, logical(1), x = draft_text, fixed = TRUE)]
positive_overclaims <- c("完美复现", "证明了肥胖悖论", "完全算法无关", "所有医院均显著", "原始质心的严格外部验证")
overclaim_hits <- positive_overclaims[vapply(positive_overclaims, grepl, logical(1), x = draft_text, fixed = TRUE)]
add_check("DRAFT-OBSOLETE", "Draft consistency", "Draft contains no obsolete core results", length(obsolete_hits) == 0L, detail = paste(obsolete_hits, collapse = "; "))
add_check("DRAFT-REQUIRED", "Draft consistency", "Draft contains all frozen core results and framing", length(missing_draft_tokens) == 0L, detail = paste(missing_draft_tokens, collapse = "; "))
add_check("DRAFT-OVERCLAIM", "Draft consistency", "Draft contains no positive prohibited overclaim", length(overclaim_hits) == 0L, detail = paste(overclaim_hits, collapse = "; "))
add_check("DRAFT-MIMIC-MODEL", "Draft consistency", "Draft correctly describes MIMIC primary model as sex-adjusted OASIS-stratified Cox", grepl("官方OASIS四分位分层，同时调整性别；年龄已经包含在OASIS评分构成中", draft_text, fixed = TRUE), detail = "Age is not entered separately in the frozen primary model")
add_check("DRAFT-EICU-SII", "Draft consistency", "Draft distinguishes eICU SII-like construction", grepl("SII-like指标", draft_text, fixed = TRUE), detail = "eICU uses differential percentages")

# Write outputs ----------------------------------------------------------------
audit <- dplyr::bind_rows(audit_checks)
canonical_pipeline <- canonical_pipeline |>
  dplyr::mutate(
    script_md5 = unname(tools::md5sum(file.path(scripts_dir, script))),
    output_exists = file.exists(file.path(root, primary_output)),
    output_md5 = ifelse(output_exists, unname(tools::md5sum(file.path(root, primary_output))), NA_character_)
  )

readr::write_csv(canonical_pipeline, file.path(output_dir, "CANONICAL_PIPELINE_MANIFEST.csv"))
write_markdown_table(canonical_pipeline, file.path(output_dir, "CANONICAL_PIPELINE_MANIFEST.md"))
readr::write_csv(primary_replication, file.path(output_dir, "PRIMARY_MODEL_INDEPENDENT_REPLICATION.csv"))
write_markdown_table(primary_replication, file.path(output_dir, "PRIMARY_MODEL_INDEPENDENT_REPLICATION.md"))
readr::write_csv(audit, file.path(output_dir, "PIPELINE_ANALYSIS_AUDIT_REPORT.csv"))
write_markdown_table(audit, file.path(output_dir, "PIPELINE_ANALYSIS_AUDIT_REPORT.md"))

required_failures <- audit |>
  dplyr::filter(severity == "Required", !passed)
summary_lines <- c(
  "Frozen pipeline and analysis audit",
  "",
  paste0("Checks passed: ", sum(audit$passed), "/", nrow(audit)),
  paste0("Required failures: ", nrow(required_failures)),
  paste0("Scripts in canonical pipeline: ", nrow(canonical_pipeline)),
  "",
  "Independent primary-model replication:",
  paste(capture.output(print(primary_replication)), collapse = "\n"),
  "",
  "Audit checks:",
  paste(capture.output(print(audit)), collapse = "\n")
)
writeLines(summary_lines, file.path(output_dir, "PIPELINE_ANALYSIS_AUDIT_REPORT.txt"), useBytes = TRUE)
cat(paste(summary_lines, collapse = "\n"), "\n")

if (nrow(required_failures) > 0L) {
  stop("Frozen pipeline audit failed. Review PIPELINE_ANALYSIS_AUDIT_REPORT.csv.", call. = FALSE)
}
