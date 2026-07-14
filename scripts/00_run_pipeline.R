# ==============================================================================
# Stage-specific runner for the Freeze V2 release candidate.
# Run from the repository root. Each analysis runs in a separate R process.
# ==============================================================================

root <- normalizePath(
  Sys.getenv("PROJECT_ROOT", unset = getwd()), winslash = "/", mustWork = TRUE
)
script_dir <- file.path(root, "scripts")
rscript <- file.path(R.home("bin"), "Rscript")
if (.Platform$OS.type == "windows" && file.exists(paste0(rscript, ".exe"))) {
  rscript <- paste0(rscript, ".exe")
}

stages <- list(
  nhanes = c(
    "27_rebuild_NHANES_2011_2018.R",
    "27B_augment_NHANES_covariates.R",
    "28_NHANES_robust_phenotype_reanalysis.R",
    "28B_NHANES_cluster_number_validation.R",
    "29_NHANES_albumin_and_conventional_benchmarks.R",
    "32_NHANES_cluster_aware_bootstrap_and_PH.R",
    "32B_NHANES_time_varying_score.R",
    "33_NHANES_selection_bias_sensitivity.R",
    "34_NHANES_bootstrap_Cindex_increment.R",
    "35_NHANES_leave_one_cycle_out_validation.R"
  ),
  mimic = c(
    "09A_build_NHANES_projection_reference.R",
    "09_MIMIC_IV_external_validation_formal.R",
    "09B_MIMIC_IV_de_novo_reproducibility.R",
    "25_MIMIC_OASISlike_severity_adjustment.R",
    "30_MIMIC_IV_24h_robust_revalidation.R",
    "41_MIMIC_official_OASIS_v301.R",
    "42_MIMIC_official_first_day_SOFA_v301.R"
  ),
  eicu = c(
    "10_eICU_multicenter_validation.R",
    "14_eICU_APACHE_adjusted_sensitivity.R",
    "31_eICU_24h_robust_revalidation.R"
  ),
  cross_database = c(
    "39_ICU_complete_case_selection_audit.R",
    "40_harmonised_albumin_cross_database_sensitivity.R",
    "43_eICU_hospital_heterogeneity.R",
    "44_age_sex_effect_modification.R",
    "38_NHANES_adjusted_absolute_risk_RMST.R",
    "36_build_corrected_publication_tables.R",
    "37_verify_corrected_publication_tables.R",
    "45_build_final_analysis_freeze_v2.R"
  ),
  postfreeze = c(
    "47_alternative_clustering_freeze_v2.R",
    "48_NHANES_nonlinearity_threshold_freeze_v2.R",
    "70_NHANES_amendment14_cluster_sensitivity.R",
    "71_audit_amendment14.R"
  ),
  exploratory_tool = c(
    "50_NHANES_locked_36m_risk_model.R",
    "51_NHANES_2005_2010_historical_transport_validation.R"
  ),
  internal_audit = "46_audit_frozen_pipeline_and_analysis.R"
)

args <- commandArgs(trailingOnly = TRUE)
stage <- if (length(args) > 0L) args[[1L]] else Sys.getenv("PIPELINE_STAGE", unset = "nhanes")
if (!stage %in% names(stages)) {
  stop("Unknown stage. Choose one of: ", paste(names(stages), collapse = ", "), call. = FALSE)
}

for (script in stages[[stage]]) {
  path <- file.path(script_dir, script)
  if (!file.exists(path)) stop("Missing script: ", path, call. = FALSE)
  message("Running ", script, " ...")
  status <- system2(rscript, path, stdout = "", stderr = "")
  if (!identical(status, 0L)) {
    stop("Pipeline stopped because ", script, " returned status ", status, call. = FALSE)
  }
}

message("Stage completed: ", stage)
