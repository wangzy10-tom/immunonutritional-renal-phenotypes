# ==============================================================================
# Independent audit of Freeze Amendment 14
# ==============================================================================

required_packages <- c("dplyr", "readr", "tibble", "survival", "survey")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0L) {
  stop("Missing required package(s): ", paste(missing_packages, collapse = ", "), call. = FALSE)
}
options(survey.lonely.psu = "adjust")

root <- normalizePath(
  Sys.getenv("PROJECT_ROOT", unset = getwd()), winslash = "/", mustWork = TRUE
)
output_dir <- file.path(root, "output", "amendment14_nhanes_cluster_sensitivity")
source_manifest_path <- file.path(output_dir, "AMENDMENT14_SOURCE_MANIFEST.csv")
results_path <- file.path(output_dir, "AMENDMENT14_RESULTS.rds")
qa_path <- file.path(output_dir, "AMENDMENT14_QA_REPORT.csv")
model_k_path <- file.path(output_dir, "TableA14_4_K2_K6_survey_cox.csv")
model_single_path <- file.path(output_dir, "TableA14_7_NLR_SII_single_index_survey_cox.csv")

source_manifest <- readr::read_csv(source_manifest_path, show_col_types = FALSE)
manifest_paths <- file.path(root, source_manifest$source_file)
manifest_exists <- file.exists(manifest_paths)
manifest_current_md5 <- rep(NA_character_, length(manifest_paths))
manifest_current_md5[manifest_exists] <- unname(tools::md5sum(manifest_paths[manifest_exists]))
manifest_match <- manifest_exists & manifest_current_md5 == source_manifest$md5

cohort <- readRDS(file.path(root, "output", "nhanes_covariate_upgrade", "NHANES_2011_2018_covariate_augmented.rds"))
model_object <- readRDS(file.path(root, "output", "nhanes_albumin_benchmarks", "NHANES_albumin_benchmark_results.rds"))
model_data <- model_object$model_data
frozen_assignments <- readr::read_csv(
  file.path(root, "output", "nhanes_robust_reanalysis", "NHANES_2011_2018_variant_assignments.csv"),
  show_col_types = FALSE
)
dictionary <- readr::read_csv(
  file.path(root, "output", "final_analysis_freeze_v2", "UNIQUE_RESULT_DICTIONARY.csv"),
  show_col_types = FALSE
)
saved <- readRDS(results_path)
internal_qa <- readr::read_csv(qa_path, show_col_types = FALSE)
reported_k <- readr::read_csv(model_k_path, show_col_types = FALSE)
reported_single <- readr::read_csv(model_single_path, show_col_types = FALSE)

clip_1_99 <- function(x) {
  cut_points <- stats::quantile(x, c(0.01, 0.99), na.rm = TRUE, names = FALSE)
  pmax(cut_points[1], pmin(cut_points[2], x))
}

standardise <- function(x) {
  spread <- stats::sd(x)
  if (!is.finite(spread) || spread == 0) return(rep(0, length(x)))
  (x - mean(x)) / spread
}

x_all <- data.frame(
  NLR = log(clip_1_99(cohort$NLR)),
  SII = log(clip_1_99(cohort$SII)),
  LBXHGB = clip_1_99(cohort$LBXHGB),
  LBXSTP = clip_1_99(cohort$LBXSTP),
  BMXBMI = clip_1_99(cohort$BMXBMI),
  LBXSCR = log(clip_1_99(cohort$LBXSCR))
)
allowed_feature_names <- c("NLR", "SII", "LBXHGB", "LBXSTP", "BMXBMI", "LBXSCR")

cluster_profile <- function(cluster) {
  data.frame(
    cluster_raw = sort(unique(cluster)),
    n = as.integer(table(factor(cluster, levels = sort(unique(cluster)))))
  ) |>
    dplyr::left_join(
      cohort |>
        dplyr::mutate(cluster_raw = cluster) |>
        dplyr::group_by(cluster_raw) |>
        dplyr::summarise(
          nlr = stats::median(NLR),
          sii = stats::median(SII),
          haemoglobin = stats::median(LBXHGB),
          protein = stats::median(LBXSTP),
          bmi = stats::median(BMXBMI),
          creatinine = stats::median(LBXSCR),
          .groups = "drop"
        ),
      by = "cluster_raw"
    )
}

assign_three <- function(cluster, variant) {
  profile <- cluster_profile(cluster)
  score_parts <- list()
  if (variant != "sii_only") score_parts[[length(score_parts) + 1L]] <- standardise(log(profile$nlr))
  if (variant != "nlr_only") score_parts[[length(score_parts) + 1L]] <- standardise(log(profile$sii))
  score_parts[[length(score_parts) + 1L]] <- standardise(log(profile$creatinine))
  p1_cluster <- profile$cluster_raw[which.max(Reduce(`+`, score_parts))]
  rest <- profile[profile$cluster_raw != p1_cluster, , drop = FALSE]
  reserve_score <- -standardise(rest$haemoglobin) - standardise(rest$protein) - standardise(rest$bmi)
  p2_cluster <- rest$cluster_raw[which.max(reserve_score)]
  p3_cluster <- setdiff(profile$cluster_raw, c(p1_cluster, p2_cluster))
  factor(
    dplyr::case_when(
      cluster == p1_cluster ~ "P1",
      cluster == p2_cluster ~ "P2",
      cluster == p3_cluster ~ "P3",
      TRUE ~ NA_character_
    ),
    levels = c("P3", "P2", "P1")
  )
}

assign_ranked <- function(cluster, k) {
  profile <- cluster_profile(cluster)
  vulnerability <- rowSums(data.frame(
    high_nlr = standardise(log(profile$nlr)),
    high_sii = standardise(log(profile$sii)),
    low_haemoglobin = standardise(-profile$haemoglobin),
    low_protein = standardise(-profile$protein),
    low_bmi = standardise(-profile$bmi),
    high_creatinine = standardise(log(profile$creatinine))
  ))
  ordered_clusters <- profile$cluster_raw[order(vulnerability, decreasing = TRUE)]
  labels <- if (k == 2L) c("Higher vulnerability", "Lower vulnerability") else paste0("V", seq_len(k))
  assigned <- labels[match(cluster, ordered_clusters)]
  levels_out <- if (k == 2L) c("Lower vulnerability", "Higher vulnerability") else paste0("V", k:1)
  factor(assigned, levels = levels_out)
}

independent_cluster <- function(features, k, seed, labeller, variant = NULL) {
  matrix_z <- scale(x_all[features])
  set.seed(seed)
  fit <- stats::kmeans(matrix_z, centers = k, nstart = 100, iter.max = 500, algorithm = "Lloyd")
  if (labeller == "three") assign_three(fit$cluster, variant) else assign_ranked(fit$cluster, k)
}

independent_assignments <- tibble::tibble(
  SEQN = cohort$SEQN,
  primary_k3 = as.character(independent_cluster(allowed_feature_names, 3L, 20260710L, "three", "six")),
  k2 = as.character(independent_cluster(allowed_feature_names, 2L, 20260712L, "ranked")),
  k6 = as.character(independent_cluster(allowed_feature_names, 6L, 20260716L, "ranked")),
  nlr_only_k3 = as.character(independent_cluster(c("NLR", "LBXHGB", "LBXSTP", "BMXBMI", "LBXSCR"), 3L, 20260710L, "three", "nlr_only")),
  sii_only_k3 = as.character(independent_cluster(c("SII", "LBXHGB", "LBXSTP", "BMXBMI", "LBXSCR"), 3L, 20260710L, "three", "sii_only"))
)

fit_model <- function(assignments, analysis, k, seed) {
  joined <- model_data |>
    dplyr::left_join(
      tibble::tibble(SEQN = cohort$SEQN, phenotype_audit = assignments),
      by = "SEQN"
    )
  levels_model <- if (k == 2L) {
    c("Lower vulnerability", "Higher vulnerability")
  } else if (k == 6L) {
    paste0("V", 6:1)
  } else {
    c("P3", "P2", "P1")
  }
  reference <- levels_model[1]
  joined$phenotype_audit <- factor(joined$phenotype_audit, levels = levels_model)
  design <- survey::svydesign(
    ids = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~WTMEC8YR,
    nest = TRUE, data = joined
  )
  fit <- survey::svycoxph(
    survival::Surv(PERMTH_INT, MORTSTAT) ~ phenotype_audit + RIDAGEYR + male + race +
      INDFMPIR + Comorbidity_Score_Extended + cycle + smoking + hypertension + diabetes,
    design = design
  )
  co <- summary(fit)$coefficients
  ci <- suppressMessages(stats::confint(fit))
  terms <- rownames(co)
  keep <- startsWith(terms, "phenotype_audit")
  p_col <- grep("Pr\\(", colnames(co), value = TRUE)[1]
  level <- sub("^phenotype_audit", "", terms[keep])
  tibble::tibble(
    analysis = analysis,
    k = k,
    seed = seed,
    comparison = paste0(level, " vs ", reference),
    audit_HR = exp(co[keep, "coef"]),
    audit_lower_95 = exp(ci[keep, 1]),
    audit_upper_95 = exp(ci[keep, 2]),
    audit_p_value = co[keep, p_col]
  )
}

audit_models <- dplyr::bind_rows(
  fit_model(independent_assignments$k2, "K=2 six-feature sensitivity", 2L, 20260712L),
  fit_model(independent_assignments$k6, "K=6 six-feature sensitivity", 6L, 20260716L),
  fit_model(independent_assignments$nlr_only_k3, "K=3 NLR-only inflammation sensitivity", 3L, 20260710L),
  fit_model(independent_assignments$sii_only_k3, "K=3 SII-only inflammation sensitivity", 3L, 20260710L)
)
reported_models <- dplyr::bind_rows(
  reported_k |>
    dplyr::select(analysis, k, seed, comparison, HR, lower_95, upper_95, p_value),
  reported_single |>
    dplyr::select(analysis, k, seed, comparison, HR, lower_95, upper_95, p_value)
) |>
  dplyr::rename(
    reported_HR = HR,
    reported_lower_95 = lower_95,
    reported_upper_95 = upper_95,
    reported_p_value = p_value
  )
replication <- audit_models |>
  dplyr::left_join(reported_models, by = c("analysis", "k", "seed", "comparison")) |>
  dplyr::mutate(
    HR_absolute_difference = abs(audit_HR - reported_HR),
    lower_absolute_difference = abs(audit_lower_95 - reported_lower_95),
    upper_absolute_difference = abs(audit_upper_95 - reported_upper_95),
    p_absolute_difference = abs(audit_p_value - reported_p_value)
  )

frozen_primary <- frozen_assignments$phenotype_log6_winsor[match(cohort$SEQN, frozen_assignments$SEQN)]
dictionary_p1 <- as.numeric(dictionary$estimate[dictionary$result_id == "NH-PRIMARY-P1"])
dictionary_p2 <- as.numeric(dictionary$estimate[dictionary$result_id == "NH-SECONDARY-P2"])
model_data$phenotype_audit <- factor(
  independent_assignments$primary_k3[match(model_data$SEQN, independent_assignments$SEQN)],
  levels = c("P3", "P2", "P1")
)
primary_design <- survey::svydesign(
  ids = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~WTMEC8YR,
  nest = TRUE, data = model_data
)
primary_fit <- survey::svycoxph(
  survival::Surv(PERMTH_INT, MORTSTAT) ~ phenotype_audit + RIDAGEYR + male + race +
    INDFMPIR + Comorbidity_Score_Extended + cycle + smoking + hypertension + diabetes,
  design = primary_design
)
primary_co <- stats::coef(primary_fit)
primary_p1 <- exp(primary_co[["phenotype_auditP1"]])
primary_p2 <- exp(primary_co[["phenotype_auditP2"]])

checks <- tibble::tribble(
  ~check_id, ~description, ~passed, ~detail,
  "A14A-MANIFEST", "Every locked source and output matches its recorded MD5", all(manifest_match), paste0(sum(manifest_match), "/", length(manifest_match)),
  "A14A-INTERNAL-QA", "All analysis-program QA checks passed", nrow(internal_qa) > 0L && all(internal_qa$passed), paste0(sum(internal_qa$passed), "/", nrow(internal_qa)),
  "A14A-FEATURE-BOUNDARY", "Independent clustering matrix contains only the six locked pre-outcome biomarkers", identical(names(x_all), allowed_feature_names), paste(names(x_all), collapse = "; "),
  "A14A-PRIMARY-FROZEN", "Independent K=3 assignment equals the frozen assignment for all participants", identical(independent_assignments$primary_k3, as.character(frozen_primary)), paste0(sum(independent_assignments$primary_k3 == frozen_primary), "/4636"),
  "A14A-PRIMARY-SAVED", "Independent K=3 assignment equals the Amendment 14 saved assignment", identical(independent_assignments$primary_k3, saved$assignments$primary_k3), paste0(sum(independent_assignments$primary_k3 == saved$assignments$primary_k3), "/4636"),
  "A14A-K2-SAVED", "Independent K=2 assignment equals the Amendment 14 saved assignment", identical(independent_assignments$k2, saved$assignments$k2), paste0(sum(independent_assignments$k2 == saved$assignments$k2), "/4636"),
  "A14A-K6-SAVED", "Independent K=6 assignment equals the Amendment 14 saved assignment", identical(independent_assignments$k6, saved$assignments$k6), paste0(sum(independent_assignments$k6 == saved$assignments$k6), "/4636"),
  "A14A-NLR-SAVED", "Independent NLR-only assignment equals the Amendment 14 saved assignment", identical(independent_assignments$nlr_only_k3, saved$assignments$nlr_only_k3), paste0(sum(independent_assignments$nlr_only_k3 == saved$assignments$nlr_only_k3), "/4636"),
  "A14A-SII-SAVED", "Independent SII-only assignment equals the Amendment 14 saved assignment", identical(independent_assignments$sii_only_k3, saved$assignments$sii_only_k3), paste0(sum(independent_assignments$sii_only_k3 == saved$assignments$sii_only_k3), "/4636"),
  "A14A-PRIMARY-MODEL", "Independent primary Cox model reproduces both frozen estimates within 1e-10", abs(primary_p1 - dictionary_p1) < 1e-10 && abs(primary_p2 - dictionary_p2) < 1e-10, paste(format(c(abs(primary_p1 - dictionary_p1), abs(primary_p2 - dictionary_p2)), scientific = TRUE), collapse = "/"),
  "A14A-SENSITIVITY-MODELS", "Independent sensitivity Cox models reproduce every reported number within 1e-10", nrow(replication) == 10L && all(replication$HR_absolute_difference < 1e-10) && all(replication$lower_absolute_difference < 1e-10) && all(replication$upper_absolute_difference < 1e-10) && all(replication$p_absolute_difference < 1e-10), sprintf("max difference %.3e", max(c(replication$HR_absolute_difference, replication$lower_absolute_difference, replication$upper_absolute_difference, replication$p_absolute_difference))),
  "A14A-MODEL-N", "Independent models retain the locked 3,979 participants and 720 deaths", nrow(model_data) == 3979L && sum(model_data$MORTSTAT == 1L) == 720L, paste0(nrow(model_data), "/", sum(model_data$MORTSTAT == 1L))
)

audit_report_path <- file.path(output_dir, "AMENDMENT14_INDEPENDENT_AUDIT.csv")
replication_path <- file.path(output_dir, "AMENDMENT14_INDEPENDENT_MODEL_REPLICATION.csv")
readr::write_csv(checks, audit_report_path)
readr::write_csv(replication, replication_path)

audit_manifest <- tibble::tibble(
  source_file = c(
    "71_audit_amendment14.R",
    "output/amendment14_nhanes_cluster_sensitivity/AMENDMENT14_SOURCE_MANIFEST.csv",
    "output/amendment14_nhanes_cluster_sensitivity/AMENDMENT14_INDEPENDENT_AUDIT.csv",
    "output/amendment14_nhanes_cluster_sensitivity/AMENDMENT14_INDEPENDENT_MODEL_REPLICATION.csv"
  ),
  role = c("independent_audit_script", "analysis_source_manifest", "independent_audit_report", "independent_model_replication")
) |>
  dplyr::mutate(
    bytes = as.numeric(file.info(file.path(root, source_file))$size),
    md5 = unname(tools::md5sum(file.path(root, source_file)))
  )
readr::write_csv(audit_manifest, file.path(output_dir, "AMENDMENT14_AUDIT_MANIFEST.csv"))

print(checks)
if (!all(checks$passed)) {
  stop("Independent Amendment 14 audit failed.", call. = FALSE)
}
message("Independent Amendment 14 audit passed: ", nrow(checks), "/", nrow(checks), " checks.")
