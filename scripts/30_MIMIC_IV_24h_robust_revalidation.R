# ==============================================================================
# MIMIC-IV first-24h robust de novo phenotype revalidation
# ==============================================================================

required_packages <- c("dplyr", "readr", "tibble", "survival")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop("Missing required package(s): ", paste(missing_packages, collapse = ", "), call. = FALSE)
}

project_root <- normalizePath(
  Sys.getenv("PROJECT_ROOT", unset = getwd()), winslash = "/", mustWork = TRUE
)
input_dir <- Sys.getenv(
  "MIMIC_24H_INPUT",
  unset = file.path(project_root, "output", "mimic_24h_validation")
)
input_path <- file.path(input_dir, "MIMIC_projected_albumin_proxy_dataset.csv")
oasis_path <- Sys.getenv(
  "MIMIC_SEVERITY_FILE",
  unset = file.path(
    project_root, "output", "mimic_oasislike_severity_corrected",
    "MIMIC_denovo_with_OASISlike_score.csv"
  )
)
old_path <- Sys.getenv("LEGACY_MIMIC_DATASET", unset = "")
output_dir <- Sys.getenv(
  "MIMIC_ROBUST_OUTPUT",
  unset = file.path(project_root, "output", "mimic_24h_robust_corrected_severity")
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(input_path)) {
  stop("First-24h MIMIC extraction not found. Run script 09 with the 0/+24h window first.", call. = FALSE)
}

winsorise <- function(x, lower = 0.01, upper = 0.99) {
  limits <- stats::quantile(x, c(lower, upper), na.rm = TRUE, names = FALSE)
  pmin(pmax(x, limits[1]), limits[2])
}

safe_z <- function(x) {
  current_sd <- stats::sd(x, na.rm = TRUE)
  if (!is.finite(current_sd) || current_sd == 0) return(rep(0, length(x)))
  (x - mean(x, na.rm = TRUE)) / current_sd
}

stable_ntile <- function(value, groups, tie_breaker) {
  keep <- is.finite(value) & !is.na(tie_breaker)
  output <- rep(NA_integer_, length(value))
  ordered_index <- which(keep)[order(value[keep], tie_breaker[keep])]
  output[ordered_index] <- dplyr::ntile(seq_along(ordered_index), groups)
  output
}

adjusted_rand_index <- function(x, y) {
  tab <- table(x, y)
  choose2 <- function(v) v * (v - 1) / 2
  n <- sum(tab)
  sum_cells <- sum(choose2(tab))
  sum_rows <- sum(choose2(rowSums(tab)))
  sum_cols <- sum(choose2(colSums(tab)))
  total_pairs <- choose2(n)
  expected <- sum_rows * sum_cols / total_pairs
  maximum <- 0.5 * (sum_rows + sum_cols)
  if (maximum == expected) return(1)
  (sum_cells - expected) / (maximum - expected)
}

extract_cox <- function(fit, model) {
  co <- summary(fit)$coefficients
  ci <- suppressMessages(stats::confint(fit))
  terms <- rownames(co)
  keep <- grepl("^phenotype", terms)
  p_col <- grep("Pr\\(", colnames(co), value = TRUE)[1]
  tibble::tibble(
    model = model,
    comparison = dplyr::recode(
      terms[keep],
      phenotypeP1 = "P1 vs P3",
      phenotypeP2 = "P2 vs P3"
    ),
    HR = exp(co[keep, "coef"]),
    lower_95 = exp(ci[keep, 1]),
    upper_95 = exp(ci[keep, 2]),
    p_value = co[keep, p_col],
    hazard_ratio_95ci = sprintf(
      "%.2f (%.2f-%.2f)", exp(co[keep, "coef"]), exp(ci[keep, 1]), exp(ci[keep, 2])
    )
  )
}

extract_logistic <- function(fit, model) {
  co <- summary(fit)$coefficients
  ci <- suppressMessages(stats::confint.default(fit))
  terms <- rownames(co)
  keep <- grepl("^phenotype", terms)
  tibble::tibble(
    model = model,
    comparison = dplyr::recode(
      terms[keep],
      phenotypeP1 = "P1 vs P3",
      phenotypeP2 = "P2 vs P3"
    ),
    OR = exp(co[keep, "Estimate"]),
    lower_95 = exp(ci[keep, 1]),
    upper_95 = exp(ci[keep, 2]),
    p_value = co[keep, "Pr(>|z|)"],
    odds_ratio_95ci = sprintf(
      "%.2f (%.2f-%.2f)", exp(co[keep, "Estimate"]), exp(ci[keep, 1]), exp(ci[keep, 2])
    )
  )
}

raw <- readr::read_csv(input_path, show_col_types = FALSE)
feature_vars <- c("nlr", "sii", "haemoglobin", "protein_proxy", "bmi", "creatinine")
required <- c(
  "stay_id", "anchor_age", "gender", "mortality_365d", "survival_days_365",
  "hospital_expire_flag", "intime", "dischtime", feature_vars
)
missing <- setdiff(required, names(raw))
if (length(missing) > 0L) {
  stop("Missing required column(s): ", paste(missing, collapse = ", "), call. = FALSE)
}

analysis <- raw |>
  dplyr::filter(dplyr::if_all(dplyr::all_of(required), ~ !is.na(.x))) |>
  dplyr::mutate(
    hospital_los_days = as.numeric(
      difftime(as.POSIXct(dischtime), as.POSIXct(intime), units = "days")
    )
  )

z <- analysis |>
  dplyr::transmute(
    log_nlr = log(winsorise(nlr)),
    log_sii = log(winsorise(sii)),
    haemoglobin = winsorise(haemoglobin),
    protein_proxy = winsorise(protein_proxy),
    bmi = winsorise(bmi),
    log_creatinine = log(winsorise(creatinine))
  ) |>
  as.data.frame() |>
  scale()

set.seed(20260710)
fit <- stats::kmeans(z, centers = 3, nstart = 100, iter.max = 500, algorithm = "Lloyd")
profiles <- analysis |>
  dplyr::mutate(cluster_raw = fit$cluster) |>
  dplyr::group_by(cluster_raw) |>
  dplyr::summarise(
    n = dplyr::n(),
    nlr = stats::median(nlr),
    sii = stats::median(sii),
    haemoglobin = stats::median(haemoglobin),
    protein_proxy = stats::median(protein_proxy),
    bmi = stats::median(bmi),
    creatinine = stats::median(creatinine),
    .groups = "drop"
  )

p1_score <- safe_z(log(profiles$nlr)) + safe_z(log(profiles$sii)) +
  safe_z(log(profiles$creatinine))
p1_raw <- profiles$cluster_raw[which.max(p1_score)]
remaining <- profiles |>
  dplyr::filter(cluster_raw != p1_raw)
p2_score <- -safe_z(remaining$haemoglobin) - safe_z(remaining$protein_proxy) -
  safe_z(remaining$bmi)
p2_raw <- remaining$cluster_raw[which.max(p2_score)]
p3_raw <- setdiff(profiles$cluster_raw, c(p1_raw, p2_raw))

analysis <- analysis |>
  dplyr::mutate(
    phenotype = dplyr::case_when(
      fit$cluster == p1_raw ~ "P1",
      fit$cluster == p2_raw ~ "P2",
      fit$cluster == p3_raw ~ "P3",
      TRUE ~ NA_character_
    ),
    phenotype = factor(phenotype, levels = c("P3", "P2", "P1")),
    male = as.integer(gender == "M")
  )

counts <- analysis |>
  dplyr::group_by(phenotype) |>
  dplyr::summarise(
    n = dplyr::n(),
    percent = 100 * n / nrow(analysis),
    deaths_365d = sum(mortality_365d == 1),
    mortality_365d_percent = 100 * mean(mortality_365d == 1),
    hospital_deaths = sum(hospital_expire_flag == 1),
    hospital_mortality_percent = 100 * mean(hospital_expire_flag == 1),
    .groups = "drop"
  )

fit_age_sex <- survival::coxph(
  survival::Surv(survival_days_365, mortality_365d) ~ phenotype + anchor_age + male,
  data = analysis,
  ties = "efron"
)

if (!file.exists(oasis_path)) {
  stop("OASIS-like score file not found.", call. = FALSE)
}
oasis <- readr::read_csv(oasis_path, show_col_types = FALSE) |>
  dplyr::select(stay_id, oasis_like_score, oasis_like_complete)
analysis_oasis <- analysis |>
  dplyr::left_join(oasis, by = "stay_id") |>
  dplyr::filter(oasis_like_complete == 1, is.finite(oasis_like_score)) |>
  dplyr::mutate(
    oasis_quartile = factor(
      stable_ntile(oasis_like_score, 4, stay_id), levels = 1:4,
      labels = c("Q1", "Q2", "Q3", "Q4")
    )
  )

fit_oasis_strata <- survival::coxph(
  survival::Surv(survival_days_365, mortality_365d) ~
    phenotype + male + strata(oasis_quartile),
  data = analysis_oasis,
  ties = "efron"
)
fit_hospital <- stats::glm(
  hospital_expire_flag ~ phenotype + male + oasis_like_score,
  family = stats::binomial(),
  data = analysis_oasis
)

landmark <- analysis_oasis |>
  dplyr::filter(survival_days_365 > 1) |>
  dplyr::mutate(
    landmark_time = survival_days_365 - 1,
    landmark_event = as.integer(mortality_365d == 1 & survival_days_365 > 1)
  )
fit_landmark <- survival::coxph(
  survival::Surv(landmark_time, landmark_event) ~
    phenotype + male + strata(oasis_quartile),
  data = landmark,
  ties = "efron"
)
fit_hospital_landmark <- stats::glm(
  hospital_expire_flag ~ phenotype + male + oasis_like_score,
  family = stats::binomial(),
  data = analysis_oasis |> dplyr::filter(hospital_los_days > 1)
)

cox_results <- dplyr::bind_rows(
  extract_cox(fit_age_sex, "365-day mortality, age-sex adjusted"),
  extract_cox(fit_oasis_strata, "365-day mortality, OASIS-like quartile-stratified"),
  extract_cox(fit_landmark, "365-day mortality, 24-hour landmark and OASIS-like stratified")
)
logistic_results <- dplyr::bind_rows(
  extract_logistic(fit_hospital, "Hospital mortality, OASIS-like adjusted"),
  extract_logistic(fit_hospital_landmark, "Hospital mortality, 24-hour landmark and OASIS-like adjusted")
)

ph <- survival::cox.zph(fit_oasis_strata)
ph_table <- tibble::tibble(
  variable = rownames(ph$table),
  chisq = ph$table[, "chisq"],
  df = ph$table[, "df"],
  p_value = ph$table[, "p"]
)

agreement <- NULL
if (file.exists(old_path)) {
  old <- readr::read_csv(old_path, show_col_types = FALSE) |>
    dplyr::select(stay_id, mimic_phenotype_raw)
  common <- analysis |>
    dplyr::select(stay_id, phenotype) |>
    dplyr::inner_join(old, by = "stay_id")
  agreement <- tibble::tibble(
    common_n = nrow(common),
    adjusted_rand_index = adjusted_rand_index(common$phenotype, common$mimic_phenotype_raw)
  )
}

readr::write_csv(profiles, file.path(output_dir, "Table30A_MIMIC_24h_profiles.csv"))
readr::write_csv(counts, file.path(output_dir, "Table30B_MIMIC_24h_counts.csv"))
readr::write_csv(cox_results, file.path(output_dir, "Table30C_MIMIC_24h_Cox_models.csv"))
readr::write_csv(logistic_results, file.path(output_dir, "Table30D_MIMIC_24h_hospital_models.csv"))
readr::write_csv(ph_table, file.path(output_dir, "Table30E_MIMIC_24h_PH_check.csv"))
if (!is.null(agreement)) {
  readr::write_csv(agreement, file.path(output_dir, "Table30F_old72h_vs_new24h_agreement.csv"))
}
readr::write_csv(analysis, file.path(output_dir, "MIMIC_24h_robust_dataset.csv"))
saveRDS(
  list(
    analysis = analysis,
    analysis_oasis = analysis_oasis,
    profiles = profiles,
    counts = counts,
    cox = cox_results,
    logistic = logistic_results,
    ph = ph_table,
    agreement = agreement
  ),
  file.path(output_dir, "MIMIC_24h_robust_results.rds")
)

summary_lines <- c(
  "MIMIC-IV first-24h robust phenotype revalidation",
  paste0("Complete albumin-proxy cohort: n = ", nrow(analysis)),
  paste0("OASIS-like complete cohort: n = ", nrow(analysis_oasis)),
  "",
  "Phenotype counts:",
  paste(capture.output(print(counts)), collapse = "\n"),
  "",
  "Cox models:",
  paste(capture.output(print(cox_results)), collapse = "\n"),
  "",
  "Hospital mortality models:",
  paste(capture.output(print(logistic_results)), collapse = "\n")
)
writeLines(summary_lines, file.path(output_dir, "MIMIC_24h_robust_summary.txt"))
message("MIMIC-IV first-24h robust revalidation completed.")
cat(paste(summary_lines, collapse = "\n"), "\n")
