# ==============================================================================
# eICU first-24h robust multicentre phenotype revalidation
# ==============================================================================

required_packages <- c("data.table", "dplyr", "readr", "tibble", "lme4")
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
  "EICU_24H_INPUT",
  unset = file.path(project_root, "output", "eicu_24h_extraction")
)
input_path <- file.path(input_dir, "eICU_denovo_strict_total_protein_dataset.csv")
apache_path <- Sys.getenv(
  "EICU_APACHE_FILE",
  unset = file.path(
    project_root, "output", "eicu_apache_support",
    "eICU_APACHE_adjusted_sensitivity.rds"
  )
)
old_path <- Sys.getenv("LEGACY_EICU_DATASET", unset = "")
output_dir <- Sys.getenv(
  "EICU_ROBUST_OUTPUT",
  unset = file.path(project_root, "output", "eicu_24h_robust")
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(input_path) || !file.exists(apache_path)) {
  stop("Run the first-24h eICU extraction and the APACHE audit before this analysis.", call. = FALSE)
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

extract_mixed_effects <- function(fit, model) {
  co <- summary(fit)$coefficients
  beta <- co[, "Estimate"]
  se <- co[, "Std. Error"]
  terms <- rownames(co)
  keep <- grepl("^phenotype", terms)
  z_value <- beta / se
  p_value <- 2 * stats::pnorm(abs(z_value), lower.tail = FALSE)
  tibble::tibble(
    model = model,
    comparison = dplyr::recode(
      terms[keep],
      phenotypeP1 = "P1 vs P3",
      phenotypeP2 = "P2 vs P3"
    ),
    OR = exp(beta[keep]),
    lower_95 = exp(beta[keep] - 1.96 * se[keep]),
    upper_95 = exp(beta[keep] + 1.96 * se[keep]),
    p_value = p_value[keep],
    odds_ratio_95ci = sprintf(
      "%.2f (%.2f-%.2f)",
      exp(beta[keep]), exp(beta[keep] - 1.96 * se[keep]), exp(beta[keep] + 1.96 * se[keep])
    )
  )
}

raw <- data.table::fread(input_path, showProgress = FALSE)
required <- c(
  "patientunitstayid", "hospitalid", "age_num", "gender_model",
  "hospital_mortality", "icu_mortality", "survival_days_hosp", "unitdischargeoffset",
  "nlr", "sii_like", "haemoglobin", "total_protein", "bmi", "creatinine"
)
missing <- setdiff(required, names(raw))
if (length(missing) > 0L) {
  stop("Missing required column(s): ", paste(missing, collapse = ", "), call. = FALSE)
}

analysis <- as.data.frame(raw) |>
  dplyr::filter(dplyr::if_all(dplyr::all_of(required), ~ !is.na(.x))) |>
  dplyr::mutate(
    icu_los_days = as.numeric(unitdischargeoffset) / 1440,
    male = as.integer(gender_model == "Male")
  )

z <- analysis |>
  dplyr::transmute(
    log_nlr = log(winsorise(nlr)),
    log_sii = log(winsorise(sii_like)),
    haemoglobin = winsorise(haemoglobin),
    total_protein = winsorise(total_protein),
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
    sii = stats::median(sii_like),
    haemoglobin = stats::median(haemoglobin),
    total_protein = stats::median(total_protein),
    bmi = stats::median(bmi),
    creatinine = stats::median(creatinine),
    .groups = "drop"
  )

p1_score <- safe_z(log(profiles$nlr)) + safe_z(log(profiles$sii)) +
  safe_z(log(profiles$creatinine))
p1_raw <- profiles$cluster_raw[which.max(p1_score)]
remaining <- profiles |>
  dplyr::filter(cluster_raw != p1_raw)
p2_score <- -safe_z(remaining$haemoglobin) - safe_z(remaining$total_protein) -
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
    phenotype = factor(phenotype, levels = c("P3", "P2", "P1"))
  )

apache_results <- readRDS(apache_path)
apache <- apache_results$model_data |>
  dplyr::select(patientunitstayid, apachescore) |>
  dplyr::distinct(patientunitstayid, .keep_all = TRUE)
model_data <- analysis |>
  dplyr::left_join(apache, by = "patientunitstayid") |>
  dplyr::filter(is.finite(apachescore))

fit_hospital <- lme4::glmer(
  hospital_mortality ~ phenotype + age_num + male + apachescore + (1 | hospitalid),
  family = stats::binomial(), data = model_data, nAGQ = 1,
  control = lme4::glmerControl(optimizer = "bobyqa")
)
fit_icu <- lme4::glmer(
  icu_mortality ~ phenotype + age_num + male + apachescore + (1 | hospitalid),
  family = stats::binomial(), data = model_data, nAGQ = 1,
  control = lme4::glmerControl(optimizer = "bobyqa")
)
hospital_landmark <- model_data |>
  dplyr::filter(survival_days_hosp > 1)
icu_landmark <- model_data |>
  dplyr::filter(icu_los_days > 1)
fit_hospital_landmark <- lme4::glmer(
  hospital_mortality ~ phenotype + age_num + male + apachescore + (1 | hospitalid),
  family = stats::binomial(), data = hospital_landmark, nAGQ = 1,
  control = lme4::glmerControl(optimizer = "bobyqa")
)
fit_icu_landmark <- lme4::glmer(
  icu_mortality ~ phenotype + age_num + male + apachescore + (1 | hospitalid),
  family = stats::binomial(), data = icu_landmark, nAGQ = 1,
  control = lme4::glmerControl(optimizer = "bobyqa")
)

models <- dplyr::bind_rows(
  extract_mixed_effects(fit_hospital, "Hospital mortality, APACHE-adjusted mixed effects"),
  extract_mixed_effects(fit_icu, "ICU mortality, APACHE-adjusted mixed effects"),
  extract_mixed_effects(
    fit_hospital_landmark,
    "Hospital mortality, 24-hour landmark and APACHE-adjusted"
  ),
  extract_mixed_effects(
    fit_icu_landmark,
    "ICU mortality, 24-hour landmark and APACHE-adjusted"
  )
)

counts <- analysis |>
  dplyr::group_by(phenotype) |>
  dplyr::summarise(
    n = dplyr::n(),
    percent = 100 * n / nrow(analysis),
    hospital_deaths = sum(hospital_mortality == 1),
    hospital_mortality_percent = 100 * mean(hospital_mortality == 1),
    icu_deaths = sum(icu_mortality == 1),
    icu_mortality_percent = 100 * mean(icu_mortality == 1),
    .groups = "drop"
  )

agreement <- NULL
if (file.exists(old_path)) {
  old <- data.table::fread(old_path, select = c("patientunitstayid", "eicu_phenotype_raw"))
  common <- analysis |>
    dplyr::select(patientunitstayid, phenotype) |>
    dplyr::inner_join(as.data.frame(old), by = "patientunitstayid")
  agreement <- tibble::tibble(
    common_n = nrow(common),
    adjusted_rand_index = adjusted_rand_index(common$phenotype, common$eicu_phenotype_raw)
  )
}

readr::write_csv(profiles, file.path(output_dir, "Table31A_eICU_24h_profiles.csv"))
readr::write_csv(counts, file.path(output_dir, "Table31B_eICU_24h_counts.csv"))
readr::write_csv(models, file.path(output_dir, "Table31C_eICU_24h_APACHE_models.csv"))
if (!is.null(agreement)) {
  readr::write_csv(agreement, file.path(output_dir, "Table31D_old72h_vs_new24h_agreement.csv"))
}
readr::write_csv(analysis, file.path(output_dir, "eICU_24h_robust_dataset.csv"))
saveRDS(
  list(
    analysis = analysis,
    model_data = model_data,
    profiles = profiles,
    counts = counts,
    models = models,
    agreement = agreement
  ),
  file.path(output_dir, "eICU_24h_robust_results.rds")
)

summary_lines <- c(
  "eICU first-24h robust multicentre phenotype revalidation",
  paste0("Complete strict total-protein cohort: n = ", nrow(analysis)),
  paste0("APACHE-complete cohort: n = ", nrow(model_data)),
  "",
  "Phenotype counts:",
  paste(capture.output(print(counts)), collapse = "\n"),
  "",
  "APACHE-adjusted hospital-random-intercept models:",
  paste(capture.output(print(models)), collapse = "\n")
)
writeLines(summary_lines, file.path(output_dir, "eICU_24h_robust_summary.txt"))
message("eICU first-24h robust revalidation completed.")
cat(paste(summary_lines, collapse = "\n"), "\n")
