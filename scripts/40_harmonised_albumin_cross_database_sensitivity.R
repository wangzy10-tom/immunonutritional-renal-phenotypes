# ==============================================================================
# Post-freeze harmonised albumin sensitivity across NHANES, MIMIC-IV, and eICU
# ==============================================================================

required_packages <- c(
  "data.table", "dplyr", "readr", "tibble", "survival", "survey", "lme4"
)
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
eicu_dir <- Sys.getenv(
  "EICU_DIR", unset = file.path(root, "data", "eicu-crd-2.0")
)
output_dir <- file.path(root, "output", "harmonised_albumin_sensitivity")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

winsorise <- function(x) {
  limits <- stats::quantile(x, c(0.01, 0.99), na.rm = TRUE, names = FALSE)
  pmin(pmax(as.numeric(x), limits[1]), limits[2])
}

safe_z <- function(x) {
  spread <- stats::sd(x, na.rm = TRUE)
  if (!is.finite(spread) || spread == 0) return(rep(0, length(x)))
  (x - mean(x, na.rm = TRUE)) / spread
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

fit_albumin_clusters <- function(data, nlr, sii, haemoglobin, albumin, bmi, creatinine) {
  matrix_data <- data |>
    dplyr::transmute(
      log_nlr = log(pmax(winsorise(.data[[nlr]]), .Machine$double.eps)),
      log_sii = log(pmax(winsorise(.data[[sii]]), .Machine$double.eps)),
      haemoglobin = winsorise(.data[[haemoglobin]]),
      albumin = winsorise(.data[[albumin]]),
      bmi = winsorise(.data[[bmi]]),
      log_creatinine = log(pmax(winsorise(.data[[creatinine]]), .Machine$double.eps))
    ) |>
    as.data.frame() |>
    scale()

  fit <- stats::kmeans(
    matrix_data, centers = 3, nstart = 100, iter.max = 500, algorithm = "Lloyd"
  )
  profiles <- data |>
    dplyr::mutate(cluster_raw = fit$cluster) |>
    dplyr::group_by(cluster_raw) |>
    dplyr::summarise(
      n = dplyr::n(),
      nlr = stats::median(.data[[nlr]]),
      sii = stats::median(.data[[sii]]),
      haemoglobin = stats::median(.data[[haemoglobin]]),
      albumin = stats::median(.data[[albumin]]),
      bmi = stats::median(.data[[bmi]]),
      creatinine = stats::median(.data[[creatinine]]),
      .groups = "drop"
    )

  p1_score <- safe_z(log(profiles$nlr)) + safe_z(log(profiles$sii)) +
    safe_z(log(profiles$creatinine))
  p1_raw <- profiles$cluster_raw[which.max(p1_score)]
  remaining <- profiles |>
    dplyr::filter(cluster_raw != p1_raw)
  p2_score <- -safe_z(remaining$haemoglobin) - safe_z(remaining$albumin) -
    safe_z(remaining$bmi)
  p2_raw <- remaining$cluster_raw[which.max(p2_score)]
  p3_raw <- setdiff(profiles$cluster_raw, c(p1_raw, p2_raw))

  phenotype <- dplyr::case_when(
    fit$cluster == p1_raw ~ "P1",
    fit$cluster == p2_raw ~ "P2",
    fit$cluster == p3_raw ~ "P3",
    TRUE ~ NA_character_
  )
  mapping <- tibble::tibble(
    cluster_raw = c(p1_raw, p2_raw, p3_raw),
    phenotype = c("P1", "P2", "P3")
  )
  profiles <- profiles |>
    dplyr::left_join(mapping, by = "cluster_raw") |>
    dplyr::arrange(factor(phenotype, levels = c("P1", "P2", "P3")))

  list(
    phenotype = factor(phenotype, levels = c("P3", "P2", "P1")),
    profiles = profiles
  )
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

make_selection_weights <- function(data, model, selected_variable, dataset, feature_n) {
  if (!isTRUE(model$converged)) stop(dataset, " selection model did not converge.")
  probability <- pmin(pmax(stats::predict(model, type = "response"), 0.001), 0.999)
  selected <- data[[selected_variable]]
  rate <- mean(selected == 1)
  weight <- ifelse(selected == 1, rate / probability, NA_real_)
  limits <- stats::quantile(weight[selected == 1], c(0.01, 0.99), na.rm = TRUE)
  weight <- ifelse(selected == 1, pmin(pmax(weight, limits[1]), limits[2]), NA_real_)
  list(
    weight = weight,
    diagnostics = tibble::tibble(
      dataset = dataset,
      denominator_n = nrow(data),
      albumin_feature_complete_n = feature_n,
      primary_model_n = sum(selected),
      primary_model_percent = 100 * rate,
      weight_1st_percentile = limits[1],
      weight_99th_percentile = limits[2],
      weight_max = max(weight, na.rm = TRUE),
      effective_sample_size = sum(weight, na.rm = TRUE)^2 / sum(weight^2, na.rm = TRUE)
    )
  )
}

extract_cox <- function(fit, dataset, model) {
  beta <- stats::coef(fit)
  se <- sqrt(diag(stats::vcov(fit)))
  terms <- names(beta)
  keep <- grepl("^phenotype", terms)
  tibble::tibble(
    dataset = dataset,
    model = model,
    effect_measure = "HR",
    comparison = dplyr::recode(
      terms[keep], phenotypeP2 = "P2 vs P3", phenotypeP1 = "P1 vs P3"
    ),
    estimate = exp(beta[keep]),
    lower_95 = exp(beta[keep] - 1.96 * se[keep]),
    upper_95 = exp(beta[keep] + 1.96 * se[keep]),
    p_value = 2 * stats::pnorm(abs(beta[keep] / se[keep]), lower.tail = FALSE)
  ) |>
    dplyr::mutate(effect_95ci = sprintf("%.2f (%.2f-%.2f)", estimate, lower_95, upper_95))
}

extract_glmer <- function(fit, dataset, model) {
  coefficients <- summary(fit)$coefficients
  beta <- coefficients[, "Estimate"]
  se <- coefficients[, "Std. Error"]
  terms <- rownames(coefficients)
  keep <- grepl("^phenotype", terms)
  tibble::tibble(
    dataset = dataset,
    model = model,
    effect_measure = "OR",
    comparison = dplyr::recode(
      terms[keep], phenotypeP2 = "P2 vs P3", phenotypeP1 = "P1 vs P3"
    ),
    estimate = exp(beta[keep]),
    lower_95 = exp(beta[keep] - 1.96 * se[keep]),
    upper_95 = exp(beta[keep] + 1.96 * se[keep]),
    p_value = 2 * stats::pnorm(abs(beta[keep] / se[keep]), lower.tail = FALSE)
  ) |>
    dplyr::mutate(effect_95ci = sprintf("%.2f (%.2f-%.2f)", estimate, lower_95, upper_95))
}

# NHANES -----------------------------------------------------------------------
nhanes_results <- readRDS(file.path(
  root, "output", "nhanes_albumin_benchmarks", "NHANES_albumin_benchmark_results.rds"
))
nhanes <- nhanes_results$model_data |>
  dplyr::mutate(
    phenotype = factor(phenotype_albumin, levels = c("P3", "P2", "P1"))
  )
nhanes_design <- survey::svydesign(
  ids = ~SDMVPSU,
  strata = ~SDMVSTRA,
  weights = ~WTMEC8YR,
  nest = TRUE,
  data = nhanes
)
nhanes_fit <- survey::svycoxph(
  survival::Surv(PERMTH_INT, MORTSTAT) ~
    phenotype + RIDAGEYR + male + race + INDFMPIR + Comorbidity_Score_Extended +
    cycle + smoking + hypertension + diabetes,
  design = nhanes_design
)
nhanes_effects <- extract_cox(nhanes_fit, "NHANES", "Complex-survey albumin phenotype Cox") |>
  dplyr::mutate(n = nrow(nhanes), events = sum(nhanes$MORTSTAT == 1))
nhanes_counts <- nhanes |>
  dplyr::group_by(phenotype) |>
  dplyr::summarise(
    n = dplyr::n(),
    events = sum(MORTSTAT == 1),
    mortality_percent = 100 * mean(MORTSTAT == 1),
    .groups = "drop"
  ) |>
  dplyr::mutate(dataset = "NHANES", .before = 1)
# MIMIC-IV ---------------------------------------------------------------------
mimic_denominator <- readr::read_csv(
  file.path(
    root, "output", "mimic_24h_validation",
    "MIMIC_first_ICU_feature_availability_dataset.csv"
  ),
  show_col_types = FALSE
)
mimic_required <- c(
  "nlr", "sii", "haemoglobin", "albumin", "bmi", "creatinine",
  "anchor_age", "gender", "survival_days_365", "mortality_365d",
  "hospital_expire_flag", "intime", "dischtime"
)
mimic_albumin <- mimic_denominator |>
  dplyr::filter(dplyr::if_all(dplyr::all_of(mimic_required), ~ !is.na(.x))) |>
  dplyr::filter(
    is.finite(nlr), is.finite(sii), is.finite(haemoglobin), is.finite(albumin),
    is.finite(bmi), is.finite(creatinine), nlr > 0, sii > 0, albumin > 0,
    creatinine > 0
  )
set.seed(20260710)
mimic_cluster <- fit_albumin_clusters(
  mimic_albumin, "nlr", "sii", "haemoglobin", "albumin", "bmi", "creatinine"
)
mimic_albumin$phenotype <- mimic_cluster$phenotype

oasis <- readr::read_csv(
  file.path(
    root, "output", "mimic_official_oasis_v301",
    "MIMIC_official_OASIS_v301_scores.csv"
  ),
  show_col_types = FALSE
) |>
  dplyr::select(stay_id, oasis) |>
  dplyr::distinct(stay_id, .keep_all = TRUE)
mimic_model <- mimic_albumin |>
  dplyr::left_join(oasis, by = "stay_id") |>
  dplyr::filter(is.finite(oasis)) |>
  dplyr::mutate(
    male = as.integer(gender == "M"),
    oasis_quartile = factor(
      stable_ntile(oasis, 4, stay_id),
      levels = 1:4, labels = paste0("Q", 1:4)
    )
  )

mimic_denominator <- mimic_denominator |>
  dplyr::mutate(
    selected_feature = as.integer(stay_id %in% mimic_albumin$stay_id),
    selected_primary = as.integer(stay_id %in% mimic_model$stay_id),
    gender_selection = collapse_for_selection(gender, selected_primary),
    race_selection = collapse_for_selection(race, selected_primary),
    careunit_selection = collapse_for_selection(first_careunit, selected_primary),
    year_selection = collapse_for_selection(anchor_year_group, selected_primary)
  )
mimic_selection_fit <- stats::glm(
  selected_primary ~ splines::ns(anchor_age, df = 3) + gender_selection +
    race_selection + careunit_selection + year_selection,
  data = mimic_denominator,
  family = stats::binomial(),
  control = stats::glm.control(maxit = 100)
)
mimic_weight_object <- make_selection_weights(
  mimic_denominator, mimic_selection_fit, "selected_primary", "MIMIC-IV", nrow(mimic_albumin)
)
mimic_denominator$selection_ipw <- mimic_weight_object$weight
mimic_model <- mimic_model |>
  dplyr::left_join(
    mimic_denominator |>
      dplyr::select(stay_id, selection_ipw),
    by = "stay_id"
  ) |>
  dplyr::mutate(selection_ipw_scaled = selection_ipw / mean(selection_ipw))

mimic_fit <- survival::coxph(
  survival::Surv(survival_days_365, mortality_365d) ~
    phenotype + male + strata(oasis_quartile),
  data = mimic_model,
  ties = "efron"
)
mimic_ipw_fit <- survival::coxph(
  survival::Surv(survival_days_365, mortality_365d) ~
    phenotype + male + strata(oasis_quartile),
  data = mimic_model,
  weights = selection_ipw_scaled,
  robust = TRUE,
  cluster = stay_id,
  ties = "efron"
)
mimic_effects <- dplyr::bind_rows(
  extract_cox(mimic_fit, "MIMIC-IV", "Official OASIS stratified albumin phenotype Cox"),
  extract_cox(mimic_ipw_fit, "MIMIC-IV", "Selection-IPW official OASIS albumin phenotype Cox")
) |>
  dplyr::mutate(n = nrow(mimic_model), events = sum(mimic_model$mortality_365d == 1))
mimic_counts <- mimic_albumin |>
  dplyr::group_by(phenotype) |>
  dplyr::summarise(
    n = dplyr::n(), events = sum(mortality_365d == 1),
    mortality_percent = 100 * mean(mortality_365d == 1), .groups = "drop"
  ) |>
  dplyr::mutate(dataset = "MIMIC-IV", .before = 1)
mimic_profiles <- mimic_cluster$profiles |>
  dplyr::mutate(dataset = "MIMIC-IV", .before = 1)

# eICU -------------------------------------------------------------------------
eicu_denominator <- data.table::fread(
  file.path(
    root, "output", "eicu_24h_extraction",
    "eICU_first_ICU_feature_availability_dataset.csv"
  ),
  data.table = FALSE
) |>
  tibble::as_tibble()
eicu_required <- c(
  "nlr", "sii_like", "haemoglobin", "albumin", "bmi", "creatinine",
  "age_num", "gender_model", "hospitalid", "hospital_mortality", "icu_mortality"
)
eicu_albumin <- eicu_denominator |>
  dplyr::filter(dplyr::if_all(dplyr::all_of(eicu_required), ~ !is.na(.x))) |>
  dplyr::filter(
    is.finite(nlr), is.finite(sii_like), is.finite(haemoglobin), is.finite(albumin),
    is.finite(bmi), is.finite(creatinine), nlr > 0, nlr <= 100,
    sii_like > 0, sii_like <= 200000, albumin > 0, creatinine > 0
  )
set.seed(20260710)
eicu_cluster <- fit_albumin_clusters(
  eicu_albumin, "nlr", "sii_like", "haemoglobin", "albumin", "bmi", "creatinine"
)
eicu_albumin$phenotype <- eicu_cluster$phenotype
eicu_albumin$male <- as.integer(eicu_albumin$gender_model == "Male")

apache <- readr::read_csv(
  file.path(eicu_dir, "apachePatientResult.csv.gz"),
  col_select = c(patientunitstayid, apacheversion, apachescore),
  show_col_types = FALSE,
  progress = FALSE
) |>
  dplyr::filter(apacheversion == "IVa") |>
  dplyr::transmute(
    patientunitstayid,
    apachescore = as.numeric(apachescore)
  ) |>
  dplyr::mutate(apachescore = ifelse(apachescore < 0, NA_real_, apachescore)) |>
  dplyr::distinct(patientunitstayid, .keep_all = TRUE)
eicu_model <- eicu_albumin |>
  dplyr::left_join(apache, by = "patientunitstayid") |>
  dplyr::filter(is.finite(apachescore))

eicu_denominator <- eicu_denominator |>
  dplyr::mutate(
    selected_feature = as.integer(patientunitstayid %in% eicu_albumin$patientunitstayid),
    selected_primary = as.integer(patientunitstayid %in% eicu_model$patientunitstayid),
    gender_selection = collapse_for_selection(gender, selected_primary),
    ethnicity_selection = collapse_for_selection(ethnicity, selected_primary),
    unittype_selection = collapse_for_selection(unittype, selected_primary),
    hospital_key = ifelse(is.na(hospitalid), "Missing", as.character(hospitalid))
  )
eicu_overall_rate <- mean(eicu_denominator$selected_primary)
hospital_rates <- eicu_denominator |>
  dplyr::group_by(hospital_key) |>
  dplyr::summarise(n = dplyr::n(), selected = sum(selected_primary), .groups = "drop") |>
  dplyr::mutate(
    smoothed_rate = (selected + 20 * eicu_overall_rate) / (n + 20),
    hospital_selection_logit = stats::qlogis(pmin(pmax(smoothed_rate, 0.001), 0.999))
  )
eicu_denominator <- eicu_denominator |>
  dplyr::left_join(hospital_rates, by = "hospital_key")
eicu_selection_fit <- stats::glm(
  selected_primary ~ splines::ns(age_num, df = 3) + gender_selection +
    ethnicity_selection + unittype_selection + hospital_selection_logit,
  data = eicu_denominator,
  family = stats::binomial(),
  control = stats::glm.control(maxit = 100)
)
eicu_weight_object <- make_selection_weights(
  eicu_denominator, eicu_selection_fit, "selected_primary", "eICU", nrow(eicu_albumin)
)
eicu_denominator$selection_ipw <- eicu_weight_object$weight
eicu_model <- eicu_model |>
  dplyr::left_join(
    eicu_denominator |>
      dplyr::select(patientunitstayid, selection_ipw),
    by = "patientunitstayid"
  ) |>
  dplyr::mutate(selection_ipw_scaled = selection_ipw / mean(selection_ipw))

eicu_fit <- lme4::glmer(
  hospital_mortality ~ phenotype + age_num + male + apachescore + (1 | hospitalid),
  family = stats::binomial(), data = eicu_model, nAGQ = 1,
  control = lme4::glmerControl(optimizer = "bobyqa")
)
eicu_ipw_fit <- suppressWarnings(lme4::glmer(
  hospital_mortality ~ phenotype + age_num + male + apachescore + (1 | hospitalid),
  family = stats::binomial(), data = eicu_model, weights = selection_ipw_scaled,
  nAGQ = 1, control = lme4::glmerControl(optimizer = "bobyqa")
))
eicu_effects <- dplyr::bind_rows(
  extract_glmer(eicu_fit, "eICU", "APACHE-adjusted albumin phenotype mixed model"),
  extract_glmer(eicu_ipw_fit, "eICU", "Selection-IPW APACHE-adjusted albumin phenotype mixed model")
) |>
  dplyr::mutate(n = nrow(eicu_model), events = sum(eicu_model$hospital_mortality == 1))
eicu_counts <- eicu_albumin |>
  dplyr::group_by(phenotype) |>
  dplyr::summarise(
    n = dplyr::n(), events = sum(hospital_mortality == 1),
    mortality_percent = 100 * mean(hospital_mortality == 1), .groups = "drop"
  ) |>
  dplyr::mutate(dataset = "eICU", .before = 1)
eicu_profiles <- eicu_cluster$profiles |>
  dplyr::mutate(dataset = "eICU", .before = 1)

# Agreement with current frozen variants --------------------------------------
mimic_current <- readRDS(file.path(
  root, "output", "mimic_24h_robust_corrected_severity", "MIMIC_24h_robust_results.rds"
))$analysis |>
  dplyr::select(stay_id, current_phenotype = phenotype)
mimic_agreement_data <- mimic_albumin |>
  dplyr::select(stay_id, albumin_phenotype = phenotype) |>
  dplyr::inner_join(mimic_current, by = "stay_id")

eicu_current <- readRDS(file.path(
  root, "output", "eicu_24h_robust", "eICU_24h_robust_results.rds"
))$analysis |>
  dplyr::select(patientunitstayid, current_phenotype = phenotype)
eicu_agreement_data <- eicu_albumin |>
  dplyr::select(patientunitstayid, albumin_phenotype = phenotype) |>
  dplyr::inner_join(eicu_current, by = "patientunitstayid")

agreement <- tibble::tibble(
  dataset = c("NHANES", "MIMIC-IV", "eICU"),
  common_n = c(
    nhanes_results$agreement$common_n,
    nrow(mimic_agreement_data),
    nrow(eicu_agreement_data)
  ),
  exact_label_agreement = c(
    nhanes_results$agreement$exact_label_agreement,
    mean(mimic_agreement_data$albumin_phenotype == mimic_agreement_data$current_phenotype),
    mean(eicu_agreement_data$albumin_phenotype == eicu_agreement_data$current_phenotype)
  ),
  adjusted_rand_index = c(
    nhanes_results$agreement$adjusted_rand_index,
    adjusted_rand_index(
      mimic_agreement_data$albumin_phenotype, mimic_agreement_data$current_phenotype
    ),
    adjusted_rand_index(
      eicu_agreement_data$albumin_phenotype, eicu_agreement_data$current_phenotype
    )
  )
)

effects <- dplyr::bind_rows(nhanes_effects, mimic_effects, eicu_effects)
counts <- dplyr::bind_rows(nhanes_counts, mimic_counts, eicu_counts)
profiles <- dplyr::bind_rows(
  mimic_profiles,
  eicu_profiles
)
selection_diagnostics <- dplyr::bind_rows(
  mimic_weight_object$diagnostics,
  eicu_weight_object$diagnostics
)

readr::write_csv(counts, file.path(output_dir, "Table40A_harmonised_albumin_counts.csv"))
readr::write_csv(profiles, file.path(output_dir, "Table40B_harmonised_albumin_profiles_ICU.csv"))
readr::write_csv(effects, file.path(output_dir, "Table40C_harmonised_albumin_outcome_models.csv"))
readr::write_csv(agreement, file.path(output_dir, "Table40D_albumin_vs_current_agreement.csv"))
readr::write_csv(selection_diagnostics, file.path(output_dir, "Table40E_albumin_selection_weights.csv"))
saveRDS(
  list(
    effects = effects,
    counts = counts,
    profiles = profiles,
    agreement = agreement,
    selection_diagnostics = selection_diagnostics
  ),
  file.path(output_dir, "harmonised_albumin_results.rds")
)

summary_lines <- c(
  "Harmonised albumin cross-database sensitivity",
  "",
  "Cohort counts:",
  paste(capture.output(print(counts)), collapse = "\n"),
  "",
  "Outcome models:",
  paste(capture.output(print(effects)), collapse = "\n"),
  "",
  "Agreement with current variants:",
  paste(capture.output(print(agreement)), collapse = "\n"),
  "",
  "Selection diagnostics:",
  paste(capture.output(print(selection_diagnostics)), collapse = "\n"),
  "",
  "Interpretation: harmonised albumin results are a post-freeze sensitivity analysis and do not replace the total-protein primary phenotype."
)
writeLines(summary_lines, file.path(output_dir, "harmonised_albumin_summary.txt"))
cat(paste(summary_lines, collapse = "\n"), "\n")
