# ==============================================================================
# Amendment 17: eICU albumin-based multicentre robustness and transportability
# ==============================================================================

required_packages <- c("data.table", "dplyr", "readr", "tibble", "lme4", "metafor", "digest")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop("Missing required package(s): ", paste(missing_packages, collapse = ", "), call. = FALSE)
}

root <- normalizePath(Sys.getenv("PROJECT_ROOT", unset = getwd()), winslash = "/", mustWork = TRUE)
denominator_path <- Sys.getenv("EICU_DENOMINATOR_CSV", unset = file.path(root, "data", "derived", "eICU_first_ICU_feature_availability_dataset.csv"))
audited_apache_path <- Sys.getenv("EICU_APACHE_AUDIT_RDS", unset = file.path(root, "data", "derived", "eICU_APACHE_adjusted_sensitivity.rds"))
eicu_dir <- Sys.getenv("EICU_DIR", unset = "")
if (!nzchar(eicu_dir) && !file.exists(audited_apache_path)) stop("Set EICU_DIR or EICU_APACHE_AUDIT_RDS.", call. = FALSE)
raw_apache_path <- file.path(eicu_dir, "apachePatientResult.csv.gz")
output_dir <- file.path(root, "output", "eicu_albumin_amendment17_2026-08-30")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(denominator_path)) stop("Missing frozen eICU denominator.", call. = FALSE)
if (!file.exists(raw_apache_path) && !file.exists(audited_apache_path)) {
  stop("Neither the raw nor audited APACHE source is available.", call. = FALSE)
}

winsorise <- function(x) {
  limits <- stats::quantile(x, c(0.01, 0.99), na.rm = TRUE, names = FALSE)
  pmin(pmax(as.numeric(x), limits[1]), limits[2])
}

safe_z <- function(x) {
  current_sd <- stats::sd(x, na.rm = TRUE)
  if (!is.finite(current_sd) || current_sd == 0) return(rep(0, length(x)))
  (x - mean(x, na.rm = TRUE)) / current_sd
}

collapse_for_selection <- function(x, selected, min_total = 100L, min_selected = 10L) {
  value <- as.character(x)
  value[is.na(value) | trimws(value) == ""] <- "Missing"
  total_table <- table(value)
  selected_table <- table(value[selected == 1])
  selected_count <- as.numeric(selected_table[names(total_table)])
  selected_count[is.na(selected_count)] <- 0
  excluded_count <- as.numeric(total_table) - selected_count
  keep <- names(total_table)[
    as.numeric(total_table) >= min_total & selected_count >= min_selected & excluded_count >= min_selected
  ]
  factor(ifelse(value %in% keep, value, "Other"))
}

fit_albumin_clusters <- function(data) {
  z <- data |>
    dplyr::transmute(
      log_nlr = log(pmax(winsorise(nlr), .Machine$double.eps)),
      log_sii = log(pmax(winsorise(sii_like), .Machine$double.eps)),
      haemoglobin = winsorise(haemoglobin),
      albumin = winsorise(albumin),
      bmi = winsorise(bmi),
      log_creatinine = log(pmax(winsorise(creatinine), .Machine$double.eps))
    ) |>
    as.data.frame() |>
    scale()

  set.seed(20260710)
  fit <- stats::kmeans(z, centers = 3, nstart = 100, iter.max = 500, algorithm = "Lloyd")
  profiles <- data |>
    dplyr::mutate(cluster_raw = fit$cluster) |>
    dplyr::group_by(cluster_raw) |>
    dplyr::summarise(
      n = dplyr::n(),
      nlr = stats::median(nlr),
      sii_like = stats::median(sii_like),
      haemoglobin = stats::median(haemoglobin),
      albumin = stats::median(albumin),
      bmi = stats::median(bmi),
      creatinine = stats::median(creatinine),
      .groups = "drop"
    )

  p1_score <- safe_z(log(profiles$nlr)) + safe_z(log(profiles$sii_like)) +
    safe_z(log(profiles$creatinine))
  p1_raw <- profiles$cluster_raw[which.max(p1_score)]
  remaining <- profiles |> dplyr::filter(cluster_raw != p1_raw)
  p2_score <- -safe_z(remaining$haemoglobin) - safe_z(remaining$albumin) - safe_z(remaining$bmi)
  p2_raw <- remaining$cluster_raw[which.max(p2_score)]
  p3_raw <- setdiff(profiles$cluster_raw, c(p1_raw, p2_raw))
  mapping <- tibble::tibble(
    cluster_raw = c(p1_raw, p2_raw, p3_raw),
    phenotype = c("P1", "P2", "P3")
  )
  phenotype <- dplyr::case_when(
    fit$cluster == p1_raw ~ "P1",
    fit$cluster == p2_raw ~ "P2",
    fit$cluster == p3_raw ~ "P3",
    TRUE ~ NA_character_
  )
  list(
    phenotype = factor(phenotype, levels = c("P3", "P2", "P1")),
    profiles = profiles |>
      dplyr::left_join(mapping, by = "cluster_raw") |>
      dplyr::arrange(factor(phenotype, levels = c("P1", "P2", "P3")))
  )
}

capture_glmer <- function(formula, data, weights = NULL, optimizer = "bobyqa") {
  warnings <- character()
  fit <- tryCatch(
    withCallingHandlers(
      lme4::glmer(
        formula, data = data, weights = weights, family = stats::binomial(), nAGQ = 1,
        control = lme4::glmerControl(
          optimizer = optimizer, optCtrl = list(maxfun = 300000), calc.derivs = TRUE
        )
      ),
      warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) e
  )
  list(fit = fit, warnings = unique(warnings))
}

convergence_messages <- function(fit) {
  messages <- fit@optinfo$conv$lme4$messages
  if (is.null(messages)) character() else as.character(messages)
}

extract_phenotype_effects <- function(captured, model, data, outcome) {
  if (inherits(captured$fit, "error")) {
    return(tibble::tibble(
      model = model, outcome = outcome, n = nrow(data), events = sum(data[[outcome]] == 1),
      hospitals = dplyr::n_distinct(data$hospitalid), comparison = c("P1 vs P3", "P2 vs P3"),
      estimate = NA_real_, lower_95 = NA_real_, upper_95 = NA_real_, p_value = NA_real_,
      effect_95ci = NA_character_
    ))
  }
  coefficients <- summary(captured$fit)$coefficients
  terms <- if (all(c("phenotypeP1", "phenotypeP2") %in% rownames(coefficients))) {
    c("phenotypeP1", "phenotypeP2")
  } else if (all(c("p1", "p2") %in% rownames(coefficients))) {
    c("p1", "p2")
  } else {
    stop("Neither phenotype nor P1/P2 indicator terms are available in the fitted model.", call. = FALSE)
  }
  beta <- coefficients[terms, "Estimate"]
  se <- coefficients[terms, "Std. Error"]
  tibble::tibble(
    model = model,
    outcome = outcome,
    n = nrow(data),
    events = sum(data[[outcome]] == 1),
    hospitals = dplyr::n_distinct(data$hospitalid),
    comparison = c("P1 vs P3", "P2 vs P3"),
    estimate = exp(beta),
    lower_95 = exp(beta - 1.96 * se),
    upper_95 = exp(beta + 1.96 * se),
    p_value = 2 * stats::pnorm(abs(beta / se), lower.tail = FALSE),
    effect_95ci = sprintf("%.3f (%.3f-%.3f)", exp(beta), exp(beta - 1.96 * se), exp(beta + 1.96 * se))
  )
}

extract_p1_slope_sd <- function(fit) {
  vc <- as.data.frame(lme4::VarCorr(fit))
  diagonal <- vc[is.na(vc$var2), , drop = FALSE]
  direct <- diagonal[diagonal$var1 == "p1", , drop = FALSE]
  if (nrow(direct) >= 1L) return(as.numeric(direct$sdcor[1]))
  secondary <- diagonal[
    grepl("^hospitalid\\.", diagonal$grp) & diagonal$var1 == "(Intercept)", , drop = FALSE
  ]
  if (nrow(secondary) == 1L) return(as.numeric(secondary$sdcor[1]))
  NA_real_
}

extract_model_diagnostics <- function(captured, model, data, outcome, random_slope = FALSE) {
  if (inherits(captured$fit, "error")) {
    return(tibble::tibble(
      model = model, n = nrow(data), events = sum(data[[outcome]] == 1),
      hospitals = dplyr::n_distinct(data$hospitalid), fit_ok = FALSE,
      optimizer_code = NA_integer_, convergence_message = conditionMessage(captured$fit),
      captured_warnings = paste(captured$warnings, collapse = " | "), singular = NA,
      random_intercept_sd = NA_real_, p1_random_slope_sd = NA_real_, aic = NA_real_, bic = NA_real_
    ))
  }
  fit <- captured$fit
  opt_code <- fit@optinfo$conv$opt
  if (is.null(opt_code)) opt_code <- 0L
  vc <- as.data.frame(lme4::VarCorr(fit))
  intercept_sd <- vc$sdcor[vc$var1 == "(Intercept)" & is.na(vc$var2)][1]
  tibble::tibble(
    model = model, n = nrow(data), events = sum(data[[outcome]] == 1),
    hospitals = dplyr::n_distinct(data$hospitalid), fit_ok = TRUE,
    optimizer_code = as.integer(opt_code),
    convergence_message = paste(convergence_messages(fit), collapse = " | "),
    captured_warnings = paste(captured$warnings, collapse = " | "),
    singular = lme4::isSingular(fit, tol = 1e-5),
    random_intercept_sd = as.numeric(intercept_sd),
    p1_random_slope_sd = if (random_slope) extract_p1_slope_sd(fit) else NA_real_,
    aic = stats::AIC(fit), bic = stats::BIC(fit)
  )
}

sha256_file <- function(path) {
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}

cat("Loading the frozen eICU candidate denominator...\n")
denominator <- data.table::fread(denominator_path, data.table = FALSE, showProgress = FALSE) |>
  tibble::as_tibble()
if (nrow(denominator) != 68798L || anyDuplicated(denominator$patientunitstayid)) {
  stop("Frozen denominator dimensions or patient identifiers are unexpected.", call. = FALSE)
}

feature_names <- c("nlr", "sii_like", "haemoglobin", "albumin", "bmi", "creatinine")
required_names <- c(
  "patientunitstayid", "hospitalid", "age_num", "gender_model", "gender", "ethnicity",
  "unittype", "hospital_mortality", "icu_mortality", "survival_days_hosp",
  "unitdischargeoffset", feature_names
)
missing_names <- setdiff(required_names, names(denominator))
if (length(missing_names) > 0L) {
  stop("Missing required column(s): ", paste(missing_names, collapse = ", "), call. = FALSE)
}

feature_complete <- denominator |>
  dplyr::filter(dplyr::if_all(dplyr::all_of(feature_names), ~ !is.na(.x) & is.finite(.x))) |>
  dplyr::filter(
    nlr > 0, nlr <= 100, sii_like > 0, sii_like <= 200000,
    albumin > 0, creatinine > 0
  ) |>
  dplyr::mutate(male = as.integer(gender_model == "Male"))
if (nrow(feature_complete) != 15242L) {
  stop("Albumin six-feature cohort does not contain the expected 15,242 stays.", call. = FALSE)
}

cat("Clustering the outcome-blind six-feature cohort...\n")
cluster <- fit_albumin_clusters(feature_complete)
feature_complete$phenotype <- cluster$phenotype
feature_complete <- feature_complete |>
  dplyr::mutate(
    p1 = as.integer(phenotype == "P1"),
    p2 = as.integer(phenotype == "P2"),
    icu_los_days = as.numeric(unitdischargeoffset) / 1440
  )

if (file.exists(raw_apache_path)) {
  apache <- readr::read_csv(
    raw_apache_path,
    col_select = c(patientunitstayid, apacheversion, apachescore),
    show_col_types = FALSE, progress = FALSE
  ) |>
    dplyr::filter(apacheversion == "IVa") |>
    dplyr::transmute(patientunitstayid, apachescore = as.numeric(apachescore)) |>
    dplyr::mutate(apachescore = ifelse(apachescore < 0, NA_real_, apachescore)) |>
    dplyr::distinct(patientunitstayid, .keep_all = TRUE)
  apache_source_used <- raw_apache_path
} else {
  apache <- readRDS(audited_apache_path)$model_data |>
    dplyr::select(patientunitstayid, apachescore) |>
    dplyr::distinct(patientunitstayid, .keep_all = TRUE)
  apache_source_used <- audited_apache_path
}

apache_audit <- tibble::tibble(
  source = apache_source_used,
  rows = nrow(apache),
  unique_patientunitstayid = dplyr::n_distinct(apache$patientunitstayid),
  invalid_negative_scores_after_cleaning = sum(apache$apachescore < 0, na.rm = TRUE),
  raw_vs_audited_common_n = NA_integer_,
  raw_vs_audited_exact_score_agreement = NA_real_
)
if (file.exists(raw_apache_path) && file.exists(audited_apache_path)) {
  audited_apache <- readRDS(audited_apache_path)$model_data |>
    dplyr::select(patientunitstayid, audited_apachescore = apachescore) |>
    dplyr::distinct(patientunitstayid, .keep_all = TRUE)
  comparison <- apache |>
    dplyr::inner_join(audited_apache, by = "patientunitstayid")
  apache_audit$raw_vs_audited_common_n <- nrow(comparison)
  apache_audit$raw_vs_audited_exact_score_agreement <- mean(
    comparison$apachescore == comparison$audited_apachescore,
    na.rm = TRUE
  )
}

analysis <- feature_complete |>
  dplyr::left_join(apache, by = "patientunitstayid")
hospital_data <- analysis |>
  dplyr::filter(is.finite(apachescore), !is.na(hospital_mortality)) |>
  dplyr::mutate(
    hospitalid = factor(hospitalid),
    phenotype = droplevels(phenotype),
    age_z = as.numeric(scale(age_num)),
    apache_z = as.numeric(scale(apachescore))
  )
icu_data <- analysis |>
  dplyr::filter(is.finite(apachescore), !is.na(icu_mortality)) |>
  dplyr::mutate(hospitalid = factor(hospitalid), phenotype = droplevels(phenotype))

denominator <- denominator |>
  dplyr::mutate(
    selected_feature = as.integer(patientunitstayid %in% feature_complete$patientunitstayid),
    selected_primary = as.integer(patientunitstayid %in% hospital_data$patientunitstayid),
    apache_available = as.integer(patientunitstayid %in% apache$patientunitstayid),
    gender_selection = collapse_for_selection(gender, selected_primary),
    ethnicity_selection = collapse_for_selection(ethnicity, selected_primary),
    unittype_selection = collapse_for_selection(unittype, selected_primary),
    hospital_key = ifelse(is.na(hospitalid), "Missing", as.character(hospitalid))
  )
overall_selection_rate <- mean(denominator$selected_primary)
hospital_rates <- denominator |>
  dplyr::group_by(hospital_key) |>
  dplyr::summarise(hospital_n = dplyr::n(), selected_n = sum(selected_primary), .groups = "drop") |>
  dplyr::mutate(
    smoothed_selection_rate = (selected_n + 20 * overall_selection_rate) / (hospital_n + 20),
    hospital_selection_logit = stats::qlogis(pmin(pmax(smoothed_selection_rate, 0.001), 0.999))
  )
denominator <- denominator |>
  dplyr::left_join(hospital_rates, by = "hospital_key")

cat("Fitting the selection model and constructing locked IPW...\n")
selection_fit <- stats::glm(
  selected_primary ~ splines::ns(age_num, df = 3) + gender_selection +
    ethnicity_selection + unittype_selection + hospital_selection_logit,
  data = denominator, family = stats::binomial(), control = stats::glm.control(maxit = 100)
)
if (!isTRUE(selection_fit$converged)) stop("eICU selection model did not converge.", call. = FALSE)
denominator$selection_probability <- pmin(
  pmax(stats::predict(selection_fit, type = "response"), 0.001), 0.999
)
raw_weight <- ifelse(
  denominator$selected_primary == 1,
  overall_selection_rate / denominator$selection_probability,
  NA_real_
)
weight_limits <- stats::quantile(
  raw_weight[denominator$selected_primary == 1], c(0.01, 0.99), na.rm = TRUE, names = FALSE
)
denominator$selection_ipw <- ifelse(
  denominator$selected_primary == 1,
  pmin(pmax(raw_weight, weight_limits[1]), weight_limits[2]),
  NA_real_
)
hospital_data <- hospital_data |>
  dplyr::left_join(
    denominator |> dplyr::select(patientunitstayid, selection_probability, selection_ipw),
    by = "patientunitstayid"
  ) |>
  dplyr::mutate(selection_ipw_scaled = selection_ipw / mean(selection_ipw))

selection_diagnostics <- tibble::tibble(
  denominator_n = nrow(denominator),
  albumin_feature_complete_n = nrow(feature_complete),
  primary_model_n = nrow(hospital_data),
  primary_model_percent = 100 * overall_selection_rate,
  probability_min = min(denominator$selection_probability),
  probability_1st_percentile = stats::quantile(denominator$selection_probability, 0.01),
  probability_median = stats::median(denominator$selection_probability),
  probability_99th_percentile = stats::quantile(denominator$selection_probability, 0.99),
  probability_max = max(denominator$selection_probability),
  weight_1st_percentile = weight_limits[1],
  weight_99th_percentile = weight_limits[2],
  weight_min = min(hospital_data$selection_ipw),
  weight_max = max(hospital_data$selection_ipw),
  weight_mean_before_scaling = mean(hospital_data$selection_ipw),
  weight_mean_after_scaling = mean(hospital_data$selection_ipw_scaled),
  effective_sample_size = sum(hospital_data$selection_ipw)^2 / sum(hospital_data$selection_ipw^2),
  selection_model_converged = selection_fit$converged
)

feature_availability <- tibble::tibble(
  feature = feature_names,
  denominator_n = nrow(denominator),
  available_n = vapply(feature_names, function(v) sum(!is.na(denominator[[v]]) & is.finite(denominator[[v]])), integer(1)),
  available_percent = 100 * available_n / denominator_n
)

summarise_population <- function(data, label) {
  tibble::tibble(
    population = label,
    n = nrow(data),
    hospitals = dplyr::n_distinct(data$hospitalid),
    age_mean = mean(data$age_num),
    age_sd = stats::sd(data$age_num),
    age_median = stats::median(data$age_num),
    male_percent = 100 * mean(data$gender_model == "Male"),
    hospital_mortality_available_n = sum(!is.na(data$hospital_mortality)),
    hospital_mortality_percent_among_available = 100 * mean(data$hospital_mortality == 1, na.rm = TRUE),
    icu_mortality_available_n = sum(!is.na(data$icu_mortality)),
    icu_mortality_percent_among_available = 100 * mean(data$icu_mortality == 1, na.rm = TRUE)
  )
}
population_comparison <- dplyr::bind_rows(
  summarise_population(denominator, "Candidate denominator"),
  summarise_population(
    denominator |> dplyr::filter(selected_primary == 1),
    "Selected APACHE-complete hospital-mortality cohort"
  )
)

phenotype_counts <- feature_complete |>
  dplyr::group_by(phenotype) |>
  dplyr::summarise(
    n = dplyr::n(), percent = 100 * n / nrow(feature_complete),
    hospital_outcome_available_n = sum(!is.na(hospital_mortality)),
    hospital_deaths = sum(hospital_mortality == 1, na.rm = TRUE),
    hospital_mortality_percent = 100 * mean(hospital_mortality == 1, na.rm = TRUE),
    icu_outcome_available_n = sum(!is.na(icu_mortality)),
    icu_deaths = sum(icu_mortality == 1, na.rm = TRUE),
    icu_mortality_percent = 100 * mean(icu_mortality == 1, na.rm = TRUE),
    .groups = "drop"
  )

cat("Fitting primary, IPW, secondary, and 24-hour landmark models...\n")
formula_hospital <- hospital_mortality ~ phenotype + age_num + male + apachescore + (1 | hospitalid)
formula_icu <- icu_mortality ~ phenotype + age_num + male + apachescore + (1 | hospitalid)
fit_primary <- capture_glmer(formula_hospital, hospital_data)
fit_ipw <- capture_glmer(formula_hospital, hospital_data, hospital_data$selection_ipw_scaled)
fit_icu <- capture_glmer(formula_icu, icu_data)
hospital_landmark <- hospital_data |> dplyr::filter(survival_days_hosp > 1) |> droplevels()
icu_landmark <- icu_data |> dplyr::filter(icu_los_days > 1) |> droplevels()
fit_hospital_landmark <- capture_glmer(formula_hospital, hospital_landmark)
fit_icu_landmark <- capture_glmer(formula_icu, icu_landmark)

outcome_models <- dplyr::bind_rows(
  extract_phenotype_effects(fit_primary, "Primary APACHE-adjusted hospital random-intercept", hospital_data, "hospital_mortality"),
  extract_phenotype_effects(fit_ipw, "Selection-IPW APACHE-adjusted hospital random-intercept", hospital_data, "hospital_mortality"),
  extract_phenotype_effects(fit_icu, "APACHE-adjusted ICU random-intercept", icu_data, "icu_mortality"),
  extract_phenotype_effects(fit_hospital_landmark, "Hospital mortality, 24-hour landmark, APACHE-adjusted", hospital_landmark, "hospital_mortality"),
  extract_phenotype_effects(fit_icu_landmark, "ICU mortality, 24-hour landmark, APACHE-adjusted", icu_landmark, "icu_mortality")
)
model_diagnostics <- dplyr::bind_rows(
  extract_model_diagnostics(fit_primary, "Primary hospital mortality", hospital_data, "hospital_mortality"),
  extract_model_diagnostics(fit_ipw, "Selection-IPW hospital mortality", hospital_data, "hospital_mortality"),
  extract_model_diagnostics(fit_icu, "ICU mortality", icu_data, "icu_mortality"),
  extract_model_diagnostics(fit_hospital_landmark, "Hospital mortality 24-hour landmark", hospital_landmark, "hospital_mortality"),
  extract_model_diagnostics(fit_icu_landmark, "ICU mortality 24-hour landmark", icu_landmark, "icu_mortality")
)

cat("Auditing hospital information and fitting hospital-specific effects...\n")
hospital_info <- hospital_data |>
  dplyr::group_by(hospitalid) |>
  dplyr::summarise(
    total_n = dplyr::n(), p1_n = sum(phenotype == "P1"), p2_n = sum(phenotype == "P2"),
    p3_n = sum(phenotype == "P3"), p1_p3_n = sum(phenotype %in% c("P1", "P3")),
    p1_p3_deaths = sum(hospital_mortality == 1 & phenotype %in% c("P1", "P3")),
    p1_p3_survivors = sum(hospital_mortality == 0 & phenotype %in% c("P1", "P3")),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    meets_total_n = total_n >= 100, meets_p1_n = p1_n >= 20, meets_p3_n = p3_n >= 20,
    meets_deaths = p1_p3_deaths >= 20, meets_survivors = p1_p3_survivors >= 20,
    high_information = meets_total_n & meets_p1_n & meets_p3_n & meets_deaths & meets_survivors,
    exclusion_reason = dplyr::case_when(
      high_information ~ "Included", !meets_total_n ~ "Total n < 100", !meets_p1_n ~ "P1 n < 20",
      !meets_p3_n ~ "P3 n < 20", !meets_deaths ~ "P1/P3 deaths < 20",
      !meets_survivors ~ "P1/P3 survivors < 20", TRUE ~ "Other"
    )
  ) |>
  dplyr::arrange(dplyr::desc(high_information), dplyr::desc(total_n), hospitalid)
eligible_hospitals <- hospital_info |> dplyr::filter(high_information) |> dplyr::pull(hospitalid)
if (length(eligible_hospitals) < 10L) {
  stop("Fewer than 10 hospitals meet the locked information criteria.", call. = FALSE)
}

hospital_effects <- dplyr::bind_rows(lapply(eligible_hospitals, function(hospital) {
  current <- hospital_data |>
    dplyr::filter(hospitalid == hospital, phenotype %in% c("P1", "P3")) |>
    dplyr::mutate(p1_local = as.integer(phenotype == "P1"))
  fit <- tryCatch(
    suppressWarnings(stats::glm(
      hospital_mortality ~ p1_local + age_num + male + apachescore,
      data = current, family = stats::binomial(), control = stats::glm.control(maxit = 100)
    )),
    error = function(e) e
  )
  if (inherits(fit, "error")) {
    return(tibble::tibble(
      hospitalid = as.character(hospital), model_n = nrow(current), deaths = sum(current$hospital_mortality),
      converged = FALSE, estimable = FALSE, log_or = NA_real_, standard_error = NA_real_,
      model_note = conditionMessage(fit)
    ))
  }
  coefficients <- summary(fit)$coefficients
  estimable <- "p1_local" %in% rownames(coefficients) &&
    all(is.finite(coefficients["p1_local", c("Estimate", "Std. Error")])) &&
    coefficients["p1_local", "Std. Error"] > 0
  tibble::tibble(
    hospitalid = as.character(hospital), model_n = nrow(current), deaths = sum(current$hospital_mortality),
    converged = isTRUE(fit$converged), estimable = estimable,
    log_or = if (estimable) coefficients["p1_local", "Estimate"] else NA_real_,
    standard_error = if (estimable) coefficients["p1_local", "Std. Error"] else NA_real_,
    model_note = if (isTRUE(fit$converged) && estimable) "Included" else "Non-converged or non-estimable"
  )
})) |>
  dplyr::left_join(
    hospital_info |>
      dplyr::mutate(hospitalid = as.character(hospitalid)) |>
      dplyr::select(hospitalid, total_n, p1_n, p3_n, p1_p3_deaths, p1_p3_survivors),
    by = "hospitalid"
  ) |>
  dplyr::mutate(
    estimate = exp(log_or), lower_95 = exp(log_or - 1.96 * standard_error),
    upper_95 = exp(log_or + 1.96 * standard_error),
    effect_95ci = ifelse(estimable, sprintf("%.3f (%.3f-%.3f)", estimate, lower_95, upper_95), NA_character_)
  ) |>
  dplyr::arrange(as.numeric(hospitalid))

meta_data <- hospital_effects |>
  dplyr::filter(converged, estimable, is.finite(log_or), is.finite(standard_error))
if (nrow(meta_data) < 10L) stop("Fewer than 10 hospital effects are estimable.", call. = FALSE)
meta_fit <- metafor::rma.uni(
  yi = log_or, sei = standard_error, data = meta_data, method = "REML", test = "knha", slab = hospitalid
)
meta_prediction <- predict(meta_fit)
meta_summary <- tibble::tibble(
  model = "Hospital-specific APACHE-adjusted P1 vs P3 random-effects meta-analysis",
  hospitals = meta_fit$k, participants_p1_p3 = sum(meta_data$model_n), deaths_p1_p3 = sum(meta_data$deaths),
  pooled_or = exp(as.numeric(meta_fit$b)), lower_95 = exp(meta_fit$ci.lb), upper_95 = exp(meta_fit$ci.ub),
  prediction_lower_95 = exp(meta_prediction$pi.lb), prediction_upper_95 = exp(meta_prediction$pi.ub),
  p_value = meta_fit$pval, q_statistic = meta_fit$QE, q_df = meta_fit$k - 1,
  q_p_value = meta_fit$QEp, tau_squared = meta_fit$tau2, i_squared_percent = meta_fit$I2,
  effect_95ci = sprintf("%.3f (%.3f-%.3f)", exp(as.numeric(meta_fit$b)), exp(meta_fit$ci.lb), exp(meta_fit$ci.ub)),
  prediction_interval = sprintf("%.3f-%.3f", exp(meta_prediction$pi.lb), exp(meta_prediction$pi.ub))
)

leave_one_out <- dplyr::bind_rows(lapply(meta_data$hospitalid, function(omitted_hospital) {
  current <- meta_data |> dplyr::filter(hospitalid != omitted_hospital)
  fit <- metafor::rma.uni(yi = log_or, sei = standard_error, data = current, method = "REML", test = "knha")
  prediction <- predict(fit)
  tibble::tibble(
    omitted_hospitalid = omitted_hospital, hospitals_remaining = fit$k,
    pooled_or = exp(as.numeric(fit$b)), lower_95 = exp(fit$ci.lb), upper_95 = exp(fit$ci.ub),
    prediction_lower_95 = exp(prediction$pi.lb), prediction_upper_95 = exp(prediction$pi.ub),
    p_value = fit$pval, i_squared_percent = fit$I2, tau_squared = fit$tau2
  )
}))
leave_one_out_summary <- tibble::tibble(
  omitted_models = nrow(leave_one_out), pooled_or_min = min(leave_one_out$pooled_or),
  pooled_or_max = max(leave_one_out$pooled_or), lower_95_min = min(leave_one_out$lower_95),
  lower_95_max = max(leave_one_out$lower_95), upper_95_min = min(leave_one_out$upper_95),
  upper_95_max = max(leave_one_out$upper_95), models_with_ci_above_one = sum(leave_one_out$lower_95 > 1),
  models_with_ci_crossing_one = sum(leave_one_out$lower_95 <= 1 & leave_one_out$upper_95 >= 1),
  i_squared_min = min(leave_one_out$i_squared_percent), i_squared_max = max(leave_one_out$i_squared_percent)
)
direction_summary <- tibble::tibble(
  estimable_hospitals = nrow(meta_data), hospitals_or_above_one = sum(meta_data$estimate > 1),
  hospitals_or_below_one = sum(meta_data$estimate < 1), hospitals_ci_above_one = sum(meta_data$lower_95 > 1),
  hospitals_ci_crossing_one = sum(meta_data$lower_95 <= 1 & meta_data$upper_95 >= 1),
  hospitals_ci_below_one = sum(meta_data$upper_95 < 1)
)

cat("Fitting all-hospital and high-information random-effect diagnostics...\n")
formula_ri <- hospital_mortality ~ p1 + p2 + age_z + male + apache_z + (1 | hospitalid)
formula_rs <- hospital_mortality ~ p1 + p2 + age_z + male + apache_z + (1 + p1 || hospitalid)
high_information_data <- hospital_data |>
  dplyr::filter(hospitalid %in% eligible_hospitals) |>
  droplevels()
fit_ri_all <- capture_glmer(formula_ri, hospital_data)
fit_ri_high <- capture_glmer(formula_ri, high_information_data)
fit_rs_all <- capture_glmer(formula_rs, hospital_data)

random_effect_diagnostics <- dplyr::bind_rows(
  extract_model_diagnostics(fit_ri_all, "All hospitals: random intercept", hospital_data, "hospital_mortality"),
  extract_model_diagnostics(fit_ri_high, "High-information hospitals: random intercept", high_information_data, "hospital_mortality"),
  extract_model_diagnostics(fit_rs_all, "All hospitals: uncorrelated P1 random slope", hospital_data, "hospital_mortality", TRUE)
)
random_effect_models <- dplyr::bind_rows(
  extract_phenotype_effects(fit_ri_all, "All hospitals: random intercept", hospital_data, "hospital_mortality"),
  extract_phenotype_effects(fit_ri_high, "High-information hospitals: random intercept", high_information_data, "hospital_mortality"),
  extract_phenotype_effects(fit_rs_all, "All hospitals: uncorrelated P1 random slope", hospital_data, "hospital_mortality")
)
if (!inherits(fit_ri_all$fit, "error") && !inherits(fit_rs_all$fit, "error")) {
  reduced_ll <- stats::logLik(fit_ri_all$fit)
  expanded_ll <- stats::logLik(fit_rs_all$fit)
  lrt_df <- attr(expanded_ll, "df") - attr(reduced_ll, "df")
  lrt_chisq <- max(0, 2 * (as.numeric(expanded_ll) - as.numeric(reduced_ll)))
  random_slope_lrt <- tibble::tibble(
    comparison = "All hospitals: random intercept vs uncorrelated P1 random slope",
    chisq = lrt_chisq, df = lrt_df,
    p_value_descriptive = stats::pchisq(lrt_chisq, df = lrt_df, lower.tail = FALSE),
    caution = "Descriptive only: the null variance is on the parameter boundary; singularity is not evidence of no heterogeneity."
  )
} else {
  random_slope_lrt <- tibble::tibble(
    comparison = "All hospitals: random intercept vs uncorrelated P1 random slope",
    chisq = NA_real_, df = NA_integer_, p_value_descriptive = NA_real_,
    caution = "Model failure prevented the descriptive comparison; singularity is not evidence of no heterogeneity."
  )
}

tables <- list(
  Table82A_albumin_cluster_profiles = cluster$profiles,
  Table82B_phenotype_counts = phenotype_counts,
  Table82C_outcome_models = outcome_models,
  Table82D_model_diagnostics = model_diagnostics,
  Table82E_selection_diagnostics = selection_diagnostics,
  Table82F_candidate_vs_selected = population_comparison,
  Table82G_feature_availability = feature_availability,
  Table82H_hospital_information = hospital_info,
  Table82I_hospital_specific_P1_effects = hospital_effects,
  Table82J_REML_HK_meta_analysis = meta_summary,
  Table82K_leave_one_hospital_out = leave_one_out,
  Table82L_leave_one_out_summary = leave_one_out_summary,
  Table82M_hospital_direction_summary = direction_summary,
  Table82N_random_effect_models = random_effect_models,
  Table82O_random_effect_diagnostics = random_effect_diagnostics,
  Table82P_random_slope_LRT = random_slope_lrt,
  Table82Q_APACHE_source_audit = apache_audit
)
for (name in names(tables)) {
  readr::write_csv(tables[[name]], file.path(output_dir, paste0(name, ".csv")))
}

saveRDS(
  tables,
  file.path(output_dir, "eICU_albumin_amendment17_summary_results.rds")
)

source_paths <- c(denominator_path, apache_source_used, file.path(root, "scripts", "03_eicu_albumin_analysis.R"))
if (file.exists(audited_apache_path) && normalizePath(audited_apache_path, winslash = "/") !=
    normalizePath(apache_source_used, winslash = "/")) {
  source_paths <- c(source_paths, audited_apache_path)
}
source_manifest <- tibble::tibble(
  file = basename(source_paths),
  bytes = file.info(source_paths)$size,
  sha256 = vapply(source_paths, sha256_file, character(1))
)
readr::write_csv(source_manifest, file.path(output_dir, "SOURCE_SHA256_MANIFEST.csv"))

summary_lines <- c(
  "Amendment 17 eICU albumin-based multicentre robustness analysis",
  "",
  paste0("Candidate denominator: n=", nrow(denominator), "."),
  paste0("Outcome-blind six-feature cohort: n=", nrow(feature_complete), "."),
  paste0("APACHE-complete hospital-mortality cohort: n=", nrow(hospital_data), ", deaths=", sum(hospital_data$hospital_mortality), "."),
  paste0("Hospitals in primary model: ", dplyr::n_distinct(hospital_data$hospitalid), "."),
  paste0("High-information hospitals: ", length(eligible_hospitals), "."),
  "",
  "Phenotype counts:",
  paste(capture.output(print(phenotype_counts, width = Inf)), collapse = "\n"),
  "",
  "Outcome models:",
  paste(capture.output(print(outcome_models, width = Inf)), collapse = "\n"),
  "",
  "Selection diagnostics:",
  paste(capture.output(print(selection_diagnostics, width = Inf)), collapse = "\n"),
  "",
  "Hospital meta-analysis:",
  paste(capture.output(print(meta_summary, width = Inf)), collapse = "\n"),
  "",
  "Random-effect diagnostics:",
  paste(capture.output(print(random_effect_diagnostics, width = Inf)), collapse = "\n"),
  "",
  "Interpretation boundary: outcome data were not used for feature completeness, clustering, K selection, or phenotype labelling. Random-slope boundary or singular solutions are not evidence of absent heterogeneity."
)
writeLines(summary_lines, file.path(output_dir, "eICU_albumin_amendment17_summary.txt"))

output_files <- list.files(output_dir, full.names = TRUE, recursive = FALSE)
output_files <- output_files[basename(output_files) != "OUTPUT_SHA256_MANIFEST.csv"]
output_manifest <- tibble::tibble(
  file = basename(output_files), bytes = file.info(output_files)$size,
  sha256 = vapply(output_files, sha256_file, character(1))
)
readr::write_csv(output_manifest, file.path(output_dir, "OUTPUT_SHA256_MANIFEST.csv"))

cat(paste(summary_lines, collapse = "\n"), "\n")
message("Amendment 17 eICU albumin analysis completed.")
