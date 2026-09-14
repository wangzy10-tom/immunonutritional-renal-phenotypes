# ==============================================================================
# MIMIC-IV albumin-only primary conceptual replication (Amendment 16)
# ==============================================================================
# Purpose:
#   Reconstruct the outcome-blind albumin-only phenotype and place the already
#   approved analyte-consistent MIMIC-IV analysis in a separate, auditable
#   output directory. Existing freeze-v2 outputs are read-only inputs and are
#   not overwritten.

required_packages <- c("dplyr", "readr", "tibble", "survival")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop("Missing required package(s): ", paste(missing_packages, collapse = ", "), call. = FALSE)
}

root <- normalizePath(
  Sys.getenv("PROJECT_ROOT", unset = getwd()), winslash = "/", mustWork = TRUE
)
output_dir <- file.path(root, "output", "mimic_albumin_primary_2026-08-30")
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

collapse_for_selection <- function(x, selected, min_total = 100L, min_selected = 10L) {
  value <- as.character(x)
  value[is.na(value) | trimws(value) == ""] <- "Missing"
  tab_total <- table(value)
  tab_selected <- table(value[selected == 1])
  selected_count <- as.numeric(tab_selected[names(tab_total)])
  selected_count[is.na(selected_count)] <- 0
  excluded_count <- as.numeric(tab_total) - selected_count
  keep <- names(tab_total)[
    as.numeric(tab_total) >= min_total &
      selected_count >= min_selected & excluded_count >= min_selected
  ]
  factor(ifelse(value %in% keep, value, "Other"))
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
    dplyr::mutate(
      effect_95ci = sprintf("%.2f (%.2f-%.2f)", estimate, lower_95, upper_95)
    )
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
    dplyr::mutate(
      effect_95ci = sprintf("%.2f (%.2f-%.2f)", estimate, lower_95, upper_95)
    )
}

denominator_path <- Sys.getenv("MIMIC_DENOMINATOR_CSV", unset = file.path(root, "data", "derived", "MIMIC_first_ICU_feature_availability_dataset.csv"))
oasis_path <- Sys.getenv("MIMIC_OASIS_CSV", unset = file.path(root, "data", "derived", "MIMIC_official_OASIS_v301_scores.csv"))
sofa_path <- Sys.getenv("MIMIC_SOFA_CSV", unset = file.path(root, "data", "derived", "MIMIC_official_first_day_SOFA_v301_scores.csv"))
required_inputs <- c(denominator_path, oasis_path, sofa_path)
if (any(!file.exists(required_inputs))) {
  stop(
    "Missing required input(s): ",
    paste(required_inputs[!file.exists(required_inputs)], collapse = "; "),
    call. = FALSE
  )
}

denominator <- readr::read_csv(denominator_path, show_col_types = FALSE)
required_columns <- c(
  "stay_id", "nlr", "sii", "haemoglobin", "albumin", "bmi", "creatinine",
  "anchor_age", "gender", "race", "first_careunit", "anchor_year_group",
  "survival_days_365", "mortality_365d", "hospital_expire_flag", "intime", "dischtime"
)
if (any(!required_columns %in% names(denominator))) {
  stop(
    "Missing denominator column(s): ",
    paste(setdiff(required_columns, names(denominator)), collapse = ", "),
    call. = FALSE
  )
}

albumin_cohort <- denominator |>
  dplyr::filter(dplyr::if_all(dplyr::all_of(required_columns), ~ !is.na(.x))) |>
  dplyr::filter(
    is.finite(nlr), is.finite(sii), is.finite(haemoglobin),
    is.finite(albumin), is.finite(bmi), is.finite(creatinine),
    nlr > 0, sii > 0, albumin > 0, creatinine > 0
  )

matrix_data <- albumin_cohort |>
  dplyr::transmute(
    log_nlr = log(pmax(winsorise(nlr), .Machine$double.eps)),
    log_sii = log(pmax(winsorise(sii), .Machine$double.eps)),
    haemoglobin = winsorise(haemoglobin),
    albumin = winsorise(albumin),
    bmi = winsorise(bmi),
    log_creatinine = log(pmax(winsorise(creatinine), .Machine$double.eps))
  ) |>
  as.data.frame() |>
  scale()

set.seed(20260710)
cluster_fit <- stats::kmeans(
  matrix_data, centers = 3, nstart = 100, iter.max = 500, algorithm = "Lloyd"
)
profiles <- albumin_cohort |>
  dplyr::mutate(cluster_raw = cluster_fit$cluster) |>
  dplyr::group_by(cluster_raw) |>
  dplyr::summarise(
    n = dplyr::n(),
    nlr = stats::median(nlr),
    sii = stats::median(sii),
    haemoglobin = stats::median(haemoglobin),
    albumin = stats::median(albumin),
    bmi = stats::median(bmi),
    creatinine = stats::median(creatinine),
    .groups = "drop"
  )
p1_score <- safe_z(log(profiles$nlr)) + safe_z(log(profiles$sii)) +
  safe_z(log(profiles$creatinine))
p1_raw <- profiles$cluster_raw[which.max(p1_score)]
remaining <- profiles |> dplyr::filter(cluster_raw != p1_raw)
p2_score <- -safe_z(remaining$haemoglobin) - safe_z(remaining$albumin) -
  safe_z(remaining$bmi)
p2_raw <- remaining$cluster_raw[which.max(p2_score)]
p3_raw <- setdiff(profiles$cluster_raw, c(p1_raw, p2_raw))
mapping <- tibble::tibble(
  cluster_raw = c(p1_raw, p2_raw, p3_raw),
  phenotype_label = c("P1", "P2", "P3")
)
albumin_cohort <- albumin_cohort |>
  dplyr::mutate(cluster_raw = cluster_fit$cluster) |>
  dplyr::left_join(mapping, by = "cluster_raw") |>
  dplyr::mutate(
    phenotype = factor(phenotype_label, levels = c("P3", "P2", "P1")),
    male = as.integer(gender == "M")
  )
profiles <- profiles |>
  dplyr::left_join(mapping, by = "cluster_raw") |>
  dplyr::rename(phenotype = phenotype_label) |>
  dplyr::arrange(factor(phenotype, levels = c("P1", "P2", "P3")))

oasis <- readr::read_csv(oasis_path, show_col_types = FALSE) |>
  dplyr::select(stay_id, oasis) |>
  dplyr::distinct(stay_id, .keep_all = TRUE)
sofa <- readr::read_csv(sofa_path, show_col_types = FALSE) |>
  dplyr::select(stay_id, sofa) |>
  dplyr::distinct(stay_id, .keep_all = TRUE)

analysis <- albumin_cohort |>
  dplyr::left_join(oasis, by = "stay_id") |>
  dplyr::left_join(sofa, by = "stay_id") |>
  dplyr::filter(is.finite(oasis), is.finite(sofa)) |>
  dplyr::mutate(
    oasis_quartile = factor(
      stable_ntile(oasis, 4, stay_id),
      levels = 1:4, labels = paste0("Q", 1:4)
    ),
    sofa_quartile = factor(
      stable_ntile(sofa, 4, stay_id),
      levels = 1:4, labels = paste0("Q", 1:4)
    )
  )

if (nrow(analysis) != 1100L || sum(analysis$mortality_365d == 1) != 550L) {
  stop("Albumin cohort did not reproduce the locked n=1100 and 550 deaths.", call. = FALSE)
}
locked_counts <- c(P3 = 401L, P2 = 257L, P1 = 442L)
observed_counts <- table(analysis$phenotype)
if (!identical(as.integer(observed_counts[names(locked_counts)]), unname(locked_counts))) {
  stop(
    "Albumin phenotype counts did not reproduce the locked P3/P2/P1 counts. Observed: ",
    paste(names(observed_counts), as.integer(observed_counts), collapse = ", "),
    call. = FALSE
  )
}

reversed <- analysis[order(analysis$stay_id, decreasing = TRUE), ]
oasis_reversed <- stable_ntile(reversed$oasis, 4, reversed$stay_id)
sofa_reversed <- stable_ntile(reversed$sofa, 4, reversed$stay_id)
oasis_order_invariant <- identical(
  as.integer(analysis$oasis_quartile),
  oasis_reversed[match(analysis$stay_id, reversed$stay_id)]
)
sofa_order_invariant <- identical(
  as.integer(analysis$sofa_quartile),
  sofa_reversed[match(analysis$stay_id, reversed$stay_id)]
)
if (!oasis_order_invariant || !sofa_order_invariant) {
  stop("Deterministic severity-score quartile assignment failed.", call. = FALSE)
}

denominator <- denominator |>
  dplyr::mutate(
    selected_primary = as.integer(stay_id %in% analysis$stay_id),
    gender_selection = collapse_for_selection(gender, selected_primary),
    race_selection = collapse_for_selection(race, selected_primary),
    careunit_selection = collapse_for_selection(first_careunit, selected_primary),
    year_selection = collapse_for_selection(anchor_year_group, selected_primary)
  )
selection_fit <- stats::glm(
  selected_primary ~ splines::ns(anchor_age, df = 3) + gender_selection +
    race_selection + careunit_selection + year_selection,
  data = denominator,
  family = stats::binomial(),
  control = stats::glm.control(maxit = 100)
)
if (!isTRUE(selection_fit$converged)) {
  stop("Albumin selection model did not converge.", call. = FALSE)
}
probability <- pmin(pmax(stats::predict(selection_fit, type = "response"), 0.001), 0.999)
selection_rate <- mean(denominator$selected_primary == 1)
raw_weight <- ifelse(
  denominator$selected_primary == 1, selection_rate / probability, NA_real_
)
trim_limits <- stats::quantile(
  raw_weight[denominator$selected_primary == 1], c(0.01, 0.99), na.rm = TRUE
)
denominator$selection_ipw <- ifelse(
  denominator$selected_primary == 1,
  pmin(pmax(raw_weight, trim_limits[1]), trim_limits[2]),
  NA_real_
)
analysis <- analysis |>
  dplyr::left_join(
    denominator |> dplyr::select(stay_id, selection_ipw), by = "stay_id"
  ) |>
  dplyr::mutate(selection_ipw_scaled = selection_ipw / mean(selection_ipw))

oasis_continuous <- survival::coxph(
  survival::Surv(survival_days_365, mortality_365d) ~ phenotype + male + oasis,
  data = analysis, ties = "efron"
)
oasis_stratified <- survival::coxph(
  survival::Surv(survival_days_365, mortality_365d) ~
    phenotype + male + strata(oasis_quartile),
  data = analysis, ties = "efron"
)
oasis_ipw <- survival::coxph(
  survival::Surv(survival_days_365, mortality_365d) ~
    phenotype + male + strata(oasis_quartile),
  data = analysis,
  weights = selection_ipw_scaled,
  robust = TRUE,
  cluster = stay_id,
  ties = "efron"
)
landmark <- analysis |>
  dplyr::filter(survival_days_365 > 1) |>
  dplyr::mutate(
    landmark_time = survival_days_365 - 1,
    landmark_death = mortality_365d
  )
oasis_landmark <- survival::coxph(
  survival::Surv(landmark_time, landmark_death) ~
    phenotype + male + strata(oasis_quartile),
  data = landmark, ties = "efron"
)
oasis_hospital <- stats::glm(
  hospital_expire_flag ~ phenotype + male + oasis,
  data = analysis, family = stats::binomial()
)

sofa_continuous <- survival::coxph(
  survival::Surv(survival_days_365, mortality_365d) ~
    phenotype + anchor_age + male + sofa,
  data = analysis, ties = "efron"
)
sofa_stratified <- survival::coxph(
  survival::Surv(survival_days_365, mortality_365d) ~
    phenotype + anchor_age + male + strata(sofa_quartile),
  data = analysis, ties = "efron"
)
sofa_ipw <- survival::coxph(
  survival::Surv(survival_days_365, mortality_365d) ~
    phenotype + anchor_age + male + strata(sofa_quartile),
  data = analysis,
  weights = selection_ipw_scaled,
  robust = TRUE,
  cluster = stay_id,
  ties = "efron"
)
sofa_landmark <- survival::coxph(
  survival::Surv(landmark_time, landmark_death) ~
    phenotype + anchor_age + male + strata(sofa_quartile),
  data = landmark, ties = "efron"
)
sofa_hospital <- stats::glm(
  hospital_expire_flag ~ phenotype + anchor_age + male + sofa,
  data = analysis, family = stats::binomial()
)

counts <- analysis |>
  dplyr::group_by(phenotype) |>
  dplyr::summarise(
    n = dplyr::n(),
    deaths_365d = sum(mortality_365d == 1),
    mortality_365d_percent = 100 * mean(mortality_365d == 1),
    hospital_deaths = sum(hospital_expire_flag == 1),
    hospital_mortality_percent = 100 * mean(hospital_expire_flag == 1),
    .groups = "drop"
  ) |>
  dplyr::arrange(phenotype)

oasis_models <- dplyr::bind_rows(
  extract_cox(oasis_continuous, "Official OASIS continuous adjusted Cox", 1100L, 550L),
  extract_cox(oasis_stratified, "Official OASIS quartile-stratified Cox", 1100L, 550L),
  extract_cox(oasis_ipw, "Selection-IPW official OASIS quartile-stratified Cox", 1100L, 550L),
  extract_cox(
    oasis_landmark, "24-hour landmark official OASIS quartile-stratified Cox",
    nrow(landmark), sum(landmark$landmark_death == 1)
  ),
  extract_logistic(
    oasis_hospital, "Hospital mortality official OASIS adjusted logistic model",
    1100L, sum(analysis$hospital_expire_flag == 1)
  )
)
sofa_models <- dplyr::bind_rows(
  extract_cox(sofa_continuous, "Official first-day SOFA continuous adjusted Cox", 1100L, 550L),
  extract_cox(sofa_stratified, "Official first-day SOFA quartile-stratified Cox", 1100L, 550L),
  extract_cox(
    sofa_ipw, "Selection-IPW official first-day SOFA quartile-stratified Cox",
    1100L, 550L
  ),
  extract_cox(
    sofa_landmark, "24-hour landmark official first-day SOFA quartile-stratified Cox",
    nrow(landmark), sum(landmark$landmark_death == 1)
  ),
  extract_logistic(
    sofa_hospital, "Hospital mortality official first-day SOFA adjusted logistic model",
    1100L, sum(analysis$hospital_expire_flag == 1)
  )
)

ph_checks <- dplyr::bind_rows(
  survival::cox.zph(oasis_stratified)$table |>
    as.data.frame() |>
    tibble::rownames_to_column("term") |>
    dplyr::mutate(model = "Official OASIS quartile-stratified Cox", .before = 1),
  survival::cox.zph(sofa_stratified)$table |>
    as.data.frame() |>
    tibble::rownames_to_column("term") |>
    dplyr::mutate(model = "Official first-day SOFA quartile-stratified Cox", .before = 1)
)

selection_diagnostics <- tibble::tibble(
  denominator_source = basename(denominator_path),
  denominator_n = nrow(denominator),
  albumin_primary_n = nrow(analysis),
  albumin_primary_percent = 100 * selection_rate,
  weight_1st_percentile = trim_limits[1],
  weight_99th_percentile = trim_limits[2],
  weight_max = max(analysis$selection_ipw),
  effective_sample_size = sum(analysis$selection_ipw)^2 / sum(analysis$selection_ipw^2),
  oasis_quartile_row_order_invariant = oasis_order_invariant,
  sofa_quartile_row_order_invariant = sofa_order_invariant
)

readr::write_csv(counts, file.path(output_dir, "Table80A_MIMIC_albumin_counts.csv"))
readr::write_csv(profiles, file.path(output_dir, "Table80B_MIMIC_albumin_profiles.csv"))
readr::write_csv(oasis_models, file.path(output_dir, "Table80C_MIMIC_albumin_OASIS_models.csv"))
readr::write_csv(sofa_models, file.path(output_dir, "Table80D_MIMIC_albumin_SOFA_models.csv"))
readr::write_csv(selection_diagnostics, file.path(output_dir, "Table80E_MIMIC_albumin_selection.csv"))
readr::write_csv(ph_checks, file.path(output_dir, "Table80F_MIMIC_albumin_PH_checks.csv"))
saveRDS(
  list(
    analysis = analysis,
    counts = counts,
    profiles = profiles,
    oasis_models = oasis_models,
    sofa_models = sofa_models,
    selection_diagnostics = selection_diagnostics,
    ph_checks = ph_checks
  ),
  file.path(output_dir, "MIMIC_albumin_primary_results.rds")
)

primary_p1 <- oasis_models |>
  dplyr::filter(model == "Official OASIS quartile-stratified Cox", comparison == "P1 vs P3")
ipw_p1 <- oasis_models |>
  dplyr::filter(
    model == "Selection-IPW official OASIS quartile-stratified Cox",
    comparison == "P1 vs P3"
  )
summary_lines <- c(
  "MIMIC-IV albumin-only primary conceptual replication — Amendment 16",
  "",
  paste0("Cohort: n=", nrow(analysis), "; 365-day deaths=", sum(analysis$mortality_365d == 1), "."),
  paste0(
    "Primary official-OASIS quartile-stratified P1 vs P3: ",
    primary_p1$effect_95ci, "; P=", format(primary_p1$p_value, scientific = TRUE, digits = 4), "."
  ),
  paste0(
    "Selection-IPW official-OASIS P1 vs P3: ",
    ipw_p1$effect_95ci, "; P=", format(ipw_p1$p_value, scientific = TRUE, digits = 4), "."
  ),
  "",
  "All submission analyses in this release use albumin as the shared analyte.",
  "The NLR/SII feature set is unchanged in Amendment 16."
)
writeLines(summary_lines, file.path(output_dir, "MIMIC_albumin_primary_summary.txt"))

cat(paste(summary_lines, collapse = "\n"), "\n")
