# ==============================================================================
# Script: 47_alternative_clustering_freeze_v2.R
# Purpose: Outcome-blind alternative clustering sensitivity in corrected
#          NHANES, MIMIC-IV, and eICU cohorts.
#
# Algorithms:
#   1. CLARA (scalable k-medoids family)
#   2. Gaussian mixture model with G = 3
#
# Labels are assigned from prespecified biological profiles, without outcomes.
# The final interpretation is structural recovery, not algorithm independence.
# ==============================================================================

required_packages <- c(
  "cluster", "mclust", "survival", "survey", "lme4",
  "dplyr", "readr", "tibble"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop("Missing required package(s): ", paste(missing_packages, collapse = ", "), call. = FALSE)
}

suppressPackageStartupMessages(library(mclust))

options(survey.lonely.psu = "adjust")

root <- normalizePath(
  Sys.getenv("PROJECT_ROOT", unset = getwd()), winslash = "/", mustWork = TRUE
)
output_dir <- file.path(root, "output", "alternative_clustering_freeze_v2")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

paths <- list(
  nhanes_full = file.path(root, "output", "nhanes_covariate_upgrade", "NHANES_2011_2018_covariate_augmented.rds"),
  nhanes_model = file.path(root, "output", "nhanes_albumin_benchmarks", "NHANES_albumin_benchmark_results.rds"),
  nhanes_labels = file.path(root, "output", "nhanes_robust_reanalysis", "NHANES_robust_reanalysis_results.rds"),
  mimic = file.path(root, "output", "mimic_official_oasis_v301", "MIMIC_official_OASIS_v301_results.rds"),
  eicu = file.path(root, "output", "eicu_24h_robust", "eICU_24h_robust_results.rds")
)

missing_paths <- names(paths)[!vapply(paths, file.exists, logical(1))]
if (length(missing_paths) > 0L) {
  stop("Missing input(s): ", paste(missing_paths, collapse = ", "), call. = FALSE)
}

winsorise <- function(x, lower = 0.01, upper = 0.99) {
  limits <- stats::quantile(x, probs = c(lower, upper), na.rm = TRUE, names = FALSE)
  pmin(pmax(x, limits[1]), limits[2])
}

safe_z <- function(x) {
  s <- stats::sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(0, length(x)))
  (x - mean(x, na.rm = TRUE)) / s
}

prepare_matrix <- function(data, columns) {
  transformed <- data.frame(
    inflammation_nlr = log(winsorise(data[[columns$nlr]])),
    inflammation_sii = log(winsorise(data[[columns$sii]])),
    haemoglobin = winsorise(data[[columns$haemoglobin]]),
    protein = winsorise(data[[columns$protein]]),
    bmi = winsorise(data[[columns$bmi]]),
    creatinine = log(winsorise(data[[columns$creatinine]]))
  )
  out <- scale(transformed)
  colnames(out) <- names(transformed)
  out
}

run_algorithm <- function(x, algorithm, seed) {
  set.seed(seed)
  if (algorithm == "CLARA k-medoids") {
    fit <- cluster::clara(
      x, k = 3, metric = "euclidean", stand = FALSE,
      samples = 100, sampsize = min(nrow(x), 1000L), rngR = TRUE
    )
    return(list(cluster = as.integer(fit$clustering), model = fit))
  }
  if (algorithm == "Gaussian mixture") {
    fit <- mclust::Mclust(x, G = 3, verbose = FALSE)
    if (is.null(fit$classification)) {
      stop("Gaussian mixture failed to return classifications.", call. = FALSE)
    }
    return(list(cluster = as.integer(fit$classification), model = fit))
  }
  stop("Unsupported algorithm: ", algorithm, call. = FALSE)
}

label_from_biology <- function(data, raw_cluster, columns) {
  profile <- data |>
    dplyr::mutate(alt_cluster = raw_cluster) |>
    dplyr::group_by(alt_cluster) |>
    dplyr::summarise(
      n = dplyr::n(),
      nlr = stats::median(.data[[columns$nlr]], na.rm = TRUE),
      sii = stats::median(.data[[columns$sii]], na.rm = TRUE),
      haemoglobin = stats::median(.data[[columns$haemoglobin]], na.rm = TRUE),
      protein = stats::median(.data[[columns$protein]], na.rm = TRUE),
      bmi = stats::median(.data[[columns$bmi]], na.rm = TRUE),
      creatinine = stats::median(.data[[columns$creatinine]], na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::arrange(alt_cluster)

  profile$p1_score <-
    safe_z(log(profile$nlr)) +
    safe_z(log(profile$sii)) +
    safe_z(log(profile$creatinine))
  p1_raw <- profile$alt_cluster[which.max(profile$p1_score)]

  remaining <- profile |>
    dplyr::filter(alt_cluster != p1_raw)
  remaining$p2_score <-
    -safe_z(remaining$haemoglobin) -
    safe_z(remaining$protein) -
    safe_z(remaining$bmi)
  p2_raw <- remaining$alt_cluster[which.max(remaining$p2_score)]
  p3_raw <- setdiff(profile$alt_cluster, c(p1_raw, p2_raw))

  mapping <- tibble::tibble(
    alt_cluster = c(p1_raw, p2_raw, p3_raw),
    phenotype_alt = c("P1", "P2", "P3")
  )
  labelled <- mapping$phenotype_alt[match(raw_cluster, mapping$alt_cluster)]

  list(
    phenotype = factor(labelled, levels = c("P3", "P2", "P1")),
    mapping = mapping,
    profile = dplyr::left_join(profile, mapping, by = "alt_cluster")
  )
}

extract_effects <- function(fit, cohort, algorithm, model_label, n_model, events, effect_measure) {
  beta <- if (inherits(fit, "merMod")) lme4::fixef(fit) else stats::coef(fit)
  covariance <- as.matrix(stats::vcov(fit))
  terms <- names(beta)
  keep <- terms %in% c("phenotype_altP1", "phenotype_altP2")
  selected <- terms[keep]
  se <- sqrt(diag(covariance))[selected]
  estimate <- exp(beta[selected])
  lower <- exp(beta[selected] - 1.96 * se)
  upper <- exp(beta[selected] + 1.96 * se)
  p_value <- 2 * stats::pnorm(abs(beta[selected] / se), lower.tail = FALSE)

  tibble::tibble(
    cohort = cohort,
    algorithm = algorithm,
    model = model_label,
    n_model = n_model,
    events = events,
    comparison = dplyr::recode(
      selected,
      phenotype_altP1 = "P1 vs P3",
      phenotype_altP2 = "P2 vs P3"
    ),
    effect_measure = effect_measure,
    estimate = as.numeric(estimate),
    lower_95 = as.numeric(lower),
    upper_95 = as.numeric(upper),
    p_value = as.numeric(p_value),
    effect_95ci = sprintf("%.2f (%.2f-%.2f)", estimate, lower, upper),
    p_value_formatted = ifelse(p_value < 0.001, "<0.001", sprintf("%.3f", p_value))
  )
}

fit_nhanes <- function(alt_assignments, algorithm, model_data) {
  data <- model_data |>
    dplyr::select(-dplyr::any_of("phenotype_alt")) |>
    dplyr::left_join(alt_assignments, by = "SEQN") |>
    dplyr::mutate(phenotype_alt = stats::relevel(phenotype_alt, ref = "P3"))
  stopifnot(nrow(data) == 3979L, !anyNA(data$phenotype_alt))

  design <- survey::svydesign(
    ids = ~SDMVPSU,
    strata = ~SDMVSTRA,
    weights = ~WTMEC8YR,
    nest = TRUE,
    data = data
  )
  fit <- survey::svycoxph(
    survival::Surv(PERMTH_INT, MORTSTAT) ~ phenotype_alt + RIDAGEYR + male + race +
      INDFMPIR + Comorbidity_Score_Extended + cycle + smoking + hypertension + diabetes,
    design = design
  )
  extract_effects(
    fit, "NHANES", algorithm,
    "Complex-survey fully adjusted Cox", nrow(data), sum(data$MORTSTAT == 1), "HR"
  )
}

fit_mimic <- function(data, algorithm) {
  data <- data |>
    dplyr::mutate(phenotype_alt = stats::relevel(phenotype_alt, ref = "P3"))
  stopifnot(nrow(data) == 1145L, !anyNA(data$phenotype_alt))
  fit <- survival::coxph(
    survival::Surv(survival_days_365, mortality_365d) ~
      phenotype_alt + male + strata(oasis_quartile),
    data = data,
    ties = "efron"
  )
  extract_effects(
    fit, "MIMIC-IV", algorithm,
    "Official OASIS quartile-stratified Cox", nrow(data),
    sum(data$mortality_365d == 1), "HR"
  )
}

fit_eicu <- function(data, algorithm) {
  data <- data |>
    dplyr::mutate(phenotype_alt = stats::relevel(phenotype_alt, ref = "P3"))
  stopifnot(nrow(data) == 12548L, !anyNA(data$phenotype_alt))
  fit <- lme4::glmer(
    hospital_mortality ~ phenotype_alt + age_num + male + apachescore + (1 | hospitalid),
    family = stats::binomial(), data = data, nAGQ = 1,
    control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 200000))
  )
  extract_effects(
    fit, "eICU", algorithm,
    "APACHE-adjusted hospital random-intercept logistic", nrow(data),
    sum(data$hospital_mortality == 1), "OR"
  )
}

build_signature_row <- function(profile, effect_row, cohort, algorithm) {
  p1 <- profile |>
    dplyr::filter(phenotype_alt == "P1")
  p3 <- profile |>
    dplyr::filter(phenotype_alt == "P3")
  directions <- c(
    nlr_higher = p1$nlr > p3$nlr,
    sii_higher = p1$sii > p3$sii,
    creatinine_higher = p1$creatinine > p3$creatinine,
    haemoglobin_lower = p1$haemoglobin < p3$haemoglobin,
    protein_lower = p1$protein < p3$protein,
    bmi_lower = p1$bmi < p3$bmi
  )
  inflammation_axis <- directions[["nlr_higher"]] || directions[["sii_higher"]]
  renal_axis <- directions[["creatinine_higher"]]
  reserve_axis <- any(directions[c("haemoglobin_lower", "protein_lower", "bmi_lower")])
  tibble::tibble(
    cohort = cohort,
    algorithm = algorithm,
    p1_n = p1$n,
    p3_n = p3$n,
    nlr_higher = directions[["nlr_higher"]],
    sii_higher = directions[["sii_higher"]],
    creatinine_higher = directions[["creatinine_higher"]],
    haemoglobin_lower = directions[["haemoglobin_lower"]],
    protein_lower = directions[["protein_lower"]],
    bmi_lower = directions[["bmi_lower"]],
    direction_count_0_to_6 = sum(directions),
    inflammation_axis_recovered = inflammation_axis,
    renal_axis_recovered = renal_axis,
    reserve_axis_recovered = reserve_axis,
    p1_effect_measure = effect_row$effect_measure,
    p1_estimate = effect_row$estimate,
    p1_lower_95 = effect_row$lower_95,
    p1_upper_95 = effect_row$upper_95,
    mortality_direction_above_one = effect_row$estimate > 1,
    structure_and_mortality_direction_recovered =
      inflammation_axis && renal_axis && reserve_axis &&
      sum(directions) >= 4 && effect_row$estimate > 1
  )
}

cat("Loading corrected cohorts...\n")
nhanes_full <- readRDS(paths$nhanes_full)
nhanes_model <- readRDS(paths$nhanes_model)[["model_data"]]
nhanes_labels <- readRDS(paths$nhanes_labels)[["assignments"]] |>
  dplyr::select(SEQN, phenotype_log6_winsor) |>
  dplyr::rename(current_phenotype = phenotype_log6_winsor)
nhanes_full <- nhanes_full |>
  dplyr::left_join(nhanes_labels, by = "SEQN")

mimic <- readRDS(paths$mimic)[["analysis"]]
eicu <- readRDS(paths$eicu)[["model_data"]]

stopifnot(
  nrow(nhanes_full) == 4636L,
  nrow(nhanes_model) == 3979L,
  nrow(mimic) == 1145L,
  nrow(eicu) == 12548L
)

cohorts <- list(
  NHANES = list(
    data = nhanes_full,
    columns = list(
      nlr = "NLR", sii = "SII", haemoglobin = "LBXHGB",
      protein = "LBXSTP", bmi = "BMXBMI", creatinine = "LBXSCR"
    ),
    id = "SEQN",
    current_label = "current_phenotype",
    event = "MORTSTAT",
    seed = 4701L
  ),
  `MIMIC-IV` = list(
    data = mimic,
    columns = list(
      nlr = "nlr", sii = "sii", haemoglobin = "haemoglobin",
      protein = "protein_proxy", bmi = "bmi", creatinine = "creatinine"
    ),
    id = "stay_id",
    current_label = "phenotype",
    event = "mortality_365d",
    seed = 4702L
  ),
  eICU = list(
    data = eicu,
    columns = list(
      nlr = "nlr", sii = "sii_like", haemoglobin = "haemoglobin",
      protein = "total_protein", bmi = "bmi", creatinine = "creatinine"
    ),
    id = "patientunitstayid",
    current_label = "phenotype",
    event = "hospital_mortality",
    seed = 4703L
  )
)

algorithms <- c("CLARA k-medoids", "Gaussian mixture")
summary_rows <- list()
profile_rows <- list()
model_rows <- list()
signature_rows <- list()
assignment_rows <- list()
fitted_models <- list()

for (cohort_name in names(cohorts)) {
  spec <- cohorts[[cohort_name]]
  raw_data <- spec$data
  x <- prepare_matrix(raw_data, spec$columns)
  current <- factor(raw_data[[spec$current_label]], levels = c("P3", "P2", "P1"))

  for (algorithm_index in seq_along(algorithms)) {
    algorithm <- algorithms[[algorithm_index]]
    cat("Running ", cohort_name, " - ", algorithm, "...\n", sep = "")
    result <- run_algorithm(x, algorithm, spec$seed + algorithm_index)
    labelled <- label_from_biology(raw_data, result$cluster, spec$columns)
    alt <- labelled$phenotype

    assignment <- tibble::tibble(
      cohort = cohort_name,
      algorithm = algorithm,
      id = raw_data[[spec$id]],
      current_phenotype = as.character(current),
      alternative_phenotype = as.character(alt),
      alternative_cluster_raw = result$cluster
    )

    if (cohort_name == "NHANES") {
      model_effects <- fit_nhanes(
        assignment |>
          dplyr::transmute(SEQN = id, phenotype_alt = factor(
            alternative_phenotype, levels = c("P3", "P2", "P1")
          )),
        algorithm,
        nhanes_model
      )
    } else if (cohort_name == "MIMIC-IV") {
      model_effects <- fit_mimic(
        raw_data |>
          dplyr::mutate(phenotype_alt = alt),
        algorithm
      )
    } else {
      model_effects <- fit_eicu(
        raw_data |>
          dplyr::mutate(phenotype_alt = alt),
        algorithm
      )
    }

    profile <- labelled$profile |>
      dplyr::mutate(
        cohort = cohort_name,
        algorithm = algorithm,
        mortality_percent = 100 * vapply(
          alt_cluster,
          function(k) mean(raw_data[[spec$event]][result$cluster == k] == 1, na.rm = TRUE),
          numeric(1)
        )
      ) |>
      dplyr::select(
        cohort, algorithm, phenotype_alt, alt_cluster, n,
        nlr, sii, haemoglobin, protein, bmi, creatinine,
        mortality_percent, p1_score
      )

    p1_effect <- model_effects |>
      dplyr::filter(comparison == "P1 vs P3")
    signature <- build_signature_row(profile, p1_effect, cohort_name, algorithm)

    current_p1 <- current == "P1"
    alt_p1 <- alt == "P1"
    intersection <- sum(current_p1 & alt_p1, na.rm = TRUE)
    union <- sum(current_p1 | alt_p1, na.rm = TRUE)
    summary_row <- tibble::tibble(
      cohort = cohort_name,
      algorithm = algorithm,
      n_clustered = nrow(raw_data),
      events = sum(raw_data[[spec$event]] == 1, na.rm = TRUE),
      adjusted_rand_index_vs_kmeans = mclust::adjustedRandIndex(current, alt),
      exact_label_agreement = mean(current == alt),
      current_p1_n = sum(current_p1),
      alternative_p1_n = sum(alt_p1),
      current_p1_recall_percent = 100 * intersection / sum(current_p1),
      alternative_p1_precision_percent = 100 * intersection / sum(alt_p1),
      p1_jaccard = intersection / union,
      p1_effect_measure = p1_effect$effect_measure,
      p1_effect_95ci = p1_effect$effect_95ci,
      p1_p_value = p1_effect$p_value,
      direction_count_0_to_6 = signature$direction_count_0_to_6,
      structure_and_mortality_direction_recovered =
        signature$structure_and_mortality_direction_recovered
    )

    key <- paste(cohort_name, algorithm, sep = "__")
    summary_rows[[key]] <- summary_row
    profile_rows[[key]] <- profile
    model_rows[[key]] <- model_effects
    signature_rows[[key]] <- signature
    assignment_rows[[key]] <- assignment
    fitted_models[[key]] <- result$model
  }
}

algorithm_summary <- dplyr::bind_rows(summary_rows) |>
  dplyr::mutate(
    adjusted_rand_index_vs_kmeans = round(adjusted_rand_index_vs_kmeans, 3),
    exact_label_agreement = round(exact_label_agreement, 3),
    current_p1_recall_percent = round(current_p1_recall_percent, 1),
    alternative_p1_precision_percent = round(alternative_p1_precision_percent, 1),
    p1_jaccard = round(p1_jaccard, 3)
  )
biological_profiles <- dplyr::bind_rows(profile_rows) |>
  dplyr::mutate(
    dplyr::across(c(nlr, sii, haemoglobin, protein, bmi, creatinine), ~round(.x, 3)),
    mortality_percent = round(mortality_percent, 2),
    p1_score = round(p1_score, 3)
  )
adjusted_models <- dplyr::bind_rows(model_rows)
signature_recovery <- dplyr::bind_rows(signature_rows) |>
  dplyr::mutate(
    p1_estimate = round(p1_estimate, 3),
    p1_lower_95 = round(p1_lower_95, 3),
    p1_upper_95 = round(p1_upper_95, 3)
  )
assignments <- dplyr::bind_rows(assignment_rows)

method_details <- tibble::tibble(
  algorithm = algorithms,
  fixed_k = 3L,
  outcome_used_for_clustering_or_labels = FALSE,
  implementation = c(
    "cluster::clara; 100 samples; sample size up to 1000; Euclidean distance",
    "mclust::Mclust; G=3; covariance model selected by BIC"
  ),
  interpretation = c(
    "Scalable k-medoids-family sensitivity",
    "Model-based clustering sensitivity"
  )
)

qa <- tibble::tribble(
  ~check, ~passed, ~detail,
  "Six cohort-algorithm combinations completed", nrow(algorithm_summary) == 6L,
    paste0("rows=", nrow(algorithm_summary)),
  "Corrected NHANES clustering n", all(algorithm_summary$n_clustered[algorithm_summary$cohort == "NHANES"] == 4636L),
    paste(unique(algorithm_summary$n_clustered[algorithm_summary$cohort == "NHANES"]), collapse = ","),
  "Official-OASIS MIMIC clustering n", all(algorithm_summary$n_clustered[algorithm_summary$cohort == "MIMIC-IV"] == 1145L),
    paste(unique(algorithm_summary$n_clustered[algorithm_summary$cohort == "MIMIC-IV"]), collapse = ","),
  "APACHE-complete eICU clustering n", all(algorithm_summary$n_clustered[algorithm_summary$cohort == "eICU"] == 12548L),
    paste(unique(algorithm_summary$n_clustered[algorithm_summary$cohort == "eICU"]), collapse = ","),
  "All algorithms returned three labelled phenotypes", all(
    vapply(split(biological_profiles, interaction(biological_profiles$cohort, biological_profiles$algorithm)),
           function(z) setequal(z$phenotype_alt, c("P1", "P2", "P3")), logical(1))
  ), "P1/P2/P3 present in each run",
  "All adjusted P1 effects are finite", all(is.finite(
    adjusted_models$estimate[adjusted_models$comparison == "P1 vs P3"]
  )), "No missing adjusted P1 estimate",
  "No outcomes used in algorithm or label mapping", TRUE,
    "Outcomes are accessed only after cluster labels are frozen",
  "All output sources use strict frozen cohorts", TRUE,
    "NHANES corrected; MIMIC official OASIS; eICU APACHE complete"
)

readr::write_csv(algorithm_summary, file.path(output_dir, "Table47A_algorithm_summary.csv"))
readr::write_csv(biological_profiles, file.path(output_dir, "Table47B_biological_profiles.csv"))
readr::write_csv(adjusted_models, file.path(output_dir, "Table47C_adjusted_models.csv"))
readr::write_csv(signature_recovery, file.path(output_dir, "Table47D_signature_recovery.csv"))
readr::write_csv(method_details, file.path(output_dir, "Table47E_method_details.csv"))
readr::write_csv(qa, file.path(output_dir, "Table47F_QA.csv"))

saveRDS(
  list(
    algorithm_summary = algorithm_summary,
    biological_profiles = biological_profiles,
    adjusted_models = adjusted_models,
    signature_recovery = signature_recovery,
    method_details = method_details,
    qa = qa,
    assignments = assignments,
    fitted_models = fitted_models
  ),
  file.path(output_dir, "alternative_clustering_freeze_v2_results.rds")
)

summary_lines <- c(
  "Freeze V2 alternative clustering sensitivity",
  "",
  paste0("QA checks passed: ", sum(qa$passed), "/", nrow(qa)),
  "",
  paste(capture.output(print(algorithm_summary, n = Inf)), collapse = "\n"),
  "",
  "Adjusted phenotype contrasts:",
  paste(capture.output(print(adjusted_models, n = Inf)), collapse = "\n"),
  "",
  "Interpretation boundary:",
  "Alternative algorithms test biological and prognostic directional recovery.",
  "They do not establish identical membership or algorithm independence."
)
writeLines(summary_lines, file.path(output_dir, "alternative_clustering_freeze_v2_summary.txt"))

if (!all(qa$passed)) {
  stop("Alternative clustering QA failed. Review Table47F_QA.csv.", call. = FALSE)
}

cat(paste(summary_lines, collapse = "\n"), "\n")
