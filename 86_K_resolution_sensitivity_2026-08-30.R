# ==============================================================================
# Amendment 19: K-resolution sensitivity at K=2 and K=4
# Added after primary K=3 results were known; supportive, not preregistered.
# ==============================================================================

required_packages <- c(
  "survival", "survey", "lme4", "dplyr", "readr", "tibble",
  "data.table", "mclust", "digest"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop("Missing required package(s): ", paste(missing_packages, collapse = ", "), call. = FALSE)
}
options(survey.lonely.psu = "adjust")

root <- normalizePath(Sys.getenv("PROJECT_ROOT", unset = getwd()), winslash = "/", mustWork = TRUE)
eicu_dir <- Sys.getenv("EICU_DIR", unset = "")
protocol_path <- file.path(root, "ANALYSIS_FREEZE_AMENDMENT_19_K_RESOLUTION_PROTOCOL_2026-08-30.md")
script_path <- file.path(root, "86_K_resolution_sensitivity_2026-08-30.R")
paths <- list(
  nhanes_cluster = file.path(
    root, "output", "nhanes_albumin_amendment17_2026-08-30",
    "INTERNAL_NHANES_albumin_complete_cohort.rds"
  ),
  nhanes_model = file.path(
    root, "output", "nhanes_albumin_amendment17_2026-08-30",
    "NHANES_albumin_benchmark_results.rds"
  ),
  mimic = file.path(
    root, "output", "mimic_albumin_primary_2026-08-30",
    "MIMIC_albumin_primary_results.rds"
  ),
  eicu_denominator = file.path(
    root, "output", "eicu_24h_extraction",
    "eICU_first_ICU_feature_availability_dataset.csv"
  ),
  eicu_apache = file.path(eicu_dir, "apachePatientResult.csv.gz")
)
paths$eicu_denominator <- Sys.getenv("EICU_DENOMINATOR_CSV", unset=file.path(root,"data","derived","eICU_first_ICU_feature_availability_dataset.csv"))

required_paths <- c(unname(unlist(paths)), protocol_path, script_path)
missing_paths <- required_paths[!file.exists(required_paths)]
if (length(missing_paths) > 0L) {
  stop("Missing input(s): ", paste(missing_paths, collapse = ", "), call. = FALSE)
}
output_dir <- file.path(root, "output", "k_resolution_sensitivity_2026-08-30")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

winsorise <- function(x) {
  limits <- stats::quantile(x, c(0.01, 0.99), na.rm = TRUE, names = FALSE)
  pmin(pmax(as.numeric(x), limits[1]), limits[2])
}

safe_z <- function(x) {
  current_sd <- stats::sd(x, na.rm = TRUE)
  if (!is.finite(current_sd) || current_sd == 0) return(rep(0, length(x)))
  (x - mean(x, na.rm = TRUE)) / current_sd
}

prepare_matrix <- function(data) {
  transformed <- data.frame(
    log_nlr = log(pmax(winsorise(data$nlr), .Machine$double.eps)),
    log_sii = log(pmax(winsorise(data$sii), .Machine$double.eps)),
    haemoglobin = winsorise(data$haemoglobin),
    albumin = winsorise(data$albumin),
    bmi = winsorise(data$bmi),
    log_creatinine = log(pmax(winsorise(data$creatinine), .Machine$double.eps))
  )
  output <- scale(transformed)
  if (any(!is.finite(output))) stop("Non-finite transformed clustering value.", call. = FALSE)
  output
}

make_profile <- function(data, raw_cluster) {
  data |>
    dplyr::mutate(cluster_raw = as.integer(raw_cluster)) |>
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
    ) |>
    dplyr::arrange(cluster_raw) |>
    dplyr::mutate(
      p1_score = safe_z(log(nlr)) + safe_z(log(sii)) + safe_z(log(creatinine)),
      reserve_deficit_score = -safe_z(haemoglobin) - safe_z(albumin) - safe_z(bmi),
      favorable_reference_score = safe_z(haemoglobin) + safe_z(albumin) + safe_z(bmi) -
        safe_z(log(nlr)) - safe_z(log(sii)) - safe_z(log(creatinine))
    )
}

pick_max <- function(profile, score_name) {
  profile$cluster_raw[order(-profile[[score_name]], profile$cluster_raw)][1]
}

label_primary_k3 <- function(data, raw_cluster) {
  profile <- make_profile(data, raw_cluster)
  p1_raw <- pick_max(profile, "p1_score")
  remaining <- profile |> dplyr::filter(cluster_raw != p1_raw)
  p2_raw <- pick_max(remaining, "reserve_deficit_score")
  p3_raw <- setdiff(profile$cluster_raw, c(p1_raw, p2_raw))
  mapping <- tibble::tibble(
    cluster_raw = c(p1_raw, p2_raw, p3_raw),
    phenotype = c("P1", "P2", "P3")
  )
  factor(mapping$phenotype[match(raw_cluster, mapping$cluster_raw)], levels = c("P3", "P2", "P1"))
}

label_resolution <- function(data, raw_cluster, k) {
  profile <- make_profile(data, raw_cluster)
  p1_raw <- pick_max(profile, "p1_score")
  if (k == 2L) {
    p3_raw <- setdiff(profile$cluster_raw, p1_raw)
    mapping <- tibble::tibble(
      cluster_raw = c(p1_raw, p3_raw),
      phenotype = c("P1", "P3")
    )
    level_order <- c("P3", "P1")
  } else if (k == 4L) {
    after_p1 <- profile |> dplyr::filter(cluster_raw != p1_raw)
    p2_raw <- pick_max(after_p1, "reserve_deficit_score")
    after_p2 <- after_p1 |> dplyr::filter(cluster_raw != p2_raw)
    p3_raw <- pick_max(after_p2, "favorable_reference_score")
    p4_raw <- setdiff(profile$cluster_raw, c(p1_raw, p2_raw, p3_raw))
    mapping <- tibble::tibble(
      cluster_raw = c(p1_raw, p2_raw, p3_raw, p4_raw),
      phenotype = c("P1", "P2", "P3", "P4")
    )
    level_order <- c("P3", "P2", "P4", "P1")
  } else {
    stop("Only K=2 and K=4 are supported.", call. = FALSE)
  }
  labels <- mapping$phenotype[match(raw_cluster, mapping$cluster_raw)]
  list(
    phenotype = factor(labels, levels = level_order),
    profile = profile |>
      dplyr::left_join(mapping, by = "cluster_raw") |>
      dplyr::arrange(match(phenotype, c("P1", "P2", "P3", "P4"))),
    mapping = mapping,
    level_order = level_order
  )
}

run_resolution <- function(data, k) {
  x <- prepare_matrix(data)
  set.seed(20260710)
  fit <- stats::kmeans(
    x, centers = k, nstart = 100, iter.max = 500, algorithm = "Lloyd"
  )
  labelled <- label_resolution(data, fit$cluster, k)
  list(
    phenotype = labelled$phenotype,
    profile = labelled$profile,
    level_order = labelled$level_order,
    diagnostics = tibble::tibble(
      k = k,
      selected_model = paste0("K-means Lloyd, K=", k),
      objective = as.numeric(fit$tot.withinss),
      iterations = as.integer(fit$iter),
      converged = fit$iter < 500L,
      ifault = if (is.null(fit$ifault)) NA_integer_ else as.integer(fit$ifault)
    )
  )
}

extract_effects <- function(fit, dataset, k, model, n, events, effect_measure, comparisons) {
  beta <- if (inherits(fit, "merMod")) lme4::fixef(fit) else stats::coef(fit)
  covariance <- as.matrix(stats::vcov(fit))
  terms <- paste0("phenotype_resolution", comparisons)
  if (!all(terms %in% names(beta))) {
    stop("Resolution phenotype terms are missing: ", paste(setdiff(terms, names(beta)), collapse = ", "), call. = FALSE)
  }
  se <- sqrt(diag(covariance))[terms]
  estimate <- exp(beta[terms])
  lower <- exp(beta[terms] - 1.96 * se)
  upper <- exp(beta[terms] + 1.96 * se)
  tibble::tibble(
    dataset = dataset,
    k = k,
    model = model,
    n = n,
    events = events,
    comparison = paste(comparisons, "vs P3"),
    effect_measure = effect_measure,
    estimate = as.numeric(estimate),
    lower_95 = as.numeric(lower),
    upper_95 = as.numeric(upper),
    p_value = 2 * stats::pnorm(abs(beta[terms] / se), lower.tail = FALSE),
    effect_95ci = sprintf("%.3f (%.3f-%.3f)", estimate, lower, upper)
  )
}

fit_nhanes <- function(assignments, model_data, k, level_order) {
  current <- model_data |>
    dplyr::select(-dplyr::any_of("phenotype_resolution")) |>
    dplyr::left_join(assignments, by = "SEQN") |>
    dplyr::mutate(phenotype_resolution = factor(phenotype_resolution, levels = level_order))
  if (nrow(current) != 3979L || sum(current$MORTSTAT == 1L) != 720L || anyNA(current$phenotype_resolution)) {
    stop("NHANES resolution model cohort mismatch.", call. = FALSE)
  }
  design <- survey::svydesign(
    ids = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~WTMEC8YR,
    nest = TRUE, data = current
  )
  fit <- survey::svycoxph(
    survival::Surv(PERMTH_INT, MORTSTAT) ~ phenotype_resolution + RIDAGEYR +
      male + race + INDFMPIR + Comorbidity_Score_Extended + cycle + smoking +
      hypertension + diabetes,
    design = design
  )
  comparisons <- if (k == 2L) c("P1") else c("P1", "P2", "P4")
  list(
    effects = extract_effects(
      fit, "NHANES", k, "Complex-survey fully adjusted Cox",
      nrow(current), sum(current$MORTSTAT == 1L), "HR", comparisons
    ),
    diagnostics = tibble::tibble(
      fit_ok = TRUE, convergence_message = "", captured_warnings = "", singular = NA
    )
  )
}

fit_mimic <- function(data, k, level_order) {
  current <- data |>
    dplyr::mutate(phenotype_resolution = factor(phenotype_resolution, levels = level_order))
  if (nrow(current) != 1100L || sum(current$mortality_365d == 1L) != 550L || anyNA(current$phenotype_resolution)) {
    stop("MIMIC-IV resolution model cohort mismatch.", call. = FALSE)
  }
  fit <- survival::coxph(
    survival::Surv(survival_days_365, mortality_365d) ~
      phenotype_resolution + male + strata(oasis_quartile),
    data = current, ties = "efron"
  )
  comparisons <- if (k == 2L) c("P1") else c("P1", "P2", "P4")
  list(
    effects = extract_effects(
      fit, "MIMIC-IV", k, "Official OASIS quartile-stratified Cox",
      nrow(current), sum(current$mortality_365d == 1L), "HR", comparisons
    ),
    diagnostics = tibble::tibble(
      fit_ok = TRUE, convergence_message = "", captured_warnings = "", singular = NA
    )
  )
}

fit_eicu <- function(data, k, level_order) {
  current <- data |>
    dplyr::filter(
      !is.na(phenotype_resolution), is.finite(apachescore), !is.na(hospital_mortality)
    ) |>
    dplyr::mutate(
      phenotype_resolution = factor(phenotype_resolution, levels = level_order),
      hospitalid = factor(hospitalid)
    )
  if (nrow(current) != 13234L || sum(current$hospital_mortality == 1L) != 1990L ||
      dplyr::n_distinct(current$hospitalid) != 166L || anyNA(current$phenotype_resolution)) {
    stop("eICU resolution model cohort mismatch.", call. = FALSE)
  }
  warnings <- character()
  fit <- tryCatch(
    withCallingHandlers(
      lme4::glmer(
        hospital_mortality ~ phenotype_resolution + age_num + male + apachescore + (1 | hospitalid),
        family = stats::binomial(), data = current, nAGQ = 1,
        control = lme4::glmerControl(
          optimizer = "bobyqa", optCtrl = list(maxfun = 300000), calc.derivs = TRUE
        )
      ),
      warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) e
  )
  if (inherits(fit, "error")) stop("eICU model failed: ", conditionMessage(fit), call. = FALSE)
  messages <- fit@optinfo$conv$lme4$messages
  if (is.null(messages)) messages <- character()
  comparisons <- if (k == 2L) c("P1") else c("P1", "P2", "P4")
  list(
    effects = extract_effects(
      fit, "eICU", k, "APACHE-adjusted hospital random-intercept logistic",
      nrow(current), sum(current$hospital_mortality == 1L), "OR", comparisons
    ),
    diagnostics = tibble::tibble(
      fit_ok = TRUE,
      convergence_message = paste(messages, collapse = " | "),
      captured_warnings = paste(unique(warnings), collapse = " | "),
      singular = lme4::isSingular(fit, tol = 1e-5)
    )
  )
}

agreement_metrics <- function(primary, alternative, dataset, k) {
  primary <- as.character(primary)
  alternative <- as.character(alternative)
  primary_p1 <- primary == "P1"
  alternative_p1 <- alternative == "P1"
  intersection <- sum(primary_p1 & alternative_p1)
  union <- sum(primary_p1 | alternative_p1)
  tibble::tibble(
    dataset = dataset,
    k = k,
    n_clustered = length(primary),
    adjusted_rand_index_vs_k3 = mclust::adjustedRandIndex(primary, alternative),
    exact_mapped_label_agreement = mean(primary == alternative),
    primary_k3_p1_n = sum(primary_p1),
    resolution_p1_n = sum(alternative_p1),
    p1_recall_vs_k3 = intersection / sum(primary_p1),
    p1_precision_vs_k3 = intersection / sum(alternative_p1),
    p1_jaccard_vs_k3 = intersection / union
  )
}

biological_direction <- function(profile, dataset, k, p1_effect) {
  p1 <- profile |> dplyr::filter(phenotype == "P1")
  p3 <- profile |> dplyr::filter(phenotype == "P3")
  tibble::tibble(
    dataset = dataset,
    k = k,
    nlr_higher = p1$nlr > p3$nlr,
    sii_higher = p1$sii > p3$sii,
    haemoglobin_lower = p1$haemoglobin < p3$haemoglobin,
    albumin_lower = p1$albumin < p3$albumin,
    bmi_lower = p1$bmi < p3$bmi,
    creatinine_higher = p1$creatinine > p3$creatinine,
    core_biological_direction_recovered =
      p1$nlr > p3$nlr && p1$sii > p3$sii &&
      p1$haemoglobin < p3$haemoglobin && p1$albumin < p3$albumin &&
      p1$creatinine > p3$creatinine,
    p1_effect_measure = p1_effect$effect_measure,
    p1_estimate = p1_effect$estimate,
    p1_lower_95 = p1_effect$lower_95,
    p1_upper_95 = p1_effect$upper_95,
    p1_p_value = p1_effect$p_value,
    mortality_direction_above_one = p1_effect$estimate > 1,
    lower_95_above_one = p1_effect$lower_95 > 1
  )
}

sha256_file <- function(path) digest::digest(file = path, algo = "sha256", serialize = FALSE)

cat("Loading locked unified-albumin cohorts...\n")
nhanes_cluster <- readRDS(paths$nhanes_cluster) |>
  dplyr::transmute(
    id = SEQN, nlr = NLR, sii = SII, haemoglobin = LBXHGB, albumin = LBXSAL,
    bmi = BMXBMI, creatinine = LBXSCR,
    primary_phenotype = factor(phenotype_albumin, levels = c("P3", "P2", "P1"))
  )
nhanes_model <- readRDS(paths$nhanes_model)$model_data

mimic_source <- readRDS(paths$mimic)$analysis
mimic_cluster <- mimic_source |>
  dplyr::transmute(
    id = stay_id, nlr = nlr, sii = sii, haemoglobin = haemoglobin, albumin = albumin,
    bmi = bmi, creatinine = creatinine,
    primary_phenotype = factor(phenotype, levels = c("P3", "P2", "P1"))
  )

eicu_denominator <- data.table::fread(
  paths$eicu_denominator, data.table = FALSE, showProgress = FALSE
) |>
  tibble::as_tibble()
eicu_cluster <- eicu_denominator |>
  dplyr::filter(
    dplyr::if_all(
      dplyr::all_of(c("nlr", "sii_like", "haemoglobin", "albumin", "bmi", "creatinine")),
      ~ !is.na(.x) & is.finite(.x)
    )
  ) |>
  dplyr::filter(
    nlr > 0, nlr <= 100, sii_like > 0, sii_like <= 200000,
    albumin > 0, creatinine > 0
  ) |>
  dplyr::transmute(
    id = patientunitstayid, nlr = nlr, sii = sii_like, haemoglobin = haemoglobin,
    albumin = albumin, bmi = bmi, creatinine = creatinine
  )
set.seed(20260710)
eicu_primary_fit <- stats::kmeans(
  prepare_matrix(eicu_cluster), centers = 3, nstart = 100, iter.max = 500, algorithm = "Lloyd"
)
eicu_cluster$primary_phenotype <- label_primary_k3(eicu_cluster, eicu_primary_fit$cluster)

apache <- readr::read_csv(
  paths$eicu_apache,
  col_select = c(patientunitstayid, apacheversion, apachescore),
  show_col_types = FALSE, progress = FALSE
) |>
  dplyr::filter(apacheversion == "IVa") |>
  dplyr::transmute(
    patientunitstayid,
    apachescore = ifelse(as.numeric(apachescore) < 0, NA_real_, as.numeric(apachescore))
  ) |>
  dplyr::distinct(patientunitstayid, .keep_all = TRUE)
eicu_model_base <- eicu_denominator |>
  dplyr::transmute(
    id = patientunitstayid, hospitalid = hospitalid, age_num = age_num,
    male = as.integer(gender_model == "Male"), hospital_mortality = hospital_mortality
  ) |>
  dplyr::left_join(apache, by = c("id" = "patientunitstayid"))

cohorts <- list(
  NHANES = list(cluster = nhanes_cluster, model = nhanes_model),
  `MIMIC-IV` = list(cluster = mimic_cluster, model = mimic_source),
  eICU = list(cluster = eicu_cluster, model = eicu_model_base)
)
expected_primary_counts <- list(
  NHANES = c(P1 = 1152L, P2 = 1689L, P3 = 1796L),
  `MIMIC-IV` = c(P1 = 442L, P2 = 257L, P3 = 401L),
  eICU = c(P1 = 5737L, P2 = 4337L, P3 = 5168L)
)
for (dataset in names(cohorts)) {
  observed <- table(factor(cohorts[[dataset]]$cluster$primary_phenotype, levels = c("P1", "P2", "P3")))
  if (!identical(as.integer(observed), as.integer(expected_primary_counts[[dataset]]))) {
    stop(dataset, " primary K=3 reconstruction mismatch.", call. = FALSE)
  }
}

profile_rows <- list()
effect_rows <- list()
agreement_rows <- list()
direction_rows <- list()
algorithm_rows <- list()
model_diagnostic_rows <- list()

for (dataset in names(cohorts)) {
  cluster_data <- cohorts[[dataset]]$cluster
  for (k in c(2L, 4L)) {
    cat("Running ", dataset, " at K=", k, "...\n", sep = "")
    result <- run_resolution(cluster_data, k)
    assignment <- tibble::tibble(
      id = cluster_data$id,
      phenotype_resolution = as.character(result$phenotype)
    )
    if (dataset == "NHANES") {
      model_result <- fit_nhanes(
        assignment |> dplyr::rename(SEQN = id),
        cohorts[[dataset]]$model, k, result$level_order
      )
    } else if (dataset == "MIMIC-IV") {
      model_result <- fit_mimic(
        cohorts[[dataset]]$model |>
          dplyr::left_join(assignment, by = c("stay_id" = "id")),
        k, result$level_order
      )
    } else {
      model_result <- fit_eicu(
        cohorts[[dataset]]$model |>
          dplyr::left_join(assignment, by = "id"),
        k, result$level_order
      )
    }
    p1_effect <- model_result$effects |>
      dplyr::filter(comparison == "P1 vs P3")
    key <- paste(dataset, k, sep = "__")
    profile_rows[[key]] <- result$profile |>
      dplyr::mutate(dataset = dataset, k = k, .before = 1)
    effect_rows[[key]] <- model_result$effects
    agreement_rows[[key]] <- agreement_metrics(
      cluster_data$primary_phenotype, result$phenotype, dataset, k
    )
    direction_rows[[key]] <- biological_direction(
      result$profile, dataset, k, p1_effect
    )
    algorithm_rows[[key]] <- result$diagnostics |>
      dplyr::mutate(dataset = dataset, n_clustered = nrow(cluster_data), .before = 1)
    model_diagnostic_rows[[key]] <- model_result$diagnostics |>
      dplyr::mutate(dataset = dataset, k = k, .before = 1)
  }
}

profiles <- dplyr::bind_rows(profile_rows) |>
  dplyr::select(
    dataset, k, phenotype, cluster_raw, n, nlr, sii, haemoglobin, albumin,
    bmi, creatinine, p1_score, reserve_deficit_score, favorable_reference_score
  )
outcome_models <- dplyr::bind_rows(effect_rows)
agreement <- dplyr::bind_rows(agreement_rows)
direction <- dplyr::bind_rows(direction_rows)
algorithm_diagnostics <- dplyr::bind_rows(algorithm_rows)
model_diagnostics <- dplyr::bind_rows(model_diagnostic_rows)

method_details <- tibble::tibble(
  k = c(2L, 4L),
  features_used = "NLR, SII/SII-like, haemoglobin, albumin, BMI, creatinine",
  preprocessing = "1st/99th winsorisation; log NLR/SII/creatinine; within-database Z scores",
  clustering = c(
    "K-means Lloyd; K=2; nstart=100; iter.max=500",
    "K-means Lloyd; K=4; nstart=100; iter.max=500"
  ),
  label_rule = c(
    "P1=max inflammation+creatinine score; other cluster=P3; no P2",
    "P1=max inflammation+creatinine; P2=max reserve deficit among remaining; P3=max favorable-reference score among remaining; residual=P4"
  ),
  seed = 20260710L,
  outcome_used_for_clustering_or_labels = FALSE,
  timing = "Added after primary K=3 phenotype and outcome results were known"
)

qa <- tibble::tribble(
  ~check, ~passed, ~detail,
  "Six database-resolution combinations completed", nrow(agreement) == 6L, paste0("rows=", nrow(agreement)),
  "Expected aggregate profile rows completed", nrow(profiles) == 18L, paste0("rows=", nrow(profiles)),
  "Expected outcome contrast rows completed", nrow(outcome_models) == 12L, paste0("rows=", nrow(outcome_models)),
  "Six clustering diagnostic rows completed", nrow(algorithm_diagnostics) == 6L, paste0("rows=", nrow(algorithm_diagnostics)),
  "Primary K=3 cohort counts reproduced", TRUE, "NHANES 4637; MIMIC-IV 1100; eICU 15242",
  "Model counts and events reproduced", all(
    outcome_models$n[outcome_models$dataset == "NHANES"] == 3979L &
      outcome_models$events[outcome_models$dataset == "NHANES"] == 720L
  ) && all(
    outcome_models$n[outcome_models$dataset == "MIMIC-IV"] == 1100L &
      outcome_models$events[outcome_models$dataset == "MIMIC-IV"] == 550L
  ) && all(
    outcome_models$n[outcome_models$dataset == "eICU"] == 13234L &
      outcome_models$events[outcome_models$dataset == "eICU"] == 1990L
  ), "3979/720; 1100/550; 13234/1990",
  "Every combination contains a P1-versus-P3 estimate", nrow(
    outcome_models |> dplyr::filter(comparison == "P1 vs P3")
  ) == 6L, "six P1 core contrasts",
  "All estimates finite", all(is.finite(outcome_models$estimate)), "no missing effect estimate",
  "All outcome models fitted", all(model_diagnostics$fit_ok), "fit_ok TRUE",
  "No eICU model singular", !any(model_diagnostics$singular %in% TRUE, na.rm = TRUE),
    paste0("singular=", sum(model_diagnostics$singular %in% TRUE, na.rm = TRUE)),
  "No outcomes used for clustering or labels", TRUE, "outcomes accessed only after mapping",
  "No participant-level output written", TRUE, "aggregate tables and summary-only RDS"
)
if (!all(qa$passed)) stop("K-resolution QA failed before output write.", call. = FALSE)

tables <- list(
  Table86A_method_details = method_details,
  Table86B_phenotype_profiles = profiles,
  Table86C_outcome_models = outcome_models,
  Table86D_agreement_with_primary_K3 = agreement,
  Table86E_biological_directions = direction,
  Table86F_algorithm_diagnostics = algorithm_diagnostics,
  Table86G_outcome_model_diagnostics = model_diagnostics,
  Table86H_QA = qa
)
for (name in names(tables)) {
  readr::write_csv(tables[[name]], file.path(output_dir, paste0(name, ".csv")))
}
saveRDS(tables, file.path(output_dir, "K_resolution_sensitivity_summary_results.rds"))

source_paths <- required_paths
source_manifest <- tibble::tibble(
  file = basename(source_paths),
  path = normalizePath(source_paths, winslash = "/", mustWork = TRUE),
  bytes = file.info(source_paths)$size,
  sha256 = vapply(source_paths, sha256_file, character(1))
)
readr::write_csv(source_manifest, file.path(output_dir, "SOURCE_SHA256_MANIFEST.csv"))

p1_summary <- outcome_models |>
  dplyr::filter(comparison == "P1 vs P3") |>
  dplyr::select(dataset, k, n, events, effect_measure, estimate, lower_95, upper_95, p_value)
summary_lines <- c(
  "Amendment 19 K-resolution sensitivity results",
  "",
  "Timing: added after primary K=3 phenotype and outcome results were known.",
  "K=3 remains primary; K=2 and K=4 are supportive.",
  "",
  "P1 versus P3 outcome estimates:",
  paste(capture.output(print(p1_summary, n = Inf, width = Inf)), collapse = "\n"),
  "",
  "Biological directions:",
  paste(capture.output(print(direction, n = Inf, width = Inf)), collapse = "\n"),
  "",
  "Agreement with primary K=3:",
  paste(capture.output(print(agreement, n = Inf, width = Inf)), collapse = "\n"),
  "",
  "Interpretation boundary: these analyses assess whether the adverse P1-versus-P3 direction is visible at alternative K values. They do not identify a unique optimal K or establish invariant individual membership."
)
writeLines(summary_lines, file.path(output_dir, "K_resolution_sensitivity_summary.txt"))

output_files <- list.files(output_dir, full.names = TRUE, recursive = FALSE)
output_files <- output_files[basename(output_files) != "OUTPUT_SHA256_MANIFEST.csv"]
output_manifest <- tibble::tibble(
  file = basename(output_files),
  bytes = file.info(output_files)$size,
  sha256 = vapply(output_files, sha256_file, character(1))
)
readr::write_csv(output_manifest, file.path(output_dir, "OUTPUT_SHA256_MANIFEST.csv"))

cat(paste(summary_lines, collapse = "\n"), "\n")
message("Amendment 19 K-resolution sensitivity completed.")
