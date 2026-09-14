# ==============================================================================
# Amendment 17: cross-database structural sensitivities
# - paired five-feature K-means variants
# - CLARA K-medoids and Gaussian mixture alternatives at K/G = 3
# ==============================================================================

required_packages <- c(
  "cluster", "mclust", "survival", "survey", "lme4", "dplyr", "readr",
  "tibble", "data.table", "digest"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop("Missing required package(s): ", paste(missing_packages, collapse = ", "), call. = FALSE)
}
options(survey.lonely.psu = "adjust")
suppressPackageStartupMessages(library(mclust))

root <- normalizePath(Sys.getenv("PROJECT_ROOT", unset = getwd()), winslash = "/", mustWork = TRUE)
eicu_dir <- Sys.getenv("EICU_DIR", unset = "")
if (!nzchar(eicu_dir)) stop("Set EICU_DIR to the authorised eICU data directory.", call. = FALSE)
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
  eicu_denominator = Sys.getenv("EICU_DENOMINATOR_CSV", unset = file.path(root, "data", "derived", "eICU_first_ICU_feature_availability_dataset.csv")),
  eicu_apache = file.path(eicu_dir, "apachePatientResult.csv.gz")
)
missing_paths <- names(paths)[!vapply(paths, file.exists, logical(1))]
if (length(missing_paths) > 0L) {
  stop("Missing input(s): ", paste(missing_paths, collapse = ", "), call. = FALSE)
}
output_dir <- file.path(root, "output", "cross_database_structural_sensitivity_2026-08-30")
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

adjusted_rand_index <- function(x, y) {
  mclust::adjustedRandIndex(as.character(x), as.character(y))
}

prepare_matrix <- function(data, include_nlr = TRUE, include_sii = TRUE) {
  transformed <- list()
  if (include_nlr) transformed$log_nlr <- log(pmax(winsorise(data$nlr), .Machine$double.eps))
  if (include_sii) transformed$log_sii <- log(pmax(winsorise(data$sii), .Machine$double.eps))
  transformed$haemoglobin <- winsorise(data$haemoglobin)
  transformed$albumin <- winsorise(data$albumin)
  transformed$bmi <- winsorise(data$bmi)
  transformed$log_creatinine <- log(pmax(winsorise(data$creatinine), .Machine$double.eps))
  output <- scale(as.data.frame(transformed))
  if (any(!is.finite(output))) stop("Non-finite transformed clustering value.", call. = FALSE)
  output
}

label_from_biology <- function(data, raw_cluster, include_nlr = TRUE, include_sii = TRUE) {
  profile <- data |>
    dplyr::mutate(cluster_raw = raw_cluster) |>
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
    dplyr::arrange(cluster_raw)

  inflammation_score <- rep(0, nrow(profile))
  if (include_nlr) inflammation_score <- inflammation_score + safe_z(log(profile$nlr))
  if (include_sii) inflammation_score <- inflammation_score + safe_z(log(profile$sii))
  p1_score <- inflammation_score + safe_z(log(profile$creatinine))
  p1_raw <- profile$cluster_raw[which.max(p1_score)]
  remaining <- profile |> dplyr::filter(cluster_raw != p1_raw)
  p2_score <- -safe_z(remaining$haemoglobin) - safe_z(remaining$albumin) - safe_z(remaining$bmi)
  p2_raw <- remaining$cluster_raw[which.max(p2_score)]
  p3_raw <- setdiff(profile$cluster_raw, c(p1_raw, p2_raw))
  if (length(p3_raw) != 1L) stop("Biological label mapping did not yield one P3 cluster.", call. = FALSE)

  mapping <- tibble::tibble(
    cluster_raw = c(p1_raw, p2_raw, p3_raw),
    phenotype = c("P1", "P2", "P3")
  )
  labels <- mapping$phenotype[match(raw_cluster, mapping$cluster_raw)]
  list(
    phenotype = factor(labels, levels = c("P3", "P2", "P1")),
    profile = profile |>
      dplyr::mutate(p1_score = p1_score) |>
      dplyr::left_join(mapping, by = "cluster_raw") |>
      dplyr::arrange(factor(phenotype, levels = c("P1", "P2", "P3"))),
    mapping = mapping
  )
}

run_variant <- function(data, variant) {
  include_nlr <- variant != "Five-feature: SII only"
  include_sii <- variant != "Five-feature: NLR only"
  x <- prepare_matrix(data, include_nlr, include_sii)
  set.seed(20260710)

  if (variant %in% c(
      "Primary K-means six-feature", "Five-feature: NLR only", "Five-feature: SII only"
    )) {
    fit <- stats::kmeans(x, centers = 3, nstart = 100, iter.max = 500, algorithm = "Lloyd")
    raw_cluster <- as.integer(fit$cluster)
    diagnostics <- tibble::tibble(
      selected_model = "K-means Lloyd, K=3",
      objective = as.numeric(fit$tot.withinss),
      bic = NA_real_, log_likelihood = NA_real_, converged = fit$iter < 500L,
      iterations = fit$iter
    )
  } else if (variant == "CLARA K-medoids") {
    fit <- cluster::clara(
      x, k = 3, metric = "euclidean", stand = FALSE, samples = 100,
      sampsize = min(nrow(x), 1000L), rngR = TRUE, keep.data = FALSE
    )
    raw_cluster <- as.integer(fit$clustering)
    diagnostics <- tibble::tibble(
      selected_model = "CLARA Euclidean, K=3, 100 samples, sample size <=1000",
      objective = as.numeric(fit$objective)[1],
      bic = NA_real_, log_likelihood = NA_real_, converged = TRUE,
      iterations = NA_integer_
    )
  } else if (variant == "Gaussian mixture") {
    fit <- mclust::Mclust(x, G = 3, verbose = FALSE)
    if (is.null(fit$classification)) stop("Gaussian mixture returned no classification.", call. = FALSE)
    raw_cluster <- as.integer(fit$classification)
    diagnostics <- tibble::tibble(
      selected_model = paste0("Mclust G=3, covariance=", fit$modelName),
      objective = NA_real_, bic = as.numeric(fit$bic),
      log_likelihood = as.numeric(fit$loglik), converged = TRUE,
      iterations = NA_integer_
    )
  } else {
    stop("Unsupported variant: ", variant, call. = FALSE)
  }

  labelled <- label_from_biology(data, raw_cluster, include_nlr, include_sii)
  list(
    phenotype = labelled$phenotype,
    profile = labelled$profile,
    diagnostics = diagnostics,
    include_nlr = include_nlr,
    include_sii = include_sii
  )
}

extract_effects <- function(fit, dataset, variant, model, n, events, effect_measure) {
  beta <- if (inherits(fit, "merMod")) lme4::fixef(fit) else stats::coef(fit)
  covariance <- as.matrix(stats::vcov(fit))
  terms <- c("phenotype_sensitivityP1", "phenotype_sensitivityP2")
  if (!all(terms %in% names(beta))) stop("Sensitivity phenotype terms are missing.", call. = FALSE)
  se <- sqrt(diag(covariance))[terms]
  estimate <- exp(beta[terms])
  lower <- exp(beta[terms] - 1.96 * se)
  upper <- exp(beta[terms] + 1.96 * se)
  tibble::tibble(
    dataset = dataset, variant = variant, model = model, n = n, events = events,
    comparison = c("P1 vs P3", "P2 vs P3"), effect_measure = effect_measure,
    estimate = as.numeric(estimate), lower_95 = as.numeric(lower), upper_95 = as.numeric(upper),
    p_value = 2 * stats::pnorm(abs(beta[terms] / se), lower.tail = FALSE),
    effect_95ci = sprintf("%.3f (%.3f-%.3f)", estimate, lower, upper)
  )
}

fit_nhanes <- function(assignments, model_data, variant) {
  current <- model_data |>
    dplyr::select(-dplyr::any_of("phenotype_sensitivity")) |>
    dplyr::left_join(assignments, by = "SEQN") |>
    dplyr::mutate(phenotype_sensitivity = factor(phenotype_sensitivity, levels = c("P3", "P2", "P1")))
  if (nrow(current) != 3979L || anyNA(current$phenotype_sensitivity)) {
    stop("NHANES sensitivity model cohort mismatch.", call. = FALSE)
  }
  design <- survey::svydesign(
    ids = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~WTMEC8YR,
    nest = TRUE, data = current
  )
  fit <- survey::svycoxph(
    survival::Surv(PERMTH_INT, MORTSTAT) ~ phenotype_sensitivity + RIDAGEYR +
      male + race + INDFMPIR + Comorbidity_Score_Extended + cycle + smoking +
      hypertension + diabetes,
    design = design
  )
  list(
    effects = extract_effects(
      fit, "NHANES", variant, "Complex-survey fully adjusted Cox",
      nrow(current), sum(current$MORTSTAT == 1L), "HR"
    ),
    diagnostics = tibble::tibble(
      fit_ok = TRUE, convergence_message = "", captured_warnings = "", singular = NA
    )
  )
}

fit_mimic <- function(data, variant) {
  current <- data |>
    dplyr::mutate(phenotype_sensitivity = factor(phenotype_sensitivity, levels = c("P3", "P2", "P1")))
  if (nrow(current) != 1100L || anyNA(current$phenotype_sensitivity)) {
    stop("MIMIC-IV sensitivity model cohort mismatch.", call. = FALSE)
  }
  fit <- survival::coxph(
    survival::Surv(survival_days_365, mortality_365d) ~
      phenotype_sensitivity + male + strata(oasis_quartile),
    data = current, ties = "efron"
  )
  list(
    effects = extract_effects(
      fit, "MIMIC-IV", variant, "Official OASIS quartile-stratified Cox",
      nrow(current), sum(current$mortality_365d == 1L), "HR"
    ),
    diagnostics = tibble::tibble(
      fit_ok = TRUE, convergence_message = "", captured_warnings = "", singular = NA
    )
  )
}

fit_eicu <- function(data, variant) {
  current <- data |>
    dplyr::filter(
      !is.na(phenotype_sensitivity), is.finite(apachescore), !is.na(hospital_mortality)
    ) |>
    dplyr::mutate(
      phenotype_sensitivity = factor(phenotype_sensitivity, levels = c("P3", "P2", "P1")),
      hospitalid = factor(hospitalid)
    )
  if (nrow(current) != 13234L || sum(current$hospital_mortality == 1L) != 1990L ||
      dplyr::n_distinct(current$hospitalid) != 166L || anyNA(current$phenotype_sensitivity)) {
    stop(
      "eICU sensitivity model cohort mismatch: n=", nrow(current),
      ", events=", sum(current$hospital_mortality == 1L),
      ", hospitals=", dplyr::n_distinct(current$hospitalid),
      ", missing labels=", sum(is.na(current$phenotype_sensitivity)), ".",
      call. = FALSE
    )
  }
  warnings <- character()
  fit <- tryCatch(
    withCallingHandlers(
      lme4::glmer(
        hospital_mortality ~ phenotype_sensitivity + age_num + male + apachescore + (1 | hospitalid),
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
  list(
    effects = extract_effects(
      fit, "eICU", variant, "APACHE-adjusted hospital random-intercept logistic",
      nrow(current), sum(current$hospital_mortality == 1L), "OR"
    ),
    diagnostics = tibble::tibble(
      fit_ok = TRUE, convergence_message = paste(messages, collapse = " | "),
      captured_warnings = paste(unique(warnings), collapse = " | "),
      singular = lme4::isSingular(fit, tol = 1e-5)
    )
  )
}

agreement_metrics <- function(primary, alternative, dataset, variant) {
  primary <- factor(as.character(primary), levels = c("P3", "P2", "P1"))
  alternative <- factor(as.character(alternative), levels = c("P3", "P2", "P1"))
  primary_p1 <- primary == "P1"
  alternative_p1 <- alternative == "P1"
  intersection <- sum(primary_p1 & alternative_p1)
  union <- sum(primary_p1 | alternative_p1)
  tibble::tibble(
    dataset = dataset, variant = variant, n_clustered = length(primary),
    adjusted_rand_index = adjusted_rand_index(primary, alternative),
    exact_label_agreement = mean(primary == alternative),
    primary_p1_n = sum(primary_p1), sensitivity_p1_n = sum(alternative_p1),
    p1_recall = intersection / sum(primary_p1),
    p1_precision = intersection / sum(alternative_p1),
    p1_jaccard = intersection / union
  )
}

profile_output <- function(profile, dataset, variant, data, phenotype) {
  outcome <- if (dataset == "NHANES") data$event else data$event
  profile |>
    dplyr::mutate(
      dataset = dataset, variant = variant,
      events_available = vapply(
        cluster_raw,
        function(k) sum(outcome[as.integer(factor(phenotype, levels = c("P3", "P2", "P1"))) ==
          as.integer(factor(profile$phenotype[profile$cluster_raw == k], levels = c("P3", "P2", "P1")))] == 1, na.rm = TRUE),
        integer(1)
      )
    ) |>
    dplyr::select(
      dataset, variant, phenotype, cluster_raw, n, events_available,
      nlr, sii, haemoglobin, albumin, bmi, creatinine, p1_score
    )
}

signature_output <- function(profile, dataset, variant, include_nlr, include_sii, p1_effect) {
  p1 <- profile |> dplyr::filter(phenotype == "P1")
  p3 <- profile |> dplyr::filter(phenotype == "P3")
  nlr_higher <- p1$nlr > p3$nlr
  sii_higher <- p1$sii > p3$sii
  tibble::tibble(
    dataset = dataset, variant = variant,
    retained_inflammation = if (include_nlr && include_sii) "NLR and SII/SII-like" else if (include_nlr) "NLR" else "SII/SII-like",
    nlr_higher = nlr_higher, sii_higher = sii_higher,
    retained_inflammation_higher = (include_nlr && nlr_higher) || (include_sii && sii_higher),
    haemoglobin_lower = p1$haemoglobin < p3$haemoglobin,
    albumin_lower = p1$albumin < p3$albumin,
    bmi_lower = p1$bmi < p3$bmi,
    creatinine_higher = p1$creatinine > p3$creatinine,
    p1_effect_measure = p1_effect$effect_measure,
    p1_estimate = p1_effect$estimate,
    p1_lower_95 = p1_effect$lower_95,
    p1_upper_95 = p1_effect$upper_95,
    p1_p_value = p1_effect$p_value,
    mortality_direction_above_one = p1_effect$estimate > 1
  )
}

sha256_file <- function(path) digest::digest(file = path, algo = "sha256", serialize = FALSE)

cat("Loading locked albumin cohorts...\n")
nhanes_cluster <- readRDS(paths$nhanes_cluster) |>
  dplyr::transmute(
    id = SEQN, nlr = NLR, sii = SII, haemoglobin = LBXHGB, albumin = LBXSAL,
    bmi = BMXBMI, creatinine = LBXSCR,
    primary_phenotype = factor(phenotype_albumin, levels = c("P3", "P2", "P1")),
    event = MORTSTAT
  )
nhanes_model <- readRDS(paths$nhanes_model)$model_data

mimic_source <- readRDS(paths$mimic)$analysis
mimic_cluster <- mimic_source |>
  dplyr::transmute(
    id = stay_id, nlr = nlr, sii = sii, haemoglobin = haemoglobin, albumin = albumin,
    bmi = bmi, creatinine = creatinine,
    primary_phenotype = factor(phenotype, levels = c("P3", "P2", "P1")),
    event = mortality_365d
  )

eicu_denominator <- data.table::fread(paths$eicu_denominator, data.table = FALSE, showProgress = FALSE) |>
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
    albumin = albumin, bmi = bmi, creatinine = creatinine,
    event = hospital_mortality
  )
primary_eicu <- run_variant(eicu_cluster, "Primary K-means six-feature")
eicu_cluster$primary_phenotype <- primary_eicu$phenotype
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

expected_primary_counts <- list(
  NHANES = c(P1 = 1152L, P2 = 1689L, P3 = 1796L),
  `MIMIC-IV` = c(P1 = 442L, P2 = 257L, P3 = 401L),
  eICU = c(P1 = 5737L, P2 = 4337L, P3 = 5168L)
)
cohorts <- list(
  NHANES = list(cluster = nhanes_cluster, model = nhanes_model),
  `MIMIC-IV` = list(cluster = mimic_cluster, model = mimic_source),
  eICU = list(cluster = eicu_cluster, model = eicu_model_base)
)
for (dataset in names(cohorts)) {
  observed <- table(factor(cohorts[[dataset]]$cluster$primary_phenotype, levels = c("P1", "P2", "P3")))
  if (!identical(as.integer(observed), as.integer(expected_primary_counts[[dataset]]))) {
    stop(dataset, " primary phenotype reconstruction mismatch.", call. = FALSE)
  }
}

variants <- c(
  "Five-feature: NLR only", "Five-feature: SII only",
  "CLARA K-medoids", "Gaussian mixture"
)
agreement_rows <- list()
profile_rows <- list()
effect_rows <- list()
signature_rows <- list()
diagnostic_rows <- list()
model_diagnostic_rows <- list()

for (dataset in names(cohorts)) {
  cluster_data <- cohorts[[dataset]]$cluster
  for (variant in variants) {
    cat("Running ", dataset, " — ", variant, "...\n", sep = "")
    result <- run_variant(cluster_data, variant)
    assignment <- tibble::tibble(
      id = cluster_data$id,
      phenotype_sensitivity = as.character(result$phenotype)
    )

    if (dataset == "NHANES") {
      model_result <- fit_nhanes(
        assignment |> dplyr::rename(SEQN = id), cohorts[[dataset]]$model, variant
      )
    } else if (dataset == "MIMIC-IV") {
      model_result <- fit_mimic(
        cohorts[[dataset]]$model |>
          dplyr::left_join(assignment, by = c("stay_id" = "id")),
        variant
      )
    } else {
      model_result <- fit_eicu(
        cohorts[[dataset]]$model |>
          dplyr::left_join(assignment, by = "id"),
        variant
      )
    }

    p1_effect <- model_result$effects |> dplyr::filter(comparison == "P1 vs P3")
    key <- paste(dataset, variant, sep = "__")
    agreement_rows[[key]] <- agreement_metrics(
      cluster_data$primary_phenotype, result$phenotype, dataset, variant
    )
    profile_rows[[key]] <- result$profile |>
      dplyr::mutate(dataset = dataset, variant = variant, .before = 1)
    effect_rows[[key]] <- model_result$effects
    signature_rows[[key]] <- signature_output(
      result$profile, dataset, variant, result$include_nlr, result$include_sii, p1_effect
    )
    diagnostic_rows[[key]] <- result$diagnostics |>
      dplyr::mutate(
        dataset = dataset, variant = variant, n_clustered = nrow(cluster_data),
        .before = 1
      )
    model_diagnostic_rows[[key]] <- model_result$diagnostics |>
      dplyr::mutate(dataset = dataset, variant = variant, .before = 1)
  }
}

agreement <- dplyr::bind_rows(agreement_rows)
profiles <- dplyr::bind_rows(profile_rows) |>
  dplyr::select(
    dataset, variant, phenotype, cluster_raw, n, nlr, sii, haemoglobin,
    albumin, bmi, creatinine, p1_score
  )
outcome_models <- dplyr::bind_rows(effect_rows)
signature <- dplyr::bind_rows(signature_rows)
algorithm_diagnostics <- dplyr::bind_rows(diagnostic_rows)
model_diagnostics <- dplyr::bind_rows(model_diagnostic_rows)

method_details <- tibble::tibble(
  variant = variants,
  features_used = c(
    "NLR, haemoglobin, albumin, BMI, creatinine",
    "SII/SII-like, haemoglobin, albumin, BMI, creatinine",
    "NLR, SII/SII-like, haemoglobin, albumin, BMI, creatinine",
    "NLR, SII/SII-like, haemoglobin, albumin, BMI, creatinine"
  ),
  method = c(
    "K-means Lloyd; K=3; nstart=100; iter.max=500",
    "K-means Lloyd; K=3; nstart=100; iter.max=500",
    "CLARA; K=3; Euclidean; 100 samples; sample size up to 1000",
    "mclust::Mclust; G=3; covariance model selected by BIC"
  ),
  seed = 20260710L,
  outcome_used_for_clustering_or_labels = FALSE
)

qa <- tibble::tribble(
  ~check, ~passed, ~detail,
  "Twelve database-variant combinations completed", nrow(agreement) == 12L, paste0("rows=", nrow(agreement)),
  "Twenty-four P1/P2 outcome rows completed", nrow(outcome_models) == 24L, paste0("rows=", nrow(outcome_models)),
  "All combinations contain P1/P2/P3", all(vapply(
    split(profiles, interaction(profiles$dataset, profiles$variant, drop = TRUE)),
    function(x) setequal(x$phenotype, c("P1", "P2", "P3")), logical(1)
  )), "Three labels in every combination",
  "Primary cohort counts reproduced", TRUE, "NHANES 4637; MIMIC-IV 1100; eICU 15242",
  "Model counts reproduced", all(
    outcome_models$n[outcome_models$dataset == "NHANES"] == 3979L &
      outcome_models$events[outcome_models$dataset == "NHANES"] == 720L
  ) && all(
    outcome_models$n[outcome_models$dataset == "MIMIC-IV"] == 1100L &
      outcome_models$events[outcome_models$dataset == "MIMIC-IV"] == 550L
  ) && all(
    outcome_models$n[outcome_models$dataset == "eICU"] == 13234L &
      outcome_models$events[outcome_models$dataset == "eICU"] == 1990L
  ), "3979/720; 1100/550; 13234/1990",
  "All effect estimates finite", all(is.finite(outcome_models$estimate)), "No missing P1/P2 estimate",
  "All fitted outcome models succeeded", all(model_diagnostics$fit_ok), "fit_ok TRUE",
  "No eICU outcome model singular", !any(model_diagnostics$singular %in% TRUE, na.rm = TRUE),
    paste0("singular=", sum(model_diagnostics$singular %in% TRUE, na.rm = TRUE)),
  "No outcomes used in clustering or labels", TRUE, "Outcomes accessed after labels were frozen",
  "No participant-level output written", TRUE, "Only aggregate tables and summary RDS are saved"
)
if (!all(qa$passed)) stop("Structural-sensitivity QA failed before output write.", call. = FALSE)

tables <- list(
  Table83A_method_details = method_details,
  Table83B_assignment_agreement = agreement,
  Table83C_phenotype_profiles = profiles,
  Table83D_outcome_models = outcome_models,
  Table83E_biological_directions = signature,
  Table83F_algorithm_diagnostics = algorithm_diagnostics,
  Table83G_outcome_model_diagnostics = model_diagnostics,
  Table83H_QA = qa
)
for (name in names(tables)) {
  readr::write_csv(tables[[name]], file.path(output_dir, paste0(name, ".csv")))
}
saveRDS(tables, file.path(output_dir, "cross_database_structural_sensitivity_summary_results.rds"))

source_paths <- c(unname(unlist(paths)), file.path(root, "scripts", "04_cross_database_structural_sensitivity.R"))
source_manifest <- tibble::tibble(
  file = basename(source_paths),
  bytes = file.info(source_paths)$size,
  sha256 = vapply(source_paths, sha256_file, character(1))
)
readr::write_csv(source_manifest, file.path(output_dir, "SOURCE_SHA256_MANIFEST.csv"))

summary_lines <- c(
  "Amendment 17 cross-database structural sensitivities",
  "",
  paste0("Completed combinations: ", nrow(agreement), "."),
  paste0("Outcome rows: ", nrow(outcome_models), "."),
  "",
  "Assignment agreement:",
  paste(capture.output(print(agreement, n = Inf, width = Inf)), collapse = "\n"),
  "",
  "Outcome models:",
  paste(capture.output(print(outcome_models, n = Inf, width = Inf)), collapse = "\n"),
  "",
  "Biological directions:",
  paste(capture.output(print(signature, n = Inf, width = Inf)), collapse = "\n"),
  "",
  "Interpretation boundary: five-feature analyses test dependence on joint NLR/SII representation. Alternative algorithms test method dependence and cannot establish algorithm independence. No cross-database effect pooling is permitted."
)
writeLines(summary_lines, file.path(output_dir, "cross_database_structural_sensitivity_summary.txt"))

output_files <- list.files(output_dir, full.names = TRUE, recursive = FALSE)
output_files <- output_files[basename(output_files) != "OUTPUT_SHA256_MANIFEST.csv"]
output_manifest <- tibble::tibble(
  file = basename(output_files), bytes = file.info(output_files)$size,
  sha256 = vapply(output_files, sha256_file, character(1))
)
readr::write_csv(output_manifest, file.path(output_dir, "OUTPUT_SHA256_MANIFEST.csv"))

cat(paste(summary_lines, collapse = "\n"), "\n")
message("Amendment 17 cross-database structural sensitivities completed.")
