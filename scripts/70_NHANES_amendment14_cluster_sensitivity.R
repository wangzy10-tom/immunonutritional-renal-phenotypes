# ==============================================================================
# Freeze Amendment 14: final-cohort NHANES K and inflammation-axis sensitivities
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
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

paths <- list(
  phenotype_cohort = file.path(root, "output", "nhanes_covariate_upgrade", "NHANES_2011_2018_covariate_augmented.rds"),
  final_model_object = file.path(root, "output", "nhanes_albumin_benchmarks", "NHANES_albumin_benchmark_results.rds"),
  frozen_assignments = file.path(root, "output", "nhanes_robust_reanalysis", "NHANES_2011_2018_variant_assignments.csv"),
  result_dictionary = file.path(root, "output", "final_analysis_freeze_v2", "UNIQUE_RESULT_DICTIONARY.csv"),
  protocol = file.path(root, "documentation", "ANALYSIS_FREEZE_AMENDMENT_14_PROTOCOL_2026-07-13.md"),
  analysis_script = file.path(root, "scripts", "70_NHANES_amendment14_cluster_sensitivity.R")
)
missing_inputs <- names(paths)[!file.exists(unlist(paths))]
if (length(missing_inputs) > 0L) {
  stop("Missing locked input(s): ", paste(missing_inputs, collapse = ", "), call. = FALSE)
}

cohort <- readRDS(paths$phenotype_cohort)
model_object <- readRDS(paths$final_model_object)
model_data <- model_object$model_data
frozen_assignments <- readr::read_csv(paths$frozen_assignments, show_col_types = FALSE)
dictionary <- readr::read_csv(paths$result_dictionary, show_col_types = FALSE)

required_cohort_columns <- c(
  "SEQN", "NLR", "SII", "LBXHGB", "LBXSTP", "BMXBMI", "LBXSCR",
  "PERMTH_INT", "MORTSTAT", "SDMVPSU", "SDMVSTRA", "WTMEC8YR"
)
required_model_columns <- c(
  "SEQN", "PERMTH_INT", "MORTSTAT", "RIDAGEYR", "male", "race", "INDFMPIR",
  "Comorbidity_Score_Extended", "cycle", "smoking", "hypertension", "diabetes",
  "SDMVPSU", "SDMVSTRA", "WTMEC8YR", "phenotype_total_protein"
)
if (length(setdiff(required_cohort_columns, names(cohort))) > 0L) {
  stop("Phenotype cohort is missing required fields.", call. = FALSE)
}
if (length(setdiff(required_model_columns, names(model_data))) > 0L) {
  stop("Final model cohort is missing required fields.", call. = FALSE)
}
if (nrow(cohort) != 4636L || sum(cohort$MORTSTAT == 1L) != 831L) {
  stop("Locked phenotype cohort size/event gate failed.", call. = FALSE)
}
if (nrow(model_data) != 3979L || sum(model_data$MORTSTAT == 1L) != 720L) {
  stop("Locked final model size/event gate failed.", call. = FALSE)
}
if (dplyr::n_distinct(cohort$SEQN) != nrow(cohort) || dplyr::n_distinct(model_data$SEQN) != nrow(model_data)) {
  stop("SEQN uniqueness gate failed.", call. = FALSE)
}
if (!all(model_data$SEQN %in% cohort$SEQN)) {
  stop("Final model cohort is not nested in phenotype cohort.", call. = FALSE)
}

winsorise <- function(x, lower = 0.01, upper = 0.99) {
  limits <- stats::quantile(x, probs = c(lower, upper), na.rm = TRUE, names = FALSE)
  pmin(pmax(x, limits[1]), limits[2])
}

z_score <- function(x) {
  current_sd <- stats::sd(x)
  if (!is.finite(current_sd) || current_sd == 0) return(rep(0, length(x)))
  (x - mean(x)) / current_sd
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

format_p <- function(p) {
  ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))
}

transformed_features <- cohort |>
  dplyr::transmute(
    NLR = log(winsorise(NLR)),
    SII = log(winsorise(SII)),
    LBXHGB = winsorise(LBXHGB),
    LBXSTP = winsorise(LBXSTP),
    BMXBMI = winsorise(BMXBMI),
    LBXSCR = log(winsorise(LBXSCR))
  )

prepare_matrix <- function(selected) {
  matrix_scaled <- scale(as.data.frame(transformed_features[selected]))
  if (any(!is.finite(matrix_scaled))) stop("Non-finite value in clustering matrix.", call. = FALSE)
  matrix_scaled
}

profile_clusters <- function(raw_cluster) {
  cohort |>
    dplyr::mutate(cluster_raw = raw_cluster) |>
    dplyr::group_by(cluster_raw) |>
    dplyr::summarise(
      n = dplyr::n(),
      nlr_median = stats::median(NLR),
      sii_median = stats::median(SII),
      haemoglobin_median = stats::median(LBXHGB),
      total_protein_median = stats::median(LBXSTP),
      bmi_median = stats::median(BMXBMI),
      creatinine_median = stats::median(LBXSCR),
      .groups = "drop"
    )
}

label_k3 <- function(raw_cluster, variant) {
  profiles <- profile_clusters(raw_cluster)
  p1_components <- list()
  if (variant != "log5_sii_only") {
    p1_components[[length(p1_components) + 1L]] <- z_score(log(profiles$nlr_median))
  }
  if (variant != "log5_nlr_only") {
    p1_components[[length(p1_components) + 1L]] <- z_score(log(profiles$sii_median))
  }
  p1_components[[length(p1_components) + 1L]] <- z_score(log(profiles$creatinine_median))
  p1_score <- Reduce(`+`, p1_components)
  p1_raw <- profiles$cluster_raw[which.max(p1_score)]

  remaining <- profiles |>
    dplyr::filter(cluster_raw != p1_raw)
  p2_score <- -z_score(remaining$haemoglobin_median) -
    z_score(remaining$total_protein_median) -
    z_score(remaining$bmi_median)
  p2_raw <- remaining$cluster_raw[which.max(p2_score)]
  p3_raw <- setdiff(profiles$cluster_raw, c(p1_raw, p2_raw))

  mapping <- tibble::tibble(
    cluster_raw = c(p1_raw, p2_raw, p3_raw),
    phenotype = c("P1", "P2", "P3")
  )
  assigned <- mapping$phenotype[match(raw_cluster, mapping$cluster_raw)]
  list(
    assigned = factor(assigned, levels = c("P3", "P2", "P1")),
    mapping = mapping,
    profiles = dplyr::left_join(profiles, mapping, by = "cluster_raw")
  )
}

label_ranked_k <- function(raw_cluster, k) {
  profiles <- profile_clusters(raw_cluster)
  adverse_matrix <- data.frame(
    high_log_nlr = log(profiles$nlr_median),
    high_log_sii = log(profiles$sii_median),
    low_haemoglobin = -profiles$haemoglobin_median,
    low_total_protein = -profiles$total_protein_median,
    low_bmi = -profiles$bmi_median,
    high_log_creatinine = log(profiles$creatinine_median)
  )
  scaled_adverse <- as.data.frame(lapply(adverse_matrix, z_score))
  vulnerability_score <- rowSums(scaled_adverse)
  p1_like_score <- z_score(log(profiles$nlr_median)) +
    z_score(log(profiles$sii_median)) +
    z_score(log(profiles$creatinine_median))
  order_high_to_low <- order(vulnerability_score, decreasing = TRUE)
  labels <- if (k == 2L) {
    c("Higher vulnerability", "Lower vulnerability")
  } else {
    paste0("V", seq_len(k))
  }
  mapping <- tibble::tibble(
    cluster_raw = profiles$cluster_raw[order_high_to_low],
    phenotype = labels,
    vulnerability_rank = seq_len(k)
  )
  mapping$p1_like <- mapping$cluster_raw == profiles$cluster_raw[which.max(p1_like_score)]
  assigned_character <- mapping$phenotype[match(raw_cluster, mapping$cluster_raw)]
  factor_levels <- if (k == 2L) {
    c("Lower vulnerability", "Higher vulnerability")
  } else {
    paste0("V", k:1)
  }
  list(
    assigned = factor(assigned_character, levels = factor_levels),
    mapping = mapping,
    profiles = profiles |>
      dplyr::mutate(
        vulnerability_score = vulnerability_score,
        p1_like_score = p1_like_score
      ) |>
      dplyr::left_join(mapping, by = "cluster_raw")
  )
}

fit_kmeans <- function(selected, k, seed, label_type, variant = NULL) {
  z <- prepare_matrix(selected)
  set.seed(seed)
  fit <- stats::kmeans(z, centers = k, nstart = 100, iter.max = 500, algorithm = "Lloyd")
  labelled <- if (label_type == "k3") {
    label_k3(fit$cluster, variant)
  } else {
    label_ranked_k(fit$cluster, k)
  }
  list(fit = fit, assigned = labelled$assigned, mapping = labelled$mapping, profiles = labelled$profiles)
}

dictionary_value <- function(result_id, column) {
  value <- dictionary[dictionary$result_id == result_id, column, drop = TRUE]
  if (length(value) != 1L) stop("Dictionary lookup failed for ", result_id, call. = FALSE)
  as.numeric(value)
}

model_data <- model_data |>
  dplyr::mutate(
    phenotype_total_protein = factor(as.character(phenotype_total_protein), levels = c("P3", "P2", "P1"))
  )
model_formula_primary <- survival::Surv(PERMTH_INT, MORTSTAT) ~
  phenotype_total_protein + RIDAGEYR + male + race + INDFMPIR +
  Comorbidity_Score_Extended + cycle + smoking + hypertension + diabetes
primary_design <- survey::svydesign(
  ids = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~WTMEC8YR,
  nest = TRUE, data = model_data
)
primary_model <- survey::svycoxph(model_formula_primary, design = primary_design)

extract_terms <- function(fit, term_prefix, reference, analysis, k, seed) {
  coefficients <- summary(fit)$coefficients
  confidence <- suppressMessages(stats::confint(fit))
  terms <- rownames(coefficients)
  keep <- startsWith(terms, term_prefix)
  p_column <- grep("Pr\\(", colnames(coefficients), value = TRUE)[1]
  comparison_level <- sub(paste0("^", term_prefix), "", terms[keep])
  tibble::tibble(
    analysis = analysis,
    k = k,
    seed = seed,
    n_model = nrow(model_data),
    events = sum(model_data$MORTSTAT == 1L),
    comparison = paste0(comparison_level, " vs ", reference),
    HR = exp(coefficients[keep, "coef"]),
    lower_95 = exp(confidence[keep, 1]),
    upper_95 = exp(confidence[keep, 2]),
    p_value = coefficients[keep, p_column],
    hazard_ratio_95ci = sprintf(
      "%.2f (%.2f-%.2f)",
      exp(coefficients[keep, "coef"]), exp(confidence[keep, 1]), exp(confidence[keep, 2])
    ),
    p_value_formatted = format_p(coefficients[keep, p_column])
  )
}

primary_terms <- extract_terms(primary_model, "phenotype_total_protein", "P3", "Frozen primary K=3 reproduction", 3L, 20260710L)
frozen_p1 <- dictionary_value("NH-PRIMARY-P1", "estimate")
frozen_p2 <- dictionary_value("NH-SECONDARY-P2", "estimate")
primary_p1 <- primary_terms$HR[primary_terms$comparison == "P1 vs P3"]
primary_p2 <- primary_terms$HR[primary_terms$comparison == "P2 vs P3"]

six_features <- c("NLR", "SII", "LBXHGB", "LBXSTP", "BMXBMI", "LBXSCR")
primary_cluster <- fit_kmeans(six_features, 3L, 20260710L, "k3", "log6_winsor")
frozen_primary <- frozen_assignments$phenotype_log6_winsor[match(cohort$SEQN, frozen_assignments$SEQN)]
if (any(is.na(frozen_primary))) stop("Frozen primary assignment join failed.", call. = FALSE)
if (!identical(as.character(primary_cluster$assigned), as.character(frozen_primary))) {
  stop("Mandatory gate failed: K=3 assignments were not reproduced exactly.", call. = FALSE)
}
primary_counts <- table(primary_cluster$assigned)
if (!identical(as.integer(primary_counts[c("P1", "P2", "P3")]), c(817L, 1894L, 1925L))) {
  stop("Mandatory gate failed: frozen K=3 counts were not reproduced.", call. = FALSE)
}
if (abs(primary_p1 - frozen_p1) >= 1e-10 || abs(primary_p2 - frozen_p2) >= 1e-10) {
  stop("Mandatory gate failed: frozen primary Cox estimates were not reproduced.", call. = FALSE)
}

fit_amendment_model <- function(assignment, analysis, k, seed) {
  assignment_table <- tibble::tibble(
    SEQN = cohort$SEQN,
    phenotype_amendment = as.character(assignment)
  )
  joined <- model_data |>
    dplyr::select(-dplyr::any_of("phenotype_amendment")) |>
    dplyr::left_join(assignment_table, by = "SEQN")
  if (any(is.na(joined$phenotype_amendment)) || nrow(joined) != 3979L) {
    stop("Amendment assignment join failed for ", analysis, call. = FALSE)
  }
  reference <- if (k == 2L) "Lower vulnerability" else if (k == 6L) "V6" else "P3"
  levels_current <- if (k == 2L) {
    c("Lower vulnerability", "Higher vulnerability")
  } else if (k == 6L) {
    paste0("V", 6:1)
  } else {
    c("P3", "P2", "P1")
  }
  joined$phenotype_amendment <- factor(joined$phenotype_amendment, levels = levels_current)
  formula <- survival::Surv(PERMTH_INT, MORTSTAT) ~
    phenotype_amendment + RIDAGEYR + male + race + INDFMPIR +
    Comorbidity_Score_Extended + cycle + smoking + hypertension + diabetes
  design <- survey::svydesign(
    ids = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~WTMEC8YR,
    nest = TRUE, data = joined
  )
  fit <- survey::svycoxph(formula, design = design)
  list(
    data = joined,
    fit = fit,
    estimates = extract_terms(fit, "phenotype_amendment", reference, analysis, k, seed)
  )
}

summarise_counts <- function(assignment, analysis, k) {
  descriptive <- cohort |>
    dplyr::mutate(phenotype = factor(as.character(assignment), levels = levels(assignment)))
  design <- survey::svydesign(
    ids = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~WTMEC8YR,
    nest = TRUE, data = descriptive
  )
  weighted <- survey::svyby(~MORTSTAT, ~phenotype, design, survey::svymean, na.rm = TRUE) |>
    as.data.frame() |>
    dplyr::transmute(
      phenotype = as.character(phenotype),
      survey_weighted_mortality_percent = 100 * MORTSTAT,
      survey_weighted_mortality_se_percent = 100 * se
    )
  descriptive |>
    dplyr::group_by(phenotype) |>
    dplyr::summarise(
      n = dplyr::n(),
      percent = 100 * n / nrow(descriptive),
      deaths = sum(MORTSTAT == 1L),
      unweighted_mortality_percent = 100 * mean(MORTSTAT == 1L),
      .groups = "drop"
    ) |>
    dplyr::left_join(weighted, by = "phenotype") |>
    dplyr::mutate(analysis = analysis, k = k, .before = 1)
}

message("Mandatory reproduction gate passed; running locked Amendment 14 analyses...")
k2 <- fit_kmeans(six_features, 2L, 20260712L, "ranked")
k6 <- fit_kmeans(six_features, 6L, 20260716L, "ranked")
nlr_only <- fit_kmeans(c("NLR", "LBXHGB", "LBXSTP", "BMXBMI", "LBXSCR"), 3L, 20260710L, "k3", "log5_nlr_only")
sii_only <- fit_kmeans(c("SII", "LBXHGB", "LBXSTP", "BMXBMI", "LBXSCR"), 3L, 20260710L, "k3", "log5_sii_only")

model_k2 <- fit_amendment_model(k2$assigned, "K=2 six-feature sensitivity", 2L, 20260712L)
model_k6 <- fit_amendment_model(k6$assigned, "K=6 six-feature sensitivity", 6L, 20260716L)
model_nlr <- fit_amendment_model(nlr_only$assigned, "K=3 NLR-only inflammation sensitivity", 3L, 20260710L)
model_sii <- fit_amendment_model(sii_only$assigned, "K=3 SII-only inflammation sensitivity", 3L, 20260710L)

profiles_k <- dplyr::bind_rows(
  k2$profiles |> dplyr::mutate(analysis = "K=2 six-feature sensitivity", k = 2L, seed = 20260712L, .before = 1),
  k6$profiles |> dplyr::mutate(analysis = "K=6 six-feature sensitivity", k = 6L, seed = 20260716L, .before = 1)
)
counts_k <- dplyr::bind_rows(
  summarise_counts(k2$assigned, "K=2 six-feature sensitivity", 2L),
  summarise_counts(k6$assigned, "K=6 six-feature sensitivity", 6L)
)
models_k <- dplyr::bind_rows(model_k2$estimates, model_k6$estimates)

profiles_single <- dplyr::bind_rows(
  nlr_only$profiles |> dplyr::mutate(analysis = "K=3 NLR-only inflammation sensitivity", seed = 20260710L, .before = 1),
  sii_only$profiles |> dplyr::mutate(analysis = "K=3 SII-only inflammation sensitivity", seed = 20260710L, .before = 1)
)
counts_single <- dplyr::bind_rows(
  summarise_counts(nlr_only$assigned, "K=3 NLR-only inflammation sensitivity", 3L),
  summarise_counts(sii_only$assigned, "K=3 SII-only inflammation sensitivity", 3L)
)
models_single <- dplyr::bind_rows(model_nlr$estimates, model_sii$estimates)

assignments <- tibble::tibble(
  SEQN = cohort$SEQN,
  primary_k3 = as.character(primary_cluster$assigned),
  k2 = as.character(k2$assigned),
  k6 = as.character(k6$assigned),
  nlr_only_k3 = as.character(nlr_only$assigned),
  sii_only_k3 = as.character(sii_only$assigned)
)

agreement <- tibble::tibble(
  comparison = c(
    "K=2 six-feature vs primary K=3",
    "K=6 six-feature vs primary K=3",
    "K=3 NLR-only vs primary K=3",
    "K=3 SII-only vs primary K=3",
    "K=3 NLR-only vs K=3 SII-only"
  ),
  adjusted_rand_index = c(
    adjusted_rand_index(assignments$k2, assignments$primary_k3),
    adjusted_rand_index(assignments$k6, assignments$primary_k3),
    adjusted_rand_index(assignments$nlr_only_k3, assignments$primary_k3),
    adjusted_rand_index(assignments$sii_only_k3, assignments$primary_k3),
    adjusted_rand_index(assignments$nlr_only_k3, assignments$sii_only_k3)
  )
)

alignment_long <- dplyr::bind_rows(
  assignments |>
    dplyr::count(primary_k3, alternative = k2, name = "n") |>
    dplyr::mutate(analysis = "K=2 six-feature sensitivity", .before = 1),
  assignments |>
    dplyr::count(primary_k3, alternative = k6, name = "n") |>
    dplyr::mutate(analysis = "K=6 six-feature sensitivity", .before = 1),
  assignments |>
    dplyr::count(primary_k3, alternative = nlr_only_k3, name = "n") |>
    dplyr::mutate(analysis = "K=3 NLR-only inflammation sensitivity", .before = 1),
  assignments |>
    dplyr::count(primary_k3, alternative = sii_only_k3, name = "n") |>
    dplyr::mutate(analysis = "K=3 SII-only inflammation sensitivity", .before = 1)
) |>
  dplyr::group_by(analysis, alternative) |>
  dplyr::mutate(proportion_within_alternative = n / sum(n)) |>
  dplyr::ungroup()

primary_reproduction <- dplyr::bind_rows(
  primary_terms |>
    dplyr::mutate(
      frozen_estimate = dplyr::case_when(
        comparison == "P1 vs P3" ~ frozen_p1,
        comparison == "P2 vs P3" ~ frozen_p2,
        TRUE ~ NA_real_
      ),
      absolute_difference = abs(HR - frozen_estimate)
    )
)

method_lock <- tibble::tribble(
  ~analysis, ~features, ~k, ~seed, ~nstart, ~iter_max, ~labelling, ~model_n, ~events,
  "Frozen primary reproduction", "NLR; SII; haemoglobin; total protein; BMI; creatinine", 3L, 20260710L, 100L, 500L, "Frozen P1/P2/P3 outcome-blind rule", 3979L, 720L,
  "K=2 sensitivity", "NLR; SII; haemoglobin; total protein; BMI; creatinine", 2L, 20260712L, 100L, 500L, "Six-domain vulnerability rank", 3979L, 720L,
  "K=6 sensitivity", "NLR; SII; haemoglobin; total protein; BMI; creatinine", 6L, 20260716L, 100L, 500L, "Six-domain vulnerability rank; P1-like flag locked by inflammation plus creatinine", 3979L, 720L,
  "NLR-only sensitivity", "NLR; haemoglobin; total protein; BMI; creatinine", 3L, 20260710L, 100L, 500L, "Frozen P1/P2/P3 outcome-blind rule", 3979L, 720L,
  "SII-only sensitivity", "SII; haemoglobin; total protein; BMI; creatinine", 3L, 20260710L, 100L, 500L, "Frozen P1/P2/P3 outcome-blind rule", 3979L, 720L
)

qa <- tibble::tribble(
  ~check_id, ~description, ~passed, ~detail,
  "A14-INPUT-FULL", "Phenotype cohort is exactly n=4,636 with 831 deaths", nrow(cohort) == 4636L && sum(cohort$MORTSTAT == 1L) == 831L, paste0(nrow(cohort), "/", sum(cohort$MORTSTAT == 1L)),
  "A14-INPUT-MODEL", "Final model cohort is exactly n=3,979 with 720 deaths", nrow(model_data) == 3979L && sum(model_data$MORTSTAT == 1L) == 720L, paste0(nrow(model_data), "/", sum(model_data$MORTSTAT == 1L)),
  "A14-PRIMARY-ASSIGN", "Primary K=3 assignments reproduce all frozen assignments", identical(as.character(primary_cluster$assigned), as.character(frozen_primary)), paste0(sum(as.character(primary_cluster$assigned) == as.character(frozen_primary)), "/4636"),
  "A14-PRIMARY-COUNTS", "Primary K=3 counts reproduce 817/1,894/1,925", identical(as.integer(primary_counts[c("P1", "P2", "P3")]), c(817L, 1894L, 1925L)), paste(primary_counts[c("P1", "P2", "P3")], collapse = "/"),
  "A14-PRIMARY-P1", "Primary P1 HR reproduces frozen estimate within 1e-10", abs(primary_p1 - frozen_p1) < 1e-10, format(abs(primary_p1 - frozen_p1), scientific = TRUE),
  "A14-PRIMARY-P2", "Primary P2 HR reproduces frozen estimate within 1e-10", abs(primary_p2 - frozen_p2) < 1e-10, format(abs(primary_p2 - frozen_p2), scientific = TRUE),
  "A14-K2-COMPLETE", "K=2 assigns all participants to two nonempty clusters", length(k2$assigned) == 4636L && all(table(k2$assigned) > 0), paste(table(k2$assigned), collapse = "/"),
  "A14-K6-COMPLETE", "K=6 assigns all participants to six nonempty clusters", length(k6$assigned) == 4636L && all(table(k6$assigned) > 0), paste(table(k6$assigned), collapse = "/"),
  "A14-NLR-COMPLETE", "NLR-only K=3 assigns all participants", length(nlr_only$assigned) == 4636L && !any(is.na(nlr_only$assigned)), paste(table(nlr_only$assigned), collapse = "/"),
  "A14-SII-COMPLETE", "SII-only K=3 assigns all participants", length(sii_only$assigned) == 4636L && !any(is.na(sii_only$assigned)), paste(table(sii_only$assigned), collapse = "/"),
  "A14-MODEL-ROWS", "Every sensitivity model uses n=3,979 and 720 events", all(c(models_k$n_model, models_single$n_model) == 3979L) && all(c(models_k$events, models_single$events) == 720L), "All locked model rows/events checked",
  "A14-ESTIMATES", "All sensitivity HRs and confidence limits are finite and positive", all(is.finite(c(models_k$HR, models_k$lower_95, models_k$upper_95, models_single$HR, models_single$lower_95, models_single$upper_95))) && all(c(models_k$HR, models_k$lower_95, models_k$upper_95, models_single$HR, models_single$lower_95, models_single$upper_95) > 0), "All estimates checked",
  "A14-K6-P1LIKE", "Exactly one K=6 cluster is designated P1-like before outcome analysis", sum(k6$mapping$p1_like) == 1L, paste(k6$mapping$phenotype[k6$mapping$p1_like], collapse = "")
)
if (!all(qa$passed)) {
  print(qa)
  stop("Amendment 14 QA failed; no result is admissible.", call. = FALSE)
}

output_files <- c(
  method_lock = file.path(output_dir, "TableA14_0_method_lock.csv"),
  primary_reproduction = file.path(output_dir, "TableA14_1_primary_reproduction.csv"),
  k_profiles = file.path(output_dir, "TableA14_2_K2_K6_profiles.csv"),
  k_counts = file.path(output_dir, "TableA14_3_K2_K6_counts_mortality.csv"),
  k_models = file.path(output_dir, "TableA14_4_K2_K6_survey_cox.csv"),
  single_profiles = file.path(output_dir, "TableA14_5_NLR_SII_single_index_profiles.csv"),
  single_counts = file.path(output_dir, "TableA14_6_NLR_SII_single_index_counts_mortality.csv"),
  single_models = file.path(output_dir, "TableA14_7_NLR_SII_single_index_survey_cox.csv"),
  agreement = file.path(output_dir, "TableA14_8_assignment_agreement.csv"),
  alignment = file.path(output_dir, "TableA14_9_alignment_with_primary_K3.csv"),
  qa = file.path(output_dir, "AMENDMENT14_QA_REPORT.csv"),
  results_rds = file.path(output_dir, "AMENDMENT14_RESULTS.rds"),
  summary = file.path(output_dir, "AMENDMENT14_SUMMARY.txt")
)

readr::write_csv(method_lock, output_files[["method_lock"]])
readr::write_csv(primary_reproduction, output_files[["primary_reproduction"]])
readr::write_csv(profiles_k, output_files[["k_profiles"]])
readr::write_csv(counts_k, output_files[["k_counts"]])
readr::write_csv(models_k, output_files[["k_models"]])
readr::write_csv(profiles_single, output_files[["single_profiles"]])
readr::write_csv(counts_single, output_files[["single_counts"]])
readr::write_csv(models_single, output_files[["single_models"]])
readr::write_csv(agreement, output_files[["agreement"]])
readr::write_csv(alignment_long, output_files[["alignment"]])
readr::write_csv(qa, output_files[["qa"]])

saveRDS(
  list(
    method_lock = method_lock,
    primary_reproduction = primary_reproduction,
    k2 = k2,
    k6 = k6,
    nlr_only = nlr_only,
    sii_only = sii_only,
    models_k = models_k,
    models_single = models_single,
    counts_k = counts_k,
    counts_single = counts_single,
    agreement = agreement,
    alignment = alignment_long,
    assignments = assignments,
    qa = qa
  ),
  output_files[["results_rds"]]
)

p1_like_label <- k6$mapping$phenotype[k6$mapping$p1_like]
summary_lines <- c(
  "Freeze Amendment 14: final-cohort NHANES sensitivity analyses",
  "",
  "Mandatory reproduction gate: PASSED",
  sprintf("Primary K=3 counts: P1=%d; P2=%d; P3=%d", primary_counts[["P1"]], primary_counts[["P2"]], primary_counts[["P3"]]),
  sprintf("Primary P1 HR absolute difference from Freeze V2: %.3e", abs(primary_p1 - frozen_p1)),
  sprintf("Primary P2 HR absolute difference from Freeze V2: %.3e", abs(primary_p2 - frozen_p2)),
  "",
  "K=2 and K=6 survey-weighted Cox models:",
  paste(capture.output(print(models_k)), collapse = "\n"),
  "",
  paste0("Outcome-blind K=6 P1-like cluster: ", p1_like_label),
  "",
  "NLR-only and SII-only survey-weighted Cox models:",
  paste(capture.output(print(models_single)), collapse = "\n"),
  "",
  "Assignment agreement:",
  paste(capture.output(print(agreement)), collapse = "\n"),
  "",
  "Interpretation boundary: sensitivity analyses only; the six-feature K=3 analysis remains primary."
)
writeLines(summary_lines, output_files[["summary"]], useBytes = TRUE)

manifest_input <- tibble::tibble(
  source_file = unname(unlist(paths)),
  role = names(paths),
  bytes = as.numeric(file.info(unname(unlist(paths)))$size),
  md5 = unname(tools::md5sum(unname(unlist(paths))))
)
manifest_output <- tibble::tibble(
  source_file = unname(output_files),
  role = names(output_files),
  bytes = as.numeric(file.info(unname(output_files))$size),
  md5 = unname(tools::md5sum(unname(output_files)))
)
manifest <- dplyr::bind_rows(manifest_input, manifest_output) |>
  dplyr::mutate(source_file = sub(paste0("^", root, "/?"), "", source_file))
readr::write_csv(manifest, file.path(output_dir, "AMENDMENT14_SOURCE_MANIFEST.csv"))

cat(paste(summary_lines, collapse = "\n"), "\n")
message("Amendment 14 analysis completed; all internal QA checks passed.")
