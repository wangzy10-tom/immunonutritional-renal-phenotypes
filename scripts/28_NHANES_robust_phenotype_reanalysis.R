# ==============================================================================
# Robust NHANES 2011-2018 phenotype reanalysis
#
# Objectives:
#   1. Reproduce the original raw-z K-means approach in the corrected cohort.
#   2. Evaluate log/winsorised preprocessing for skewed variables.
#   3. Test NLR/SII redundancy and creatinine feature dependence.
#   4. Fit unweighted and complex-survey Cox models without outcome-based labels.
# ==============================================================================

required_packages <- c("dplyr", "readr", "tibble", "survival", "survey", "cluster")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0L) {
  stop("Missing required package(s): ", paste(missing_packages, collapse = ", "), call. = FALSE)
}
options(survey.lonely.psu = "adjust")

project_root <- normalizePath(
  Sys.getenv("PROJECT_ROOT", unset = getwd()), winslash = "/", mustWork = TRUE
)
input_path <- file.path(
  project_root, "output", "nhanes_2011_2018_rebuild",
  "NHANES_2011_2018_rebuilt_analytical_cohort.rds"
)
output_dir <- file.path(project_root, "output", "nhanes_robust_reanalysis")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(input_path)) {
  stop("Rebuilt NHANES cohort not found. Run script 27 first.", call. = FALSE)
}

dat <- readRDS(input_path)

required_columns <- c(
  "SEQN", "Cycle_ID", "RIDAGEYR", "RIAGENDR", "RIDRETH3", "DMDEDUC2",
  "DMDMARTL", "INDFMPIR", "WTMEC8YR", "SDMVPSU", "SDMVSTRA",
  "NLR", "SII", "LBXHGB", "LBXSTP", "BMXBMI", "LBXSCR",
  "PERMTH_INT", "MORTSTAT", "Comorbidity_Score_Extended"
)
missing_columns <- setdiff(required_columns, names(dat))
if (length(missing_columns) > 0L) {
  stop("Missing required column(s): ", paste(missing_columns, collapse = ", "), call. = FALSE)
}

winsorise <- function(x, lower = 0.01, upper = 0.99) {
  limits <- stats::quantile(x, probs = c(lower, upper), na.rm = TRUE, names = FALSE)
  pmin(pmax(x, limits[1]), limits[2])
}

z_score <- function(x) {
  as.numeric(scale(x))
}

adjusted_rand_index <- function(x, y) {
  tab <- table(x, y)
  choose2 <- function(v) v * (v - 1) / 2
  n <- sum(tab)
  if (n < 2) return(NA_real_)
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
  ifelse(is.na(p), NA_character_, ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
}

prepare_matrix <- function(data, variant) {
  transformed <- data |>
    dplyr::transmute(
      NLR = NLR,
      SII = SII,
      LBXHGB = LBXHGB,
      LBXSTP = LBXSTP,
      BMXBMI = BMXBMI,
      LBXSCR = LBXSCR
    )

  if (variant != "raw6_original") {
    transformed <- transformed |>
      dplyr::mutate(
        NLR = log(winsorise(NLR)),
        SII = log(winsorise(SII)),
        LBXHGB = winsorise(LBXHGB),
        LBXSTP = winsorise(LBXSTP),
        BMXBMI = winsorise(BMXBMI),
        LBXSCR = log(winsorise(LBXSCR))
      )
  }

  selected <- switch(
    variant,
    raw6_original = c("NLR", "SII", "LBXHGB", "LBXSTP", "BMXBMI", "LBXSCR"),
    log6_winsor = c("NLR", "SII", "LBXHGB", "LBXSTP", "BMXBMI", "LBXSCR"),
    log5_nlr_only = c("NLR", "LBXHGB", "LBXSTP", "BMXBMI", "LBXSCR"),
    log5_sii_only = c("SII", "LBXHGB", "LBXSTP", "BMXBMI", "LBXSCR"),
    log5_no_creatinine = c("NLR", "SII", "LBXHGB", "LBXSTP", "BMXBMI"),
    stop("Unknown variant: ", variant, call. = FALSE)
  )

  matrix_scaled <- scale(as.data.frame(transformed[selected]))
  attr(matrix_scaled, "selected_features") <- selected
  matrix_scaled
}

label_clusters <- function(data, raw_cluster, variant) {
  profiles <- data |>
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

  p1_components <- list()
  if (variant != "log5_sii_only") {
    p1_components[[length(p1_components) + 1L]] <- z_score(log(profiles$nlr_median))
  }
  if (variant != "log5_nlr_only") {
    p1_components[[length(p1_components) + 1L]] <- z_score(log(profiles$sii_median))
  }
  if (variant != "log5_no_creatinine") {
    p1_components[[length(p1_components) + 1L]] <- z_score(log(profiles$creatinine_median))
  }
  p1_score <- Reduce(`+`, p1_components)
  p1_raw <- profiles$cluster_raw[which.max(p1_score)]

  remaining <- profiles |>
    dplyr::filter(cluster_raw != p1_raw)
  p2_score <- -z_score(remaining$haemoglobin_median) -
    z_score(remaining$total_protein_median) -
    z_score(remaining$bmi_median)
  p2_raw <- remaining$cluster_raw[which.max(p2_score)]
  p3_raw <- setdiff(profiles$cluster_raw, c(p1_raw, p2_raw))

  phenotype <- dplyr::case_when(
    raw_cluster == p1_raw ~ "P1",
    raw_cluster == p2_raw ~ "P2",
    raw_cluster == p3_raw ~ "P3",
    TRUE ~ NA_character_
  )

  list(
    phenotype = factor(phenotype, levels = c("P3", "P2", "P1")),
    raw_profiles = profiles,
    mapping = tibble::tibble(
      cluster_raw = c(p1_raw, p2_raw, p3_raw),
      phenotype = c("P1", "P2", "P3")
    )
  )
}

extract_cox_terms <- function(fit, variant, model, weighted, n_model) {
  co <- summary(fit)$coefficients
  ci <- suppressMessages(stats::confint(fit))
  terms <- rownames(co)
  keep <- terms %in% c("phenotypeP1", "phenotypeP2")
  p_col <- grep("Pr\\(", colnames(co), value = TRUE)[1]

  tibble::tibble(
    variant = variant,
    model = model,
    weighted = weighted,
    n_model = n_model,
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
      "%.2f (%.2f-%.2f)",
      exp(co[keep, "coef"]), exp(ci[keep, 1]), exp(ci[keep, 2])
    ),
    p_value_formatted = format_p(co[keep, p_col])
  )
}

fit_models <- function(data, phenotype, variant) {
  model_data <- data |>
    dplyr::mutate(
      phenotype = stats::relevel(phenotype, ref = "P3"),
      male = as.integer(RIAGENDR == 1),
      race = factor(RIDRETH3),
      education = factor(DMDEDUC2),
      marital = factor(DMDMARTL),
      cycle = factor(Cycle_ID)
    )

  model1_data <- model_data |>
    dplyr::filter(
      !is.na(PERMTH_INT), !is.na(MORTSTAT), !is.na(RIDAGEYR), !is.na(male),
      !is.na(phenotype)
    )
  model2_data <- model1_data |>
    dplyr::filter(!is.na(Comorbidity_Score_Extended))
  model3_data <- model2_data |>
    dplyr::filter(
      !is.na(race), !is.na(education), !is.na(marital), !is.na(INDFMPIR),
      !is.na(cycle)
    )

  fit1 <- survival::coxph(
    survival::Surv(PERMTH_INT, MORTSTAT) ~ phenotype + RIDAGEYR + male,
    data = model1_data,
    ties = "efron"
  )
  fit2 <- survival::coxph(
    survival::Surv(PERMTH_INT, MORTSTAT) ~ phenotype + RIDAGEYR + male +
      Comorbidity_Score_Extended,
    data = model2_data,
    ties = "efron"
  )
  fit3 <- survival::coxph(
    survival::Surv(PERMTH_INT, MORTSTAT) ~ phenotype + RIDAGEYR + male +
      race + education + marital + INDFMPIR + Comorbidity_Score_Extended + cycle,
    data = model3_data,
    ties = "efron"
  )

  survey_design <- survey::svydesign(
    ids = ~SDMVPSU,
    strata = ~SDMVSTRA,
    weights = ~WTMEC8YR,
    nest = TRUE,
    data = model3_data
  )
  survey_candidates <- list(
    "Survey-weighted fully adjusted" =
      survival::Surv(PERMTH_INT, MORTSTAT) ~ phenotype + RIDAGEYR + male +
        race + education + marital + INDFMPIR + Comorbidity_Score_Extended,
    "Survey-weighted age, sex, race, income, comorbidity" =
      survival::Surv(PERMTH_INT, MORTSTAT) ~ phenotype + RIDAGEYR + male +
        race + INDFMPIR + Comorbidity_Score_Extended,
    "Survey-weighted age, sex, comorbidity" =
      survival::Surv(PERMTH_INT, MORTSTAT) ~ phenotype + RIDAGEYR + male +
        Comorbidity_Score_Extended
  )
  fit_survey <- NULL
  survey_model_label <- NULL
  for (candidate_label in names(survey_candidates)) {
    candidate_fit <- tryCatch(
      survey::svycoxph(survey_candidates[[candidate_label]], design = survey_design),
      error = function(e) NULL
    )
    if (!is.null(candidate_fit)) {
      fit_survey <- candidate_fit
      survey_model_label <- candidate_label
      break
    }
  }
  if (is.null(fit_survey)) {
    stop("All survey-weighted Cox candidate models failed for ", variant, call. = FALSE)
  }

  dplyr::bind_rows(
    extract_cox_terms(fit1, variant, "Age and sex", FALSE, nrow(model1_data)),
    extract_cox_terms(fit2, variant, "Age, sex, extended comorbidity", FALSE, nrow(model2_data)),
    extract_cox_terms(
      fit3, variant,
      "Fully adjusted demographics, socioeconomic factors, comorbidity, cycle",
      FALSE, nrow(model3_data)
    ),
    extract_cox_terms(
      fit_survey, variant,
      survey_model_label,
      TRUE, nrow(model3_data)
    )
  )
}

cluster_diagnostics <- function(z, variant, sample_n = 2500L) {
  set.seed(20260710)
  diagnostic_rows <- if (nrow(z) > sample_n) sample(seq_len(nrow(z)), sample_n) else seq_len(nrow(z))
  z_sample <- z[diagnostic_rows, , drop = FALSE]
  distance <- stats::dist(z_sample)

  dplyr::bind_rows(lapply(2:6, function(k) {
    set.seed(20260710 + k)
    fit <- stats::kmeans(
      z_sample, centers = k, nstart = 50, iter.max = 500, algorithm = "Lloyd"
    )
    silhouette_mean <- mean(cluster::silhouette(fit$cluster, distance)[, "sil_width"])
    n <- nrow(z_sample)
    p <- ncol(z_sample)
    between <- fit$betweenss
    within <- fit$tot.withinss
    ch <- (between / (k - 1)) / (within / (n - k))
    tibble::tibble(
      variant = variant,
      k = k,
      diagnostic_n = n,
      p = p,
      total_within_ss = within,
      average_silhouette = silhouette_mean,
      calinski_harabasz = ch
    )
  }))
}

variants <- c(
  "raw6_original",
  "log6_winsor",
  "log5_nlr_only",
  "log5_sii_only",
  "log5_no_creatinine"
)

assignment_data <- dat |>
  dplyr::select(SEQN)
all_profiles <- list()
all_mappings <- list()
all_counts <- list()
all_models <- list()
all_diagnostics <- list()

for (variant in variants) {
  message("Running clustering variant: ", variant)
  z <- prepare_matrix(dat, variant)
  set.seed(20260710)
  fit <- stats::kmeans(z, centers = 3, nstart = 100, iter.max = 500, algorithm = "Lloyd")
  labelled <- label_clusters(dat, fit$cluster, variant)

  assignment_data[[paste0("phenotype_", variant)]] <- as.character(labelled$phenotype)
  all_profiles[[variant]] <- labelled$raw_profiles |>
    dplyr::left_join(labelled$mapping, by = "cluster_raw") |>
    dplyr::mutate(variant = variant)
  all_mappings[[variant]] <- labelled$mapping |>
    dplyr::mutate(variant = variant)
  all_counts[[variant]] <- dat |>
    dplyr::mutate(phenotype = labelled$phenotype) |>
    dplyr::group_by(phenotype) |>
    dplyr::summarise(
      n = dplyr::n(),
      percent = 100 * n / nrow(dat),
      deaths = sum(MORTSTAT == 1),
      mortality_percent = 100 * mean(MORTSTAT == 1),
      median_followup_months = stats::median(PERMTH_INT),
      .groups = "drop"
    ) |>
    dplyr::mutate(variant = variant)
  all_models[[variant]] <- fit_models(dat, labelled$phenotype, variant)
  all_diagnostics[[variant]] <- cluster_diagnostics(z, variant)
}

profiles <- dplyr::bind_rows(all_profiles)
mappings <- dplyr::bind_rows(all_mappings)
counts <- dplyr::bind_rows(all_counts)
model_results <- dplyr::bind_rows(all_models)
diagnostics <- dplyr::bind_rows(all_diagnostics)

reference_assignment <- assignment_data$phenotype_raw6_original
agreement <- dplyr::bind_rows(lapply(variants, function(variant) {
  candidate_labels <- assignment_data[[paste0("phenotype_", variant)]]
  tibble::tibble(
    reference = "raw6_original",
    candidate = variant,
    exact_label_agreement = mean(reference_assignment == candidate_labels),
    adjusted_rand_index = adjusted_rand_index(reference_assignment, candidate_labels)
  )
}))

old_assignment_comparison <- NULL
old_path <- Sys.getenv("LEGACY_NHANES_CLUSTERED_RDATA", unset = "")
if (file.exists(old_path)) {
  old_environment <- new.env(parent = emptyenv())
  load(old_path, envir = old_environment)
  if (exists("nhanes_clustered", envir = old_environment, inherits = FALSE)) {
    old <- get("nhanes_clustered", envir = old_environment) |>
      dplyr::transmute(SEQN, old_phenotype = as.character(Phenotype))
    common <- assignment_data |>
      dplyr::inner_join(old, by = "SEQN")
    old_labels <- common[["old_phenotype"]]
    corrected_labels <- common[["phenotype_raw6_original"]]
    old_vs_corrected_ari <- if (
      !is.null(old_labels) && !is.null(corrected_labels) &&
        length(old_labels) == length(corrected_labels) && length(old_labels) > 1L
    ) {
      adjusted_rand_index(old_labels, corrected_labels)
    } else {
      NA_real_
    }
    old_assignment_comparison <- tibble::tibble(
      common_n = nrow(common),
      old_cycles = "G/H/I",
      new_cycles = "G/H/I/J",
      adjusted_rand_index_unmapped = old_vs_corrected_ari
    )
  }
}

readr::write_csv(diagnostics, file.path(output_dir, "Table28A_cluster_number_diagnostics.csv"))
readr::write_csv(counts, file.path(output_dir, "Table28B_phenotype_counts_and_mortality.csv"))
readr::write_csv(profiles, file.path(output_dir, "Table28C_raw_cluster_profiles.csv"))
readr::write_csv(mappings, file.path(output_dir, "Table28D_cluster_label_mapping.csv"))
readr::write_csv(model_results, file.path(output_dir, "Table28E_Cox_models_all_variants.csv"))
readr::write_csv(agreement, file.path(output_dir, "Table28F_variant_agreement.csv"))
readr::write_csv(assignment_data, file.path(output_dir, "NHANES_2011_2018_variant_assignments.csv"))
if (!is.null(old_assignment_comparison)) {
  readr::write_csv(
    old_assignment_comparison,
    file.path(output_dir, "Table28G_old_vs_corrected_raw_cluster_comparison.csv")
  )
}

saveRDS(
  list(
    assignments = assignment_data,
    diagnostics = diagnostics,
    counts = counts,
    profiles = profiles,
    mappings = mappings,
    model_results = model_results,
    agreement = agreement,
    old_assignment_comparison = old_assignment_comparison
  ),
  file.path(output_dir, "NHANES_robust_reanalysis_results.rds")
)

primary_summary <- model_results |>
  dplyr::filter(
    variant == "log6_winsor",
    model %in% c(
      "Age, sex, extended comorbidity",
      "Fully adjusted demographics, socioeconomic factors, comorbidity, cycle",
      "Survey-weighted fully adjusted"
    )
  ) |>
  dplyr::select(variant, model, weighted, n_model, comparison, hazard_ratio_95ci, p_value_formatted)

summary_lines <- c(
  "NHANES 2011-2018 robust phenotype reanalysis",
  paste0("Analytical cohort: n = ", nrow(dat), "; deaths = ", sum(dat$MORTSTAT == 1)),
  "",
  "Primary robust preprocessing: 1st/99th percentile winsorisation, log transformation",
  "of NLR/SII/creatinine, followed by z-standardisation and K-means (K=3, nstart=100).",
  "Phenotype labels were assigned from biomarker profiles without mortality information.",
  "",
  "Primary robust model summary:",
  paste(capture.output(print(primary_summary)), collapse = "\n"),
  "",
  "Agreement with corrected raw-z clustering:",
  paste(capture.output(print(agreement)), collapse = "\n")
)
writeLines(summary_lines, file.path(output_dir, "NHANES_robust_reanalysis_summary.txt"))

message("Robust NHANES phenotype reanalysis completed.")
cat(paste(summary_lines, collapse = "\n"), "\n")
