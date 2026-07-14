# ==============================================================================
# NHANES cluster-number validation after cohort and preprocessing correction
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
output_dir <- file.path(project_root, "output", "nhanes_cluster_number_validation")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

dat <- readRDS(input_path)

winsorise <- function(x, lower = 0.01, upper = 0.99) {
  limits <- stats::quantile(x, c(lower, upper), na.rm = TRUE, names = FALSE)
  pmin(pmax(x, limits[1]), limits[2])
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

z <- dat |>
  dplyr::transmute(
    log_nlr = log(winsorise(NLR)),
    log_sii = log(winsorise(SII)),
    haemoglobin = winsorise(LBXHGB),
    total_protein = winsorise(LBXSTP),
    bmi = winsorise(BMXBMI),
    log_creatinine = log(winsorise(LBXSCR))
  ) |>
  as.data.frame() |>
  scale()

cluster_function <- function(x, k) {
  stats::kmeans(x, centers = k, nstart = 20, iter.max = 500, algorithm = "Lloyd")
}

message("Running gap statistic (B=50, K.max=6)...")
set.seed(20260710)
gap <- cluster::clusGap(
  z,
  FUNcluster = cluster_function,
  K.max = 6,
  B = 50,
  verbose = interactive()
)
gap_table <- as.data.frame(gap$Tab) |>
  tibble::rownames_to_column("k") |>
  dplyr::mutate(k = as.integer(k))

gap_selected <- tibble::tibble(
  rule = c("global maximum", "Tibshirani 1-SE", "first local maximum within 1-SE"),
  selected_k = c(
    which.max(gap$Tab[, "gap"]),
    cluster::maxSE(gap$Tab[, "gap"], gap$Tab[, "SE.sim"], method = "Tibs2001SEmax"),
    cluster::maxSE(gap$Tab[, "gap"], gap$Tab[, "SE.sim"], method = "firstSEmax")
  )
)

label_by_profile <- function(kmeans_cluster, k) {
  profiles <- dat |>
    dplyr::mutate(cluster_raw = kmeans_cluster) |>
    dplyr::group_by(cluster_raw) |>
    dplyr::summarise(
      n = dplyr::n(),
      nlr = stats::median(NLR),
      sii = stats::median(SII),
      haemoglobin = stats::median(LBXHGB),
      total_protein = stats::median(LBXSTP),
      bmi = stats::median(BMXBMI),
      creatinine = stats::median(LBXSCR),
      .groups = "drop"
    )

  profile_matrix <- data.frame(
    log_nlr = log(profiles$nlr),
    log_sii = log(profiles$sii),
    negative_haemoglobin = -profiles$haemoglobin,
    negative_total_protein = -profiles$total_protein,
    negative_bmi = -profiles$bmi,
    log_creatinine = log(profiles$creatinine)
  )
  safe_scale <- function(x) {
    current_sd <- stats::sd(x)
    if (!is.finite(current_sd) || current_sd == 0) return(rep(0, length(x)))
    (x - mean(x)) / current_sd
  }
  scaled_profiles <- as.data.frame(lapply(profile_matrix, safe_scale))
  vulnerability <- rowSums(scaled_profiles)
  profiles$vulnerability_score <- vulnerability
  order_high_to_low <- order(vulnerability, decreasing = TRUE)

  labels <- if (k == 2L) c("High vulnerability", "Lower vulnerability") else paste0("V", seq_len(k))
  mapping <- tibble::tibble(
    cluster_raw = profiles$cluster_raw[order_high_to_low],
    phenotype = labels
  )
  assigned <- mapping$phenotype[match(kmeans_cluster, mapping$cluster_raw)]

  list(
    assigned = factor(assigned, levels = rev(labels)),
    mapping = mapping,
    profiles = dplyr::left_join(profiles, mapping, by = "cluster_raw")
  )
}

fit_cluster_models <- function(k) {
  set.seed(20260710 + k)
  fit <- stats::kmeans(z, centers = k, nstart = 100, iter.max = 500, algorithm = "Lloyd")
  labelled <- label_by_profile(fit$cluster, k)

  model_data <- dat |>
    dplyr::mutate(
      phenotype = labelled$assigned,
      male = as.integer(RIAGENDR == 1),
      race = factor(RIDRETH3)
    ) |>
    dplyr::filter(
      !is.na(Comorbidity_Score_Extended), !is.na(INDFMPIR), !is.na(race)
    )

  fit_unweighted <- survival::coxph(
    survival::Surv(PERMTH_INT, MORTSTAT) ~ phenotype + RIDAGEYR + male +
      race + INDFMPIR + Comorbidity_Score_Extended,
    data = model_data,
    ties = "efron"
  )
  design <- survey::svydesign(
    ids = ~SDMVPSU,
    strata = ~SDMVSTRA,
    weights = ~WTMEC8YR,
    nest = TRUE,
    data = model_data
  )
  fit_weighted <- survey::svycoxph(
    survival::Surv(PERMTH_INT, MORTSTAT) ~ phenotype + RIDAGEYR + male +
      race + INDFMPIR + Comorbidity_Score_Extended,
    design = design
  )

  extract <- function(model, weighted) {
    co <- summary(model)$coefficients
    ci <- suppressMessages(stats::confint(model))
    terms <- rownames(co)
    keep <- grepl("^phenotype", terms)
    p_col <- grep("Pr\\(", colnames(co), value = TRUE)[1]
    tibble::tibble(
      k = k,
      weighted = weighted,
      n_model = nrow(model_data),
      comparison = sub("^phenotype", "", terms[keep]),
      HR = exp(co[keep, "coef"]),
      lower_95 = exp(ci[keep, 1]),
      upper_95 = exp(ci[keep, 2]),
      p_value = co[keep, p_col],
      hazard_ratio_95ci = sprintf(
        "%.2f (%.2f-%.2f)", exp(co[keep, "coef"]), exp(ci[keep, 1]), exp(ci[keep, 2])
      )
    )
  }

  counts <- dat |>
    dplyr::mutate(phenotype = labelled$assigned) |>
    dplyr::group_by(phenotype) |>
    dplyr::summarise(
      n = dplyr::n(),
      percent = 100 * n / nrow(dat),
      deaths = sum(MORTSTAT == 1),
      mortality_percent = 100 * mean(MORTSTAT == 1),
      .groups = "drop"
    ) |>
    dplyr::mutate(k = k)

  list(
    full_cluster = fit$cluster,
    assigned = labelled$assigned,
    profiles = dplyr::mutate(labelled$profiles, k = k),
    counts = counts,
    models = dplyr::bind_rows(extract(fit_unweighted, FALSE), extract(fit_weighted, TRUE))
  )
}

message("Fitting K=2, K=3, and K=4 clinical structures...")
k_results <- lapply(2:4, fit_cluster_models)
names(k_results) <- as.character(2:4)

message("Running 200 repeated 80% subsampling stability analyses...")
set.seed(20260710)
stability_rows <- list()
for (k in 2:4) {
  full_cluster <- k_results[[as.character(k)]]$full_cluster
  for (replicate_id in seq_len(200)) {
    idx <- sample(seq_len(nrow(z)), size = floor(0.8 * nrow(z)), replace = FALSE)
    sub_fit <- stats::kmeans(
      z[idx, , drop = FALSE],
      centers = k,
      nstart = 30,
      iter.max = 500,
      algorithm = "Lloyd"
    )
    stability_rows[[length(stability_rows) + 1L]] <- tibble::tibble(
      k = k,
      replicate = replicate_id,
      adjusted_rand_index = adjusted_rand_index(full_cluster[idx], sub_fit$cluster)
    )
  }
}
stability <- dplyr::bind_rows(stability_rows)
stability_summary <- stability |>
  dplyr::group_by(k) |>
  dplyr::summarise(
    replicates = dplyr::n(),
    median_ari = stats::median(adjusted_rand_index),
    q1_ari = stats::quantile(adjusted_rand_index, 0.25),
    q3_ari = stats::quantile(adjusted_rand_index, 0.75),
    minimum_ari = min(adjusted_rand_index),
    .groups = "drop"
  )

all_profiles <- dplyr::bind_rows(lapply(k_results, `[[`, "profiles"))
all_counts <- dplyr::bind_rows(lapply(k_results, `[[`, "counts"))
all_models <- dplyr::bind_rows(lapply(k_results, `[[`, "models"))

readr::write_csv(gap_table, file.path(output_dir, "Table28B_A_gap_statistic.csv"))
readr::write_csv(gap_selected, file.path(output_dir, "Table28B_B_gap_selected_k.csv"))
readr::write_csv(stability_summary, file.path(output_dir, "Table28B_C_subsampling_stability_summary.csv"))
readr::write_csv(stability, file.path(output_dir, "Table28B_D_subsampling_stability_replicates.csv"))
readr::write_csv(all_profiles, file.path(output_dir, "Table28B_E_K2_K3_K4_profiles.csv"))
readr::write_csv(all_counts, file.path(output_dir, "Table28B_F_K2_K3_K4_counts_mortality.csv"))
readr::write_csv(all_models, file.path(output_dir, "Table28B_G_K2_K3_K4_Cox_models.csv"))
saveRDS(
  list(
    gap = gap,
    gap_table = gap_table,
    gap_selected = gap_selected,
    stability = stability,
    stability_summary = stability_summary,
    k_results = k_results,
    profiles = all_profiles,
    counts = all_counts,
    models = all_models
  ),
  file.path(output_dir, "NHANES_cluster_number_validation_results.rds")
)

summary_lines <- c(
  "NHANES corrected-cohort cluster-number validation",
  "",
  "Gap-statistic selections:",
  paste(capture.output(print(gap_selected)), collapse = "\n"),
  "",
  "Subsampling stability:",
  paste(capture.output(print(stability_summary)), collapse = "\n"),
  "",
  "K=2/K=3/K=4 mortality:",
  paste(capture.output(print(all_counts)), collapse = "\n"),
  "",
  "Adjusted Cox models:",
  paste(capture.output(print(all_models)), collapse = "\n")
)
writeLines(summary_lines, file.path(output_dir, "NHANES_cluster_number_validation_summary.txt"))
message("Cluster-number validation completed.")
cat(paste(summary_lines, collapse = "\n"), "\n")
