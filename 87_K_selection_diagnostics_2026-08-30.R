# ==============================================================================
# Amendment 21: Cross-database clusterability and K-selection audit
# Post-result, outcome-blind diagnostic module. Outcomes are not loaded.
# ==============================================================================

required_packages <- c("cluster", "dplyr", "readr", "tibble", "data.table", "digest", "mclust")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop("Missing required package(s): ", paste(missing_packages, collapse = ", "), call. = FALSE)
}

root <- normalizePath(Sys.getenv("PROJECT_ROOT", unset = getwd()), winslash = "/", mustWork = TRUE)
protocol_path <- file.path(root, "ANALYSIS_FREEZE_AMENDMENT_21_K_SELECTION_PROTOCOL_2026-08-30.md")
script_path <- file.path(root, "87_K_selection_diagnostics_2026-08-30.R")
paths <- list(
  nhanes = file.path(
    root, "output", "nhanes_albumin_amendment17_2026-08-30",
    "INTERNAL_NHANES_albumin_complete_cohort.rds"
  ),
  mimic = file.path(
    root, "output", "mimic_albumin_primary_2026-08-30",
    "MIMIC_albumin_primary_results.rds"
  ),
  eicu = file.path(
    root, "output", "eicu_24h_extraction",
    "eICU_first_ICU_feature_availability_dataset.csv"
  )
)
paths$eicu <- Sys.getenv("EICU_DENOMINATOR_CSV", unset=file.path(root,"data","derived","eICU_first_ICU_feature_availability_dataset.csv"))

required_paths <- c(unname(unlist(paths)), protocol_path, script_path)
missing_paths <- required_paths[!file.exists(required_paths)]
if (length(missing_paths) > 0L) {
  stop("Missing input(s): ", paste(missing_paths, collapse = ", "), call. = FALSE)
}

output_dir <- file.path(root, "output", "k_selection_diagnostics_amendment21_2026-08-30")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

winsorise <- function(x) {
  limits <- stats::quantile(x, c(0.01, 0.99), na.rm = TRUE, names = FALSE)
  pmin(pmax(as.numeric(x), limits[1]), limits[2])
}

prepare_matrix <- function(data) {
  required <- c("nlr", "sii", "haemoglobin", "albumin", "bmi", "creatinine")
  if (!identical(names(data), required)) {
    stop("Clustering input contains fields outside the locked six-feature matrix.", call. = FALSE)
  }
  transformed <- data.frame(
    log_nlr = log(pmax(winsorise(data$nlr), .Machine$double.eps)),
    log_sii = log(pmax(winsorise(data$sii), .Machine$double.eps)),
    haemoglobin = winsorise(data$haemoglobin),
    albumin = winsorise(data$albumin),
    bmi = winsorise(data$bmi),
    log_creatinine = log(pmax(winsorise(data$creatinine), .Machine$double.eps))
  )
  z <- scale(transformed)
  if (any(!is.finite(z))) stop("Non-finite transformed clustering value.", call. = FALSE)
  z
}

cluster_function <- function(x, k) {
  stats::kmeans(
    x, centers = k, nstart = 20, iter.max = 500, algorithm = "Lloyd"
  )
}

fit_reference_solutions <- function(z) {
  output <- vector("list", 5L)
  names(output) <- as.character(2:6)
  for (k in 2:6) {
    set.seed(20260710)
    output[[as.character(k)]] <- stats::kmeans(
      z, centers = k, nstart = 100, iter.max = 500, algorithm = "Lloyd"
    )
  }
  output
}

select_smallest_within <- function(values, tolerance, direction = c("max", "min")) {
  direction <- match.arg(direction)
  target <- if (direction == "max") max(values) else min(values)
  eligible <- if (direction == "max") {
    which(values >= target - tolerance)
  } else {
    which(values <= target + tolerance)
  }
  min(eligible)
}

run_gap <- function(z, dataset) {
  message("Running gap statistic for ", dataset, "...")
  set.seed(20260710)
  result <- cluster::clusGap(
    z, FUNcluster = cluster_function, K.max = 6, B = 100, verbose = FALSE
  )
  table <- as.data.frame(result$Tab) |>
    tibble::rownames_to_column("k") |>
    dplyr::mutate(dataset = dataset, k = as.integer(k), .before = 1)
  selected <- tibble::tibble(
    dataset = dataset,
    rule = c("global maximum", "Tibshirani 1-SE", "first local maximum within 1-SE"),
    selected_k = c(
      which.max(result$Tab[, "gap"]),
      cluster::maxSE(result$Tab[, "gap"], result$Tab[, "SE.sim"], method = "Tibs2001SEmax"),
      cluster::maxSE(result$Tab[, "gap"], result$Tab[, "SE.sim"], method = "firstSEmax")
    )
  )
  list(table = table, selected = selected)
}

run_silhouette <- function(z, references, dataset, dataset_index) {
  message("Running silhouette diagnostic for ", dataset, "...")
  if (nrow(z) > 5000L) {
    set.seed(20260711 + dataset_index)
    idx <- sort(sample(seq_len(nrow(z)), 5000L, replace = FALSE))
  } else {
    idx <- seq_len(nrow(z))
  }
  distance <- stats::dist(z[idx, , drop = FALSE])
  rows <- lapply(2:6, function(k) {
    labels <- references[[as.character(k)]]$cluster[idx]
    widths <- cluster::silhouette(labels, distance)[, "sil_width"]
    tibble::tibble(
      dataset = dataset,
      k = k,
      n_evaluated = length(idx),
      mean_silhouette = mean(widths),
      median_silhouette = stats::median(widths),
      q1_silhouette = stats::quantile(widths, 0.25),
      q3_silhouette = stats::quantile(widths, 0.75),
      negative_fraction = mean(widths < 0)
    )
  })
  table <- dplyr::bind_rows(rows)
  selected_position <- select_smallest_within(table$mean_silhouette, 0.001, "max")
  list(table = table, selected_k = table$k[selected_position], n_evaluated = length(idx))
}

run_consensus <- function(z, dataset, dataset_index) {
  message("Running consensus diagnostic for ", dataset, "...")
  if (nrow(z) > 2000L) {
    set.seed(20260712 + dataset_index)
    audit_idx <- sort(sample(seq_len(nrow(z)), 2000L, replace = FALSE))
    audit_z <- z[audit_idx, , drop = FALSE]
  } else {
    audit_z <- z
  }
  n <- nrow(audit_z)
  summary_rows <- list()
  cdf_rows <- list()
  for (k in 2:6) {
    co_sampled <- matrix(0L, nrow = n, ncol = n)
    co_clustered <- matrix(0L, nrow = n, ncol = n)
    set.seed(2026071200 + dataset_index * 100L + k)
    for (replicate_id in seq_len(100L)) {
      idx <- sort(sample(seq_len(n), floor(0.8 * n), replace = FALSE))
      fit <- stats::kmeans(
        audit_z[idx, , drop = FALSE], centers = k, nstart = 20,
        iter.max = 500, algorithm = "Lloyd"
      )
      co_sampled[idx, idx] <- co_sampled[idx, idx] + 1L
      for (cluster_id in seq_len(k)) {
        members <- idx[fit$cluster == cluster_id]
        co_clustered[members, members] <- co_clustered[members, members] + 1L
      }
    }
    upper <- upper.tri(co_sampled)
    valid <- upper & co_sampled > 0L
    consensus_values <- co_clustered[valid] / co_sampled[valid]
    grid <- seq(0, 1, by = 0.01)
    cdf <- vapply(grid, function(value) mean(consensus_values <= value), numeric(1))
    auc <- sum(diff(grid) * (head(cdf, -1L) + tail(cdf, -1L)) / 2)
    pac <- mean(consensus_values > 0.1 & consensus_values < 0.9)
    summary_rows[[as.character(k)]] <- tibble::tibble(
      dataset = dataset,
      k = k,
      n_audit_sample = n,
      repetitions = 100L,
      consensus_cdf_auc = auc,
      pac_0_1_0_9 = pac
    )
    cdf_rows[[as.character(k)]] <- tibble::tibble(
      dataset = dataset, k = k, consensus = grid, empirical_cdf = cdf
    )
    rm(co_sampled, co_clustered, consensus_values)
    invisible(gc(FALSE))
  }
  table <- dplyr::bind_rows(summary_rows) |>
    dplyr::arrange(k) |>
    dplyr::mutate(delta_cdf_auc = consensus_cdf_auc - dplyr::lag(consensus_cdf_auc))
  selected_position <- select_smallest_within(table$pac_0_1_0_9, 0.01, "min")
  list(
    table = table,
    cdf = dplyr::bind_rows(cdf_rows),
    selected_k = table$k[selected_position],
    n_audit_sample = n
  )
}

cluster_jaccard <- function(reference_labels, comparison_labels, k) {
  vapply(seq_len(k), function(reference_cluster) {
    reference_members <- reference_labels == reference_cluster
    max(vapply(seq_len(k), function(comparison_cluster) {
      comparison_members <- comparison_labels == comparison_cluster
      intersection <- sum(reference_members & comparison_members)
      union <- sum(reference_members | comparison_members)
      if (union == 0L) 0 else intersection / union
    }, numeric(1)))
  }, numeric(1))
}

run_stability <- function(z, references, dataset, dataset_index) {
  message("Running subsampling stability for ", dataset, "...")
  ari_rows <- list()
  jaccard_rows <- list()
  size_rows <- list()
  diagnostic_rows <- list()
  for (k in 2:6) {
    reference <- references[[as.character(k)]]
    counts <- table(factor(reference$cluster, levels = seq_len(k)))
    size_rows[[as.character(k)]] <- tibble::tibble(
      dataset = dataset,
      k = k,
      cluster_raw = seq_len(k),
      n = as.integer(counts),
      fraction = as.integer(counts) / nrow(z)
    )
    diagnostic_rows[[as.character(k)]] <- tibble::tibble(
      dataset = dataset,
      k = k,
      tot_withinss = reference$tot.withinss,
      betweenss = reference$betweenss,
      totalss = reference$totss,
      iterations = reference$iter,
      converged_before_itermax = reference$iter < 500L
    )
    set.seed(2026072000 + dataset_index * 100L + k)
    for (replicate_id in seq_len(200L)) {
      idx <- sort(sample(seq_len(nrow(z)), floor(0.8 * nrow(z)), replace = FALSE))
      fit <- stats::kmeans(
        z[idx, , drop = FALSE], centers = k, nstart = 20,
        iter.max = 500, algorithm = "Lloyd"
      )
      reference_labels <- reference$cluster[idx]
      ari_rows[[length(ari_rows) + 1L]] <- tibble::tibble(
        dataset = dataset,
        k = k,
        replicate = replicate_id,
        adjusted_rand_index = mclust::adjustedRandIndex(reference_labels, fit$cluster)
      )
      values <- cluster_jaccard(reference_labels, fit$cluster, k)
      jaccard_rows[[length(jaccard_rows) + 1L]] <- tibble::tibble(
        dataset = dataset,
        k = k,
        replicate = replicate_id,
        reference_cluster = seq_len(k),
        maximum_jaccard = values
      )
    }
  }
  ari <- dplyr::bind_rows(ari_rows)
  jaccard <- dplyr::bind_rows(jaccard_rows)
  sizes <- dplyr::bind_rows(size_rows)
  ari_summary <- ari |>
    dplyr::group_by(dataset, k) |>
    dplyr::summarise(
      repetitions = dplyr::n(),
      median_ari = stats::median(adjusted_rand_index),
      q1_ari = stats::quantile(adjusted_rand_index, 0.25),
      q3_ari = stats::quantile(adjusted_rand_index, 0.75),
      minimum_ari = min(adjusted_rand_index),
      .groups = "drop"
    )
  jaccard_summary <- jaccard |>
    dplyr::group_by(dataset, k, reference_cluster) |>
    dplyr::summarise(
      repetitions = dplyr::n(),
      mean_jaccard = mean(maximum_jaccard),
      median_jaccard = stats::median(maximum_jaccard),
      q1_jaccard = stats::quantile(maximum_jaccard, 0.25),
      q3_jaccard = stats::quantile(maximum_jaccard, 0.75),
      minimum_jaccard = min(maximum_jaccard),
      .groups = "drop"
    )
  k_gates <- jaccard_summary |>
    dplyr::group_by(dataset, k) |>
    dplyr::summarise(
      minimum_cluster_mean_jaccard = min(mean_jaccard),
      all_cluster_mean_jaccard_ge_0_75 = all(mean_jaccard >= 0.75),
      .groups = "drop"
    ) |>
    dplyr::left_join(
      sizes |>
        dplyr::group_by(dataset, k) |>
        dplyr::summarise(
          minimum_cluster_fraction = min(fraction),
          all_clusters_ge_5_percent = all(fraction >= 0.05),
          .groups = "drop"
        ),
      by = c("dataset", "k")
    )
  list(
    ari = ari,
    ari_summary = ari_summary,
    jaccard = jaccard,
    jaccard_summary = jaccard_summary,
    sizes = sizes,
    k_gates = k_gates,
    diagnostics = dplyr::bind_rows(diagnostic_rows)
  )
}

message("Loading locked outcome-blind clustering cohorts...")
nhanes <- readRDS(paths$nhanes) |>
  dplyr::transmute(
    nlr = NLR, sii = SII, haemoglobin = LBXHGB, albumin = LBXSAL,
    bmi = BMXBMI, creatinine = LBXSCR
  )
mimic <- readRDS(paths$mimic)$analysis |>
  dplyr::transmute(
    nlr = nlr, sii = sii, haemoglobin = haemoglobin, albumin = albumin,
    bmi = bmi, creatinine = creatinine
  )
eicu <- data.table::fread(paths$eicu, data.table = FALSE, showProgress = FALSE) |>
  tibble::as_tibble() |>
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
    nlr = nlr, sii = sii_like, haemoglobin = haemoglobin, albumin = albumin,
    bmi = bmi, creatinine = creatinine
  )

cohorts <- list(NHANES = nhanes, `MIMIC-IV` = mimic, eICU = eicu)
expected_n <- c(NHANES = 4637L, `MIMIC-IV` = 1100L, eICU = 15242L)
observed_n <- vapply(cohorts, nrow, integer(1))
if (!identical(observed_n, expected_n)) {
  stop(
    "Locked cohort size mismatch. Observed: ",
    paste(names(observed_n), observed_n, collapse = "; "), call. = FALSE
  )
}

all_gap <- list()
all_gap_selected <- list()
all_silhouette <- list()
all_consensus <- list()
all_consensus_cdf <- list()
all_ari <- list()
all_ari_summary <- list()
all_jaccard <- list()
all_jaccard_summary <- list()
all_sizes <- list()
all_gates <- list()
all_diagnostics <- list()
decision_rows <- list()

for (dataset_index in seq_along(cohorts)) {
  dataset <- names(cohorts)[dataset_index]
  message("==== ", dataset, " ====")
  z <- prepare_matrix(cohorts[[dataset]])
  references <- fit_reference_solutions(z)
  gap <- run_gap(z, dataset)
  silhouette <- run_silhouette(z, references, dataset, dataset_index)
  consensus <- run_consensus(z, dataset, dataset_index)
  stability <- run_stability(z, references, dataset, dataset_index)

  gap_k <- gap$selected |>
    dplyr::filter(rule == "Tibshirani 1-SE") |>
    dplyr::pull(selected_k)
  silhouette_k <- silhouette$selected_k
  pac_k <- consensus$selected_k
  gate <- stability$k_gates |>
    dplyr::filter(k == gap_k)
  stable_gate <- if (gap_k >= 2L && nrow(gate) == 1L) gate$all_cluster_mean_jaccard_ge_0_75 else FALSE
  size_gate <- if (gap_k >= 2L && nrow(gate) == 1L) gate$all_clusters_ge_5_percent else FALSE
  agreement_count <- as.integer(silhouette_k == gap_k) + as.integer(pac_k == gap_k)
  conclusion <- if (
    gap_k >= 2L && agreement_count == 2L && stable_gate && size_gate
  ) {
    "unique natural K"
  } else if (
    gap_k >= 2L && agreement_count == 1L && stable_gate && size_gate
  ) {
    "partially supported K"
  } else {
    "no uniquely supported natural K"
  }
  decision_rows[[dataset]] <- tibble::tibble(
    dataset = dataset,
    n = nrow(z),
    gap_tibshirani_1se_k = gap_k,
    silhouette_selected_k = silhouette_k,
    pac_selected_k = pac_k,
    supporting_selector_agreement_count = agreement_count,
    gap_selected_k_stability_gate = stable_gate,
    gap_selected_k_size_gate = size_gate,
    dataset_level_conclusion = conclusion
  )

  all_gap[[dataset]] <- gap$table
  all_gap_selected[[dataset]] <- gap$selected
  all_silhouette[[dataset]] <- silhouette$table
  all_consensus[[dataset]] <- consensus$table
  all_consensus_cdf[[dataset]] <- consensus$cdf
  all_ari[[dataset]] <- stability$ari
  all_ari_summary[[dataset]] <- stability$ari_summary
  all_jaccard[[dataset]] <- stability$jaccard
  all_jaccard_summary[[dataset]] <- stability$jaccard_summary
  all_sizes[[dataset]] <- stability$sizes
  all_gates[[dataset]] <- stability$k_gates
  all_diagnostics[[dataset]] <- stability$diagnostics
  rm(z, references)
  invisible(gc(FALSE))
}

decision_matrix <- dplyr::bind_rows(decision_rows)
unique_rows <- decision_matrix |>
  dplyr::filter(dataset_level_conclusion == "unique natural K")
common_natural_k <- if (
  nrow(unique_rows) == 3L && length(unique(unique_rows$gap_tibshirani_1se_k)) == 1L
) unique_rows$gap_tibshirani_1se_k[1] else NA_integer_
cross_database_decision <- tibble::tibble(
  all_three_support_same_unique_k = !is.na(common_natural_k),
  common_natural_k = common_natural_k,
  primary_manuscript_resolution = ifelse(is.na(common_natural_k), 3L, common_natural_k),
  statistical_conclusion = ifelse(
    is.na(common_natural_k),
    "No common unique natural K was established across all three databases.",
    paste0("All three databases independently supported the same unique natural K=", common_natural_k, ".")
  ),
  manuscript_consequence = ifelse(
    is.na(common_natural_k),
    "Retain K=3 only as the main descriptive resolution of a continuum; do not call it optimal, natural, or prospectively prespecified.",
    paste0("Reconsider the manuscript main resolution as K=", common_natural_k, " under the locked post-result decision rule.")
  )
)

tables <- list(
  Table87A_gap_statistic = dplyr::bind_rows(all_gap),
  Table87B_gap_selected_k = dplyr::bind_rows(all_gap_selected),
  Table87C_silhouette = dplyr::bind_rows(all_silhouette),
  Table87D_consensus_summary = dplyr::bind_rows(all_consensus),
  Table87E_consensus_cdf = dplyr::bind_rows(all_consensus_cdf),
  Table87F_subsampling_ARI_summary = dplyr::bind_rows(all_ari_summary),
  Table87G_clusterwise_Jaccard_summary = dplyr::bind_rows(all_jaccard_summary),
  Table87H_cluster_sizes = dplyr::bind_rows(all_sizes),
  Table87I_stability_and_size_gates = dplyr::bind_rows(all_gates),
  Table87J_kmeans_diagnostics = dplyr::bind_rows(all_diagnostics),
  Table87K_database_decision_matrix = decision_matrix,
  Table87L_cross_database_decision = cross_database_decision
)

qa <- tibble::tribble(
  ~check, ~passed, ~detail,
  "Locked cohort sizes reproduced", identical(observed_n, expected_n), paste(observed_n, collapse = "/"),
  "Gap tables contain K=1 through K=6", nrow(tables$Table87A_gap_statistic) == 18L, paste0("rows=", nrow(tables$Table87A_gap_statistic)),
  "Silhouette tables contain K=2 through K=6", nrow(tables$Table87C_silhouette) == 15L, paste0("rows=", nrow(tables$Table87C_silhouette)),
  "Consensus completed 100 repetitions per database and K", all(tables$Table87D_consensus$repetitions == 100L), "all=100",
  "Stability completed 200 repetitions per database and K", all(tables$Table87F_subsampling_ARI_summary$repetitions == 200L), "all=200",
  "Cluster-wise Jaccard completed 200 repetitions", all(tables$Table87G_clusterwise_Jaccard_summary$repetitions == 200L), "all=200",
  "All K-means reference fits converged before iter.max", all(tables$Table87J_kmeans_diagnostics$converged_before_itermax), "all TRUE",
  "Decision matrix has one row per database", nrow(decision_matrix) == 3L, paste0("rows=", nrow(decision_matrix)),
  "Only six locked features entered clustering", TRUE, "outcome-like fields rejected by prepare_matrix",
  "No participant-level assignments written", TRUE, "aggregate diagnostic outputs only"
)
if (!all(qa$passed)) {
  stop("K-selection QA failed before output write: ", paste(qa$check[!qa$passed], collapse = "; "), call. = FALSE)
}
tables$Table87M_QA <- qa

for (name in names(tables)) {
  readr::write_csv(tables[[name]], file.path(output_dir, paste0(name, ".csv")))
}
readr::write_csv(
  dplyr::bind_rows(all_ari),
  file.path(output_dir, "INTERNAL_AGGREGATE_stability_ARI_replicates.csv")
)
readr::write_csv(
  dplyr::bind_rows(all_jaccard),
  file.path(output_dir, "INTERNAL_AGGREGATE_clusterwise_Jaccard_replicates.csv")
)
saveRDS(tables, file.path(output_dir, "K_selection_diagnostics_aggregate_results.rds"))

package_versions <- tibble::tibble(
  package = required_packages,
  version = vapply(required_packages, function(x) as.character(utils::packageVersion(x)), character(1))
)
readr::write_csv(package_versions, file.path(output_dir, "PACKAGE_VERSIONS.csv"))

source_manifest <- tibble::tibble(
  file = basename(required_paths),
  path = normalizePath(required_paths, winslash = "/", mustWork = TRUE),
  bytes = file.info(required_paths)$size,
  sha256 = vapply(required_paths, digest::digest, character(1), file = TRUE, algo = "sha256")
)
readr::write_csv(source_manifest, file.path(output_dir, "SOURCE_SHA256_MANIFEST.csv"))

summary_lines <- c(
  "Amendment 21 cross-database clusterability and K-selection audit",
  "",
  "Timing: post-result, outcome-blind methodological audit; not prospectively registered.",
  "",
  "Database decision matrix:",
  paste(capture.output(print(decision_matrix, n = Inf, width = Inf)), collapse = "\n"),
  "",
  "Cross-database decision:",
  paste(capture.output(print(cross_database_decision, n = Inf, width = Inf)), collapse = "\n"),
  "",
  "Gap selections:",
  paste(capture.output(print(tables$Table87B_gap_selected_k, n = Inf, width = Inf)), collapse = "\n"),
  "",
  "Silhouette:",
  paste(capture.output(print(tables$Table87C_silhouette, n = Inf, width = Inf)), collapse = "\n"),
  "",
  "Consensus PAC:",
  paste(capture.output(print(tables$Table87D_consensus_summary, n = Inf, width = Inf)), collapse = "\n"),
  "",
  "ARI stability:",
  paste(capture.output(print(tables$Table87F_subsampling_ARI_summary, n = Inf, width = Inf)), collapse = "\n"),
  "",
  "Cluster-wise Jaccard:",
  paste(capture.output(print(tables$Table87G_clusterwise_Jaccard_summary, n = Inf, width = Inf)), collapse = "\n"),
  "",
  "Interpretation boundary: the locked rules distinguish evidence for a common natural K from use of K=3 as a descriptive partition of a continuum."
)
writeLines(summary_lines, file.path(output_dir, "K_selection_diagnostics_summary.txt"))

output_files <- list.files(output_dir, full.names = TRUE, recursive = FALSE)
output_files <- output_files[basename(output_files) != "OUTPUT_SHA256_MANIFEST.csv"]
output_files <- output_files[!file.info(output_files)$isdir]
output_manifest <- tibble::tibble(
  file = basename(output_files),
  bytes = file.info(output_files)$size,
  sha256 = vapply(output_files, digest::digest, character(1), file = TRUE, algo = "sha256")
)
readr::write_csv(output_manifest, file.path(output_dir, "OUTPUT_SHA256_MANIFEST.csv"))

cat(paste(summary_lines, collapse = "\n"), "\n")
message("Amendment 21 K-selection diagnostics completed.")
