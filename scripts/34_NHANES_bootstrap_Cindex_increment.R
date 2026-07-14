# ==============================================================================
# Cluster-aware out-of-bag bootstrap validation of incremental discrimination
# ==============================================================================

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
input_path <- file.path(
  root, "output", "nhanes_albumin_benchmarks",
  "NHANES_albumin_benchmark_results.rds"
)
output_dir <- file.path(root, "output", "nhanes_bootstrap_cindex")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

n_boot <- as.integer(Sys.getenv("CINDEX_BOOTSTRAP_REPS", unset = "1000"))
if (!is.finite(n_boot) || n_boot < 100L) {
  stop("CINDEX_BOOTSTRAP_REPS must be at least 100.", call. = FALSE)
}

model_data <- readRDS(input_path)$model_data
original_benchmarks <- readRDS(input_path)$benchmarks |>
  dplyr::select(predictor, apparent_c_index = c_index, apparent_delta_c_index = delta_c_index)

winsor_limits <- function(x) {
  stats::quantile(x, c(0.01, 0.99), na.rm = TRUE, names = FALSE)
}

apply_winsor <- function(x, limits) {
  pmin(pmax(as.numeric(x), limits[1]), limits[2])
}

safe_standardise <- function(x, centre, spread) {
  if (!is.finite(spread) || spread == 0) return(rep(0, length(x)))
  (x - centre) / spread
}

z_three <- function(x) {
  spread <- stats::sd(x)
  if (!is.finite(spread) || spread == 0) return(rep(0, length(x)))
  (x - mean(x)) / spread
}

fit_cluster_model <- function(train, new_data, protein_variable) {
  raw_variables <- c("NLR", "SII", "LBXHGB", protein_variable, "BMXBMI", "LBXSCR")
  transformed_names <- c("log_nlr", "log_sii", "haemoglobin", "protein", "bmi", "log_creatinine")

  limits <- lapply(raw_variables, function(variable) winsor_limits(train[[variable]]))
  names(limits) <- raw_variables

  transform_matrix <- function(dat) {
    transformed <- cbind(
      log(pmax(apply_winsor(dat$NLR, limits$NLR), .Machine$double.eps)),
      log(pmax(apply_winsor(dat$SII, limits$SII), .Machine$double.eps)),
      apply_winsor(dat$LBXHGB, limits$LBXHGB),
      apply_winsor(dat[[protein_variable]], limits[[protein_variable]]),
      apply_winsor(dat$BMXBMI, limits$BMXBMI),
      log(pmax(apply_winsor(dat$LBXSCR, limits$LBXSCR), .Machine$double.eps))
    )
    colnames(transformed) <- transformed_names
    transformed
  }

  train_matrix_unscaled <- transform_matrix(train)
  centres <- colMeans(train_matrix_unscaled)
  spreads <- apply(train_matrix_unscaled, 2, stats::sd)
  spreads[!is.finite(spreads) | spreads == 0] <- 1
  train_matrix <- sweep(sweep(train_matrix_unscaled, 2, centres, "-"), 2, spreads, "/")

  cluster_fit <- stats::kmeans(
    train_matrix,
    centers = 3,
    nstart = 25,
    iter.max = 300,
    algorithm = "Lloyd"
  )

  profiles <- train |>
    dplyr::mutate(cluster_raw = cluster_fit$cluster) |>
    dplyr::group_by(cluster_raw) |>
    dplyr::summarise(
      nlr = stats::median(NLR),
      sii = stats::median(SII),
      haemoglobin = stats::median(LBXHGB),
      protein = stats::median(.data[[protein_variable]]),
      bmi = stats::median(BMXBMI),
      creatinine = stats::median(LBXSCR),
      .groups = "drop"
    )

  p1_score <- z_three(log(profiles$nlr)) + z_three(log(profiles$sii)) +
    z_three(log(profiles$creatinine))
  p1_raw <- profiles$cluster_raw[which.max(p1_score)]
  remaining <- profiles |>
    dplyr::filter(cluster_raw != p1_raw)
  p2_score <- -z_three(remaining$haemoglobin) - z_three(remaining$protein) -
    z_three(remaining$bmi)
  p2_raw <- remaining$cluster_raw[which.max(p2_score)]
  p3_raw <- setdiff(profiles$cluster_raw, c(p1_raw, p2_raw))
  label_map <- stats::setNames(c("P1", "P2", "P3"), c(p1_raw, p2_raw, p3_raw))

  assign_cluster <- function(dat) {
    unscaled <- transform_matrix(dat)
    scaled <- sweep(sweep(unscaled, 2, centres, "-"), 2, spreads, "/")
    squared_distances <- vapply(
      seq_len(3),
      function(cluster) rowSums((scaled - matrix(
        cluster_fit$centers[cluster, ], nrow = nrow(scaled), ncol = ncol(scaled), byrow = TRUE
      ))^2),
      numeric(nrow(scaled))
    )
    raw_assignment <- max.col(-squared_distances, ties.method = "first")
    factor(unname(label_map[as.character(raw_assignment)]), levels = c("P3", "P2", "P1"))
  }

  list(train = assign_cluster(train), new = assign_cluster(new_data))
}

concordance_on_new_data <- function(fit, new_data) {
  linear_predictor <- stats::predict(fit, newdata = new_data, type = "lp")
  unname(
    survival::concordance(
      survival::Surv(PERMTH_INT, MORTSTAT) ~ linear_predictor,
      data = new_data,
      reverse = TRUE
    )$concordance
  )
}

base_formula <- survival::Surv(PERMTH_INT, MORTSTAT) ~
  RIDAGEYR + male + race + INDFMPIR + Comorbidity_Score_Extended + cycle +
  smoking + hypertension + diabetes

predictor_terms <- c(
  phenotype_total_protein = "phenotype_total_boot",
  phenotype_albumin = "phenotype_albumin_boot",
  domain_balanced_score = "domain_balanced_score",
  NLR = "risk_log_nlr",
  SII = "risk_log_sii",
  PNI = "risk_low_pni",
  GNRI = "risk_low_gnri",
  HALP = "risk_low_halp"
)

set.seed(20260710)
n <- nrow(model_data)
replicate_rows <- vector("list", n_boot)
failure_messages <- character()

for (b in seq_len(n_boot)) {
  current <- tryCatch(withCallingHandlers({
    sampled_index <- sample.int(n, size = n, replace = TRUE)
    out_of_bag_index <- setdiff(seq_len(n), unique(sampled_index))
    if (length(out_of_bag_index) < 200L) stop("Too few out-of-bag observations")

    train <- model_data[sampled_index, , drop = FALSE]
    out_of_bag <- model_data[out_of_bag_index, , drop = FALSE]

    total_clusters <- fit_cluster_model(train, out_of_bag, "LBXSTP")
    albumin_clusters <- fit_cluster_model(train, out_of_bag, "albumin_g_dl")
    train$phenotype_total_boot <- total_clusters$train
    out_of_bag$phenotype_total_boot <- total_clusters$new
    train$phenotype_albumin_boot <- albumin_clusters$train
    out_of_bag$phenotype_albumin_boot <- albumin_clusters$new

    base_fit <- survival::coxph(base_formula, data = train, ties = "efron")
    base_c <- concordance_on_new_data(base_fit, out_of_bag)

    dplyr::bind_rows(lapply(names(predictor_terms), function(label) {
      term <- predictor_terms[[label]]
      extended_fit <- survival::coxph(
        stats::update.formula(base_formula, paste(". ~ . +", term)),
        data = train,
        ties = "efron"
      )
      extended_c <- concordance_on_new_data(extended_fit, out_of_bag)
      tibble::tibble(
        replicate = b,
        predictor = label,
        out_of_bag_n = nrow(out_of_bag),
        out_of_bag_events = sum(out_of_bag$MORTSTAT == 1),
        base_c_index = base_c,
        extended_c_index = extended_c,
        delta_c_index = extended_c - base_c
      )
    }))
  }, warning = function(w) {
    if (grepl("coefficient may be infinite", conditionMessage(w), fixed = TRUE)) {
      stop(conditionMessage(w), call. = FALSE)
    }
    invokeRestart("muffleWarning")
  }), error = function(e) {
    failure_messages <<- c(failure_messages, paste0("Replicate ", b, ": ", conditionMessage(e)))
    NULL
  })

  replicate_rows[[b]] <- current
  if (b %% 100L == 0L) {
    message("Completed ", b, "/", n_boot, " bootstrap replicates.")
  }
}

replicates <- dplyr::bind_rows(replicate_rows)
valid_counts <- replicates |>
  dplyr::count(predictor, name = "valid_replicates")
if (nrow(valid_counts) != length(predictor_terms) || min(valid_counts$valid_replicates) < 0.90 * n_boot) {
  stop("Fewer than 90% valid bootstrap replicates for at least one predictor.", call. = FALSE)
}

bootstrap_summary <- replicates |>
  dplyr::group_by(predictor) |>
  dplyr::summarise(
    valid_replicates = dplyr::n(),
    median_out_of_bag_n = stats::median(out_of_bag_n),
    median_out_of_bag_c_index = stats::median(extended_c_index),
    median_delta_c_index = stats::median(delta_c_index),
    lower_95_delta = stats::quantile(delta_c_index, 0.025, names = FALSE),
    upper_95_delta = stats::quantile(delta_c_index, 0.975, names = FALSE),
    probability_delta_above_zero = mean(delta_c_index > 0),
    .groups = "drop"
  ) |>
  dplyr::left_join(original_benchmarks, by = "predictor") |>
  dplyr::arrange(dplyr::desc(median_delta_c_index))

method_note <- tibble::tibble(
  item = c(
    "Bootstrap repetitions requested",
    "Failed repetitions",
    "Cluster treatment",
    "Validation sample",
    "Interpretation"
  ),
  value = c(
    as.character(n_boot),
    as.character(length(failure_messages)),
    "Total-protein and albumin K-means refitted in every bootstrap sample",
    "Out-of-bag observations not sampled in the corresponding bootstrap replicate",
    "Distribution estimates internal discrimination uncertainty; it is not a new external validation"
  )
)

readr::write_csv(replicates, file.path(output_dir, "Table34A_bootstrap_Cindex_replicates.csv"))
readr::write_csv(bootstrap_summary, file.path(output_dir, "Table34B_bootstrap_Cindex_summary.csv"))
readr::write_csv(method_note, file.path(output_dir, "Table34C_bootstrap_method_note.csv"))
if (length(failure_messages) > 0L) {
  writeLines(failure_messages, file.path(output_dir, "bootstrap_failure_log.txt"))
}
saveRDS(
  list(replicates = replicates, summary = bootstrap_summary, method_note = method_note),
  file.path(output_dir, "NHANES_bootstrap_Cindex_results.rds")
)

summary_lines <- c(
  "NHANES cluster-aware out-of-bag bootstrap C-index analysis",
  "",
  paste0("Requested replicates: ", n_boot, "; failed replicates: ", length(failure_messages)),
  "",
  paste(capture.output(print(bootstrap_summary)), collapse = "\n"),
  "",
  paste0(
    "Interpretation: percentile intervals describe paired out-of-bag changes in Harrell's C-index. ",
    "K-means was refitted for both phenotype predictors in every bootstrap sample."
  )
)
writeLines(summary_lines, file.path(output_dir, "NHANES_bootstrap_Cindex_summary.txt"))
cat(paste(summary_lines, collapse = "\n"), "\n")
