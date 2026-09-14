# ==============================================================================
# Corrected NHANES cluster-aware bootstrap and proportional-hazards diagnostics
# ==============================================================================

required_packages <- c("dplyr", "readr", "tibble", "survival")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop("Missing required package(s): ", paste(missing_packages, collapse = ", "), call. = FALSE)
}

input_path <- paste0(
  "./output/nhanes_albumin_amendment17_2026-08-30/",
  "NHANES_albumin_benchmark_results.rds"
)
output_dir <- "./output/nhanes_albumin_amendment17_2026-08-30/cluster_bootstrap_PH"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(input_path)) {
  stop("Run script 29 before this analysis.", call. = FALSE)
}

winsorise <- function(x, lower = 0.01, upper = 0.99) {
  limits <- stats::quantile(x, c(lower, upper), na.rm = TRUE, names = FALSE)
  pmin(pmax(x, limits[1]), limits[2])
}

safe_z <- function(x) {
  current_sd <- stats::sd(x, na.rm = TRUE)
  if (!is.finite(current_sd) || current_sd == 0) return(rep(0, length(x)))
  (x - mean(x, na.rm = TRUE)) / current_sd
}

cluster_and_label <- function(data, seed) {
  z <- data |>
    dplyr::transmute(
      log_nlr = log(winsorise(NLR)),
      log_sii = log(winsorise(SII)),
      haemoglobin = winsorise(LBXHGB),
      albumin = winsorise(LBXSAL),
      bmi = winsorise(BMXBMI),
      log_creatinine = log(winsorise(LBXSCR))
    ) |>
    as.data.frame() |>
    scale()

  set.seed(seed)
  fit <- stats::kmeans(
    z, centers = 3, nstart = 10, iter.max = 300, algorithm = "Lloyd"
  )
  profiles <- data |>
    dplyr::mutate(cluster_raw = fit$cluster) |>
    dplyr::group_by(cluster_raw) |>
    dplyr::summarise(
      nlr = stats::median(NLR),
      sii = stats::median(SII),
      haemoglobin = stats::median(LBXHGB),
      albumin = stats::median(LBXSAL),
      bmi = stats::median(BMXBMI),
      creatinine = stats::median(LBXSCR),
      .groups = "drop"
    )

  p1_score <- safe_z(log(profiles$nlr)) + safe_z(log(profiles$sii)) +
    safe_z(log(profiles$creatinine))
  p1_raw <- profiles$cluster_raw[which.max(p1_score)]
  remaining <- profiles |>
    dplyr::filter(cluster_raw != p1_raw)
  p2_score <- -safe_z(remaining$haemoglobin) - safe_z(remaining$albumin) -
    safe_z(remaining$bmi)
  p2_raw <- remaining$cluster_raw[which.max(p2_score)]
  p3_raw <- setdiff(profiles$cluster_raw, c(p1_raw, p2_raw))

  factor(
    dplyr::case_when(
      fit$cluster == p1_raw ~ "P1",
      fit$cluster == p2_raw ~ "P2",
      fit$cluster == p3_raw ~ "P3",
      TRUE ~ NA_character_
    ),
    levels = c("P3", "P2", "P1")
  )
}

results <- readRDS(input_path)
dat <- results$model_data |>
  dplyr::mutate(
    phenotype_albumin = factor(
      phenotype_albumin, levels = c("P3", "P2", "P1")
    ),
    phenotype_albumin = factor(phenotype_albumin, levels = c("P3", "P2", "P1"))
  )

base_terms <- "RIDAGEYR + male + race + INDFMPIR + Comorbidity_Score_Extended + cycle + smoking + hypertension + diabetes"
primary_formula <- stats::as.formula(paste(
  "survival::Surv(PERMTH_INT, MORTSTAT) ~ phenotype_albumin +", base_terms
))
albumin_formula <- stats::as.formula(paste(
  "survival::Surv(PERMTH_INT, MORTSTAT) ~ phenotype_albumin +", base_terms
))
score_formula <- stats::as.formula(paste(
  "survival::Surv(PERMTH_INT, MORTSTAT) ~ domain_balanced_score +", base_terms
))

primary_fit <- survival::coxph(primary_formula, data = dat, ties = "efron", x = TRUE)
albumin_fit <- survival::coxph(albumin_formula, data = dat, ties = "efron", x = TRUE)
score_fit <- survival::coxph(score_formula, data = dat, ties = "efron", x = TRUE)

extract_ph <- function(fit, model) {
  check <- survival::cox.zph(fit)
  tibble::tibble(
    model = model,
    variable = rownames(check$table),
    chisq = check$table[, "chisq"],
    df = check$table[, "df"],
    p_value = check$table[, "p"]
  )
}

ph_table <- dplyr::bind_rows(
  extract_ph(primary_fit, "Albumin phenotype primary"),
  extract_ph(albumin_fit, "Albumin phenotype"),
  extract_ph(score_fit, "Domain-balanced continuous score")
)

original_co <- summary(primary_fit)$coefficients
original_ci <- suppressMessages(stats::confint(primary_fit))
original_p1 <- tibble::tibble(
  estimate_type = "Corrected-cohort primary model",
  HR = exp(original_co["phenotype_albuminP1", "coef"]),
  lower_95 = exp(original_ci["phenotype_albuminP1", 1]),
  upper_95 = exp(original_ci["phenotype_albuminP1", 2]),
  p_value = original_co["phenotype_albuminP1", "Pr(>|z|)"]
)

bootstrap_replicates <- 1000L
bootstrap_rows <- vector("list", bootstrap_replicates)
set.seed(20260710)
bootstrap_seeds <- sample.int(.Machine$integer.max, bootstrap_replicates)

message("Running ", bootstrap_replicates, " cluster-aware bootstrap replicates...")
for (b in seq_len(bootstrap_replicates)) {
  idx <- sample.int(nrow(dat), size = nrow(dat), replace = TRUE)
  boot <- dat[idx, , drop = FALSE]
  bootstrap_rows[[b]] <- tryCatch({
    boot$phenotype_bootstrap <- cluster_and_label(boot, bootstrap_seeds[b])
    fit <- survival::coxph(
      stats::as.formula(paste(
        "survival::Surv(PERMTH_INT, MORTSTAT) ~ phenotype_bootstrap +", base_terms
      )),
      data = boot,
      ties = "efron"
    )
    co <- stats::coef(fit)
    tibble::tibble(
      replicate = b,
      log_HR = unname(co["phenotype_bootstrapP1"]),
      HR = exp(log_HR),
      valid = is.finite(log_HR)
    )
  }, error = function(e) {
    tibble::tibble(replicate = b, log_HR = NA_real_, HR = NA_real_, valid = FALSE)
  })
  if (b %% 100L == 0L) message("  completed bootstrap replicates: ", b)
}

bootstrap_distribution <- dplyr::bind_rows(bootstrap_rows)
valid_hr <- bootstrap_distribution$HR[bootstrap_distribution$valid]
bootstrap_summary <- tibble::tibble(
  requested_replicates = bootstrap_replicates,
  valid_replicates = length(valid_hr),
  valid_percent = 100 * length(valid_hr) / bootstrap_replicates,
  bootstrap_median_HR = stats::median(valid_hr),
  bootstrap_lower_95 = stats::quantile(valid_hr, 0.025, names = FALSE),
  bootstrap_upper_95 = stats::quantile(valid_hr, 0.975, names = FALSE),
  bootstrap_q1 = stats::quantile(valid_hr, 0.25, names = FALSE),
  bootstrap_q3 = stats::quantile(valid_hr, 0.75, names = FALSE)
)

readr::write_csv(ph_table, file.path(output_dir, "Table32A_corrected_PH_checks.csv"))
readr::write_csv(original_p1, file.path(output_dir, "Table32B_corrected_primary_P1.csv"))
readr::write_csv(bootstrap_summary, file.path(output_dir, "Table32C_cluster_aware_bootstrap_summary.csv"))
readr::write_csv(
  bootstrap_distribution,
  file.path(output_dir, "Table32D_cluster_aware_bootstrap_replicates.csv")
)
saveRDS(
  list(
    primary_fit = primary_fit,
    albumin_fit = albumin_fit,
    score_fit = score_fit,
    ph = ph_table,
    original_p1 = original_p1,
    bootstrap_summary = bootstrap_summary,
    bootstrap_distribution = bootstrap_distribution
  ),
  file.path(output_dir, "NHANES_corrected_bootstrap_results.rds")
)

summary_lines <- c(
  "Corrected NHANES cluster-aware bootstrap and PH diagnostics",
  paste0("Analytical sample: n = ", nrow(dat), "; events = ", sum(dat$MORTSTAT == 1)),
  "",
  "Primary P1 estimate:",
  paste(capture.output(print(original_p1)), collapse = "\n"),
  "",
  "Cluster-aware bootstrap:",
  paste(capture.output(print(bootstrap_summary)), collapse = "\n"),
  "",
  "Proportional-hazards checks:",
  paste(capture.output(print(ph_table)), collapse = "\n")
)
writeLines(summary_lines, file.path(output_dir, "NHANES_corrected_bootstrap_summary.txt"))
message("Corrected NHANES bootstrap completed.")
cat(paste(summary_lines, collapse = "\n"), "\n")
