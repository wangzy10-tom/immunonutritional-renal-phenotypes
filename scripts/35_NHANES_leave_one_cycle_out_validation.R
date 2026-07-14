# ==============================================================================
# NHANES leave-one-cycle-out phenotype transportability analysis
# ==============================================================================

required_packages <- c("dplyr", "readr", "tibble", "survival", "survey")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop("Missing required package(s): ", paste(missing_packages, collapse = ", "), call. = FALSE)
}
options(survey.lonely.psu = "adjust")

root <- normalizePath(
  Sys.getenv("PROJECT_ROOT", unset = getwd()), winslash = "/", mustWork = TRUE
)
cohort_path <- file.path(
  root, "output", "nhanes_covariate_upgrade",
  "NHANES_2011_2018_covariate_augmented.rds"
)
assignment_path <- file.path(
  root, "output", "nhanes_robust_reanalysis",
  "NHANES_2011_2018_variant_assignments.csv"
)
output_dir <- file.path(root, "output", "nhanes_leave_one_cycle_out")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

winsor_limits <- function(x) {
  stats::quantile(x, c(0.01, 0.99), na.rm = TRUE, names = FALSE)
}

apply_winsor <- function(x, limits) {
  pmin(pmax(as.numeric(x), limits[1]), limits[2])
}

z_small <- function(x) {
  spread <- stats::sd(x)
  if (!is.finite(spread) || spread == 0) return(rep(0, length(x)))
  (x - mean(x)) / spread
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

fit_and_project <- function(train, test) {
  raw_variables <- c("NLR", "SII", "LBXHGB", "LBXSTP", "BMXBMI", "LBXSCR")
  limits <- lapply(raw_variables, function(variable) winsor_limits(train[[variable]]))
  names(limits) <- raw_variables

  transform_matrix <- function(dat) {
    cbind(
      log_nlr = log(pmax(apply_winsor(dat$NLR, limits$NLR), .Machine$double.eps)),
      log_sii = log(pmax(apply_winsor(dat$SII, limits$SII), .Machine$double.eps)),
      haemoglobin = apply_winsor(dat$LBXHGB, limits$LBXHGB),
      total_protein = apply_winsor(dat$LBXSTP, limits$LBXSTP),
      bmi = apply_winsor(dat$BMXBMI, limits$BMXBMI),
      log_creatinine = log(pmax(apply_winsor(dat$LBXSCR, limits$LBXSCR), .Machine$double.eps))
    )
  }

  train_unscaled <- transform_matrix(train)
  centres <- colMeans(train_unscaled)
  spreads <- apply(train_unscaled, 2, stats::sd)
  spreads[!is.finite(spreads) | spreads == 0] <- 1
  train_scaled <- sweep(sweep(train_unscaled, 2, centres, "-"), 2, spreads, "/")

  cluster_fit <- stats::kmeans(
    train_scaled, centers = 3, nstart = 100, iter.max = 500, algorithm = "Lloyd"
  )

  profiles <- train |>
    dplyr::mutate(cluster_raw = cluster_fit$cluster) |>
    dplyr::group_by(cluster_raw) |>
    dplyr::summarise(
      nlr = stats::median(NLR),
      sii = stats::median(SII),
      haemoglobin = stats::median(LBXHGB),
      total_protein = stats::median(LBXSTP),
      bmi = stats::median(BMXBMI),
      creatinine = stats::median(LBXSCR),
      .groups = "drop"
    )

  p1_score <- z_small(log(profiles$nlr)) + z_small(log(profiles$sii)) +
    z_small(log(profiles$creatinine))
  p1_raw <- profiles$cluster_raw[which.max(p1_score)]
  remaining <- profiles |>
    dplyr::filter(cluster_raw != p1_raw)
  p2_score <- -z_small(remaining$haemoglobin) - z_small(remaining$total_protein) -
    z_small(remaining$bmi)
  p2_raw <- remaining$cluster_raw[which.max(p2_score)]
  p3_raw <- setdiff(profiles$cluster_raw, c(p1_raw, p2_raw))
  label_map <- stats::setNames(c("P1", "P2", "P3"), c(p1_raw, p2_raw, p3_raw))

  assign_test <- function(dat) {
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

  list(projected = assign_test(test), profiles = profiles)
}

extract_cox_phenotypes <- function(fit, held_out_cycle, model_label) {
  beta <- stats::coef(fit)
  se <- sqrt(diag(stats::vcov(fit)))
  terms <- names(beta)
  keep <- grepl("^phenotype_loo", terms)
  tibble::tibble(
    held_out_cycle = held_out_cycle,
    model = model_label,
    comparison = dplyr::recode(
      terms[keep], phenotype_looP2 = "P2 vs P3", phenotype_looP1 = "P1 vs P3"
    ),
    log_HR = beta[keep],
    standard_error = se[keep],
    HR = exp(beta[keep]),
    lower_95 = exp(beta[keep] - 1.96 * se[keep]),
    upper_95 = exp(beta[keep] + 1.96 * se[keep]),
    p_value = 2 * stats::pnorm(abs(beta[keep] / se[keep]), lower.tail = FALSE)
  ) |>
    dplyr::mutate(
      hazard_ratio_95ci = sprintf("%.2f (%.2f-%.2f)", HR, lower_95, upper_95)
    )
}

cohort <- readRDS(cohort_path)
full_assignments <- readr::read_csv(assignment_path, show_col_types = FALSE) |>
  dplyr::select(SEQN, phenotype_log6_winsor)

dat <- cohort |>
  dplyr::left_join(full_assignments, by = "SEQN") |>
  dplyr::mutate(
    phenotype_full = factor(phenotype_log6_winsor, levels = c("P3", "P2", "P1")),
    male = as.integer(RIAGENDR == 1),
    race = factor(RIDRETH3),
    cycle = factor(Cycle_ID)
  )

cycles <- sort(unique(dat$Cycle_ID))
set.seed(20260710)
held_out_rows <- vector("list", length(cycles))
agreement_rows <- vector("list", length(cycles))
count_rows <- vector("list", length(cycles))
cox_rows <- vector("list", length(cycles))

for (i in seq_along(cycles)) {
  held_out_cycle <- cycles[i]
  train <- dat |>
    dplyr::filter(Cycle_ID != held_out_cycle)
  test <- dat |>
    dplyr::filter(Cycle_ID == held_out_cycle)

  projector <- fit_and_project(train, test)
  test$phenotype_loo <- projector$projected
  test$held_out_cycle <- held_out_cycle
  held_out_rows[[i]] <- test

  agreement_rows[[i]] <- tibble::tibble(
    held_out_cycle = held_out_cycle,
    held_out_n = nrow(test),
    exact_label_agreement = mean(test$phenotype_loo == test$phenotype_full),
    adjusted_rand_index = adjusted_rand_index(test$phenotype_loo, test$phenotype_full)
  )

  count_rows[[i]] <- test |>
    dplyr::group_by(phenotype_loo) |>
    dplyr::summarise(
      n = dplyr::n(),
      deaths = sum(MORTSTAT == 1),
      mortality_percent = 100 * mean(MORTSTAT == 1),
      .groups = "drop"
    ) |>
    dplyr::mutate(held_out_cycle = held_out_cycle, .before = 1)

  cycle_model_data <- test |>
    dplyr::filter(
      !is.na(PERMTH_INT), !is.na(MORTSTAT), !is.na(phenotype_loo),
      !is.na(RIDAGEYR), !is.na(male), !is.na(Comorbidity_Score_Extended)
    )
  cycle_fit <- survival::coxph(
    survival::Surv(PERMTH_INT, MORTSTAT) ~
      phenotype_loo + RIDAGEYR + male + Comorbidity_Score_Extended,
    data = cycle_model_data,
    ties = "efron"
  )
  cox_rows[[i]] <- extract_cox_phenotypes(
    cycle_fit, held_out_cycle, "Held-out cycle age-sex-comorbidity model"
  ) |>
    dplyr::mutate(n_model = nrow(cycle_model_data), events = sum(cycle_model_data$MORTSTAT == 1))
}

loo_data <- dplyr::bind_rows(held_out_rows) |>
  dplyr::mutate(phenotype_loo = factor(phenotype_loo, levels = c("P3", "P2", "P1")))
agreement <- dplyr::bind_rows(agreement_rows)
counts <- dplyr::bind_rows(count_rows)
cycle_cox <- dplyr::bind_rows(cox_rows)

pooled_data <- loo_data |>
  dplyr::filter(
    !is.na(PERMTH_INT), !is.na(MORTSTAT), !is.na(phenotype_loo),
    !is.na(RIDAGEYR), !is.na(male), !is.na(Comorbidity_Score_Extended), !is.na(cycle)
  )
pooled_fit <- survival::coxph(
  survival::Surv(PERMTH_INT, MORTSTAT) ~
    phenotype_loo + RIDAGEYR + male + Comorbidity_Score_Extended + cycle,
  data = pooled_data,
  ties = "efron"
)
pooled_results <- extract_cox_phenotypes(
  pooled_fit, "All", "Pooled leave-one-cycle-out model"
) |>
  dplyr::mutate(n_model = nrow(pooled_data), events = sum(pooled_data$MORTSTAT == 1))

survey_data <- loo_data |>
  dplyr::filter(
    !is.na(PERMTH_INT), !is.na(MORTSTAT), !is.na(phenotype_loo),
    !is.na(RIDAGEYR), !is.na(male), !is.na(race), !is.na(INDFMPIR),
    !is.na(Comorbidity_Score_Extended), !is.na(cycle), !is.na(smoking),
    !is.na(hypertension), !is.na(diabetes)
  )
survey_design <- survey::svydesign(
  ids = ~SDMVPSU,
  strata = ~SDMVSTRA,
  weights = ~WTMEC8YR,
  nest = TRUE,
  data = survey_data
)
survey_fit <- survey::svycoxph(
  survival::Surv(PERMTH_INT, MORTSTAT) ~
    phenotype_loo + RIDAGEYR + male + race + INDFMPIR +
    Comorbidity_Score_Extended + cycle + smoking + hypertension + diabetes,
  design = survey_design
)
survey_results <- extract_cox_phenotypes(
  survey_fit, "All", "Survey-weighted pooled leave-one-cycle-out model"
) |>
  dplyr::mutate(n_model = nrow(survey_data), events = sum(survey_data$MORTSTAT == 1))

interaction_data <- dat |>
  dplyr::filter(
    !is.na(PERMTH_INT), !is.na(MORTSTAT), !is.na(phenotype_full),
    !is.na(RIDAGEYR), !is.na(male), !is.na(Comorbidity_Score_Extended), !is.na(cycle)
  )
interaction_main <- survival::coxph(
  survival::Surv(PERMTH_INT, MORTSTAT) ~
    phenotype_full + RIDAGEYR + male + Comorbidity_Score_Extended + cycle,
  data = interaction_data,
  ties = "efron"
)
interaction_fit <- survival::coxph(
  survival::Surv(PERMTH_INT, MORTSTAT) ~
    phenotype_full * cycle + RIDAGEYR + male + Comorbidity_Score_Extended,
  data = interaction_data,
  ties = "efron"
)
interaction_test_raw <- stats::anova(interaction_main, interaction_fit, test = "LRT")
interaction_test <- tibble::tibble(
  comparison = "Phenotype-by-cycle interaction",
  n = nrow(interaction_data),
  events = sum(interaction_data$MORTSTAT == 1),
  likelihood_ratio_chisq = interaction_test_raw$Chisq[2],
  df = interaction_test_raw$Df[2],
  p_value = interaction_test_raw$`Pr(>|Chi|)`[2]
)

p1_cycle <- cycle_cox |>
  dplyr::filter(comparison == "P1 vs P3")
meta_weight <- 1 / p1_cycle$standard_error^2
pooled_log_hr <- sum(meta_weight * p1_cycle$log_HR) / sum(meta_weight)
q_statistic <- sum(meta_weight * (p1_cycle$log_HR - pooled_log_hr)^2)
q_df <- nrow(p1_cycle) - 1
heterogeneity <- tibble::tibble(
  cycles = nrow(p1_cycle),
  fixed_effect_HR = exp(pooled_log_hr),
  fixed_effect_lower_95 = exp(pooled_log_hr - 1.96 / sqrt(sum(meta_weight))),
  fixed_effect_upper_95 = exp(pooled_log_hr + 1.96 / sqrt(sum(meta_weight))),
  Q = q_statistic,
  df = q_df,
  p_heterogeneity = stats::pchisq(q_statistic, df = q_df, lower.tail = FALSE),
  I2_percent = max(0, 100 * (q_statistic - q_df) / q_statistic)
)

readr::write_csv(agreement, file.path(output_dir, "Table35A_cycle_projection_agreement.csv"))
readr::write_csv(counts, file.path(output_dir, "Table35B_heldout_cycle_counts.csv"))
readr::write_csv(cycle_cox, file.path(output_dir, "Table35C_heldout_cycle_Cox.csv"))
readr::write_csv(
  dplyr::bind_rows(pooled_results, survey_results),
  file.path(output_dir, "Table35D_pooled_LOCO_Cox.csv")
)
readr::write_csv(interaction_test, file.path(output_dir, "Table35E_cycle_interaction.csv"))
readr::write_csv(heterogeneity, file.path(output_dir, "Table35F_P1_cycle_heterogeneity.csv"))
readr::write_csv(
  loo_data |>
    dplyr::select(SEQN, Cycle_ID, phenotype_full, phenotype_loo, MORTSTAT, PERMTH_INT),
  file.path(output_dir, "NHANES_LOCO_assignments.csv")
)
saveRDS(
  list(
    agreement = agreement,
    counts = counts,
    cycle_cox = cycle_cox,
    pooled = pooled_results,
    survey = survey_results,
    interaction = interaction_test,
    heterogeneity = heterogeneity
  ),
  file.path(output_dir, "NHANES_LOCO_results.rds")
)

summary_lines <- c(
  "NHANES leave-one-cycle-out phenotype validation",
  "",
  "Agreement with full-data robust solution:",
  paste(capture.output(print(agreement)), collapse = "\n"),
  "",
  "Held-out cycle phenotype models:",
  paste(capture.output(print(cycle_cox)), collapse = "\n"),
  "",
  "Pooled models:",
  paste(capture.output(print(dplyr::bind_rows(pooled_results, survey_results))), collapse = "\n"),
  "",
  "Phenotype-by-cycle interaction:",
  paste(capture.output(print(interaction_test)), collapse = "\n"),
  "",
  "P1 cycle heterogeneity:",
  paste(capture.output(print(heterogeneity)), collapse = "\n")
)
writeLines(summary_lines, file.path(output_dir, "NHANES_LOCO_summary.txt"))
cat(paste(summary_lines, collapse = "\n"), "\n")
