# ==============================================================================
# NHANES albumin sensitivity and conventional immunonutritional benchmarks
# ==============================================================================

required_packages <- c("dplyr", "readr", "tibble", "survival", "survey")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop("Missing required package(s): ", paste(missing_packages, collapse = ", "), call. = FALSE)
}
options(survey.lonely.psu = "adjust")

project_root <- normalizePath(
  Sys.getenv("PROJECT_ROOT", unset = getwd()), winslash = "/", mustWork = TRUE
)
cohort_path <- file.path(
  project_root, "output", "nhanes_covariate_upgrade",
  "NHANES_2011_2018_covariate_augmented.rds"
)
assignment_path <- file.path(
  project_root, "output", "nhanes_robust_reanalysis",
  "NHANES_2011_2018_variant_assignments.csv"
)
output_dir <- file.path(project_root, "output", "nhanes_albumin_benchmarks")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(cohort_path) || !file.exists(assignment_path)) {
  stop("Run scripts 27, 27B, and 28 before this analysis.", call. = FALSE)
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
  ifelse(is.na(p), NA_character_, ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
}

cohort <- readRDS(cohort_path)
assignments <- readr::read_csv(assignment_path, show_col_types = FALSE) |>
  dplyr::select(SEQN, phenotype_log6_winsor)

dat <- cohort |>
  dplyr::left_join(assignments, by = "SEQN") |>
  dplyr::mutate(
    phenotype_total_protein = factor(
      phenotype_log6_winsor, levels = c("P3", "P2", "P1")
    )
  )

required <- c(
  "LBXSAL", "BMXWT", "BMXHT", "LBDLYMNO", "LBXPLTSI", "LBXHGB",
  "NLR", "SII", "LBXSCR", "BMXBMI", "WTMEC8YR", "SDMVPSU", "SDMVSTRA"
)
missing <- setdiff(required, names(dat))
if (length(missing) > 0L) {
  stop("Missing required column(s): ", paste(missing, collapse = ", "), call. = FALSE)
}

dat <- dat |>
  dplyr::mutate(
    albumin_g_dl = as.numeric(LBXSAL),
    height_m = as.numeric(BMXHT) / 100,
    ideal_weight_kg = 22 * height_m^2,
    weight_to_ideal_ratio = pmin(as.numeric(BMXWT) / ideal_weight_kg, 1),
    PNI = 10 * albumin_g_dl + 5 * as.numeric(LBDLYMNO),
    GNRI = 14.89 * albumin_g_dl + 41.7 * weight_to_ideal_ratio,
    HALP = (as.numeric(LBXHGB) * 10) * (albumin_g_dl * 10) *
      as.numeric(LBDLYMNO) / as.numeric(LBXPLTSI)
  )

albumin_complete <- dat |>
  dplyr::filter(
    is.finite(NLR), is.finite(SII), is.finite(LBXHGB), is.finite(albumin_g_dl),
    is.finite(BMXBMI), is.finite(LBXSCR), albumin_g_dl > 0
  )

albumin_matrix <- albumin_complete |>
  dplyr::transmute(
    log_nlr = log(winsorise(NLR)),
    log_sii = log(winsorise(SII)),
    haemoglobin = winsorise(LBXHGB),
    albumin = winsorise(albumin_g_dl),
    bmi = winsorise(BMXBMI),
    log_creatinine = log(winsorise(LBXSCR))
  ) |>
  as.data.frame() |>
  scale()

set.seed(20260710)
albumin_fit <- stats::kmeans(
  albumin_matrix, centers = 3, nstart = 100, iter.max = 500, algorithm = "Lloyd"
)

albumin_profiles <- albumin_complete |>
  dplyr::mutate(cluster_raw = albumin_fit$cluster) |>
  dplyr::group_by(cluster_raw) |>
  dplyr::summarise(
    n = dplyr::n(),
    nlr = stats::median(NLR),
    sii = stats::median(SII),
    haemoglobin = stats::median(LBXHGB),
    albumin = stats::median(albumin_g_dl),
    bmi = stats::median(BMXBMI),
    creatinine = stats::median(LBXSCR),
    .groups = "drop"
  )

p1_score <- safe_z(log(albumin_profiles$nlr)) +
  safe_z(log(albumin_profiles$sii)) + safe_z(log(albumin_profiles$creatinine))
p1_raw <- albumin_profiles$cluster_raw[which.max(p1_score)]
remaining <- albumin_profiles |>
  dplyr::filter(cluster_raw != p1_raw)
p2_score <- -safe_z(remaining$haemoglobin) - safe_z(remaining$albumin) - safe_z(remaining$bmi)
p2_raw <- remaining$cluster_raw[which.max(p2_score)]
p3_raw <- setdiff(albumin_profiles$cluster_raw, c(p1_raw, p2_raw))

albumin_complete <- albumin_complete |>
  dplyr::mutate(
    phenotype_albumin = dplyr::case_when(
      albumin_fit$cluster == p1_raw ~ "P1",
      albumin_fit$cluster == p2_raw ~ "P2",
      albumin_fit$cluster == p3_raw ~ "P3",
      TRUE ~ NA_character_
    ),
    phenotype_albumin = factor(phenotype_albumin, levels = c("P3", "P2", "P1"))
  )

albumin_assignments <- albumin_complete |>
  dplyr::select(SEQN, phenotype_albumin)
dat <- dat |>
  dplyr::left_join(albumin_assignments, by = "SEQN")

dat <- dat |>
  dplyr::mutate(
    inflammation_axis = rowMeans(
      cbind(safe_z(log(winsorise(NLR))), safe_z(log(winsorise(SII)))),
      na.rm = FALSE
    ),
    reserve_depletion_axis = rowMeans(
      cbind(
        -safe_z(winsorise(LBXHGB)),
        -safe_z(winsorise(albumin_g_dl)),
        -safe_z(winsorise(BMXBMI))
      ),
      na.rm = FALSE
    ),
    renal_stress_axis = safe_z(log(winsorise(LBXSCR))),
    domain_balanced_score = safe_z(
      rowMeans(
        cbind(inflammation_axis, reserve_depletion_axis, renal_stress_axis),
        na.rm = FALSE
      )
    ),
    risk_log_nlr = safe_z(log(winsorise(NLR))),
    risk_log_sii = safe_z(log(winsorise(SII))),
    risk_low_pni = -safe_z(PNI),
    risk_low_gnri = -safe_z(GNRI),
    risk_low_halp = -safe_z(log(winsorise(HALP))),
    male = as.integer(RIAGENDR == 1),
    race = factor(RIDRETH3),
    education = factor(DMDEDUC2),
    cycle = factor(Cycle_ID)
  )

predictor_terms <- c(
  phenotype_total_protein = "phenotype_total_protein",
  phenotype_albumin = "phenotype_albumin",
  domain_balanced_score = "domain_balanced_score",
  NLR = "risk_log_nlr",
  SII = "risk_log_sii",
  PNI = "risk_low_pni",
  GNRI = "risk_low_gnri",
  HALP = "risk_low_halp"
)

common_required <- c(
  "PERMTH_INT", "MORTSTAT", "RIDAGEYR", "male", "race", "INDFMPIR",
  "Comorbidity_Score_Extended", "cycle", "smoking", "hypertension", "diabetes",
  unname(predictor_terms)
)
model_data <- dat |>
  dplyr::filter(dplyr::if_all(dplyr::all_of(common_required), ~ !is.na(.x)))

base_formula <- survival::Surv(PERMTH_INT, MORTSTAT) ~
  RIDAGEYR + male + race + INDFMPIR + Comorbidity_Score_Extended + cycle +
  smoking + hypertension + diabetes
base_fit <- survival::coxph(base_formula, data = model_data, ties = "efron")

benchmark_rows <- lapply(names(predictor_terms), function(label) {
  term <- predictor_terms[[label]]
  fit <- survival::coxph(
    stats::update.formula(base_formula, paste(". ~ . +", term)),
    data = model_data,
    ties = "efron"
  )
  lr <- stats::anova(base_fit, fit, test = "LRT")
  tibble::tibble(
    predictor = label,
    n = nrow(model_data),
    events = sum(model_data$MORTSTAT == 1),
    c_index = unname(summary(fit)$concordance[1]),
    delta_c_index = unname(summary(fit)$concordance[1] - summary(base_fit)$concordance[1]),
    AIC = stats::AIC(fit),
    delta_AIC_vs_base = stats::AIC(fit) - stats::AIC(base_fit),
    likelihood_ratio_chisq = lr$Chisq[2],
    likelihood_ratio_df = lr$Df[2],
    likelihood_ratio_p = lr$`Pr(>|Chi|)`[2]
  )
}) |>
  dplyr::bind_rows()

incremental_rows <- lapply(c("risk_log_nlr", "risk_log_sii", "risk_low_pni", "risk_low_gnri", "risk_low_halp"), function(term) {
  comparator_fit <- survival::coxph(
    stats::update.formula(base_formula, paste(". ~ . +", term)),
    data = model_data,
    ties = "efron"
  )
  combined_fit <- survival::coxph(
    stats::update.formula(base_formula, paste(". ~ . +", term, "+ phenotype_total_protein")),
    data = model_data,
    ties = "efron"
  )
  lr <- stats::anova(comparator_fit, combined_fit, test = "LRT")
  co <- summary(combined_fit)$coefficients
  ci <- suppressMessages(stats::confint(combined_fit))
  p_col <- grep("Pr\\(", colnames(co), value = TRUE)[1]
  p1_term <- "phenotype_total_proteinP1"
  tibble::tibble(
    comparator = term,
    P1_HR = exp(co[p1_term, "coef"]),
    P1_lower_95 = exp(ci[p1_term, 1]),
    P1_upper_95 = exp(ci[p1_term, 2]),
    P1_p = co[p1_term, p_col],
    P1_hazard_ratio_95ci = sprintf(
      "%.2f (%.2f-%.2f)", exp(co[p1_term, "coef"]), exp(ci[p1_term, 1]), exp(ci[p1_term, 2])
    ),
    likelihood_ratio_p_for_adding_phenotype = lr$`Pr(>|Chi|)`[2]
  )
}) |>
  dplyr::bind_rows()

survey_data <- model_data
survey_design <- survey::svydesign(
  ids = ~SDMVPSU,
  strata = ~SDMVSTRA,
  weights = ~WTMEC8YR,
  nest = TRUE,
  data = survey_data
)
survey_predictors <- c(
  "phenotype_total_protein", "phenotype_albumin", "domain_balanced_score",
  "risk_log_nlr", "risk_log_sii", "risk_low_pni", "risk_low_gnri", "risk_low_halp"
)
survey_rows <- lapply(survey_predictors, function(term) {
  fit <- survey::svycoxph(
    stats::update.formula(
      survival::Surv(PERMTH_INT, MORTSTAT) ~
        RIDAGEYR + male + race + INDFMPIR + Comorbidity_Score_Extended +
        smoking + hypertension + diabetes,
      paste(". ~ . +", term)
    ),
    design = survey_design
  )
  co <- summary(fit)$coefficients
  ci <- suppressMessages(stats::confint(fit))
  p_col <- grep("Pr\\(", colnames(co), value = TRUE)[1]
  terms <- rownames(co)
  keep <- if (grepl("^phenotype", term)) grepl(paste0("^", term), terms) else terms == term
  tibble::tibble(
    predictor = term,
    comparison = terms[keep],
    HR = exp(co[keep, "coef"]),
    lower_95 = exp(ci[keep, 1]),
    upper_95 = exp(ci[keep, 2]),
    p_value = co[keep, p_col],
    hazard_ratio_95ci = sprintf(
      "%.2f (%.2f-%.2f)", exp(co[keep, "coef"]), exp(ci[keep, 1]), exp(ci[keep, 2])
    ),
    p_value_formatted = format_p(co[keep, p_col])
  )
}) |>
  dplyr::bind_rows()

albumin_counts <- albumin_complete |>
  dplyr::group_by(phenotype_albumin) |>
  dplyr::summarise(
    n = dplyr::n(),
    deaths = sum(MORTSTAT == 1),
    mortality_percent = 100 * mean(MORTSTAT == 1),
    .groups = "drop"
  )

agreement_data <- dat |>
  dplyr::filter(!is.na(phenotype_total_protein), !is.na(phenotype_albumin))
agreement <- tibble::tibble(
  common_n = nrow(agreement_data),
  exact_label_agreement = mean(
    agreement_data$phenotype_total_protein == agreement_data$phenotype_albumin
  ),
  adjusted_rand_index = adjusted_rand_index(
    agreement_data$phenotype_total_protein, agreement_data$phenotype_albumin
  )
)

readr::write_csv(albumin_profiles, file.path(output_dir, "Table29A_albumin_cluster_profiles.csv"))
readr::write_csv(albumin_counts, file.path(output_dir, "Table29B_albumin_cluster_counts.csv"))
readr::write_csv(agreement, file.path(output_dir, "Table29C_total_protein_vs_albumin_agreement.csv"))
readr::write_csv(benchmark_rows, file.path(output_dir, "Table29D_conventional_index_benchmarks.csv"))
readr::write_csv(incremental_rows, file.path(output_dir, "Table29E_P1_beyond_conventional_indices.csv"))
readr::write_csv(survey_rows, file.path(output_dir, "Table29F_survey_weighted_benchmarks.csv"))
readr::write_csv(
  dat |>
    dplyr::select(
      SEQN, phenotype_total_protein, phenotype_albumin, PNI, GNRI, HALP,
      domain_balanced_score
    ),
  file.path(output_dir, "NHANES_albumin_benchmark_assignments.csv")
)
saveRDS(
  list(
    model_data = model_data,
    albumin_profiles = albumin_profiles,
    albumin_counts = albumin_counts,
    agreement = agreement,
    benchmarks = benchmark_rows,
    incremental = incremental_rows,
    survey = survey_rows
  ),
  file.path(output_dir, "NHANES_albumin_benchmark_results.rds")
)

summary_lines <- c(
  "NHANES albumin sensitivity and conventional benchmark analysis",
  paste0("Common benchmark cohort: n = ", nrow(model_data), "; events = ", sum(model_data$MORTSTAT == 1)),
  "",
  "Total-protein versus albumin phenotype agreement:",
  paste(capture.output(print(agreement)), collapse = "\n"),
  "",
  "Conventional index benchmark models:",
  paste(capture.output(print(benchmark_rows)), collapse = "\n"),
  "",
  "P1 incremental value beyond conventional indices:",
  paste(capture.output(print(incremental_rows)), collapse = "\n")
)
writeLines(summary_lines, file.path(output_dir, "NHANES_albumin_benchmark_summary.txt"))
message("NHANES albumin and conventional benchmark analysis completed.")
cat(paste(summary_lines, collapse = "\n"), "\n")
