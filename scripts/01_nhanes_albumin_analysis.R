# ==============================================================================
# NHANES Amendment 17: unified albumin phenotype analysis suite
# Date: 2026-08-30
# ==============================================================================

required_packages <- c(
  "dplyr", "readr", "tibble", "tidyr", "survival", "survey", "cluster", "digest", "haven"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop("Missing required package(s): ", paste(missing_packages, collapse = ", "), call. = FALSE)
}
options(survey.lonely.psu = "adjust")

root <- normalizePath(Sys.getenv("PROJECT_ROOT", unset = getwd()), winslash = "/", mustWork = TRUE)
script_dir <- file.path(root, "scripts")
output_root <- file.path(root, "output", "nhanes_albumin_amendment17_2026-08-30")
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

eligible_path <- Sys.getenv("NHANES_ELIGIBLE_RDS", unset = file.path(root, "data", "derived", "NHANES_2011_2018_age65_mortality_eligible.rds"))
covariate_cache <- Sys.getenv("NHANES_COVARIATE_CACHE", unset = file.path(root, "data", "nhanes_xpt_cache"))
stopifnot(file.exists(eligible_path), dir.exists(covariate_cache))

winsorise <- function(x, lower = 0.01, upper = 0.99) {
  limits <- stats::quantile(x, c(lower, upper), na.rm = TRUE, names = FALSE)
  pmin(pmax(as.numeric(x), limits[1]), limits[2])
}
safe_z <- function(x) {
  spread <- stats::sd(x, na.rm = TRUE)
  if (!is.finite(spread) || spread == 0) return(rep(0, length(x)))
  (x - mean(x, na.rm = TRUE)) / spread
}
yes_no_flag <- function(x) {
  dplyr::case_when(x == 1 ~ 1L, x == 2 ~ 0L, TRUE ~ NA_integer_)
}

message("[1/9] Rebuilding the albumin-only cohort from the mortality-eligible denominator...")
mcq_variables <- c(
  "MCQ220", "MCQ160A", "MCQ160B", "MCQ160C", "MCQ160D", "MCQ160E",
  "MCQ160F", "MCQ160G", "MCQ160K", "MCQ160L", "MCQ160M", "MCQ160N"
)
complete_cohort <- readRDS(eligible_path) |>
  dplyr::filter(
    is.finite(NLR), NLR > 0, is.finite(SII), SII > 0,
    is.finite(LBXHGB), is.finite(LBXSAL), LBXSAL > 0,
    is.finite(BMXBMI), is.finite(LBXSCR), LBXSCR > 0
  ) |>
  dplyr::mutate(
    dplyr::across(dplyr::all_of(mcq_variables), yes_no_flag, .names = "flag_{.col}")
  )
flag_variables <- paste0("flag_", mcq_variables)
complete_cohort <- complete_cohort |>
  dplyr::mutate(
    Comorbidity_Unknown_Count = rowSums(is.na(dplyr::pick(dplyr::all_of(flag_variables)))),
    Comorbidity_Score_Extended = dplyr::if_else(
      Comorbidity_Unknown_Count == 0,
      rowSums(dplyr::pick(dplyr::all_of(flag_variables))),
      NA_real_
    )
  )

cycle_covariates <- dplyr::bind_rows(lapply(c("G", "H", "I", "J"), function(cycle) {
  bpq <- haven::read_xpt(file.path(covariate_cache, paste0("BPQ_", cycle, ".xpt"))) |>
    dplyr::select(SEQN, BPQ020)
  diq <- haven::read_xpt(file.path(covariate_cache, paste0("DIQ_", cycle, ".xpt"))) |>
    dplyr::select(SEQN, DIQ010)
  smq <- haven::read_xpt(file.path(covariate_cache, paste0("SMQ_", cycle, ".xpt"))) |>
    dplyr::select(SEQN, SMQ020, SMQ040)
  dplyr::full_join(bpq, diq, by = "SEQN") |>
    dplyr::full_join(smq, by = "SEQN") |>
    dplyr::mutate(Cycle_ID = cycle)
}))
complete_cohort <- complete_cohort |>
  dplyr::left_join(cycle_covariates, by = c("SEQN", "Cycle_ID")) |>
  dplyr::mutate(
    hypertension = factor(dplyr::case_when(BPQ020 == 1 ~ "Yes", BPQ020 == 2 ~ "No", TRUE ~ NA_character_), levels = c("No", "Yes")),
    diabetes = factor(dplyr::case_when(DIQ010 == 1 ~ "Yes", DIQ010 == 3 ~ "Borderline", DIQ010 == 2 ~ "No", TRUE ~ NA_character_), levels = c("No", "Borderline", "Yes")),
    smoking = factor(dplyr::case_when(
      SMQ020 == 2 ~ "Never", SMQ020 == 1 & SMQ040 %in% c(1, 2) ~ "Current",
      SMQ020 == 1 & SMQ040 == 3 ~ "Former", TRUE ~ NA_character_
    ), levels = c("Never", "Former", "Current"))
  )

cluster_matrix <- complete_cohort |>
  dplyr::transmute(
    log_nlr = log(winsorise(NLR)), log_sii = log(winsorise(SII)),
    haemoglobin = winsorise(LBXHGB), albumin = winsorise(LBXSAL),
    bmi = winsorise(BMXBMI), log_creatinine = log(winsorise(LBXSCR))
  ) |>
  as.data.frame() |>
  scale()
set.seed(20260710)
cluster_fit <- stats::kmeans(
  cluster_matrix, centers = 3, nstart = 100, iter.max = 500, algorithm = "Lloyd"
)
albumin_profiles <- complete_cohort |>
  dplyr::mutate(cluster_raw = cluster_fit$cluster) |>
  dplyr::group_by(cluster_raw) |>
  dplyr::summarise(
    n = dplyr::n(), nlr = stats::median(NLR), sii = stats::median(SII),
    haemoglobin = stats::median(LBXHGB), albumin = stats::median(LBXSAL),
    bmi = stats::median(BMXBMI), creatinine = stats::median(LBXSCR), .groups = "drop"
  )
p1_score <- safe_z(log(albumin_profiles$nlr)) + safe_z(log(albumin_profiles$sii)) +
  safe_z(log(albumin_profiles$creatinine))
p1_raw <- albumin_profiles$cluster_raw[which.max(p1_score)]
remaining <- dplyr::filter(albumin_profiles, cluster_raw != p1_raw)
p2_score <- -safe_z(remaining$haemoglobin) - safe_z(remaining$albumin) - safe_z(remaining$bmi)
p2_raw <- remaining$cluster_raw[which.max(p2_score)]
p3_raw <- setdiff(albumin_profiles$cluster_raw, c(p1_raw, p2_raw))

complete_cohort <- complete_cohort |>
  dplyr::mutate(
    phenotype_albumin = factor(dplyr::case_when(
      cluster_fit$cluster == p1_raw ~ "P1", cluster_fit$cluster == p2_raw ~ "P2",
      cluster_fit$cluster == p3_raw ~ "P3", TRUE ~ NA_character_
    ), levels = c("P3", "P2", "P1")),
    albumin_g_dl = as.numeric(LBXSAL), height_m = as.numeric(BMXHT) / 100,
    ideal_weight_kg = 22 * height_m^2,
    weight_to_ideal_ratio = pmin(as.numeric(BMXWT) / ideal_weight_kg, 1),
    PNI = 10 * albumin_g_dl + 5 * as.numeric(LBDLYMNO),
    GNRI = 14.89 * albumin_g_dl + 41.7 * weight_to_ideal_ratio,
    HALP = (as.numeric(LBXHGB) * 10) * (albumin_g_dl * 10) * as.numeric(LBDLYMNO) / as.numeric(LBXPLTSI),
    inflammation_axis = rowMeans(cbind(safe_z(log(winsorise(NLR))), safe_z(log(winsorise(SII))))),
    reserve_depletion_axis = rowMeans(cbind(-safe_z(winsorise(LBXHGB)), -safe_z(winsorise(LBXSAL)), -safe_z(winsorise(BMXBMI)))),
    renal_stress_axis = safe_z(log(winsorise(LBXSCR))),
    domain_balanced_score = safe_z(rowMeans(cbind(inflammation_axis, reserve_depletion_axis, renal_stress_axis))),
    risk_log_nlr = safe_z(log(winsorise(NLR))), risk_log_sii = safe_z(log(winsorise(SII))),
    risk_low_pni = -safe_z(PNI), risk_low_gnri = -safe_z(GNRI),
    risk_low_halp = -safe_z(log(winsorise(HALP))),
    male = as.integer(RIAGENDR == 1), race = factor(RIDRETH3),
    education = factor(DMDEDUC2), cycle = factor(Cycle_ID)
  )

common_required <- c(
  "PERMTH_INT", "MORTSTAT", "RIDAGEYR", "male", "race", "INDFMPIR",
  "Comorbidity_Score_Extended", "cycle", "smoking", "hypertension", "diabetes",
  "phenotype_albumin", "domain_balanced_score", "risk_log_nlr", "risk_log_sii",
  "risk_low_pni", "risk_low_gnri", "risk_low_halp"
)
model_data <- complete_cohort |>
  dplyr::filter(dplyr::if_all(dplyr::all_of(common_required), ~ !is.na(.x)))
albumin_results <- list(model_data = model_data, albumin_profiles = albumin_profiles)
benchmark_path <- file.path(output_root, "NHANES_albumin_benchmark_results.rds")
assignment_path <- file.path(output_root, "NHANES_albumin_benchmark_assignments.csv")
saveRDS(albumin_results, benchmark_path)
readr::write_csv(
  dplyr::select(complete_cohort, SEQN, phenotype_albumin),
  assignment_path
)

stopifnot(
  nrow(complete_cohort) == 4637L,
  sum(complete_cohort$MORTSTAT == 1L) == 831L,
  nrow(model_data) == 3979L,
  sum(model_data$MORTSTAT == 1L) == 720L
)

saveRDS(
  complete_cohort,
  file.path(output_root, "INTERNAL_NHANES_albumin_complete_cohort.rds")
)

format_effect <- function(beta, se) {
  beta <- as.numeric(beta)
  se <- as.numeric(se)
  tibble::tibble(
    HR = exp(beta),
    lower_95 = exp(beta - 1.96 * se),
    upper_95 = exp(beta + 1.96 * se),
    p_value = 2 * stats::pnorm(abs(beta / se), lower.tail = FALSE)
  ) |>
    dplyr::mutate(hazard_ratio_95ci = sprintf("%.2f (%.2f-%.2f)", HR, lower_95, upper_95))
}

extract_phenotype <- function(fit, model_label) {
  beta <- stats::coef(fit)
  se <- sqrt(diag(stats::vcov(fit)))
  keep <- names(beta) %in% c("phenotype_albuminP2", "phenotype_albuminP1")
  format_effect(beta[keep], se[keep]) |>
    dplyr::mutate(
      model = model_label,
      comparison = dplyr::recode(
        names(beta)[keep],
        phenotype_albuminP2 = "P2 vs P3",
        phenotype_albuminP1 = "P1 vs P3"
      ),
      n = nrow(model_data),
      events = sum(model_data$MORTSTAT == 1L),
      .before = 1
    )
}

message("[2/9] Primary survey Cox model and weighted descriptive tables...")
survey_design <- survey::svydesign(
  ids = ~SDMVPSU,
  strata = ~SDMVSTRA,
  weights = ~WTMEC8YR,
  nest = TRUE,
  data = model_data
)
primary_formula <- survival::Surv(PERMTH_INT, MORTSTAT) ~
  phenotype_albumin + RIDAGEYR + male + race + INDFMPIR +
  Comorbidity_Score_Extended + cycle + smoking + hypertension + diabetes
primary_fit <- survey::svycoxph(primary_formula, design = survey_design)
primary_effects <- extract_phenotype(primary_fit, "Primary complex-survey Cox")

observed <- stats::setNames(primary_effects$HR, primary_effects$comparison)
if (!all(is.finite(observed[c("P2 vs P3", "P1 vs P3")]))) {
  stop("Primary albumin survey-Cox estimate is non-finite.", call. = FALSE)
}

phenotype_counts <- complete_cohort |>
  dplyr::group_by(phenotype_albumin) |>
  dplyr::summarise(n = dplyr::n(), deaths = sum(MORTSTAT == 1L), .groups = "drop")
analytic_counts <- model_data |>
  dplyr::group_by(phenotype_albumin) |>
  dplyr::summarise(n = dplyr::n(), deaths = sum(MORTSTAT == 1L), .groups = "drop")

continuous_variables <- c(
  age = "RIDAGEYR", nlr = "NLR", sii = "SII", haemoglobin = "LBXHGB",
  albumin = "LBXSAL", bmi = "BMXBMI", creatinine = "LBXSCR",
  income_to_poverty_ratio = "INDFMPIR", comorbidity_score = "Comorbidity_Score_Extended"
)
weighted_continuous <- dplyr::bind_rows(lapply(names(continuous_variables), function(label) {
  variable <- continuous_variables[[label]]
  message("  weighted continuous: ", label)
  estimate <- survey::svyby(
    stats::as.formula(paste0("~", variable)), ~phenotype_albumin,
    survey_design, survey::svymean, na.rm = TRUE, vartype = "se"
  )
  current_mean <- estimate[[variable]]
  current_se <- estimate[[if ("se" %in% names(estimate)) "se" else paste0("se.", variable)]]
  tibble::tibble(
    phenotype = as.character(estimate$phenotype_albumin),
    variable = label,
    weighted_mean = current_mean,
    standard_error = current_se
  )
}))

binary_expressions <- c(
  male = "male == 1",
  ever_smoking = "smoking != 'Never'",
  hypertension = "hypertension == 'Yes'",
  diabetes = "diabetes == 'Yes'"
)
weighted_binary <- dplyr::bind_rows(lapply(
  names(binary_expressions),
  function(label) {
    expression <- binary_expressions[[label]]
    message("  weighted binary: ", label)
    estimate <- survey::svyby(
      stats::as.formula(paste0("~I(", expression, ")")), ~phenotype_albumin,
      survey_design, survey::svymean, na.rm = TRUE, vartype = "se"
    )
    tibble::tibble(
      phenotype = as.character(estimate$phenotype_albumin),
      variable = label,
      weighted_proportion = estimate[[grep("TRUE$", names(estimate), value = TRUE)[1]]],
      standard_error = estimate[[grep("^se\\..*TRUE$", names(estimate), value = TRUE)[1]]]
    )
  }
), .id = "measure") |>
  dplyr::select(-measure)

readr::write_csv(phenotype_counts, file.path(output_root, "Table81A_full_phenotype_counts.csv"))
readr::write_csv(analytic_counts, file.path(output_root, "Table81B_primary_model_counts.csv"))
readr::write_csv(albumin_results$albumin_profiles, file.path(output_root, "Table81C_albumin_cluster_profiles.csv"))
readr::write_csv(weighted_continuous, file.path(output_root, "Table81D_weighted_continuous_profiles.csv"))
readr::write_csv(weighted_binary, file.path(output_root, "Table81E_weighted_binary_profiles.csv"))
readr::write_csv(primary_effects, file.path(output_root, "Table81F_primary_survey_Cox.csv"))

message("[3/9] Conventional-index benchmark models...")
predictors <- c(
  albumin_phenotype = "phenotype_albumin", domain_balanced_score = "domain_balanced_score",
  NLR = "risk_log_nlr", SII = "risk_log_sii", PNI = "risk_low_pni",
  GNRI = "risk_low_gnri", HALP = "risk_low_halp"
)
base_terms <- paste(
  "RIDAGEYR + male + race + INDFMPIR + Comorbidity_Score_Extended +",
  "cycle + smoking + hypertension + diabetes"
)
benchmark_effects <- dplyr::bind_rows(lapply(names(predictors), function(label) {
  term <- predictors[[label]]
  fit <- survey::svycoxph(
    stats::as.formula(paste("survival::Surv(PERMTH_INT, MORTSTAT) ~", term, "+", base_terms)),
    design = survey_design
  )
  beta <- stats::coef(fit)
  se <- sqrt(diag(stats::vcov(fit)))
  keep <- if (term == "phenotype_albumin") grepl("^phenotype_albumin", names(beta)) else names(beta) == term
  format_effect(beta[keep], se[keep]) |>
    dplyr::mutate(predictor = label, term = names(beta)[keep], .before = 1)
}))

incremental_effects <- dplyr::bind_rows(lapply(
  c("risk_log_nlr", "risk_log_sii", "risk_low_pni", "risk_low_gnri", "risk_low_halp"),
  function(comparator) {
    fit <- survey::svycoxph(
      stats::as.formula(paste(
        "survival::Surv(PERMTH_INT, MORTSTAT) ~", comparator,
        "+ phenotype_albumin +", base_terms
      )),
      design = survey_design
    )
    beta <- stats::coef(fit)["phenotype_albuminP1"]
    se <- sqrt(stats::vcov(fit)["phenotype_albuminP1", "phenotype_albuminP1"])
    joint <- survey::regTermTest(fit, ~phenotype_albumin)
    format_effect(beta, se) |>
      dplyr::mutate(
        comparator = comparator,
        global_phenotype_wald_p = as.numeric(joint$p),
        .before = 1
      )
  }
))
message("  benchmark column classes: ", paste(vapply(benchmark_effects, function(x) paste(class(x), collapse = "/"), character(1)), collapse = ", "))
message("  incremental column classes: ", paste(vapply(incremental_effects, function(x) paste(class(x), collapse = "/"), character(1)), collapse = ", "))
readr::write_csv(as.data.frame(benchmark_effects), file.path(output_root, "Table81G_conventional_index_benchmarks.csv"))
readr::write_csv(as.data.frame(incremental_effects), file.path(output_root, "Table81H_P1_beyond_conventional_indices.csv"))

message("[4/9] Albumin-completeness selection IPW...")
source(file.path(script_dir, "nhanes_modules", "33_nhanes_albumin_selection_ipw.R"), local = FALSE)

message("[5/9] Adjusted absolute risk, RMST, and 1,000-replicate survey bootstrap...")
Sys.setenv(ABSOLUTE_RISK_BOOTSTRAP_REPS = "1000")
source(file.path(script_dir, "nhanes_modules", "38_nhanes_albumin_absolute_risk_rmst.R"), local = FALSE)

message("[6/9] PH checks and 1,000-replicate cluster-aware bootstrap...")
source(file.path(script_dir, "nhanes_modules", "32_nhanes_albumin_cluster_bootstrap_ph.R"), local = FALSE)
source(file.path(script_dir, "nhanes_modules", "32b_nhanes_time_varying_score.R"), local = FALSE)

message("[7/9] K validation and repeated subsampling...")
source(file.path(script_dir, "nhanes_modules", "28b_nhanes_albumin_cluster_number.R"), local = FALSE)

message("[8/9] Leave-one-cycle-out projection and early-death exclusions...")
source(file.path(script_dir, "nhanes_modules", "35_nhanes_albumin_leave_one_cycle_out.R"), local = FALSE)

early_death_rows <- dplyr::bind_rows(lapply(c(12, 24), function(landmark) {
  current <- model_data |>
    dplyr::filter(PERMTH_INT > landmark) |>
    dplyr::mutate(time_after_landmark = PERMTH_INT - landmark)
  design <- survey::svydesign(
    ids = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~WTMEC8YR,
    nest = TRUE, data = current
  )
  fit <- survey::svycoxph(
    survival::Surv(time_after_landmark, MORTSTAT) ~
      phenotype_albumin + RIDAGEYR + male + race + INDFMPIR +
      Comorbidity_Score_Extended + cycle + smoking + hypertension + diabetes,
    design = design
  )
  beta <- stats::coef(fit)
  se <- sqrt(diag(stats::vcov(fit)))
  keep <- names(beta) %in% c("phenotype_albuminP2", "phenotype_albuminP1")
  format_effect(beta[keep], se[keep]) |>
    dplyr::mutate(
      landmark_months = landmark,
      comparison = dplyr::recode(
        names(beta)[keep], phenotype_albuminP2 = "P2 vs P3", phenotype_albuminP1 = "P1 vs P3"
      ),
      n = nrow(current), events = sum(current$MORTSTAT == 1L),
      early_deaths_excluded = sum(model_data$MORTSTAT == 1L & model_data$PERMTH_INT <= landmark),
      short_followup_censored_excluded = sum(model_data$MORTSTAT == 0L & model_data$PERMTH_INT <= landmark),
      .before = 1
    )
}))
readr::write_csv(early_death_rows, file.path(output_root, "Table81I_early_death_landmark_sensitivity.csv"))

# The unweighted diagnostic identifies a borderline phenotype-level PH signal.
# Characterise it directly rather than relying only on the continuous-score model.
time_varying_data <- model_data |>
  dplyr::mutate(
    p1_indicator = as.integer(phenotype_albumin == "P1"),
    p2_indicator = as.integer(phenotype_albumin == "P2"),
    survey_cluster = interaction(SDMVSTRA, SDMVPSU, drop = TRUE),
    analysis_weight = WTMEC8YR / mean(WTMEC8YR)
  )
time_varying_fit <- survival::coxph(
  survival::Surv(PERMTH_INT, MORTSTAT) ~
    p1_indicator + p2_indicator + tt(p1_indicator) + tt(p2_indicator) +
    RIDAGEYR + male + race + INDFMPIR + Comorbidity_Score_Extended +
    cycle + smoking + hypertension + diabetes,
  data = time_varying_data,
  weights = analysis_weight,
  cluster = survey_cluster,
  robust = TRUE,
  ties = "efron",
  tt = function(x, t, ...) x * log(pmax(t, 1))
)
tv_beta <- stats::coef(time_varying_fit)
tv_vcov <- stats::vcov(time_varying_fit)
time_varying_phenotype <- dplyr::bind_rows(lapply(c("P2 vs P3", "P1 vs P3"), function(comparison) {
  main_term <- if (comparison == "P1 vs P3") "p1_indicator" else "p2_indicator"
  time_term <- paste0("tt(", main_term, ")")
  dplyr::bind_rows(lapply(c(12, 36, 60, 96), function(month) {
    log_month <- log(month)
    log_hr <- tv_beta[main_term] + tv_beta[time_term] * log_month
    variance <- tv_vcov[main_term, main_term] + log_month^2 * tv_vcov[time_term, time_term] +
      2 * log_month * tv_vcov[main_term, time_term]
    se <- sqrt(variance)
    format_effect(log_hr, se) |>
      dplyr::mutate(comparison = comparison, followup_month = month, .before = 1)
  }))
}))
time_varying_terms <- tibble::tibble(
  comparison = c("P2 vs P3", "P1 vs P3"),
  log_time_interaction = c(tv_beta["tt(p2_indicator)"], tv_beta["tt(p1_indicator)"]),
  standard_error = sqrt(c(
    tv_vcov["tt(p2_indicator)", "tt(p2_indicator)"],
    tv_vcov["tt(p1_indicator)", "tt(p1_indicator)"]
  ))
) |>
  dplyr::mutate(
    p_value = 2 * stats::pnorm(abs(log_time_interaction / standard_error), lower.tail = FALSE)
  )
readr::write_csv(time_varying_phenotype, file.path(output_root, "Table81J_time_varying_phenotype_effects.csv"))
readr::write_csv(time_varying_terms, file.path(output_root, "Table81K_time_varying_phenotype_tests.csv"))

message("[9/9] Continuous-score nonlinearity checks and output manifest...")
source(file.path(script_dir, "nhanes_modules", "48_nhanes_albumin_nonlinearity.R"), local = FALSE)

saveRDS(
  list(
    primary_fit = primary_fit,
    primary_effects = primary_effects,
    phenotype_counts = phenotype_counts,
    analytic_counts = analytic_counts,
    weighted_continuous = weighted_continuous,
    weighted_binary = weighted_binary,
    benchmark_effects = benchmark_effects,
    incremental_effects = incremental_effects,
    early_death = early_death_rows,
    time_varying_phenotype = time_varying_phenotype,
    time_varying_terms = time_varying_terms
  ),
  file.path(output_root, "NHANES_Amendment17_core_results.rds")
)

summary_lines <- c(
  "NHANES Amendment 17 unified-albumin analysis suite",
  paste0("Full phenotype cohort: n=", nrow(complete_cohort), "; deaths=", sum(complete_cohort$MORTSTAT == 1L)),
  paste0("Primary model: n=", nrow(model_data), "; deaths=", sum(model_data$MORTSTAT == 1L)),
  "",
  "Primary survey Cox:",
  paste(capture.output(print(primary_effects)), collapse = "\n"),
  "",
  "Early-death landmark sensitivity:",
  paste(capture.output(print(early_death_rows)), collapse = "\n"),
  "",
  "Time-varying phenotype sensitivity:",
  paste(capture.output(print(time_varying_phenotype)), collapse = "\n")
)
writeLines(summary_lines, file.path(output_root, "NHANES_AMENDMENT17_SUMMARY.txt"))
cat(paste(summary_lines, collapse = "\n"), "\n")

source_files <- c(eligible_path, file.path(script_dir, "01_nhanes_albumin_analysis.R"), list.files(file.path(script_dir, "nhanes_modules"), pattern = "\\.R$", full.names = TRUE), list.files(covariate_cache, pattern = "^(BPQ|DIQ|SMQ)_[G-J]\\.xpt$", full.names = TRUE))
source_manifest <- tibble::tibble(
  source = vapply(source_files, function(path) {
    current <- normalizePath(path, winslash = "/", mustWork = TRUE)
    root_prefix <- paste0(normalizePath(root, winslash = "/", mustWork = TRUE), "/")
    if (startsWith(current, root_prefix)) substring(current, nchar(root_prefix) + 1L) else basename(current)
  }, character(1)),
  bytes = file.info(source_files)$size,
  sha256 = vapply(
    source_files,
    function(path) digest::digest(file = path, algo = "sha256", serialize = FALSE),
    character(1)
  )
)
readr::write_csv(source_manifest, file.path(output_root, "AMENDMENT17_NHANES_SOURCE_MANIFEST_SHA256.csv"))

output_manifest_path <- file.path(output_root, "AMENDMENT17_NHANES_OUTPUT_MANIFEST_SHA256.csv")
output_files <- list.files(output_root, recursive = TRUE, full.names = TRUE)
output_files <- output_files[
  !grepl("INTERNAL_NHANES_albumin_complete_cohort\\.rds$", output_files) &
    normalizePath(output_files, winslash = "/", mustWork = FALSE) !=
      normalizePath(output_manifest_path, winslash = "/", mustWork = FALSE)
]
manifest <- tibble::tibble(
  relative_path = substring(output_files, nchar(output_root) + 2L),
  bytes = file.info(output_files)$size,
  sha256 = vapply(
    output_files,
    function(path) digest::digest(file = path, algo = "sha256", serialize = FALSE),
    character(1)
  )
) |>
  dplyr::arrange(relative_path)
readr::write_csv(manifest, output_manifest_path)
message("Final non-row-level output files hashed: ", nrow(manifest))
