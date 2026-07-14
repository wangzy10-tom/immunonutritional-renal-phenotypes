# ==============================================================================
# Limited age and sex effect modification across NHANES, MIMIC-IV, and eICU
# ==============================================================================

required_packages <- c("dplyr", "readr", "tibble", "survival", "survey", "lme4")
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
output_dir <- file.path(root, "output", "age_sex_effect_modification")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

nhanes_path <- file.path(
  root, "output", "nhanes_albumin_benchmarks",
  "NHANES_albumin_benchmark_results.rds"
)
mimic_path <- file.path(
  root, "output", "mimic_official_oasis_v301",
  "MIMIC_official_OASIS_v301_results.rds"
)
eicu_path <- file.path(
  root, "output", "eicu_24h_robust",
  "eICU_24h_robust_results.rds"
)

coefficient_parts <- function(fit, modifier) {
  beta <- if (inherits(fit, "merMod")) lme4::fixef(fit) else stats::coef(fit)
  covariance <- as.matrix(stats::vcov(fit))
  terms <- names(beta)
  main_index <- which(
    grepl("^phenotype.*P1$", terms) & !grepl(":", terms, fixed = TRUE)
  )
  interaction_index <- which(
    grepl("phenotype", terms, fixed = TRUE) &
      grepl("P1", terms, fixed = TRUE) &
      grepl(":", terms, fixed = TRUE) &
      grepl(modifier, terms, fixed = TRUE)
  )
  global_index <- which(
    grepl("phenotype", terms, fixed = TRUE) &
      grepl(":", terms, fixed = TRUE) &
      grepl(modifier, terms, fixed = TRUE)
  )
  if (length(main_index) != 1L || length(interaction_index) != 1L || length(global_index) != 2L) {
    stop(
      "Could not identify the expected phenotype interaction coefficients for ", modifier,
      ". Terms were: ", paste(terms, collapse = ", "), call. = FALSE
    )
  }
  list(
    beta = beta,
    covariance = covariance,
    terms = terms,
    main_index = main_index,
    interaction_index = interaction_index,
    global_index = global_index
  )
}

global_interaction_test <- function(fit, dataset, modifier, model, effect_measure, n, events) {
  parts <- coefficient_parts(fit, modifier)
  interaction_beta <- parts$beta[parts$global_index]
  interaction_covariance <- parts$covariance[parts$global_index, parts$global_index, drop = FALSE]
  statistic <- as.numeric(t(interaction_beta) %*% qr.solve(interaction_covariance, interaction_beta))
  if (inherits(fit, "svycoxph")) {
    test_formula <- stats::as.formula(paste0("~ phenotype:", modifier))
    survey_test <- survey::regTermTest(fit, test_formula, method = "Wald")
    p_value <- as.numeric(survey_test$p)
    test_method <- "Design-based Wald F test"
  } else {
    p_value <- stats::pchisq(statistic, df = length(parts$global_index), lower.tail = FALSE)
    test_method <- "Wald chi-square test"
  }
  tibble::tibble(
    dataset = dataset,
    modifier = modifier,
    model = model,
    effect_measure = effect_measure,
    n = n,
    events = events,
    test_method = test_method,
    degrees_of_freedom = length(parts$global_index),
    wald_chisq = statistic,
    p_value = p_value
  )
}

p1_interaction_ratio <- function(fit, dataset, modifier, contrast, effect_measure, n, events) {
  parts <- coefficient_parts(fit, modifier)
  log_estimate <- unname(parts$beta[parts$interaction_index])
  standard_error <- sqrt(parts$covariance[parts$interaction_index, parts$interaction_index])
  tibble::tibble(
    dataset = dataset,
    modifier = modifier,
    contrast = contrast,
    effect_measure = paste0(effect_measure, " ratio"),
    n = n,
    events = events,
    estimate = exp(log_estimate),
    lower_95 = exp(log_estimate - 1.96 * standard_error),
    upper_95 = exp(log_estimate + 1.96 * standard_error),
    p_value = 2 * stats::pnorm(abs(log_estimate / standard_error), lower.tail = FALSE),
    effect_95ci = sprintf(
      "%.2f (%.2f-%.2f)", exp(log_estimate),
      exp(log_estimate - 1.96 * standard_error), exp(log_estimate + 1.96 * standard_error)
    )
  )
}

p1_conditional_effect <- function(
  fit, dataset, modifier, level, modifier_value, effect_measure, n, events
) {
  parts <- coefficient_parts(fit, modifier)
  contrast <- numeric(length(parts$beta))
  contrast[parts$main_index] <- 1
  contrast[parts$interaction_index] <- modifier_value
  log_estimate <- as.numeric(sum(contrast * parts$beta))
  variance <- as.numeric(t(contrast) %*% parts$covariance %*% contrast)
  standard_error <- sqrt(max(variance, 0))
  tibble::tibble(
    dataset = dataset,
    modifier = modifier,
    level = level,
    effect_measure = effect_measure,
    comparison = "P1 vs P3",
    n = n,
    events = events,
    estimate = exp(log_estimate),
    lower_95 = exp(log_estimate - 1.96 * standard_error),
    upper_95 = exp(log_estimate + 1.96 * standard_error),
    effect_95ci = sprintf(
      "%.2f (%.2f-%.2f)", exp(log_estimate),
      exp(log_estimate - 1.96 * standard_error), exp(log_estimate + 1.96 * standard_error)
    )
  )
}

cat("Fitting NHANES complex-survey interaction models...\n")
nhanes_results <- readRDS(nhanes_path)
nhanes <- nhanes_results$model_data |>
  dplyr::mutate(
    phenotype = factor(as.character(phenotype_total_protein), levels = c("P3", "P2", "P1")),
    age10 = (RIDAGEYR - 75) / 10,
    male = as.integer(male)
  )
if (nrow(nhanes) != 3979L || sum(nhanes$MORTSTAT == 1) != 720L) {
  stop("Unexpected NHANES interaction-model cohort dimensions.", call. = FALSE)
}
nhanes_design <- survey::svydesign(
  ids = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~WTMEC8YR,
  nest = TRUE, data = nhanes
)
nhanes_age_fit <- survey::svycoxph(
  survival::Surv(PERMTH_INT, MORTSTAT) ~ phenotype * age10 + male + race +
    INDFMPIR + Comorbidity_Score_Extended + cycle + smoking + hypertension + diabetes,
  design = nhanes_design
)
nhanes_sex_fit <- survey::svycoxph(
  survival::Surv(PERMTH_INT, MORTSTAT) ~ phenotype * male + age10 + race +
    INDFMPIR + Comorbidity_Score_Extended + cycle + smoking + hypertension + diabetes,
  design = nhanes_design
)

cat("Fitting MIMIC-IV official OASIS-stratified interaction models...\n")
mimic_results <- readRDS(mimic_path)
mimic <- mimic_results$analysis |>
  dplyr::mutate(
    phenotype = factor(as.character(phenotype), levels = c("P3", "P2", "P1")),
    age10 = (anchor_age - 75) / 10,
    male = as.integer(male),
    oasis_quartile = factor(oasis_quartile)
  )
if (nrow(mimic) != 1145L || sum(mimic$mortality_365d == 1) != 574L) {
  stop("Unexpected MIMIC interaction-model cohort dimensions.", call. = FALSE)
}
mimic_age_fit <- survival::coxph(
  survival::Surv(survival_days_365, mortality_365d) ~ phenotype * age10 + male +
    strata(oasis_quartile), data = mimic, ties = "efron"
)
mimic_sex_fit <- survival::coxph(
  survival::Surv(survival_days_365, mortality_365d) ~ phenotype * male + age10 +
    strata(oasis_quartile), data = mimic, ties = "efron"
)

cat("Fitting eICU APACHE-adjusted hospital mixed interaction models...\n")
eicu_results <- readRDS(eicu_path)
eicu <- eicu_results$model_data |>
  dplyr::mutate(
    hospitalid = factor(hospitalid),
    phenotype = factor(as.character(phenotype), levels = c("P3", "P2", "P1")),
    age10 = (age_num - 75) / 10,
    male = as.integer(male)
  )
if (nrow(eicu) != 12548L || sum(eicu$hospital_mortality == 1) != 1884L) {
  stop("Unexpected eICU interaction-model cohort dimensions.", call. = FALSE)
}
eicu_control <- lme4::glmerControl(
  optimizer = "bobyqa", optCtrl = list(maxfun = 200000)
)
eicu_age_fit <- lme4::glmer(
  hospital_mortality ~ phenotype * age10 + male + apachescore + (1 | hospitalid),
  family = stats::binomial(), data = eicu, nAGQ = 1, control = eicu_control
)
eicu_sex_fit <- lme4::glmer(
  hospital_mortality ~ phenotype * male + age10 + apachescore + (1 | hospitalid),
  family = stats::binomial(), data = eicu, nAGQ = 1, control = eicu_control
)

fit_specifications <- list(
  list("NHANES 2011-2018", "age10", nhanes_age_fit, "Complex-survey fully adjusted Cox", "HR", nrow(nhanes), sum(nhanes$MORTSTAT == 1)),
  list("NHANES 2011-2018", "male", nhanes_sex_fit, "Complex-survey fully adjusted Cox", "HR", nrow(nhanes), sum(nhanes$MORTSTAT == 1)),
  list("MIMIC-IV v3.1", "age10", mimic_age_fit, "Official OASIS quartile-stratified Cox", "HR", nrow(mimic), sum(mimic$mortality_365d == 1)),
  list("MIMIC-IV v3.1", "male", mimic_sex_fit, "Official OASIS quartile-stratified Cox", "HR", nrow(mimic), sum(mimic$mortality_365d == 1)),
  list("eICU v2.0", "age10", eicu_age_fit, "APACHE-adjusted hospital random-intercept logistic mixed model", "OR", nrow(eicu), sum(eicu$hospital_mortality == 1)),
  list("eICU v2.0", "male", eicu_sex_fit, "APACHE-adjusted hospital random-intercept logistic mixed model", "OR", nrow(eicu), sum(eicu$hospital_mortality == 1))
)

global_tests <- dplyr::bind_rows(lapply(fit_specifications, function(specification) {
  global_interaction_test(
    fit = specification[[3]], dataset = specification[[1]], modifier = specification[[2]],
    model = specification[[4]], effect_measure = specification[[5]],
    n = specification[[6]], events = specification[[7]]
  )
})) |>
  dplyr::mutate(
    modifier = dplyr::recode(modifier, age10 = "Age", male = "Sex"),
    p_value_bh = stats::p.adjust(p_value, method = "BH"),
    significant_bh_0_05 = p_value_bh < 0.05
  )

interaction_ratios <- dplyr::bind_rows(lapply(fit_specifications, function(specification) {
  modifier <- specification[[2]]
  p1_interaction_ratio(
    fit = specification[[3]], dataset = specification[[1]], modifier = modifier,
    contrast = if (modifier == "age10") "Per 10-year increase" else "Male vs female",
    effect_measure = specification[[5]], n = specification[[6]], events = specification[[7]]
  )
})) |>
  dplyr::mutate(modifier = dplyr::recode(modifier, age10 = "Age", male = "Sex"))

conditional_effects <- dplyr::bind_rows(lapply(fit_specifications, function(specification) {
  modifier <- specification[[2]]
  if (modifier == "age10") {
    levels <- tibble::tibble(level = c("Age 70 years", "Age 80 years"), value = c(-0.5, 0.5))
  } else {
    levels <- tibble::tibble(level = c("Female", "Male"), value = c(0, 1))
  }
  dplyr::bind_rows(lapply(seq_len(nrow(levels)), function(index) {
    p1_conditional_effect(
      fit = specification[[3]], dataset = specification[[1]], modifier = modifier,
      level = levels$level[index], modifier_value = levels$value[index],
      effect_measure = specification[[5]], n = specification[[6]], events = specification[[7]]
    )
  }))
})) |>
  dplyr::mutate(modifier = dplyr::recode(modifier, age10 = "Age", male = "Sex"))

model_diagnostics <- tibble::tibble(
  dataset = c("eICU v2.0", "eICU v2.0"),
  modifier = c("Age", "Sex"),
  converged = c(
    is.null(eicu_age_fit@optinfo$conv$lme4$messages),
    is.null(eicu_sex_fit@optinfo$conv$lme4$messages)
  ),
  singular = c(
    lme4::isSingular(eicu_age_fit, tol = 1e-5),
    lme4::isSingular(eicu_sex_fit, tol = 1e-5)
  )
)

readr::write_csv(global_tests, file.path(output_dir, "Table44A_global_interaction_tests.csv"))
readr::write_csv(interaction_ratios, file.path(output_dir, "Table44B_P1_interaction_ratios.csv"))
readr::write_csv(conditional_effects, file.path(output_dir, "Table44C_P1_prespecified_conditional_effects.csv"))
readr::write_csv(model_diagnostics, file.path(output_dir, "Table44D_model_diagnostics.csv"))

saveRDS(
  list(
    global_tests = global_tests,
    interaction_ratios = interaction_ratios,
    conditional_effects = conditional_effects,
    model_diagnostics = model_diagnostics,
    models = list(
      nhanes_age = nhanes_age_fit, nhanes_sex = nhanes_sex_fit,
      mimic_age = mimic_age_fit, mimic_sex = mimic_sex_fit,
      eicu_age = eicu_age_fit, eicu_sex = eicu_sex_fit
    )
  ),
  file.path(output_dir, "age_sex_effect_modification_results.rds")
)

summary_lines <- c(
  "Age and sex effect-modification analysis",
  "",
  paste0("Global interaction tests: ", nrow(global_tests)),
  paste0("BH-significant interactions: ", sum(global_tests$significant_bh_0_05)),
  "",
  paste(capture.output(print(global_tests)), collapse = "\n"),
  "",
  "Pre-specified P1 versus P3 conditional effects:",
  paste(capture.output(print(conditional_effects)), collapse = "\n"),
  "",
  "Interpretation: interaction inference uses the two-degree-of-freedom global test and BH correction; conditional estimates are descriptive."
)
writeLines(summary_lines, file.path(output_dir, "SUMMARY.txt"))
cat(paste(summary_lines, collapse = "\n"), "\n")
