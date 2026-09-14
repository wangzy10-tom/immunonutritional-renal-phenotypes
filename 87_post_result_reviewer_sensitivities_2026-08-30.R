# ==============================================================================
# Post-result, reviewer-driven exploratory sensitivity analyses (2026-08-30)
#
# These analyses were requested after the primary and earlier sensitivity
# results were known. They must therefore be reported as outcome-aware,
# exploratory analyses rather than prospectively specified analyses.
#
# Modules:
#   1. NHANES primary model additionally adjusted for the exact standardized
#      log-creatinine representation used in clustering.
#   2. K=3 clustering after omitting BMI, within the same six-feature-complete
#      cohort in each database (sample membership held fixed).
#   3. One-, three-, and seven-day landmark analyses in MIMIC-IV and eICU.
#      All patients discharged or otherwise no longer under observation before
#      each eICU landmark are excluded together with early deaths.
#
# No row-level outputs are written.
# ==============================================================================

required_packages <- c(
  'dplyr', 'readr', 'tibble', 'survival', 'survey', 'data.table',
  'lme4', 'mclust', 'digest'
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop('Missing required package(s): ', paste(missing_packages, collapse = ', '), call. = FALSE)
}
options(survey.lonely.psu = 'adjust')

root <- normalizePath(Sys.getenv('PROJECT_ROOT', unset = getwd()), winslash = '/', mustWork = TRUE)
# The 30 August directory contains the preliminary pre-Amendment-22 execution.
# Formal executions after Amendment 22 use the separate default directory below.
output_dir <- file.path(root, 'output', 'amendment22_reviewer_sensitivities_2026-08-31')
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

paths <- list(
  nhanes_cluster = file.path(
    root, 'output', 'nhanes_albumin_amendment17_2026-08-30',
    'INTERNAL_NHANES_albumin_complete_cohort.rds'
  ),
  nhanes_model = file.path(
    root, 'output', 'nhanes_albumin_amendment17_2026-08-30',
    'NHANES_albumin_benchmark_results.rds'
  ),
  mimic = file.path(
    root, 'output', 'mimic_albumin_primary_2026-08-30',
    'MIMIC_albumin_primary_results.rds'
  ),
  eicu_denominator = file.path(
    root, 'output', 'eicu_24h_extraction',
    'eICU_first_ICU_feature_availability_dataset.csv'
  ),
  eicu_apache_raw = file.path(Sys.getenv('EICU_DIR', unset=''), 'apachePatientResult.csv.gz'),
  eicu_apache_audited = file.path(
    root, 'output', 'eicu_validation', 'eICU_APACHE_adjusted_sensitivity.rds'
  )
)
paths$eicu_denominator <- Sys.getenv("EICU_DENOMINATOR_CSV", unset=file.path(root,"data","derived","eICU_first_ICU_feature_availability_dataset.csv"))
paths$eicu_apache_audited <- Sys.getenv("EICU_APACHE_AUDIT_RDS", unset=file.path(root,"data","derived","eICU_APACHE_adjusted_sensitivity.rds"))

required_inputs <- paths[c('nhanes_cluster', 'nhanes_model', 'mimic', 'eicu_denominator')]
missing_inputs <- names(required_inputs)[!vapply(required_inputs, file.exists, logical(1))]
if (length(missing_inputs) > 0L) {
  stop('Missing input(s): ', paste(missing_inputs, collapse = ', '), call. = FALSE)
}
if (!file.exists(paths$eicu_apache_raw) && !file.exists(paths$eicu_apache_audited)) {
  stop('No eICU APACHE source is available.', call. = FALSE)
}

winsorise <- function(x) {
  limits <- stats::quantile(x, c(0.01, 0.99), na.rm = TRUE, names = FALSE)
  pmin(pmax(as.numeric(x), limits[1]), limits[2])
}

safe_z <- function(x) {
  current_sd <- stats::sd(x, na.rm = TRUE)
  if (!is.finite(current_sd) || current_sd == 0) return(rep(0, length(x)))
  (x - mean(x, na.rm = TRUE)) / current_sd
}

sha256_file <- function(path) {
  digest::digest(file = path, algo = 'sha256', serialize = FALSE)
}

prepare_no_bmi_matrix <- function(data) {
  output <- data |>
    dplyr::transmute(
      log_nlr = log(pmax(winsorise(nlr), .Machine$double.eps)),
      log_sii = log(pmax(winsorise(sii), .Machine$double.eps)),
      haemoglobin = winsorise(haemoglobin),
      albumin = winsorise(albumin),
      log_creatinine = log(pmax(winsorise(creatinine), .Machine$double.eps))
    ) |>
    as.data.frame() |>
    scale()
  if (any(!is.finite(output))) stop('Non-finite value in no-BMI clustering matrix.', call. = FALSE)
  output
}

label_no_bmi <- function(data, raw_cluster) {
  profile <- data |>
    dplyr::mutate(cluster_raw = raw_cluster) |>
    dplyr::group_by(cluster_raw) |>
    dplyr::summarise(
      n = dplyr::n(),
      nlr = stats::median(nlr),
      sii = stats::median(sii),
      haemoglobin = stats::median(haemoglobin),
      albumin = stats::median(albumin),
      bmi = stats::median(bmi),
      creatinine = stats::median(creatinine),
      .groups = 'drop'
    ) |>
    dplyr::arrange(cluster_raw)
  p1_score <- safe_z(log(profile$nlr)) + safe_z(log(profile$sii)) +
    safe_z(log(profile$creatinine))
  p1_raw <- profile$cluster_raw[which.max(p1_score)]
  remaining <- profile |> dplyr::filter(cluster_raw != p1_raw)
  p2_score <- -safe_z(remaining$haemoglobin) - safe_z(remaining$albumin)
  p2_raw <- remaining$cluster_raw[which.max(p2_score)]
  p3_raw <- setdiff(profile$cluster_raw, c(p1_raw, p2_raw))
  if (length(p3_raw) != 1L) stop('No-BMI mapping failed to yield one P3 cluster.', call. = FALSE)
  mapping <- tibble::tibble(
    cluster_raw = c(p1_raw, p2_raw, p3_raw),
    phenotype_no_bmi = c('P1', 'P2', 'P3')
  )
  labels <- mapping$phenotype_no_bmi[match(raw_cluster, mapping$cluster_raw)]
  list(
    phenotype = factor(labels, levels = c('P3', 'P2', 'P1')),
    profile = profile |>
      dplyr::mutate(p1_score = p1_score) |>
      dplyr::left_join(mapping, by = 'cluster_raw') |>
      dplyr::arrange(factor(phenotype_no_bmi, levels = c('P1', 'P2', 'P3'))),
    mapping = mapping
  )
}

run_no_bmi <- function(data) {
  x <- prepare_no_bmi_matrix(data)
  set.seed(20260710)
  fit <- stats::kmeans(x, centers = 3, nstart = 100, iter.max = 500, algorithm = 'Lloyd')
  labelled <- label_no_bmi(data, fit$cluster)
  list(
    phenotype = labelled$phenotype,
    profile = labelled$profile,
    tot_withinss = as.numeric(fit$tot.withinss),
    iterations = fit$iter
  )
}

assignment_agreement <- function(dataset, primary, alternative) {
  primary <- factor(as.character(primary), levels = c('P3', 'P2', 'P1'))
  alternative <- factor(as.character(alternative), levels = c('P3', 'P2', 'P1'))
  primary_p1 <- primary == 'P1'
  alternative_p1 <- alternative == 'P1'
  intersection <- sum(primary_p1 & alternative_p1)
  union <- sum(primary_p1 | alternative_p1)
  tibble::tibble(
    dataset = dataset,
    n_clustered = length(primary),
    adjusted_rand_index = mclust::adjustedRandIndex(primary, alternative),
    exact_label_agreement = mean(primary == alternative),
    primary_p1_n = sum(primary_p1),
    no_bmi_p1_n = sum(alternative_p1),
    p1_recall = intersection / sum(primary_p1),
    p1_precision = intersection / sum(alternative_p1),
    p1_jaccard = intersection / union
  )
}

extract_effects <- function(fit, term_prefix, dataset, model, measure, n, events) {
  beta <- if (inherits(fit, 'merMod')) lme4::fixef(fit) else stats::coef(fit)
  covariance <- as.matrix(stats::vcov(fit))
  terms <- paste0(term_prefix, c('P1', 'P2'))
  if (!all(terms %in% names(beta))) {
    stop('Missing phenotype coefficient(s) in ', dataset, ' model: ', paste(setdiff(terms, names(beta)), collapse = ', '), call. = FALSE)
  }
  se <- sqrt(diag(covariance))[terms]
  estimate <- exp(beta[terms])
  lower <- exp(beta[terms] - 1.96 * se)
  upper <- exp(beta[terms] + 1.96 * se)
  tibble::tibble(
    dataset = dataset, model = model, effect_measure = measure,
    comparison = c('P1 vs P3', 'P2 vs P3'), n = n, events = events,
    estimate = as.numeric(estimate), lower_95 = as.numeric(lower), upper_95 = as.numeric(upper),
    p_value = 2 * stats::pnorm(abs(beta[terms] / se), lower.tail = FALSE),
    effect_95ci = sprintf('%.3f (%.3f-%.3f)', estimate, lower, upper)
  )
}

capture_glmer <- function(formula, data) {
  warnings <- character()
  fit <- tryCatch(
    withCallingHandlers(
      lme4::glmer(
        formula, data = data, family = stats::binomial(), nAGQ = 1,
        control = lme4::glmerControl(
          optimizer = 'bobyqa', optCtrl = list(maxfun = 300000), calc.derivs = TRUE
        )
      ),
      warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart('muffleWarning')
      }
    ),
    error = function(e) e
  )
  if (inherits(fit, 'error')) stop('glmer failed: ', conditionMessage(fit), call. = FALSE)
  list(fit = fit, warnings = unique(warnings))
}

glmer_diagnostics <- function(captured, dataset, model, data, outcome) {
  fit <- captured$fit
  messages <- fit@optinfo$conv$lme4$messages
  if (is.null(messages)) messages <- character()
  tibble::tibble(
    dataset = dataset, model = model, n = nrow(data), events = sum(data[[outcome]] == 1),
    hospitals = dplyr::n_distinct(data$hospitalid),
    convergence_message = paste(messages, collapse = ' | '),
    captured_warnings = paste(captured$warnings, collapse = ' | '),
    singular = lme4::isSingular(fit, tol = 1e-5)
  )
}

cat('[1/5] Loading frozen cohorts and reconstructing current eICU albumin labels...\n')
nhanes_full <- readRDS(paths$nhanes_cluster)
nhanes_model <- readRDS(paths$nhanes_model)$model_data
mimic <- readRDS(paths$mimic)$analysis
eicu_denominator <- data.table::fread(
  paths$eicu_denominator, data.table = FALSE, showProgress = FALSE
) |>
  tibble::as_tibble()

if (file.exists(paths$eicu_apache_raw)) {
  apache <- readr::read_csv(
    paths$eicu_apache_raw,
    col_select = c(patientunitstayid, apacheversion, apachescore),
    show_col_types = FALSE, progress = FALSE
  ) |>
    dplyr::filter(apacheversion == 'IVa') |>
    dplyr::transmute(
      patientunitstayid,
      apachescore = ifelse(as.numeric(apachescore) < 0, NA_real_, as.numeric(apachescore))
    ) |>
    dplyr::distinct(patientunitstayid, .keep_all = TRUE)
  apache_source <- paths$eicu_apache_raw
} else {
  apache <- readRDS(paths$eicu_apache_audited)$model_data |>
    dplyr::select(patientunitstayid, apachescore) |>
    dplyr::distinct(patientunitstayid, .keep_all = TRUE)
  apache_source <- paths$eicu_apache_audited
}

eicu_feature_complete <- eicu_denominator |>
  dplyr::filter(
    dplyr::if_all(
      dplyr::all_of(c('nlr', 'sii_like', 'haemoglobin', 'albumin', 'bmi', 'creatinine')),
      ~ !is.na(.x) & is.finite(.x)
    )
  ) |>
  dplyr::filter(
    nlr > 0, nlr <= 100, sii_like > 0, sii_like <= 200000,
    albumin > 0, creatinine > 0
  ) |>
  dplyr::mutate(
    male = as.integer(gender_model == 'Male'),
    icu_los_days = as.numeric(unitdischargeoffset) / 1440
  )
if (nrow(eicu_feature_complete) != 15242L) {
  stop('eICU albumin feature-complete cohort did not reproduce n=15,242.', call. = FALSE)
}

eicu_primary_input <- eicu_feature_complete |>
  dplyr::transmute(
    id = patientunitstayid, nlr = nlr, sii = sii_like,
    haemoglobin = haemoglobin, albumin = albumin, bmi = bmi, creatinine = creatinine
  )
primary_x <- eicu_primary_input |>
  dplyr::transmute(
    log_nlr = log(pmax(winsorise(nlr), .Machine$double.eps)),
    log_sii = log(pmax(winsorise(sii), .Machine$double.eps)),
    haemoglobin = winsorise(haemoglobin), albumin = winsorise(albumin),
    bmi = winsorise(bmi),
    log_creatinine = log(pmax(winsorise(creatinine), .Machine$double.eps))
  ) |>
  as.data.frame() |>
  scale()
set.seed(20260710)
primary_fit <- stats::kmeans(primary_x, centers = 3, nstart = 100, iter.max = 500, algorithm = 'Lloyd')
primary_profile <- eicu_primary_input |>
  dplyr::mutate(cluster_raw = primary_fit$cluster) |>
  dplyr::group_by(cluster_raw) |>
  dplyr::summarise(
    n = dplyr::n(), nlr = stats::median(nlr), sii = stats::median(sii),
    haemoglobin = stats::median(haemoglobin), albumin = stats::median(albumin),
    bmi = stats::median(bmi), creatinine = stats::median(creatinine), .groups = 'drop'
  )
primary_p1_score <- safe_z(log(primary_profile$nlr)) + safe_z(log(primary_profile$sii)) +
  safe_z(log(primary_profile$creatinine))
primary_p1_raw <- primary_profile$cluster_raw[which.max(primary_p1_score)]
primary_remaining <- primary_profile |> dplyr::filter(cluster_raw != primary_p1_raw)
primary_p2_score <- -safe_z(primary_remaining$haemoglobin) - safe_z(primary_remaining$albumin) -
  safe_z(primary_remaining$bmi)
primary_p2_raw <- primary_remaining$cluster_raw[which.max(primary_p2_score)]
primary_p3_raw <- setdiff(primary_profile$cluster_raw, c(primary_p1_raw, primary_p2_raw))
eicu_feature_complete$phenotype_primary <- factor(
  dplyr::case_when(
    primary_fit$cluster == primary_p1_raw ~ 'P1',
    primary_fit$cluster == primary_p2_raw ~ 'P2',
    primary_fit$cluster == primary_p3_raw ~ 'P3',
    TRUE ~ NA_character_
  ),
  levels = c('P3', 'P2', 'P1')
)
if (!identical(
    as.integer(table(factor(eicu_feature_complete$phenotype_primary, levels = c('P1', 'P2', 'P3')))),
    c(5737L, 4337L, 5168L)
  )) {
  stop('Current eICU primary labels did not reproduce locked counts.', call. = FALSE)
}

cat('[2/5] NHANES standardized log-creatinine component adjustment...\n')
renal_assignment <- nhanes_full |>
  dplyr::transmute(
    SEQN,
    reviewer_log_creatinine_z = safe_z(log(pmax(winsorise(LBXSCR), .Machine$double.eps)))
  )
nhanes_renal <- nhanes_model |>
  dplyr::select(-dplyr::any_of('reviewer_log_creatinine_z')) |>
  dplyr::left_join(renal_assignment, by = 'SEQN')
if (nrow(nhanes_renal) != 3979L || anyNA(nhanes_renal$reviewer_log_creatinine_z)) {
  stop('NHANES renal-component model cohort mismatch.', call. = FALSE)
}
nhanes_design <- survey::svydesign(
  ids = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~WTMEC8YR,
  nest = TRUE, data = nhanes_renal
)
nhanes_base_formula <- survival::Surv(PERMTH_INT, MORTSTAT) ~
  phenotype_albumin + RIDAGEYR + male + race + INDFMPIR +
  Comorbidity_Score_Extended + cycle + smoking + hypertension + diabetes
nhanes_renal_formula <- update(nhanes_base_formula, . ~ . + reviewer_log_creatinine_z)
nhanes_primary_fit <- survey::svycoxph(nhanes_base_formula, design = nhanes_design)
nhanes_renal_fit <- survey::svycoxph(nhanes_renal_formula, design = nhanes_design)
nhanes_creatinine_models <- dplyr::bind_rows(
  extract_effects(
    nhanes_primary_fit, 'phenotype_albumin', 'NHANES',
    'Primary fully adjusted complex-survey Cox', 'HR', 3979L, 720L
  ),
  extract_effects(
    nhanes_renal_fit, 'phenotype_albumin', 'NHANES',
    'Primary model plus standardized log-creatinine component', 'HR', 3979L, 720L
  )
)

cat('[3/5] Same-cohort K=3 clustering after omitting BMI...\n')
nhanes_cluster_input <- nhanes_full |>
  dplyr::transmute(
    id = SEQN, nlr = NLR, sii = SII, haemoglobin = LBXHGB, albumin = LBXSAL,
    bmi = BMXBMI, creatinine = LBXSCR,
    primary_phenotype = factor(phenotype_albumin, levels = c('P3', 'P2', 'P1'))
  )
mimic_cluster_input <- mimic |>
  dplyr::transmute(
    id = stay_id, nlr = nlr, sii = sii, haemoglobin = haemoglobin, albumin = albumin,
    bmi = bmi, creatinine = creatinine,
    primary_phenotype = factor(phenotype, levels = c('P3', 'P2', 'P1'))
  )
eicu_cluster_input <- eicu_feature_complete |>
  dplyr::transmute(
    id = patientunitstayid, nlr = nlr, sii = sii_like, haemoglobin = haemoglobin,
    albumin = albumin, bmi = bmi, creatinine = creatinine,
    primary_phenotype = phenotype_primary
  )

cluster_inputs <- list(
  NHANES = nhanes_cluster_input,
  'MIMIC-IV' = mimic_cluster_input,
  eICU = eicu_cluster_input
)
no_bmi_results <- lapply(cluster_inputs, run_no_bmi)

no_bmi_agreement <- dplyr::bind_rows(lapply(names(cluster_inputs), function(dataset) {
  assignment_agreement(
    dataset,
    cluster_inputs[[dataset]]$primary_phenotype,
    no_bmi_results[[dataset]]$phenotype
  )
}))
no_bmi_profiles <- dplyr::bind_rows(lapply(names(cluster_inputs), function(dataset) {
  no_bmi_results[[dataset]]$profile |>
    dplyr::mutate(dataset = dataset, .before = 1)
}))
no_bmi_method <- dplyr::bind_rows(lapply(names(cluster_inputs), function(dataset) {
  tibble::tibble(
    dataset = dataset, n_clustered = nrow(cluster_inputs[[dataset]]),
    features_used = 'NLR, SII/SII-like, haemoglobin, albumin, creatinine',
    deliberately_omitted_feature = 'BMI',
    sample_membership = 'Same six-feature-complete cohort as the primary analysis',
    clustering = 'K-means Lloyd; K=3; nstart=100; iter.max=500',
    p1_mapping = 'Maximum standardized log(NLR) + log(SII/SII-like) + log(creatinine)',
    p2_mapping = 'Among remaining clusters, minimum haemoglobin + albumin reserve',
    seed = 20260710L,
    tot_withinss = no_bmi_results[[dataset]]$tot_withinss,
    iterations = no_bmi_results[[dataset]]$iterations,
    outcome_used_for_clustering_or_mapping = FALSE
  )
}))

nhanes_no_bmi_assignment <- tibble::tibble(
  SEQN = nhanes_cluster_input$id,
  phenotype_no_bmi = as.character(no_bmi_results$NHANES$phenotype)
)
nhanes_no_bmi_model <- nhanes_model |>
  dplyr::left_join(nhanes_no_bmi_assignment, by = 'SEQN') |>
  dplyr::mutate(phenotype_no_bmi = factor(phenotype_no_bmi, levels = c('P3', 'P2', 'P1')))
nhanes_no_bmi_design <- survey::svydesign(
  ids = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~WTMEC8YR,
  nest = TRUE, data = nhanes_no_bmi_model
)
nhanes_no_bmi_fit <- survey::svycoxph(
  survival::Surv(PERMTH_INT, MORTSTAT) ~
    phenotype_no_bmi + RIDAGEYR + male + race + INDFMPIR +
    Comorbidity_Score_Extended + cycle + smoking + hypertension + diabetes,
  design = nhanes_no_bmi_design
)

mimic_no_bmi_assignment <- tibble::tibble(
  stay_id = mimic_cluster_input$id,
  phenotype_no_bmi = as.character(no_bmi_results[['MIMIC-IV']]$phenotype)
)
mimic_no_bmi <- mimic |>
  dplyr::left_join(mimic_no_bmi_assignment, by = 'stay_id') |>
  dplyr::mutate(phenotype_no_bmi = factor(phenotype_no_bmi, levels = c('P3', 'P2', 'P1')))
mimic_no_bmi_fit <- survival::coxph(
  survival::Surv(survival_days_365, mortality_365d) ~
    phenotype_no_bmi + male + strata(oasis_quartile),
  data = mimic_no_bmi, ties = 'efron'
)

eicu_no_bmi_assignment <- tibble::tibble(
  patientunitstayid = eicu_cluster_input$id,
  phenotype_no_bmi = as.character(no_bmi_results$eICU$phenotype)
)
eicu_no_bmi <- eicu_feature_complete |>
  dplyr::left_join(eicu_no_bmi_assignment, by = 'patientunitstayid') |>
  dplyr::left_join(apache, by = 'patientunitstayid') |>
  dplyr::filter(is.finite(apachescore), !is.na(hospital_mortality)) |>
  dplyr::mutate(
    phenotype_no_bmi = factor(phenotype_no_bmi, levels = c('P3', 'P2', 'P1')),
    hospitalid = factor(hospitalid),
    age_z = as.numeric(scale(age_num)),
    apache_z = as.numeric(scale(apachescore))
  ) |>
  droplevels()
eicu_no_bmi_capture <- capture_glmer(
  hospital_mortality ~ phenotype_no_bmi + age_z + male + apache_z + (1 | hospitalid),
  eicu_no_bmi
)

no_bmi_outcome_models <- dplyr::bind_rows(
  extract_effects(
    nhanes_no_bmi_fit, 'phenotype_no_bmi', 'NHANES',
    'No-BMI K=3 fully adjusted complex-survey Cox', 'HR',
    nrow(nhanes_no_bmi_model), sum(nhanes_no_bmi_model$MORTSTAT == 1)
  ),
  extract_effects(
    mimic_no_bmi_fit, 'phenotype_no_bmi', 'MIMIC-IV',
    'No-BMI K=3 official OASIS quartile-stratified Cox', 'HR',
    nrow(mimic_no_bmi), sum(mimic_no_bmi$mortality_365d == 1)
  ),
  extract_effects(
    eicu_no_bmi_capture$fit, 'phenotype_no_bmi', 'eICU',
    'No-BMI K=3 APACHE-adjusted hospital random-intercept logistic', 'OR',
    nrow(eicu_no_bmi), sum(eicu_no_bmi$hospital_mortality == 1)
  )
)
no_bmi_diagnostics <- glmer_diagnostics(
  eicu_no_bmi_capture, 'eICU', 'No-BMI K=3 hospital mortality',
  eicu_no_bmi, 'hospital_mortality'
)

cat('[4/5] One-, three-, and seven-day ICU landmark analyses...\n')
landmark_effect_rows <- list()
landmark_count_rows <- list()
landmark_diagnostic_rows <- list()

for (landmark_day in c(1, 3, 7)) {
  mimic_landmark <- mimic |>
    dplyr::filter(survival_days_365 > landmark_day) |>
    dplyr::mutate(time_after_landmark = survival_days_365 - landmark_day)
  oasis_fit <- survival::coxph(
    survival::Surv(time_after_landmark, mortality_365d) ~
      phenotype + male + strata(oasis_quartile),
    data = mimic_landmark, ties = 'efron'
  )
  sofa_fit <- survival::coxph(
    survival::Surv(time_after_landmark, mortality_365d) ~
      phenotype + anchor_age + male + strata(sofa_quartile),
    data = mimic_landmark, ties = 'efron'
  )
  landmark_effect_rows[[paste0('mimic_oasis_', landmark_day)]] <- extract_effects(
    oasis_fit, 'phenotype', 'MIMIC-IV',
    paste0(landmark_day, '-day landmark, official OASIS quartile-stratified Cox'),
    'HR', nrow(mimic_landmark), sum(mimic_landmark$mortality_365d == 1)
  ) |>
    dplyr::mutate(landmark_day = landmark_day, endpoint = '365-day all-cause mortality', .before = 3)
  landmark_effect_rows[[paste0('mimic_sofa_', landmark_day)]] <- extract_effects(
    sofa_fit, 'phenotype', 'MIMIC-IV',
    paste0(landmark_day, '-day landmark, official SOFA quartile-stratified Cox'),
    'HR', nrow(mimic_landmark), sum(mimic_landmark$mortality_365d == 1)
  ) |>
    dplyr::mutate(landmark_day = landmark_day, endpoint = '365-day all-cause mortality', .before = 3)
  landmark_count_rows[[paste0('mimic_', landmark_day)]] <- tibble::tibble(
    dataset = 'MIMIC-IV', endpoint = '365-day all-cause mortality',
    landmark_day = landmark_day, original_n = nrow(mimic),
    original_events = sum(mimic$mortality_365d == 1),
    landmark_n = nrow(mimic_landmark),
    landmark_events = sum(mimic_landmark$mortality_365d == 1),
    early_events_excluded = sum(mimic$mortality_365d == 1 & mimic$survival_days_365 <= landmark_day),
    early_non_events_or_no_longer_observed_excluded = sum(
      mimic$mortality_365d == 0 & mimic$survival_days_365 <= landmark_day
    )
  )
}

eicu_primary_analysis <- eicu_feature_complete |>
  dplyr::left_join(apache, by = 'patientunitstayid') |>
  dplyr::filter(is.finite(apachescore)) |>
  dplyr::mutate(
    phenotype = factor(phenotype_primary, levels = c('P3', 'P2', 'P1')),
    hospitalid = factor(hospitalid),
    age_z = as.numeric(scale(age_num)),
    apache_z = as.numeric(scale(apachescore))
  )
eicu_hospital_base <- eicu_primary_analysis |>
  dplyr::filter(!is.na(hospital_mortality), is.finite(survival_days_hosp))
eicu_icu_base <- eicu_primary_analysis |>
  dplyr::filter(!is.na(icu_mortality), is.finite(icu_los_days))

for (landmark_day in c(1, 3, 7)) {
  hospital_landmark <- eicu_hospital_base |>
    dplyr::filter(survival_days_hosp > landmark_day) |>
    droplevels()
  icu_landmark <- eicu_icu_base |>
    dplyr::filter(icu_los_days > landmark_day) |>
    droplevels()
  hospital_capture <- capture_glmer(
    hospital_mortality ~ phenotype + age_z + male + apache_z + (1 | hospitalid),
    hospital_landmark
  )
  icu_capture <- capture_glmer(
    icu_mortality ~ phenotype + age_z + male + apache_z + (1 | hospitalid),
    icu_landmark
  )
  landmark_effect_rows[[paste0('eicu_hospital_', landmark_day)]] <- extract_effects(
    hospital_capture$fit, 'phenotype', 'eICU',
    paste0(landmark_day, '-day landmark, APACHE-adjusted hospital random-intercept logistic'),
    'OR', nrow(hospital_landmark), sum(hospital_landmark$hospital_mortality == 1)
  ) |>
    dplyr::mutate(landmark_day = landmark_day, endpoint = 'In-hospital mortality', .before = 3)
  landmark_effect_rows[[paste0('eicu_icu_', landmark_day)]] <- extract_effects(
    icu_capture$fit, 'phenotype', 'eICU',
    paste0(landmark_day, '-day landmark, APACHE-adjusted ICU random-intercept logistic'),
    'OR', nrow(icu_landmark), sum(icu_landmark$icu_mortality == 1)
  ) |>
    dplyr::mutate(landmark_day = landmark_day, endpoint = 'ICU mortality', .before = 3)
  landmark_count_rows[[paste0('eicu_hospital_', landmark_day)]] <- tibble::tibble(
    dataset = 'eICU', endpoint = 'In-hospital mortality', landmark_day = landmark_day,
    original_n = nrow(eicu_hospital_base),
    original_events = sum(eicu_hospital_base$hospital_mortality == 1),
    landmark_n = nrow(hospital_landmark),
    landmark_events = sum(hospital_landmark$hospital_mortality == 1),
    early_events_excluded = sum(
      eicu_hospital_base$hospital_mortality == 1 &
        eicu_hospital_base$survival_days_hosp <= landmark_day
    ),
    early_non_events_or_no_longer_observed_excluded = sum(
      eicu_hospital_base$hospital_mortality == 0 &
        eicu_hospital_base$survival_days_hosp <= landmark_day
    )
  )
  landmark_count_rows[[paste0('eicu_icu_', landmark_day)]] <- tibble::tibble(
    dataset = 'eICU', endpoint = 'ICU mortality', landmark_day = landmark_day,
    original_n = nrow(eicu_icu_base),
    original_events = sum(eicu_icu_base$icu_mortality == 1),
    landmark_n = nrow(icu_landmark),
    landmark_events = sum(icu_landmark$icu_mortality == 1),
    early_events_excluded = sum(
      eicu_icu_base$icu_mortality == 1 & eicu_icu_base$icu_los_days <= landmark_day
    ),
    early_non_events_or_no_longer_observed_excluded = sum(
      eicu_icu_base$icu_mortality == 0 & eicu_icu_base$icu_los_days <= landmark_day
    )
  )
  landmark_diagnostic_rows[[paste0('eicu_hospital_', landmark_day)]] <- glmer_diagnostics(
    hospital_capture, 'eICU', paste0(landmark_day, '-day hospital landmark'),
    hospital_landmark, 'hospital_mortality'
  )
  landmark_diagnostic_rows[[paste0('eicu_icu_', landmark_day)]] <- glmer_diagnostics(
    icu_capture, 'eICU', paste0(landmark_day, '-day ICU landmark'),
    icu_landmark, 'icu_mortality'
  )
}

landmark_models <- dplyr::bind_rows(landmark_effect_rows)
landmark_counts <- dplyr::bind_rows(landmark_count_rows)
landmark_diagnostics <- dplyr::bind_rows(landmark_diagnostic_rows)

cat('[5/5] QA, aggregate outputs, and result summary...\n')
qa <- tibble::tribble(
  ~check, ~passed, ~detail,
  'NHANES model cohort reproduced', nrow(nhanes_renal) == 3979L && sum(nhanes_renal$MORTSTAT == 1) == 720L, 'Expected 3979/720',
  'No-BMI clustering retained fixed cohorts',
    nrow(nhanes_cluster_input) == 4637L && nrow(mimic_cluster_input) == 1100L && nrow(eicu_cluster_input) == 15242L,
    'Expected 4637; 1100; 15242',
  'No-BMI outcome model cohorts reproduced',
    all(no_bmi_outcome_models$n[no_bmi_outcome_models$dataset == 'NHANES'] == 3979L) &&
      all(no_bmi_outcome_models$n[no_bmi_outcome_models$dataset == 'MIMIC-IV'] == 1100L) &&
      all(no_bmi_outcome_models$n[no_bmi_outcome_models$dataset == 'eICU'] == 13234L),
    'Expected 3979; 1100; 13234',
  'All no-BMI estimates finite', all(is.finite(no_bmi_outcome_models$estimate)), 'No missing estimate',
  'All landmark estimates finite', all(is.finite(landmark_models$estimate)), 'No missing estimate',
  'Landmark samples decline monotonically',
    all(vapply(split(landmark_counts, interaction(landmark_counts$dataset, landmark_counts$endpoint)), function(x) {
      all(diff(x$landmark_n[order(x$landmark_day)]) <= 0)
    }, logical(1))),
    'n at day 1 >= day 3 >= day 7',
  'No eICU model singular',
    !no_bmi_diagnostics$singular && !any(landmark_diagnostics$singular),
    'All random-intercept fits non-singular',
  'No eICU convergence messages or captured warnings',
    no_bmi_diagnostics$convergence_message == '' && no_bmi_diagnostics$captured_warnings == '' &&
      all(landmark_diagnostics$convergence_message == '') &&
      all(landmark_diagnostics$captured_warnings == ''),
    'All scaled continuous-covariate fits completed cleanly',
  'No outcomes used in no-BMI clustering or mapping', TRUE, 'Outcome models fitted only after assignments were produced',
  'No row-level outputs written', TRUE, 'Only aggregate tables and summary RDS'
)
if (!all(qa$passed)) {
  print(qa, n = Inf, width = Inf)
  stop('Reviewer-sensitivity QA failed before output write.', call. = FALSE)
}

tables <- list(
  Table87A_NHANES_creatinine_component_adjustment = nhanes_creatinine_models,
  Table87B_no_BMI_method = no_bmi_method,
  Table87C_no_BMI_assignment_agreement = no_bmi_agreement,
  Table87D_no_BMI_profiles = no_bmi_profiles,
  Table87E_no_BMI_outcome_models = no_bmi_outcome_models,
  Table87F_no_BMI_model_diagnostics = no_bmi_diagnostics,
  Table87G_landmark_counts = landmark_counts,
  Table87H_landmark_models = landmark_models,
  Table87I_landmark_model_diagnostics = landmark_diagnostics,
  Table87J_QA = qa
)
for (name in names(tables)) {
  readr::write_csv(tables[[name]], file.path(output_dir, paste0(name, '.csv')))
}
saveRDS(tables, file.path(output_dir, 'post_result_reviewer_sensitivity_results.rds'))

source_paths <- c(
  unname(unlist(required_inputs)), apache_source,
  file.path(root, '87_post_result_reviewer_sensitivities_2026-08-30.R')
)
source_manifest <- tibble::tibble(
  file = basename(source_paths),
  path = normalizePath(source_paths, winslash = '/', mustWork = TRUE),
  bytes = file.info(source_paths)$size,
  sha256 = vapply(source_paths, sha256_file, character(1))
)
readr::write_csv(source_manifest, file.path(output_dir, 'SOURCE_SHA256_MANIFEST.csv'))

p1_creatinine <- nhanes_creatinine_models |> dplyr::filter(comparison == 'P1 vs P3')
p1_no_bmi <- no_bmi_outcome_models |> dplyr::filter(comparison == 'P1 vs P3')
p1_landmark <- landmark_models |> dplyr::filter(comparison == 'P1 vs P3')
summary_lines <- c(
  'Post-result reviewer-driven exploratory sensitivity analyses',
  '',
  'Status: all modules were specified after primary results were known and must be labelled exploratory.',
  '',
  'NHANES creatinine-component adjustment (P1 vs P3):',
  paste(capture.output(print(p1_creatinine, n = Inf, width = Inf)), collapse = '\n'),
  '',
  'No-BMI same-cohort K=3 assignment agreement:',
  paste(capture.output(print(no_bmi_agreement, n = Inf, width = Inf)), collapse = '\n'),
  '',
  'No-BMI same-cohort K=3 outcome models (P1 vs P3):',
  paste(capture.output(print(p1_no_bmi, n = Inf, width = Inf)), collapse = '\n'),
  '',
  'ICU landmark outcome models (P1 vs P3):',
  paste(capture.output(print(p1_landmark, n = Inf, width = Inf)), collapse = '\n'),
  '',
  'Interpretation boundaries:',
  '- Creatinine adjustment is a part-whole/non-reducibility analysis because creatinine helped define P1; it is not ordinary confounder control.',
  '- The no-BMI analysis tests dependence on BMI within fixed cohorts; it does not validate a new five-variable classifier.',
  '- Landmark restriction reduces sensitivity to very early deaths but cannot eliminate reverse causation and introduces conditioning on remaining under observation.',
  '- All estimates must be reported irrespective of direction or statistical significance.'
)
writeLines(summary_lines, file.path(output_dir, 'POST_RESULT_REVIEWER_SENSITIVITY_SUMMARY.md'))

output_files <- list.files(output_dir, full.names = TRUE, recursive = FALSE)
output_files <- output_files[basename(output_files) != 'OUTPUT_SHA256_MANIFEST.csv']
output_manifest <- tibble::tibble(
  file = basename(output_files), bytes = file.info(output_files)$size,
  sha256 = vapply(output_files, sha256_file, character(1))
)
readr::write_csv(output_manifest, file.path(output_dir, 'OUTPUT_SHA256_MANIFEST.csv'))

cat(paste(summary_lines, collapse = '\n'), '\n')
message('Post-result reviewer sensitivity analyses completed.')
