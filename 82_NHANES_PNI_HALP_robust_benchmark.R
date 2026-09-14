# ==============================================================================
# NHANES Amendment 20: robust PNI/HALP benchmark scaling
# Protocol locked before corrected outcome estimates were generated.
# Date: 2026-08-30
# ==============================================================================

required_packages <- c("dplyr", "readr", "tibble", "survey", "survival", "digest")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop("Missing required package(s): ", paste(missing_packages, collapse = ", "), call. = FALSE)
}
options(survey.lonely.psu = "adjust")

root <- normalizePath(Sys.getenv("PROJECT_ROOT", unset=getwd()), winslash="/", mustWork=TRUE)
source_root <- file.path(root, "output", "nhanes_albumin_amendment17_2026-08-30")
output_root <- file.path(root, "output", "nhanes_pni_halp_amendment20_2026-08-30")
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

source_path <- file.path(source_root, "NHANES_albumin_benchmark_results.rds")
stopifnot(file.exists(source_path))
source_results <- readRDS(source_path)
model_data <- source_results$model_data

required_columns <- c(
  "SEQN", "PERMTH_INT", "MORTSTAT", "SDMVPSU", "SDMVSTRA", "WTMEC8YR",
  "RIDAGEYR", "male", "race", "INDFMPIR", "Comorbidity_Score_Extended",
  "cycle", "smoking", "hypertension", "diabetes", "phenotype_albumin",
  "LBXSAL", "LBDLYMNO", "LBXHGB", "LBXPLTSI"
)
stopifnot(all(required_columns %in% names(model_data)))
stopifnot(nrow(model_data) == 3979L, sum(model_data$MORTSTAT == 1L) == 720L)

winsorise_with_limits <- function(x, lower = 0.01, upper = 0.99) {
  limits <- stats::quantile(x, c(lower, upper), na.rm = TRUE, names = FALSE)
  list(value = pmin(pmax(as.numeric(x), limits[1]), limits[2]), limits = limits)
}
safe_z <- function(x) {
  spread <- stats::sd(x, na.rm = TRUE)
  if (!is.finite(spread) || spread == 0) stop("Cannot standardize a constant variable.", call. = FALSE)
  (x - mean(x, na.rm = TRUE)) / spread
}

model_data <- model_data |>
  dplyr::mutate(
    albumin_g_dl_a20 = as.numeric(LBXSAL),
    lymphocytes_10e9_l_a20 = as.numeric(LBDLYMNO),
    platelets_10e9_l_a20 = as.numeric(LBXPLTSI),
    haemoglobin_g_l_a20 = as.numeric(LBXHGB) * 10,
    PNI_a20_raw = 10 * albumin_g_dl_a20 + 5 * lymphocytes_10e9_l_a20,
    HALP_a20_raw = haemoglobin_g_l_a20 * (albumin_g_dl_a20 * 10) *
      lymphocytes_10e9_l_a20 / platelets_10e9_l_a20
  )

stopifnot(
  all(is.finite(model_data$PNI_a20_raw)), all(model_data$PNI_a20_raw > 0),
  all(is.finite(model_data$HALP_a20_raw)), all(model_data$HALP_a20_raw > 0)
)

pni_w <- winsorise_with_limits(model_data$PNI_a20_raw)
halp_w <- winsorise_with_limits(model_data$HALP_a20_raw)
model_data <- model_data |>
  dplyr::mutate(
    PNI_a20_winsor = pni_w$value,
    HALP_a20_winsor = halp_w$value,
    risk_low_pni_a20 = -safe_z(PNI_a20_winsor),
    risk_low_halp_a20 = -safe_z(log(HALP_a20_winsor))
  )

distribution_row <- function(index, raw, transformed, limits, risk_score) {
  tibble::tibble(
    index = index,
    n = length(raw),
    raw_min = min(raw),
    raw_p01 = stats::quantile(raw, 0.01, names = FALSE),
    raw_median = stats::median(raw),
    raw_p99 = stats::quantile(raw, 0.99, names = FALSE),
    raw_max = max(raw),
    lower_winsor_limit = limits[1],
    upper_winsor_limit = limits[2],
    n_below_lower_limit = sum(raw < limits[1]),
    n_above_upper_limit = sum(raw > limits[2]),
    transformed_min = min(transformed),
    transformed_max = max(transformed),
    risk_score_min = min(risk_score),
    risk_score_max = max(risk_score)
  )
}

distribution_audit <- dplyr::bind_rows(
  distribution_row(
    "PNI", model_data$PNI_a20_raw, model_data$PNI_a20_winsor,
    pni_w$limits, model_data$risk_low_pni_a20
  ),
  distribution_row(
    "HALP", model_data$HALP_a20_raw, log(model_data$HALP_a20_winsor),
    halp_w$limits, model_data$risk_low_halp_a20
  )
)

survey_design <- survey::svydesign(
  ids = ~SDMVPSU,
  strata = ~SDMVSTRA,
  weights = ~WTMEC8YR,
  nest = TRUE,
  data = model_data
)
base_terms <- paste(
  "RIDAGEYR + male + race + INDFMPIR + Comorbidity_Score_Extended +",
  "cycle + smoking + hypertension + diabetes"
)

format_effect <- function(beta, se) {
  tibble::tibble(
    HR = exp(as.numeric(beta)),
    lower_95 = exp(as.numeric(beta) - 1.96 * as.numeric(se)),
    upper_95 = exp(as.numeric(beta) + 1.96 * as.numeric(se)),
    p_value = 2 * stats::pnorm(abs(as.numeric(beta) / as.numeric(se)), lower.tail = FALSE)
  ) |>
    dplyr::mutate(hazard_ratio_95ci = sprintf("%.2f (%.2f-%.2f)", HR, lower_95, upper_95))
}

score_map <- c(PNI = "risk_low_pni_a20", HALP = "risk_low_halp_a20")
benchmark_effects <- dplyr::bind_rows(lapply(names(score_map), function(index) {
  term <- score_map[[index]]
  fit <- survey::svycoxph(
    stats::as.formula(paste("survival::Surv(PERMTH_INT, MORTSTAT) ~", term, "+", base_terms)),
    design = survey_design
  )
  beta <- stats::coef(fit)[term]
  se <- sqrt(stats::vcov(fit)[term, term])
  format_effect(beta, se) |>
    dplyr::mutate(
      index = index,
      contrast = "Per 1-SD higher lower-value risk direction",
      n = nrow(model_data),
      events = sum(model_data$MORTSTAT == 1L),
      .before = 1
    )
}))

incremental_effects <- dplyr::bind_rows(lapply(names(score_map), function(index) {
  comparator <- score_map[[index]]
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
      comparator = index,
      n = nrow(model_data),
      events = sum(model_data$MORTSTAT == 1L),
      global_phenotype_wald_p = as.numeric(joint$p),
      .before = 1
    )
}))

readr::write_csv(distribution_audit, file.path(output_root, "Table82A_PNI_HALP_distribution_audit.csv"))
readr::write_csv(benchmark_effects, file.path(output_root, "Table82B_corrected_PNI_HALP_benchmarks.csv"))
readr::write_csv(incremental_effects, file.path(output_root, "Table82C_P1_beyond_corrected_PNI_HALP.csv"))
saveRDS(
  list(
    distribution_audit = distribution_audit,
    benchmark_effects = benchmark_effects,
    incremental_effects = incremental_effects
  ),
  file.path(output_root, "NHANES_Amendment20_PNI_HALP_results.rds")
)

summary_lines <- c(
  "# Amendment 20 results — robust PNI/HALP benchmark scaling",
  "",
  paste0("Analysis sample: n=", nrow(model_data), "; deaths=", sum(model_data$MORTSTAT == 1L), "."),
  "",
  "## Distribution audit",
  "",
  paste(capture.output(print(distribution_audit)), collapse = "\n"),
  "",
  "## Corrected standalone benchmarks",
  "",
  paste(capture.output(print(benchmark_effects)), collapse = "\n"),
  "",
  "## P1 after additional comparator adjustment",
  "",
  paste(capture.output(print(incremental_effects)), collapse = "\n"),
  "",
  "The corrected results supersede only the PNI/HALP rows in the current manuscript and supplement. All primary phenotype analyses remain unchanged."
)
writeLines(summary_lines, file.path(output_root, "AMENDMENT20_RESULTS_SUMMARY.md"))

source_files <- c(
  source_path,
  file.path(root, "82_NHANES_PNI_HALP_robust_benchmark.R"),
  file.path(root, "ANALYSIS_FREEZE_AMENDMENT_20_PNI_HALP_ROBUST_SCALING_PROTOCOL_2026-08-30.md")
)
source_manifest <- tibble::tibble(
  source = source_files,
  bytes = file.info(source_files)$size,
  sha256 = vapply(
    source_files,
    function(path) digest::digest(file = path, algo = "sha256", serialize = FALSE),
    character(1)
  )
)
readr::write_csv(source_manifest, file.path(output_root, "AMENDMENT20_SOURCE_MANIFEST_SHA256.csv"))

output_manifest_path <- file.path(output_root, "AMENDMENT20_OUTPUT_MANIFEST_SHA256.csv")
output_files <- list.files(output_root, full.names = TRUE)
output_files <- output_files[
  normalizePath(output_files, winslash = "/", mustWork = FALSE) !=
    normalizePath(output_manifest_path, winslash = "/", mustWork = FALSE)
]
output_manifest <- tibble::tibble(
  file = basename(output_files),
  bytes = file.info(output_files)$size,
  sha256 = vapply(
    output_files,
    function(path) digest::digest(file = path, algo = "sha256", serialize = FALSE),
    character(1)
  )
) |>
  dplyr::arrange(file)
readr::write_csv(output_manifest, output_manifest_path)

print(distribution_audit)
print(benchmark_effects)
print(incremental_effects)
