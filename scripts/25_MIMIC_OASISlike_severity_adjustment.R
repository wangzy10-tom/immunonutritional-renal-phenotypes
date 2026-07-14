# ==============================================================================
# Script: 25_MIMIC_OASISlike_severity_adjustment.R
# Purpose: Add a first-24h OASIS-like acute severity adjustment to the MIMIC-IV
#          de novo validation analysis.
#
# Important interpretation:
#   This reconstructs an OASIS-like score from local MIMIC-IV v3.1 raw tables
#   using the official OASIS component thresholds. Local raw-table mappings are
#   retained, so it is still labelled OASIS-like rather than official OASIS.
#   It should be reported as "OASIS-like" or "first-24h acute severity score",
#   not as an official validated OASIS derived table.
# ==============================================================================

options(stringsAsFactors = FALSE)

required_pkgs <- c("DBI", "duckdb", "readr", "dplyr", "tibble", "survival", "ggplot2")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "))
}

suppressPackageStartupMessages({
  library(DBI)
  library(duckdb)
  library(readr)
  library(dplyr)
  library(tibble)
  library(survival)
  library(ggplot2)
})

project_root <- normalizePath(
  Sys.getenv("PROJECT_ROOT", unset = getwd()), winslash = "/", mustWork = TRUE
)
mimic_root <- Sys.getenv(
  "MIMIC_ROOT", unset = file.path(project_root, "data", "mimic-iv-3.1")
)
validation_dir <- Sys.getenv(
  "MIMIC_SEVERITY_SOURCE_DIR",
  unset = file.path(project_root, "output", "mimic_24h_validation")
)
root_out <- Sys.getenv(
  "MIMIC_SEVERITY_OUTPUT",
  unset = file.path(project_root, "output", "mimic_oasislike_severity_corrected")
)
dir.create(root_out, recursive = TRUE, showWarnings = FALSE)

source_dataset <- file.path(validation_dir, "MIMIC_denovo_albumin_proxy_dataset.csv")
if (!file.exists(source_dataset)) {
  stop("MIMIC de novo dataset not found: ", source_dataset, call. = FALSE)
}

paths <- list(
  admissions = file.path(mimic_root, "hosp", "admissions.csv.gz"),
  services = file.path(mimic_root, "hosp", "services.csv.gz"),
  chartevents = file.path(mimic_root, "icu", "chartevents.csv.gz"),
  outputevents = file.path(mimic_root, "icu", "outputevents.csv.gz"),
  procedureevents = file.path(mimic_root, "icu", "procedureevents.csv.gz")
)
missing_files <- paths[!vapply(paths, file.exists, logical(1))]
if (length(missing_files) > 0) {
  stop("Missing MIMIC files: ", paste(unlist(missing_files), collapse = "; "), call. = FALSE)
}

fmt_p <- function(p) ifelse(is.na(p), NA_character_, ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
fmt_est <- function(est, lo, hi) sprintf("%.2f (%.2f-%.2f)", est, lo, hi)
z <- function(x) as.numeric(scale(as.numeric(x)))

stable_ntile <- function(value, groups, tie_breaker) {
  keep <- is.finite(value) & !is.na(tie_breaker)
  output <- rep(NA_integer_, length(value))
  ordered_index <- which(keep)[order(value[keep], tie_breaker[keep])]
  output[ordered_index] <- dplyr::ntile(seq_along(ordered_index), groups)
  output
}

score_hr_one <- function(x) {
  case_when(
    is.na(x) ~ NA_real_,
    x < 33 ~ 4,
    x < 89 ~ 0,
    x < 107 ~ 1,
    x < 126 ~ 3,
    TRUE ~ 6
  )
}

score_map_one <- function(x) {
  case_when(
    is.na(x) ~ NA_real_,
    x < 20.65 ~ 4,
    x < 51 ~ 3,
    x < 61.33 ~ 2,
    x < 143.44 ~ 0,
    TRUE ~ 3
  )
}

score_rr_one <- function(x) {
  case_when(
    is.na(x) ~ NA_real_,
    x < 6 ~ 10,
    x < 14 ~ 1,
    x < 23 ~ 0,
    x < 31 ~ 1,
    x < 45 ~ 6,
    TRUE ~ 9
  )
}

score_temp_one <- function(x) {
  case_when(
    is.na(x) ~ NA_real_,
    x < 33.22 ~ 3,
    x < 35.93 ~ 4,
    x < 36.39 ~ 2,
    x < 36.89 ~ 0,
    x < 39.88 ~ 2,
    TRUE ~ 6
  )
}

score_urine <- function(x) {
  case_when(
    is.na(x) ~ NA_real_,
    x < 671.1 ~ 10,
    x < 1427 ~ 5,
    x < 2544.1 ~ 1,
    x < 6896.8 ~ 0,
    TRUE ~ 8
  )
}

score_gcs <- function(x) {
  case_when(
    is.na(x) ~ NA_real_,
    x <= 7 ~ 10,
    x < 14 ~ 4,
    x == 14 ~ 3,
    TRUE ~ 0
  )
}

score_age <- function(x) {
  case_when(
    is.na(x) ~ NA_real_,
    x < 24 ~ 0,
    x <= 53 ~ 3,
    x <= 77 ~ 6,
    x <= 89 ~ 9,
    TRUE ~ 7
  )
}

score_preicu <- function(hours) {
  minutes <- hours * 60
  case_when(
    is.na(minutes) ~ NA_real_,
    minutes < 10.2 ~ 5,
    minutes < 297 ~ 3,
    minutes < 1440 ~ 0,
    minutes < 18708 ~ 2,
    TRUE ~ 1
  )
}

extract_cox <- function(fit) {
  s <- summary(fit)
  ci <- as.data.frame(s$conf.int)
  co <- as.data.frame(s$coefficients)
  tibble(
    variable = rownames(ci),
    hazard_ratio_95ci = fmt_est(ci[, "exp(coef)"], ci[, "lower .95"], ci[, "upper .95"]),
    HR = ci[, "exp(coef)"],
    lower_95 = ci[, "lower .95"],
    upper_95 = ci[, "upper .95"],
    p_value = co[, "Pr(>|z|)"],
    p_value_formatted = fmt_p(p_value)
  )
}

cat("Loading MIMIC de novo validation dataset...\n")
mimic <- read_csv(source_dataset, show_col_types = FALSE) %>%
  mutate(
    stay_id = as.integer(stay_id),
    hadm_id = as.integer(hadm_id),
    subject_id = as.integer(subject_id),
    intime = as.POSIXct(intime, tz = "UTC"),
    admittime = as.POSIXct(admittime, tz = "UTC"),
    anchor_age = as.numeric(anchor_age),
    male = as.integer(gender == "M"),
    death365 = as.integer(mortality_365d == 1),
    time365 = as.numeric(survival_days_365),
    mimic_phenotype = relevel(factor(as.character(mimic_phenotype_raw), levels = c("1", "2", "3")), ref = "3"),
    p1_binary = as.integer(mimic_phenotype_raw == 1),
    vulnerability_score_raw = z(nlr) + z(sii) + z(creatinine) - z(haemoglobin) - z(protein_proxy) - z(bmi),
    vulnerability_score = z(vulnerability_score_raw)
  ) %>%
  filter(complete.cases(stay_id, hadm_id, subject_id, intime, admittime, anchor_age, male, death365, time365))

cohort_for_duckdb <- mimic %>%
  distinct(stay_id, hadm_id, subject_id, intime, admittime) %>%
  mutate(
    intime = as.POSIXct(intime, tz = "UTC"),
    admittime = as.POSIXct(admittime, tz = "UTC")
  )

cat("Connecting to DuckDB and extracting first-24h severity components...\n")
con <- DBI::dbConnect(duckdb::duckdb(), dbdir = ":memory:")
on.exit(try(DBI::dbDisconnect(con, shutdown = TRUE), silent = TRUE), add = TRUE)
DBI::dbWriteTable(con, "cohort", cohort_for_duckdb, overwrite = TRUE)

chartevents_path <- normalizePath(paths$chartevents, winslash = "/", mustWork = TRUE)
outputevents_path <- normalizePath(paths$outputevents, winslash = "/", mustWork = TRUE)
procedureevents_path <- normalizePath(paths$procedureevents, winslash = "/", mustWork = TRUE)
admissions_path <- normalizePath(paths$admissions, winslash = "/", mustWork = TRUE)
services_path <- normalizePath(paths$services, winslash = "/", mustWork = TRUE)

vitals_sql <- sprintf(
  "
  SELECT
    ce.stay_id,
    MIN(CASE WHEN ce.itemid = 220045 AND ce.valuenum BETWEEN 1 AND 300 THEN ce.valuenum END) AS hr_min,
    MAX(CASE WHEN ce.itemid = 220045 AND ce.valuenum BETWEEN 1 AND 300 THEN ce.valuenum END) AS hr_max,
    MIN(CASE WHEN ce.itemid IN (220052, 220181) AND ce.valuenum BETWEEN 1 AND 300 THEN ce.valuenum END) AS map_min,
    MAX(CASE WHEN ce.itemid IN (220052, 220181) AND ce.valuenum BETWEEN 1 AND 300 THEN ce.valuenum END) AS map_max,
    MIN(CASE WHEN ce.itemid IN (220210, 224689, 224690) AND ce.valuenum BETWEEN 1 AND 80 THEN ce.valuenum END) AS rr_min,
    MAX(CASE WHEN ce.itemid IN (220210, 224689, 224690) AND ce.valuenum BETWEEN 1 AND 80 THEN ce.valuenum END) AS rr_max,
    MIN(
      CASE
        WHEN ce.itemid IN (223762, 226329) AND ce.valuenum BETWEEN 20 AND 45 THEN ce.valuenum
        WHEN ce.itemid = 223761 AND ce.valuenum BETWEEN 60 AND 115 THEN (ce.valuenum - 32) / 1.8
      END
    ) AS temp_c_min,
    MAX(
      CASE
        WHEN ce.itemid IN (223762, 226329) AND ce.valuenum BETWEEN 20 AND 45 THEN ce.valuenum
        WHEN ce.itemid = 223761 AND ce.valuenum BETWEEN 60 AND 115 THEN (ce.valuenum - 32) / 1.8
      END
    ) AS temp_c_max
  FROM read_csv_auto('%s', union_by_name = true) ce
  INNER JOIN cohort c ON ce.stay_id = c.stay_id
  WHERE ce.itemid IN (220045, 220052, 220181, 220210, 224689, 224690, 223761, 223762, 226329)
    AND ce.valuenum IS NOT NULL
    AND ce.charttime >= c.intime
    AND ce.charttime < c.intime + INTERVAL 24 HOUR
  GROUP BY ce.stay_id
  ",
  chartevents_path
)
vitals <- DBI::dbGetQuery(con, vitals_sql)

gcs_sql <- sprintf(
  "
  WITH gcs_by_time AS (
    SELECT
      ce.stay_id,
      ce.charttime,
      MAX(CASE WHEN ce.itemid IN (220739, 226756, 227011) AND ce.valuenum BETWEEN 1 AND 4 THEN ce.valuenum END) AS eye,
      MAX(CASE WHEN ce.itemid IN (223901, 226757, 227012) AND ce.valuenum BETWEEN 1 AND 6 THEN ce.valuenum END) AS motor,
      MAX(CASE WHEN ce.itemid IN (223900, 226758, 227014, 228112) AND ce.valuenum BETWEEN 1 AND 5 THEN ce.valuenum END) AS verbal,
      MAX(CASE WHEN ce.itemid IN (226755, 227013) AND ce.valuenum BETWEEN 3 AND 15 THEN ce.valuenum END) AS direct_gcs
    FROM read_csv_auto('%s', union_by_name = true) ce
    INNER JOIN cohort c ON ce.stay_id = c.stay_id
    WHERE ce.itemid IN (220739, 226756, 227011, 223901, 226757, 227012, 223900, 226758, 227014, 228112, 226755, 227013)
      AND ce.valuenum IS NOT NULL
      AND ce.charttime >= c.intime
      AND ce.charttime < c.intime + INTERVAL 24 HOUR
    GROUP BY ce.stay_id, ce.charttime
  ),
  gcs_scored AS (
    SELECT
      stay_id,
      CASE
        WHEN direct_gcs BETWEEN 3 AND 15 THEN direct_gcs
        WHEN eye IS NOT NULL AND motor IS NOT NULL AND verbal IS NOT NULL
          AND eye + motor + verbal BETWEEN 3 AND 15 THEN eye + motor + verbal
      END AS gcs_total
    FROM gcs_by_time
  )
  SELECT stay_id, MIN(gcs_total) AS gcs_min
  FROM gcs_scored
  WHERE gcs_total IS NOT NULL
  GROUP BY stay_id
  ",
  chartevents_path
)
gcs <- DBI::dbGetQuery(con, gcs_sql)

urine_itemids <- c(226557, 226558, 226559, 226560, 226561, 226563, 226564, 226565, 226566, 226567, 226627, 226631, 227489)
urine_sql <- sprintf(
  "
  SELECT
    oe.stay_id,
    SUM(TRY_CAST(oe.value AS DOUBLE)) AS urine_output_24h
  FROM read_csv_auto('%s', union_by_name = true) oe
  INNER JOIN cohort c ON oe.stay_id = c.stay_id
  WHERE oe.itemid IN (%s)
    AND TRY_CAST(oe.value AS DOUBLE) IS NOT NULL
    AND TRY_CAST(oe.value AS DOUBLE) BETWEEN 0 AND 10000
    AND oe.charttime >= c.intime
    AND oe.charttime < c.intime + INTERVAL 24 HOUR
  GROUP BY oe.stay_id
  ",
  outputevents_path,
  paste(urine_itemids, collapse = ",")
)
urine <- DBI::dbGetQuery(con, urine_sql)

vent_sql <- sprintf(
  "
  SELECT
    pe.stay_id,
    1 AS mechanical_ventilation_24h
  FROM read_csv_auto('%s', union_by_name = true) pe
  INNER JOIN cohort c ON pe.stay_id = c.stay_id
  WHERE pe.itemid IN (225792, 225794)
    AND pe.starttime < c.intime + INTERVAL 24 HOUR
    AND COALESCE(pe.endtime, pe.starttime) >= c.intime
  GROUP BY pe.stay_id
  ",
  procedureevents_path
)
vent <- DBI::dbGetQuery(con, vent_sql)

admission_sql <- sprintf(
  "
  SELECT hadm_id, admission_type
  FROM read_csv_auto('%s', union_by_name = true)
  ",
  admissions_path
)
admission_types <- DBI::dbGetQuery(con, admission_sql)

services_sql <- sprintf(
  "
  SELECT
    c.stay_id,
    MAX(
      CASE
        WHEN LOWER(s.curr_service) LIKE '%%surg%%' OR s.curr_service = 'ORTHO' THEN 1
        ELSE 0
      END
    ) AS surgical_service
  FROM cohort c
  LEFT JOIN read_csv_auto('%s', union_by_name = true) s
    ON c.hadm_id = s.hadm_id
   AND s.transfertime < c.intime + INTERVAL 24 HOUR
  GROUP BY c.stay_id
  ",
  services_path
)
surgical_services <- DBI::dbGetQuery(con, services_sql)

cat("Scoring OASIS-like components...\n")
severity <- cohort_for_duckdb %>%
  left_join(vitals, by = "stay_id") %>%
  left_join(gcs, by = "stay_id") %>%
  left_join(urine, by = "stay_id") %>%
  left_join(vent, by = "stay_id") %>%
  left_join(admission_types, by = "hadm_id") %>%
  left_join(surgical_services, by = "stay_id") %>%
  mutate(
    mechanical_ventilation_24h = if_else(is.na(mechanical_ventilation_24h), 0, mechanical_ventilation_24h),
    admission_type = if_else(is.na(admission_type), "", admission_type),
    preicu_los_hours = pmax(as.numeric(difftime(intime, admittime, units = "hours")), 0),
    elective_surgery = as.integer(
      grepl("ELECTIVE", admission_type, ignore.case = TRUE) & surgical_service == 1
    ),
    age_score = score_age(mimic$anchor_age[match(stay_id, mimic$stay_id)]),
    preicu_score = score_preicu(preicu_los_hours),
    hr_score = pmax(score_hr_one(hr_min), score_hr_one(hr_max), na.rm = FALSE),
    map_score = pmax(score_map_one(map_min), score_map_one(map_max), na.rm = FALSE),
    rr_score = pmax(score_rr_one(rr_min), score_rr_one(rr_max), na.rm = FALSE),
    temp_score = pmax(score_temp_one(temp_c_min), score_temp_one(temp_c_max), na.rm = FALSE),
    gcs_score = score_gcs(gcs_min),
    urine_score = score_urine(urine_output_24h),
    vent_score = if_else(mechanical_ventilation_24h == 1, 9, 0),
    elective_score = if_else(elective_surgery == 1, 0, 6)
  )

component_scores <- c(
  "age_score", "preicu_score", "hr_score", "map_score", "rr_score",
  "temp_score", "gcs_score", "urine_score", "vent_score", "elective_score"
)
severity <- severity %>%
  mutate(
    oasis_missing_components = rowSums(is.na(across(all_of(component_scores)))),
    oasis_like_complete = oasis_missing_components == 0,
    oasis_like_score = if_else(
      oasis_like_complete,
      rowSums(across(all_of(component_scores)), na.rm = FALSE),
      NA_real_
    )
  )

component_availability <- tibble(
  component = c("heart_rate", "mean_arterial_pressure", "respiratory_rate", "temperature", "gcs", "urine_output", "mechanical_ventilation", "admission_type/elective", "complete_oasis_like_score"),
  available_n = c(
    sum(!is.na(severity$hr_score)),
    sum(!is.na(severity$map_score)),
    sum(!is.na(severity$rr_score)),
    sum(!is.na(severity$temp_score)),
    sum(!is.na(severity$gcs_score)),
    sum(!is.na(severity$urine_score)),
    sum(!is.na(severity$vent_score)),
    sum(!is.na(severity$elective_score)),
    sum(severity$oasis_like_complete)
  ),
  denominator = nrow(severity)
) %>%
  mutate(available_percent = round(100 * available_n / denominator, 2))

mimic_severity <- mimic %>%
  left_join(
    severity %>%
      select(
        stay_id, preicu_los_hours, admission_type, elective_surgery,
        hr_min, hr_max, map_min, map_max, rr_min, rr_max, temp_c_min, temp_c_max,
        gcs_min, urine_output_24h, mechanical_ventilation_24h,
        all_of(component_scores), oasis_missing_components,
        oasis_like_complete, oasis_like_score
      ),
    by = "stay_id"
  )

analysis_complete <- mimic_severity %>%
  filter(complete.cases(time365, death365, mimic_phenotype, male, oasis_like_score)) %>%
  mutate(
    oasis_quartile = factor(
      stable_ntile(oasis_like_score, 4, stay_id),
      levels = 1:4, labels = paste0("Q", 1:4)
    )
  )

cat("Fitting MIMIC Cox models with OASIS-like adjustment...\n")
fit_pheno_age_sex <- coxph(Surv(time365, death365) ~ mimic_phenotype + anchor_age + gender, data = mimic_severity)
fit_pheno_oasis <- coxph(Surv(time365, death365) ~ mimic_phenotype + male + oasis_like_score, data = analysis_complete)
fit_p1_oasis <- coxph(Surv(time365, death365) ~ p1_binary + male + oasis_like_score, data = analysis_complete)
fit_score_oasis <- coxph(Surv(time365, death365) ~ vulnerability_score + male + oasis_like_score, data = analysis_complete)
fit_pheno_oasis_strata <- coxph(Surv(time365, death365) ~ mimic_phenotype + male + strata(oasis_quartile), data = analysis_complete)
fit_p1_oasis_strata <- coxph(Surv(time365, death365) ~ p1_binary + male + strata(oasis_quartile), data = analysis_complete)
fit_score_oasis_strata <- coxph(Surv(time365, death365) ~ vulnerability_score + male + strata(oasis_quartile), data = analysis_complete)

pheno_oasis <- extract_cox(fit_pheno_oasis) %>%
  mutate(
    variable = recode(
      variable,
      "mimic_phenotype1" = "Phenotype 1 vs Phenotype 3",
      "mimic_phenotype2" = "Phenotype 2 vs Phenotype 3",
      "male" = "Male vs female",
      "oasis_like_score" = "OASIS-like score, per point"
    )
  )

p1_oasis <- extract_cox(fit_p1_oasis) %>%
  mutate(
    variable = recode(
      variable,
      "p1_binary" = "Phenotype 1 vs non-Phenotype 1",
      "male" = "Male vs female",
      "oasis_like_score" = "OASIS-like score, per point"
    )
  )

score_oasis <- extract_cox(fit_score_oasis) %>%
  mutate(
    variable = recode(
      variable,
      "vulnerability_score" = "Vulnerability score, per 1 SD",
      "male" = "Male vs female",
      "oasis_like_score" = "OASIS-like score, per point"
    )
  )

pheno_oasis_strata <- extract_cox(fit_pheno_oasis_strata) %>%
  mutate(
    variable = recode(
      variable,
      "mimic_phenotype1" = "Phenotype 1 vs Phenotype 3",
      "mimic_phenotype2" = "Phenotype 2 vs Phenotype 3",
      "male" = "Male vs female"
    )
  )

p1_oasis_strata <- extract_cox(fit_p1_oasis_strata) %>%
  mutate(
    variable = recode(
      variable,
      "p1_binary" = "Phenotype 1 vs non-Phenotype 1",
      "male" = "Male vs female"
    )
  )

score_oasis_strata <- extract_cox(fit_score_oasis_strata) %>%
  mutate(
    variable = recode(
      variable,
      "vulnerability_score" = "Vulnerability score, per 1 SD",
      "male" = "Male vs female"
    )
  )

oasis_distribution <- analysis_complete %>%
  group_by(mimic_phenotype_raw) %>%
  summarise(
    n = n(),
    deaths_365d = sum(death365 == 1, na.rm = TRUE),
    mortality_365d_percent = round(100 * mean(death365 == 1, na.rm = TRUE), 2),
    oasis_like_median = median(oasis_like_score, na.rm = TRUE),
    oasis_like_q1 = quantile(oasis_like_score, 0.25, na.rm = TRUE),
    oasis_like_q3 = quantile(oasis_like_score, 0.75, na.rm = TRUE),
    oasis_like_iqr = sprintf("%.0f (%.0f-%.0f)", oasis_like_median, oasis_like_q1, oasis_like_q3),
    .groups = "drop"
  )

model_comparison <- bind_rows(
  extract_cox(fit_pheno_age_sex) %>%
    filter(variable %in% c("mimic_phenotype1", "mimic_phenotype2")) %>%
    mutate(model = "Age-sex adjusted original MIMIC model"),
  extract_cox(fit_pheno_oasis) %>%
    filter(variable %in% c("mimic_phenotype1", "mimic_phenotype2")) %>%
    mutate(model = "OASIS-like continuous adjusted model"),
  extract_cox(fit_pheno_oasis_strata) %>%
    filter(variable %in% c("mimic_phenotype1", "mimic_phenotype2")) %>%
    mutate(model = "OASIS-like quartile-stratified model")
) %>%
  mutate(
    variable = recode(
      variable,
      "mimic_phenotype1" = "Phenotype 1 vs Phenotype 3",
      "mimic_phenotype2" = "Phenotype 2 vs Phenotype 3"
    )
  )

plot_data <- model_comparison %>%
  mutate(
    label = paste(model, variable, sep = " | "),
    label = factor(label, levels = rev(label)),
    display = paste0("HR ", hazard_ratio_95ci)
  )

forest_plot <- ggplot(plot_data, aes(x = HR, y = label, xmin = lower_95, xmax = upper_95, color = model)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey45") +
  geom_segment(aes(x = lower_95, xend = upper_95, y = label, yend = label), linewidth = 0.7) +
  geom_point(size = 2.4) +
  geom_text(aes(x = 4.7, label = display), hjust = 0, size = 3.1, color = "grey20", show.legend = FALSE) +
  scale_x_log10(breaks = c(1, 1.5, 2, 3, 4), limits = c(0.85, 6.2)) +
  labs(
    x = "Hazard ratio on log scale",
    y = NULL,
    color = NULL,
    title = "MIMIC-IV phenotype effects before and after OASIS-like adjustment",
    subtitle = "OASIS-like score reconstructed from first-24h ICU physiology; HRs are not pooled."
  ) +
  theme_minimal(base_size = 10.5) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    plot.margin = margin(8, 18, 8, 8)
  )

ph_check <- survival::cox.zph(fit_pheno_oasis)
ph_table <- tibble(
  model = "OASIS-like continuous adjusted model",
  variable = rownames(ph_check$table),
  chisq = ph_check$table[, "chisq"],
  df = ph_check$table[, "df"],
  p = ph_check$table[, "p"]
)

ph_check_strata <- survival::cox.zph(fit_pheno_oasis_strata)
ph_table_strata <- tibble(
  model = "OASIS-like quartile-stratified model",
  variable = rownames(ph_check_strata$table),
  chisq = ph_check_strata$table[, "chisq"],
  df = ph_check_strata$table[, "df"],
  p = ph_check_strata$table[, "p"]
)

readr::write_csv(component_availability, file.path(root_out, "Table16A_MIMIC_OASISlike_component_availability.csv"))
readr::write_csv(oasis_distribution, file.path(root_out, "Table16B_MIMIC_OASISlike_distribution_by_phenotype.csv"))
readr::write_csv(pheno_oasis, file.path(root_out, "Table16C_MIMIC_OASISlike_adjusted_Cox_phenotypes.csv"))
readr::write_csv(p1_oasis, file.path(root_out, "Table16D_MIMIC_OASISlike_adjusted_Cox_P1_binary.csv"))
readr::write_csv(score_oasis, file.path(root_out, "Table16E_MIMIC_OASISlike_adjusted_Cox_vulnerability_score.csv"))
readr::write_csv(model_comparison, file.path(root_out, "Table16F_MIMIC_age_sex_vs_OASISlike_model_comparison.csv"))
readr::write_csv(ph_table, file.path(root_out, "Table16G_MIMIC_OASISlike_PH_check.csv"))
readr::write_csv(pheno_oasis_strata, file.path(root_out, "Table16H_MIMIC_OASISlike_stratified_Cox_phenotypes.csv"))
readr::write_csv(p1_oasis_strata, file.path(root_out, "Table16I_MIMIC_OASISlike_stratified_Cox_P1_binary.csv"))
readr::write_csv(score_oasis_strata, file.path(root_out, "Table16J_MIMIC_OASISlike_stratified_Cox_vulnerability_score.csv"))
readr::write_csv(ph_table_strata, file.path(root_out, "Table16K_MIMIC_OASISlike_stratified_PH_check.csv"))
readr::write_csv(mimic_severity, file.path(root_out, "MIMIC_denovo_with_OASISlike_score.csv"))
saveRDS(
  list(
    component_availability = component_availability,
    oasis_distribution = oasis_distribution,
    pheno_oasis = pheno_oasis,
    p1_oasis = p1_oasis,
    score_oasis = score_oasis,
    pheno_oasis_strata = pheno_oasis_strata,
    p1_oasis_strata = p1_oasis_strata,
    score_oasis_strata = score_oasis_strata,
    model_comparison = model_comparison,
    ph_check = ph_table,
    ph_check_strata = ph_table_strata
  ),
  file.path(root_out, "MIMIC_OASISlike_severity_adjustment_results.rds")
)
if (tolower(Sys.getenv("WRITE_FIGURES", unset = "false")) %in% c("true", "1", "yes", "y")) {
  ggsave(
    file.path(root_out, "Figure12_MIMIC_OASISlike_adjusted_forest.png"),
    forest_plot, width = 10.5, height = 4.8, dpi = 300
  )
}

summary_lines <- c(
  "MIMIC-IV OASIS-like first-24h acute severity adjustment",
  "",
  paste0("Source de novo MIMIC cohort: n = ", nrow(mimic)),
  paste0("Complete OASIS-like score available: n = ", nrow(analysis_complete),
         " (", round(100 * nrow(analysis_complete) / nrow(mimic), 2), "%)"),
  "",
  "Component availability:",
  paste(capture.output(print(component_availability)), collapse = "\n"),
  "",
  "OASIS-like score by phenotype:",
  paste(capture.output(print(oasis_distribution)), collapse = "\n"),
  "",
  "OASIS-like adjusted phenotype Cox model:",
  paste(capture.output(print(pheno_oasis %>% select(variable, hazard_ratio_95ci, p_value_formatted))), collapse = "\n"),
  "",
  "OASIS-like quartile-stratified phenotype Cox model:",
  paste(capture.output(print(pheno_oasis_strata %>% select(variable, hazard_ratio_95ci, p_value_formatted))), collapse = "\n"),
  "",
  "OASIS-like adjusted binary P1 model:",
  paste(capture.output(print(p1_oasis %>% select(variable, hazard_ratio_95ci, p_value_formatted))), collapse = "\n"),
  "",
  "OASIS-like quartile-stratified binary P1 model:",
  paste(capture.output(print(p1_oasis_strata %>% select(variable, hazard_ratio_95ci, p_value_formatted))), collapse = "\n"),
  "",
  "OASIS-like adjusted continuous vulnerability score model:",
  paste(capture.output(print(score_oasis %>% select(variable, hazard_ratio_95ci, p_value_formatted))), collapse = "\n"),
  "",
  "OASIS-like quartile-stratified continuous vulnerability score model:",
  paste(capture.output(print(score_oasis_strata %>% select(variable, hazard_ratio_95ci, p_value_formatted))), collapse = "\n"),
  "",
  "Interpretation note: report this as an OASIS-like score reconstructed from MIMIC-IV raw first-24h physiology, not as an official derived OASIS table. Because continuous OASIS-like score can violate proportional hazards, the quartile-stratified Cox model should be treated as the PH-robust severity-adjusted sensitivity."
)
writeLines(summary_lines, file.path(root_out, "MIMIC_OASISlike_severity_adjustment_summary.txt"))

cat("\nCompleted. Outputs saved to: ", root_out, "\n", sep = "")
cat(paste(summary_lines, collapse = "\n"))
