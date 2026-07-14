# ==============================================================================
# Script: 10_eICU_multicenter_validation.R
# Purpose:
#   Perform a multicentre eICU-CRD 2.0 de novo reproducibility analysis for the
#   NHANES-derived immunonutritional phenotype concept.
#
# Design:
#   - Older first ICU stays, age >= 65 years
#   - Early ICU laboratory window: 0 h to +24 h after ICU admission
#   - Primary feature matrix:
#       NLR-like ratio, SII-like index, haemoglobin, serum total protein, BMI,
#       and creatinine
#   - Unsupervised K-means with K = 3
#   - Outcome: in-hospital mortality / hospital-discharge survival
#
# Important interpretation:
#   eICU does not provide long-term mortality follow-up. This script therefore
#   tests multicentre ICU transportability using in-hospital outcomes, not
#   365-day mortality.
# ==============================================================================

options(stringsAsFactors = FALSE)

project_root <- normalizePath(
  Sys.getenv("PROJECT_ROOT", unset = getwd()), winslash = "/", mustWork = TRUE
)
eicu_dir <- Sys.getenv(
  "EICU_DIR",
  unset = file.path(project_root, "data", "eicu-crd-2.0")
)

output_dir <- Sys.getenv(
  "EICU_OUTPUT",
  unset = file.path(project_root, "output", "eicu_24h_extraction")
)

lab_window_before_minutes <- as.integer(Sys.getenv("EICU_LAB_WINDOW_BEFORE_MINUTES", unset = "0"))
lab_window_after_minutes <- as.integer(Sys.getenv("EICU_LAB_WINDOW_AFTER_MINUTES", unset = "1440"))

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

required_packages <- c("data.table", "survival", "ggplot2")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop(
    "Missing required R packages: ", paste(missing_packages, collapse = ", "),
    "\nPlease install them before rerunning this script.",
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(data.table)
  library(survival)
  library(ggplot2)
})

format_p <- function(p) {
  ifelse(is.na(p), NA_character_, ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
}

iqr_text <- function(x, digits = 2) {
  x <- as.numeric(x)
  if (all(is.na(x))) {
    return(NA_character_)
  }
  sprintf(
    paste0("%.", digits, "f (%.", digits, "f-%.", digits, "f)"),
    stats::median(x, na.rm = TRUE),
    stats::quantile(x, 0.25, na.rm = TRUE, names = FALSE),
    stats::quantile(x, 0.75, na.rm = TRUE, names = FALSE)
  )
}

write_csv_utf8 <- function(x, filename) {
  data.table::fwrite(x, file.path(output_dir, filename), bom = TRUE)
}

message("Output directory: ", output_dir)
message("Reading eICU patient table...")

patient_path <- file.path(eicu_dir, "patient.csv.gz")
lab_path <- file.path(eicu_dir, "lab.csv.gz")

if (!file.exists(patient_path)) {
  stop("Missing patient table: ", patient_path, call. = FALSE)
}
if (!file.exists(lab_path)) {
  stop("Missing lab table: ", lab_path, call. = FALSE)
}

patient_cols <- c(
  "patientunitstayid", "patienthealthsystemstayid", "uniquepid",
  "gender", "age", "ethnicity", "hospitalid", "wardid", "unittype",
  "unitvisitnumber", "admissionheight", "admissionweight",
  "hospitaldischargeoffset", "unitdischargeoffset",
  "hospitaldischargestatus", "unitdischargestatus"
)
patient <- data.table::as.data.table(utils::read.csv(
  gzfile(patient_path, open = "rt"),
  stringsAsFactors = FALSE,
  na.strings = c("", "NA", "NaN")
))
patient <- patient[, ..patient_cols]

patient[, age_num := fifelse(grepl("^>", age), 90, suppressWarnings(as.numeric(age)))]
patient[, unitvisitnumber_num := suppressWarnings(as.numeric(unitvisitnumber))]
patient[, first_icu_stay := is.na(unitvisitnumber_num) | unitvisitnumber_num == 1]
patient[, bmi := suppressWarnings(as.numeric(admissionweight)) / ((suppressWarnings(as.numeric(admissionheight)) / 100)^2)]
patient[!is.finite(bmi) | bmi < 10 | bmi > 80, bmi := NA_real_]

patient[, hospital_mortality := fifelse(
  hospitaldischargestatus == "Expired", 1L,
  fifelse(hospitaldischargestatus == "Alive", 0L, NA_integer_)
)]
patient[, gender_model := fifelse(gender %in% c("Female", "Male"), gender, NA_character_)]
patient[, icu_mortality := fifelse(
  unitdischargestatus == "Expired", 1L,
  fifelse(unitdischargestatus == "Alive", 0L, NA_integer_)
)]

patient[, survival_days_hosp := suppressWarnings(as.numeric(hospitaldischargeoffset)) / 1440]
patient[!is.finite(survival_days_hosp) | survival_days_hosp <= 0,
        survival_days_hosp := suppressWarnings(as.numeric(unitdischargeoffset)) / 1440]
patient[!is.finite(survival_days_hosp) | survival_days_hosp <= 0,
        survival_days_hosp := 0.01]

older <- patient[age_num >= 65 & first_icu_stay == TRUE]

# Keep one ICU stay per unique patient identifier where possible.
older[, person_id := fifelse(is.na(uniquepid) | uniquepid == "", as.character(patientunitstayid), as.character(uniquepid))]
data.table::setorder(older, person_id, patientunitstayid)
older <- older[, .SD[1], by = person_id]

message("Older first ICU patient-level cohort: ", nrow(older))

cohort_ids <- older$patientunitstayid
target_labnames <- c("-polys", "-lymphs", "platelets x 1000", "Hgb", "total protein", "albumin", "creatinine")
feature_map <- c(
  "-polys" = "polys_percent",
  "-lymphs" = "lymphs_percent",
  "platelets x 1000" = "platelets",
  "Hgb" = "haemoglobin",
  "total protein" = "total_protein",
  "albumin" = "albumin",
  "creatinine" = "creatinine"
)

lab_cache <- file.path(output_dir, "eICU_early_labs_wide_strict_window.rds")

if (file.exists(lab_cache)) {
  message("Using cached early laboratory matrix: ", lab_cache)
  labs_wide <- readRDS(lab_cache)
} else {
  message("Streaming eICU lab table in chunks. This can take several minutes...")
  con <- gzfile(lab_path, open = "rt")
  on.exit(close(con), add = TRUE)
  header <- readLines(con, n = 1)
  selected_cols <- c("patientunitstayid", "labresultoffset", "labname", "labresult")
  chunk_size <- 500000
  chunk_id <- 0L
  lab_chunks <- list()
  keep_ids <- data.table(patientunitstayid = cohort_ids)
  data.table::setkey(keep_ids, patientunitstayid)

  repeat {
    lines <- readLines(con, n = chunk_size)
    if (length(lines) == 0) {
      break
    }
    chunk_id <- chunk_id + 1L
    dt <- data.table::fread(
      text = paste(c(header, lines), collapse = "\n"),
      select = selected_cols,
      showProgress = FALSE,
      na.strings = c("", "NA", "NaN")
    )
    dt <- dt[
      labname %chin% target_labnames &
        labresultoffset >= -lab_window_before_minutes &
        labresultoffset <= lab_window_after_minutes
    ]
    if (nrow(dt) > 0) {
      dt <- keep_ids[dt, on = "patientunitstayid", nomatch = 0]
      if (nrow(dt) > 0) {
        dt[, labresult := suppressWarnings(as.numeric(labresult))]
        dt <- dt[is.finite(labresult)]
        lab_chunks[[length(lab_chunks) + 1L]] <- dt
      }
    }
    if (chunk_id %% 10 == 0) {
      message("  processed lab chunks: ", chunk_id)
    }
  }

  if (length(lab_chunks) == 0) {
    stop("No target laboratory values found in the early ICU window.", call. = FALSE)
  }

  labs <- data.table::rbindlist(lab_chunks, use.names = TRUE, fill = TRUE)
  labs[, feature := unname(feature_map[labname])]

  # Plausibility filters to reduce obvious data-entry artefacts before K-means.
  labs <- labs[
    (feature %in% c("polys_percent", "lymphs_percent") & labresult > 0 & labresult <= 100) |
      (feature == "platelets" & labresult > 0 & labresult <= 2000) |
      (feature == "haemoglobin" & labresult >= 3 & labresult <= 25) |
      (feature == "total_protein" & labresult >= 1 & labresult <= 15) |
      (feature == "albumin" & labresult >= 0.5 & labresult <= 8) |
      (feature == "creatinine" & labresult >= 0.1 & labresult <= 25)
  ]

  labs_baseline <- labs[
    ,
    .(value = stats::median(labresult, na.rm = TRUE)),
    by = .(patientunitstayid, feature)
  ]
  labs_wide <- data.table::dcast(
    labs_baseline,
    patientunitstayid ~ feature,
    value.var = "value"
  )
  saveRDS(labs_wide, lab_cache)
  message("Saved early laboratory matrix cache: ", lab_cache)
}

analysis <- merge(older, labs_wide, by = "patientunitstayid", all.x = TRUE)
analysis[, nlr := polys_percent / lymphs_percent]
analysis[, sii_like := platelets * polys_percent / lymphs_percent]

analysis[, protein_proxy := fifelse(!is.na(total_protein), total_protein, albumin)]

write_csv_utf8(
  analysis,
  "eICU_first_ICU_feature_availability_dataset.csv"
)

availability <- data.table(
  feature = c("nlr", "sii_like", "haemoglobin", "total_protein", "albumin", "protein_proxy", "bmi", "creatinine"),
  n_nonmissing = c(
    sum(is.finite(analysis$nlr)),
    sum(is.finite(analysis$sii_like)),
    sum(is.finite(analysis$haemoglobin)),
    sum(is.finite(analysis$total_protein)),
    sum(is.finite(analysis$albumin)),
    sum(is.finite(analysis$protein_proxy)),
    sum(is.finite(analysis$bmi)),
    sum(is.finite(analysis$creatinine))
  )
)
availability[, percent_nonmissing := round(100 * n_nonmissing / nrow(analysis), 2)]
write_csv_utf8(availability, "Table6S_feature_availability_eICU.csv")

feature_vars <- c("nlr", "sii_like", "haemoglobin", "total_protein", "bmi", "creatinine")
required_cols <- c(feature_vars, "age_num", "gender_model", "hospitalid", "hospital_mortality", "survival_days_hosp")
analysis_data <- analysis[stats::complete.cases(analysis[, ..required_cols])]
analysis_data <- analysis_data[
  is.finite(nlr) & is.finite(sii_like) &
    nlr > 0 & nlr <= 100 &
    sii_like > 0 & sii_like <= 200000
]

if (nrow(analysis_data) < 1000) {
  stop("Too few complete strict total-protein cases for eICU validation.", call. = FALSE)
}

message("Strict total-protein complete eICU analytical cohort: ", nrow(analysis_data))

set.seed(20260701)
z <- scale(analysis_data[, ..feature_vars])
kmeans_fit <- stats::kmeans(z, centers = 3, nstart = 100, iter.max = 100, algorithm = "Lloyd")
analysis_data[, cluster_raw := kmeans_fit$cluster]

raw_profiles <- analysis_data[
  ,
  .(
    n = .N,
    hospital_deaths = sum(hospital_mortality == 1, na.rm = TRUE),
    hospital_mortality_percent = round(100 * mean(hospital_mortality == 1, na.rm = TRUE), 2),
    nlr_median = median(nlr, na.rm = TRUE),
    sii_like_median = median(sii_like, na.rm = TRUE),
    haemoglobin_median = median(haemoglobin, na.rm = TRUE),
    total_protein_median = median(total_protein, na.rm = TRUE),
    bmi_median = median(bmi, na.rm = TRUE),
    creatinine_median = median(creatinine, na.rm = TRUE)
  ),
  by = cluster_raw
]

inflammation_score <- as.numeric(scale(raw_profiles$nlr_median)) +
  as.numeric(scale(raw_profiles$sii_like_median)) +
  as.numeric(scale(raw_profiles$creatinine_median))
phenotype1_raw <- raw_profiles$cluster_raw[which.max(inflammation_score)]

remaining_profiles <- raw_profiles[cluster_raw != phenotype1_raw]
nutrition_depletion_score <- -as.numeric(scale(remaining_profiles$haemoglobin_median)) -
  as.numeric(scale(remaining_profiles$total_protein_median)) -
  as.numeric(scale(remaining_profiles$bmi_median))
phenotype2_raw <- remaining_profiles$cluster_raw[which.max(nutrition_depletion_score)]
phenotype3_raw <- setdiff(raw_profiles$cluster_raw, c(phenotype1_raw, phenotype2_raw))

analysis_data[
  ,
  eicu_phenotype_raw := fifelse(
    cluster_raw == phenotype1_raw, "1",
    fifelse(cluster_raw == phenotype2_raw, "2",
            fifelse(cluster_raw == phenotype3_raw, "3", NA_character_))
  )
]
analysis_data[, eicu_phenotype := relevel(factor(eicu_phenotype_raw, levels = c("1", "2", "3")), ref = "3")]
analysis_data[, plot_phenotype := factor(
  eicu_phenotype_raw,
  levels = c("1", "2", "3"),
  labels = c("Phenotype 1", "Phenotype 2", "Phenotype 3 (Ref)")
)]

write_csv_utf8(raw_profiles, "Table6_raw_eICU_denovo_cluster_profiles.csv")
write_csv_utf8(analysis_data, "eICU_denovo_strict_total_protein_dataset.csv")

counts <- analysis_data[
  ,
  .(
    n = .N,
    percent = round(100 * .N / nrow(analysis_data), 2),
    hospital_deaths = sum(hospital_mortality == 1, na.rm = TRUE),
    hospital_mortality_percent = round(100 * mean(hospital_mortality == 1, na.rm = TRUE), 2),
    icu_deaths = sum(icu_mortality == 1, na.rm = TRUE),
    icu_mortality_percent = round(100 * mean(icu_mortality == 1, na.rm = TRUE), 2)
  ),
  by = eicu_phenotype_raw
][order(eicu_phenotype_raw)]
write_csv_utf8(counts, "Table6A_eICU_denovo_counts.csv")

profiles <- analysis_data[
  ,
  .(
    n = .N,
    nlr = iqr_text(nlr),
    sii_like = iqr_text(sii_like),
    haemoglobin = iqr_text(haemoglobin),
    total_protein = iqr_text(total_protein),
    bmi = iqr_text(bmi),
    creatinine = iqr_text(creatinine)
  ),
  by = eicu_phenotype_raw
][order(eicu_phenotype_raw)]
write_csv_utf8(profiles, "Table6B_eICU_denovo_feature_profiles.csv")

cox_fit <- survival::coxph(
  survival::Surv(survival_days_hosp, hospital_mortality) ~
    eicu_phenotype + age_num + gender_model + cluster(hospitalid),
  data = analysis_data
)
cox_summary <- summary(cox_fit)
ci <- as.data.frame(cox_summary$conf.int)
co <- as.data.frame(cox_summary$coefficients)
pretty <- rownames(ci)
pretty <- sub("^eicu_phenotype1$", "Phenotype 1 vs Phenotype 3", pretty)
pretty <- sub("^eicu_phenotype2$", "Phenotype 2 vs Phenotype 3", pretty)
pretty <- sub("^age_num$", "Age, per year", pretty)
pretty <- sub("^gender_modelMale$", "Male vs female", pretty)

p_col <- grep("Pr\\(", names(co), value = TRUE)[1]
cox_table <- data.table(
  variable = pretty,
  hazard_ratio_95ci = sprintf("%.2f (%.2f-%.2f)", ci[, "exp(coef)"], ci[, "lower .95"], ci[, "upper .95"]),
  p_value = format_p(co[, p_col])
)
write_csv_utf8(cox_table, "Table6C_eICU_Cox_hospital_survival.csv")

ph <- survival::cox.zph(cox_fit)
ph_table <- data.table(
  variable = rownames(ph$table),
  chisq = ph$table[, "chisq"],
  df = ph$table[, "df"],
  p = ph$table[, "p"]
)
write_csv_utf8(ph_table, "Cox_PH_check_eICU_denovo.csv")

cluster_vcov_glm <- function(fit, cluster) {
  x <- stats::model.matrix(fit)
  y <- stats::model.response(stats::model.frame(fit))
  mu <- stats::fitted(fit)
  w <- as.numeric(mu * (1 - mu))
  bread <- solve(crossprod(x, x * w))
  score_resid <- as.numeric(y - mu)
  cluster <- as.factor(cluster)
  meat <- matrix(0, nrow = ncol(x), ncol = ncol(x))
  for (lev in levels(cluster)) {
    idx <- which(cluster == lev)
    u <- crossprod(x[idx, , drop = FALSE], score_resid[idx])
    meat <- meat + tcrossprod(u)
  }
  bread %*% meat %*% bread
}

logit_fit <- stats::glm(
  hospital_mortality ~ eicu_phenotype + age_num + gender_model,
  data = analysis_data,
  family = stats::binomial()
)
logit_vcov <- cluster_vcov_glm(logit_fit, analysis_data$hospitalid)
logit_beta <- stats::coef(logit_fit)
logit_se <- sqrt(diag(logit_vcov))
logit_p <- 2 * stats::pnorm(abs(logit_beta / logit_se), lower.tail = FALSE)
logit_keep <- setdiff(names(logit_beta), "(Intercept)")
logit_pretty <- logit_keep
logit_pretty <- sub("^eicu_phenotype1$", "Phenotype 1 vs Phenotype 3", logit_pretty)
logit_pretty <- sub("^eicu_phenotype2$", "Phenotype 2 vs Phenotype 3", logit_pretty)
logit_pretty <- sub("^age_num$", "Age, per year", logit_pretty)
logit_pretty <- sub("^gender_modelMale$", "Male vs female", logit_pretty)
logit_table <- data.table(
  variable = logit_pretty,
  odds_ratio_95ci = sprintf(
    "%.2f (%.2f-%.2f)",
    exp(logit_beta[logit_keep]),
    exp(logit_beta[logit_keep] - 1.96 * logit_se[logit_keep]),
    exp(logit_beta[logit_keep] + 1.96 * logit_se[logit_keep])
  ),
  p_value = format_p(logit_p[logit_keep])
)
write_csv_utf8(logit_table, "Table6D_eICU_logistic_hospital_mortality.csv")

logrank <- survival::survdiff(
  survival::Surv(survival_days_hosp, hospital_mortality) ~ plot_phenotype,
  data = analysis_data
)
logrank_p <- stats::pchisq(logrank$chisq, df = length(logrank$n) - 1, lower.tail = FALSE)

km_fit <- survival::survfit(
  survival::Surv(survival_days_hosp, hospital_mortality) ~ plot_phenotype,
  data = analysis_data
)
km_summary <- summary(km_fit)
km_df <- data.table(
  time = km_summary$time,
  survival = km_summary$surv,
  lower = km_summary$lower,
  upper = km_summary$upper,
  plot_phenotype = sub("^plot_phenotype=", "", km_summary$strata)
)
km_df <- rbind(
  data.table(
    time = 0,
    survival = 1,
    lower = 1,
    upper = 1,
    plot_phenotype = levels(analysis_data$plot_phenotype)
  ),
  km_df,
  fill = TRUE
)

palette <- c(
  "Phenotype 1" = "#B2182B",
  "Phenotype 2" = "#EF8A62",
  "Phenotype 3 (Ref)" = "#2166AC"
)
km_plot <- ggplot(km_df, aes(x = time, y = survival, colour = plot_phenotype)) +
  geom_step(linewidth = 0.95) +
  scale_colour_manual(values = palette, name = "eICU phenotype") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  coord_cartesian(
    xlim = c(0, stats::quantile(analysis_data$survival_days_hosp, 0.99, na.rm = TRUE)),
    ylim = c(0.75, 1.00)
  ) +
  labs(
    x = "Days from ICU admission to hospital discharge",
    y = "Hospital survival probability",
    title = "eICU multicentre de novo phenotype reproducibility",
    subtitle = paste0("Log-rank P ", format_p(logrank_p))
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "top"
  )

if (tolower(Sys.getenv("WRITE_FIGURES", unset = "false")) %in% c("true", "1", "yes", "y")) {
  ggsave(
    filename = file.path(output_dir, "Figure6_eICU_KM_hospital_survival.png"),
    plot = km_plot,
    width = 7.2,
    height = 5.2,
    dpi = 320
  )
}

summary_lines <- c(
  "eICU-CRD 2.0 multicentre validation summary",
  paste0("Output directory: ", output_dir),
  paste0("Older first ICU patient-level cohort: ", nrow(older)),
  paste0("Strict total-protein complete analytical cohort: ", nrow(analysis_data)),
  "",
  "Phenotype counts and hospital mortality:",
  paste(capture.output(print(counts)), collapse = "\n"),
  "",
  "Cox model with hospital-cluster robust standard errors:",
  paste(capture.output(print(cox_table)), collapse = "\n"),
  "",
  "Logistic model for in-hospital mortality with hospital-cluster robust standard errors:",
  paste(capture.output(print(logit_table)), collapse = "\n"),
  "",
  paste0("Log-rank P: ", format_p(logrank_p)),
  "",
  "Interpretation: eICU supports multicentre ICU transportability using in-hospital outcomes, not long-term mortality."
)
writeLines(summary_lines, file.path(output_dir, "eICU_validation_run_summary.txt"), useBytes = TRUE)

message("\nFinished eICU multicentre validation.")
message("Analytical cohort: ", nrow(analysis_data))
message("Key outputs written to: ", output_dir)
print(counts)
print(cox_table)
print(logit_table)
