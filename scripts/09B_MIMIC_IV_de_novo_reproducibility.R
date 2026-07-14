# ==============================================================================
# Script: 09B_MIMIC_IV_de_novo_reproducibility.R
# Purpose:
#   Reproduce the three nutrition-inflammation phenotypes de novo in the
#   MIMIC-IV v3.1 older ICU cohort after the strict NHANES centroid projection
#   showed an extremely small healthy-reference phenotype.
#
# Input/output directory: output/mimic_24h_validation under PROJECT_ROOT.
#
# Interpretation:
#   This is an external reproducibility / transportability analysis, not a
#   strict one-to-one centroid projection.
# ==============================================================================

options(stringsAsFactors = FALSE)

project_root <- normalizePath(
  Sys.getenv("PROJECT_ROOT", unset = getwd()), winslash = "/", mustWork = TRUE
)
output_dir <- Sys.getenv(
  "MIMIC_VALIDATION_OUTPUT",
  unset = file.path(project_root, "output", "mimic_24h_validation")
)

input_path <- file.path(output_dir, "MIMIC_projected_albumin_proxy_dataset.csv")
if (!file.exists(input_path)) {
  stop(
    "Input file not found: ", input_path, "\n",
    "Run 09_MIMIC_IV_external_validation_formal.R first.",
    call. = FALSE
  )
}

required_packages <- c("dplyr", "ggplot2", "survival")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  install.packages(missing_packages, dependencies = TRUE, repos = "https://cloud.r-project.org")
}

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(survival)
})

write_csv_utf8 <- function(x, filename) {
  utils::write.csv(
    x,
    file = file.path(output_dir, filename),
    row.names = FALSE,
    fileEncoding = "UTF-8"
  )
}

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

message("Loading MIMIC albumin-proxy validation matrix...")
mimic <- read.csv(input_path, stringsAsFactors = FALSE)

feature_vars <- c("nlr", "sii", "haemoglobin", "protein_proxy", "bmi", "creatinine")
required_cols <- c(feature_vars, "survival_days_365", "mortality_365d", "anchor_age", "gender")
missing_cols <- setdiff(required_cols, names(mimic))
if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
}

analysis_data <- mimic |>
  dplyr::filter(stats::complete.cases(dplyr::across(dplyr::all_of(required_cols)))) |>
  dplyr::mutate(gender = factor(.data$gender))

if (nrow(analysis_data) < 500) {
  stop("Too few complete cases for de novo reproducibility analysis.", call. = FALSE)
}

set.seed(20260630)
z <- scale(analysis_data[feature_vars])
kmeans_fit <- stats::kmeans(z, centers = 3, nstart = 50, iter.max = 100)
analysis_data$cluster_raw <- kmeans_fit$cluster

raw_profiles <- analysis_data |>
  dplyr::group_by(.data$cluster_raw) |>
  dplyr::summarise(
    n = dplyr::n(),
    deaths_365d = sum(.data$mortality_365d == 1, na.rm = TRUE),
    mortality_365d_percent = round(100 * mean(.data$mortality_365d == 1, na.rm = TRUE), 2),
    nlr_median = median(.data$nlr, na.rm = TRUE),
    sii_median = median(.data$sii, na.rm = TRUE),
    haemoglobin_median = median(.data$haemoglobin, na.rm = TRUE),
    protein_proxy_median = median(.data$protein_proxy, na.rm = TRUE),
    bmi_median = median(.data$bmi, na.rm = TRUE),
    creatinine_median = median(.data$creatinine, na.rm = TRUE),
    .groups = "drop"
  )

# Clinical labelling without using outcomes:
# Phenotype 1 = highest inflammation/renal-stress profile.
# Phenotype 2 = among the remaining clusters, poorer nutrition profile.
# Phenotype 3 = remaining relatively preserved nutrition/resilience profile.
inflammation_score <- as.numeric(scale(raw_profiles$nlr_median)) +
  as.numeric(scale(raw_profiles$sii_median)) +
  as.numeric(scale(raw_profiles$creatinine_median))
phenotype1_raw <- raw_profiles$cluster_raw[which.max(inflammation_score)]

remaining_profiles <- raw_profiles |>
  dplyr::filter(.data$cluster_raw != phenotype1_raw)
nutrition_depletion_score <- -as.numeric(scale(remaining_profiles$haemoglobin_median)) -
  as.numeric(scale(remaining_profiles$protein_proxy_median)) -
  as.numeric(scale(remaining_profiles$bmi_median))
phenotype2_raw <- remaining_profiles$cluster_raw[which.max(nutrition_depletion_score)]
phenotype3_raw <- setdiff(raw_profiles$cluster_raw, c(phenotype1_raw, phenotype2_raw))

analysis_data <- analysis_data |>
  dplyr::mutate(
    mimic_phenotype_raw = dplyr::case_when(
      .data$cluster_raw == phenotype1_raw ~ "1",
      .data$cluster_raw == phenotype2_raw ~ "2",
      .data$cluster_raw == phenotype3_raw ~ "3",
      TRUE ~ NA_character_
    ),
    mimic_phenotype = stats::relevel(factor(.data$mimic_phenotype_raw, levels = c("1", "2", "3")), ref = "3"),
    plot_phenotype = factor(
      .data$mimic_phenotype_raw,
      levels = c("1", "2", "3"),
      labels = c("Phenotype 1", "Phenotype 2", "Phenotype 3 (Ref)")
    )
  )

write_csv_utf8(raw_profiles, "Table5_raw_MIMIC_denovo_cluster_profiles.csv")
write_csv_utf8(analysis_data, "MIMIC_denovo_albumin_proxy_dataset.csv")

counts <- analysis_data |>
  dplyr::group_by(.data$mimic_phenotype_raw) |>
  dplyr::summarise(
    n = dplyr::n(),
    percent = round(100 * dplyr::n() / nrow(analysis_data), 2),
    deaths_365d = sum(.data$mortality_365d == 1, na.rm = TRUE),
    mortality_365d_percent = round(100 * mean(.data$mortality_365d == 1, na.rm = TRUE), 2),
    .groups = "drop"
  )
write_csv_utf8(counts, "Table5A_MIMIC_denovo_phenotype_counts.csv")

profiles <- analysis_data |>
  dplyr::group_by(.data$mimic_phenotype_raw) |>
  dplyr::summarise(
    n = dplyr::n(),
    nlr = iqr_text(.data$nlr),
    sii = iqr_text(.data$sii),
    haemoglobin = iqr_text(.data$haemoglobin),
    protein_proxy = iqr_text(.data$protein_proxy),
    bmi = iqr_text(.data$bmi),
    creatinine = iqr_text(.data$creatinine),
    .groups = "drop"
  )
write_csv_utf8(profiles, "Table5B_MIMIC_denovo_feature_profiles.csv")

cox_fit <- survival::coxph(
  survival::Surv(survival_days_365, mortality_365d) ~ mimic_phenotype + anchor_age + gender,
  data = analysis_data
)
cox_summary <- summary(cox_fit)
ci <- as.data.frame(cox_summary$conf.int)
co <- as.data.frame(cox_summary$coefficients)
pretty <- rownames(ci)
pretty <- sub("^mimic_phenotype1$", "Phenotype 1 vs Phenotype 3", pretty)
pretty <- sub("^mimic_phenotype2$", "Phenotype 2 vs Phenotype 3", pretty)
pretty <- sub("^anchor_age$", "Age, per year", pretty)
pretty <- sub("^genderM$", "Male vs female", pretty)

cox_table <- data.frame(
  variable = pretty,
  hazard_ratio_95ci = sprintf("%.2f (%.2f-%.2f)", ci[, "exp(coef)"], ci[, "lower .95"], ci[, "upper .95"]),
  p_value = format_p(co[, "Pr(>|z|)"]),
  row.names = NULL
)
write_csv_utf8(cox_table, "Table5C_MIMIC_denovo_Cox_365d.csv")

ph <- survival::cox.zph(cox_fit)
ph_table <- data.frame(
  variable = rownames(ph$table),
  chisq = ph$table[, "chisq"],
  df = ph$table[, "df"],
  p = ph$table[, "p"],
  row.names = NULL
)
write_csv_utf8(ph_table, "Cox_PH_check_MIMIC_denovo_albumin_proxy.csv")

km_fit <- survival::survfit(
  survival::Surv(survival_days_365, mortality_365d) ~ plot_phenotype,
  data = analysis_data
)
km_summary <- summary(km_fit)
km_df <- data.frame(
  time = km_summary$time,
  survival = km_summary$surv,
  lower = km_summary$lower,
  upper = km_summary$upper,
  plot_phenotype = sub("^plot_phenotype=", "", km_summary$strata)
)
km_df <- dplyr::bind_rows(
  data.frame(
    time = 0,
    survival = 1,
    lower = 1,
    upper = 1,
    plot_phenotype = levels(analysis_data$plot_phenotype)
  ),
  km_df
)
km_df$plot_phenotype <- factor(km_df$plot_phenotype, levels = levels(analysis_data$plot_phenotype))

logrank <- survival::survdiff(
  survival::Surv(survival_days_365, mortality_365d) ~ plot_phenotype,
  data = analysis_data
)
logrank_p <- stats::pchisq(logrank$chisq, df = length(logrank$n) - 1, lower.tail = FALSE)
p_label <- paste0("Log-rank P ", ifelse(logrank_p < 0.001, "<0.001", paste0("= ", sprintf("%.3f", logrank_p))))

km_plot <- ggplot(km_df, aes(x = .data$time, y = .data$survival, color = .data$plot_phenotype)) +
  geom_step(linewidth = 0.9) +
  scale_color_manual(values = c("Phenotype 1" = "#C73E3A", "Phenotype 2" = "#2F6F9F", "Phenotype 3 (Ref)" = "#3E8B4E")) +
  scale_x_continuous(limits = c(0, 365), breaks = c(0, 90, 180, 270, 365)) +
  scale_y_continuous(limits = c(0, 1), labels = function(x) paste0(round(100 * x), "%")) +
  labs(
    title = "MIMIC-IV v3.1 De Novo Reproducibility: Albumin-Proxy ICU Phenotypes",
    x = "Days after ICU admission",
    y = "Survival probability",
    color = NULL,
    caption = p_label
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "bottom",
    plot.caption = element_text(hjust = 0, face = "bold")
  )

if (tolower(Sys.getenv("WRITE_FIGURES", unset = "false")) %in% c("true", "1", "yes", "y")) {
  ggplot2::ggsave(
    filename = file.path(output_dir, "Figure5_MIMIC_denovo_KM_albumin_proxy.png"),
    plot = km_plot,
    width = 8,
    height = 6,
    dpi = 300
  )
}

method_note <- c(
  "# MIMIC-IV De Novo Reproducibility Analysis",
  "",
  "- The strict NHANES centroid projection yielded an extremely small Phenotype 3 reference group in MIMIC-IV, reflecting major case-mix differences between a community cohort and an ICU cohort.",
  "- Therefore, this analysis performs de novo K-means clustering within MIMIC-IV complete albumin-proxy cases using the same conceptual feature axes: NLR, SII, haemoglobin, protein proxy, BMI, and creatinine.",
  "- Cluster labels were assigned using clinical profiles only, not outcomes.",
  "- Phenotype 1 denotes high inflammation/renal-stress burden.",
  "- Phenotype 2 denotes poorer nutrition among the remaining patients.",
  "- Phenotype 3 denotes relatively preserved nutrition/resilience and is used as the reference group.",
  "- Endpoint: 365-day all-cause mortality after ICU admission.",
  "",
  "## Key Output Files",
  "",
  "- `MIMIC_denovo_albumin_proxy_dataset.csv`",
  "- `Table5A_MIMIC_denovo_phenotype_counts.csv`",
  "- `Table5B_MIMIC_denovo_feature_profiles.csv`",
  "- `Table5C_MIMIC_denovo_Cox_365d.csv`",
  "- `Cox_PH_check_MIMIC_denovo_albumin_proxy.csv`",
  "- `Figure5_MIMIC_denovo_KM_albumin_proxy.png`"
)
writeLines(method_note, con = file.path(output_dir, "MIMIC_denovo_method_note.md"), useBytes = TRUE)

message("Done. MIMIC de novo reproducibility outputs saved to: ", output_dir)
message("Complete cases: ", nrow(analysis_data))
message("Cox table: ", file.path(output_dir, "Table5C_MIMIC_denovo_Cox_365d.csv"))
