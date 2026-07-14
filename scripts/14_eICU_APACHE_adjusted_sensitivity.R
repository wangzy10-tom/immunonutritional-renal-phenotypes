# ==============================================================================
# Script: 14_eICU_APACHE_adjusted_sensitivity.R
# Purpose: Add APACHE IVa severity adjustment to the eICU multicentre validation.
# ==============================================================================

required_pkgs <- c("readr", "dplyr", "tibble")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "))
}

library(readr)
library(dplyr)
library(tibble)

project_root <- normalizePath(
  Sys.getenv("PROJECT_ROOT", unset = getwd()), winslash = "/", mustWork = TRUE
)
eicu_dir <- Sys.getenv(
  "EICU_DIR", unset = file.path(project_root, "data", "eicu-crd-2.0")
)
input_dir <- Sys.getenv(
  "EICU_24H_INPUT", unset = file.path(project_root, "output", "eicu_24h_extraction")
)
root_out <- Sys.getenv(
  "EICU_APACHE_OUTPUT", unset = file.path(project_root, "output", "eicu_apache_support")
)
dir.create(root_out, recursive = TRUE, showWarnings = FALSE)

cluster_robust_glm_table <- function(model, cluster) {
  X <- model.matrix(model)
  y <- model$y
  mu <- fitted(model)
  keep <- complete.cases(X, y, mu, cluster)
  X <- X[keep, , drop = FALSE]
  y <- y[keep]
  mu <- mu[keep]
  cluster <- cluster[keep]

  score_i <- X * as.numeric(y - mu)
  score_g <- rowsum(score_i, group = cluster, reorder = FALSE)

  W <- as.numeric(mu * (1 - mu))
  bread <- solve(t(X) %*% (X * W))
  meat <- t(score_g) %*% score_g

  n <- nrow(X)
  p <- ncol(X)
  g <- length(unique(cluster))
  correction <- (g / (g - 1)) * ((n - 1) / (n - p))
  vcov_cr <- correction * bread %*% meat %*% bread

  beta <- coef(model)
  se <- sqrt(diag(vcov_cr))
  z <- beta / se
  pval <- 2 * pnorm(abs(z), lower.tail = FALSE)
  ci_low <- beta - 1.96 * se
  ci_high <- beta + 1.96 * se

  tibble(
    term = names(beta),
    OR = exp(beta),
    lower_95 = exp(ci_low),
    upper_95 = exp(ci_high),
    p_value = pval,
    odds_ratio_95ci = sprintf("%.2f (%.2f-%.2f)", OR, lower_95, upper_95),
    p_value_formatted = ifelse(p_value < 0.001, "<0.001", sprintf("%.3f", p_value))
  )
}

cat("Loading eICU strict total-protein phenotype dataset...\n")
eicu <- read_csv(
  file.path(input_dir, "eICU_denovo_strict_total_protein_dataset.csv"),
  show_col_types = FALSE
)

cat("Loading APACHE patient result table...\n")
apache <- read_csv(
  file.path(eicu_dir, "apachePatientResult.csv.gz"),
  show_col_types = FALSE
) %>%
  filter(apacheversion == "IVa") %>%
  transmute(
    patientunitstayid,
    acutephysiologyscore = as.numeric(acutephysiologyscore),
    apachescore = as.numeric(apachescore),
    predictedicumortality = as.numeric(predictedicumortality),
    predictedhospitalmortality = as.numeric(predictedhospitalmortality)
  ) %>%
  mutate(
    acutephysiologyscore = ifelse(acutephysiologyscore < 0, NA_real_, acutephysiologyscore),
    apachescore = ifelse(apachescore < 0, NA_real_, apachescore),
    predictedicumortality = ifelse(predictedicumortality < 0, NA_real_, predictedicumortality),
    predictedhospitalmortality = ifelse(predictedhospitalmortality < 0, NA_real_, predictedhospitalmortality)
  )

analysis <- eicu %>%
  left_join(apache, by = "patientunitstayid") %>%
  mutate(
    phenotype = factor(eicu_phenotype_raw, levels = c(3, 1, 2),
                       labels = c("Phenotype 3", "Phenotype 1", "Phenotype 2")),
    male = ifelse(gender_model == "Male", 1, 0),
    hospital_mortality = as.integer(hospital_mortality),
    icu_mortality = as.integer(icu_mortality)
  )

availability <- analysis %>%
  summarise(
    n = n(),
    apache_available = sum(!is.na(apachescore)),
    apache_available_percent = round(100 * mean(!is.na(apachescore)), 2),
    acute_physiology_available = sum(!is.na(acutephysiologyscore)),
    predicted_hospital_mortality_available = sum(!is.na(predictedhospitalmortality))
  )

model_data <- analysis %>%
  filter(
    !is.na(phenotype),
    !is.na(age_num),
    !is.na(male),
    !is.na(hospitalid),
    !is.na(apachescore),
    !is.na(hospital_mortality),
    !is.na(icu_mortality)
  )

cat("Fitting APACHE-adjusted hospital mortality model...\n")
fit_hosp <- glm(
  hospital_mortality ~ phenotype + age_num + male + apachescore,
  family = binomial(),
  data = model_data,
  y = TRUE
)

hosp_table <- cluster_robust_glm_table(fit_hosp, model_data$hospitalid) %>%
  mutate(model = "Hospital mortality, APACHE-adjusted")

cat("Fitting APACHE-adjusted ICU mortality model...\n")
fit_icu <- glm(
  icu_mortality ~ phenotype + age_num + male + apachescore,
  family = binomial(),
  data = model_data,
  y = TRUE
)

icu_table <- cluster_robust_glm_table(fit_icu, model_data$hospitalid) %>%
  mutate(model = "ICU mortality, APACHE-adjusted")

main_terms <- c("phenotypePhenotype 1", "phenotypePhenotype 2", "age_num", "male", "apachescore")
combined_table <- bind_rows(hosp_table, icu_table) %>%
  filter(term %in% main_terms) %>%
  mutate(
    variable = recode(
      term,
      "phenotypePhenotype 1" = "Phenotype 1 vs Phenotype 3",
      "phenotypePhenotype 2" = "Phenotype 2 vs Phenotype 3",
      "age_num" = "Age, per year",
      "male" = "Male vs female",
      "apachescore" = "APACHE IVa score, per point"
    )
  ) %>%
  select(model, variable, odds_ratio_95ci, p_value_formatted, OR, lower_95, upper_95, p_value)

write_csv(availability, file.path(root_out, "Table6M_eICU_APACHE_availability.csv"))
write_csv(combined_table, file.path(root_out, "Table6N_eICU_APACHE_adjusted_logistic_models.csv"))
saveRDS(
  list(
    availability = availability,
    model_data = model_data,
    fit_hosp = fit_hosp,
    fit_icu = fit_icu,
    combined_table = combined_table
  ),
  file.path(root_out, "eICU_APACHE_adjusted_sensitivity.rds")
)

summary_lines <- c(
  "eICU APACHE-adjusted sensitivity analysis",
  paste0("Strict total-protein cohort: ", nrow(analysis)),
  paste0("APACHE IVa score available: ", availability$apache_available,
         " (", availability$apache_available_percent, "%)"),
  paste0("Complete APACHE-adjusted model cohort: ", nrow(model_data)),
  "",
  "APACHE-adjusted models:",
  paste(capture.output(print(combined_table %>% select(model, variable, odds_ratio_95ci, p_value_formatted))), collapse = "\n")
)

writeLines(summary_lines, file.path(root_out, "eICU_APACHE_adjusted_sensitivity_summary.txt"))

cat("\nCompleted. Outputs saved to: ", root_out, "\n", sep = "")
cat(paste(summary_lines, collapse = "\n"))
