# Verify aggregate rerun outputs against the frozen result anchors.
root <- normalizePath(Sys.getenv("PROJECT_ROOT", unset = getwd()), winslash = "/", mustWork = TRUE)
expected <- read.csv(file.path(root, "aggregate_results", "EXPECTED_KEY_RESULTS.csv"), check.names = FALSE)
observed <- list()
add <- function(id, estimate, lower, upper, n, events) observed[[length(observed)+1L]] <<- data.frame(result_id=id, estimate=estimate, lower_95=lower, upper_95=upper, n=n, events=events)
nh <- read.csv(file.path(root, "output", "nhanes_albumin_amendment17_2026-08-30", "Table81F_primary_survey_Cox.csv"))
for(i in seq_len(nrow(nh))) add(paste0("NHANES_PRIMARY_", sub(" vs P3", "", nh$comparison[i])), nh$HR[i], nh$lower_95[i], nh$upper_95[i], nh$n[i], nh$events[i])
mi <- read.csv(file.path(root, "output", "mimic_albumin_primary_2026-08-30", "Table80C_MIMIC_albumin_OASIS_models.csv"))
mi <- mi[mi$model %in% c("Official OASIS quartile-stratified Cox", "Selection-IPW official OASIS quartile-stratified Cox"),]
for(i in seq_len(nrow(mi))) { prefix <- if(mi$model[i] == "Official OASIS quartile-stratified Cox") "MIMIC_PRIMARY_" else "MIMIC_IPW_"; add(paste0(prefix, sub(" vs P3", "", mi$comparison[i])), mi$estimate[i], mi$lower_95[i], mi$upper_95[i], mi$n[i], mi$events[i]) }
ei <- read.csv(file.path(root, "output", "eicu_albumin_amendment17_2026-08-30", "Table82C_outcome_models.csv"))
ei <- ei[ei$model %in% c("Primary APACHE-adjusted hospital random-intercept", "Selection-IPW APACHE-adjusted hospital random-intercept"),]
for(i in seq_len(nrow(ei))) { prefix <- if(ei$model[i] == "Primary APACHE-adjusted hospital random-intercept") "EICU_PRIMARY_" else "EICU_IPW_"; add(paste0(prefix, sub(" vs P3", "", ei$comparison[i])), ei$estimate[i], ei$lower_95[i], ei$upper_95[i], ei$n[i], ei$events[i]) }
variant_code <- c("Five-feature: NLR only"="NLR5", "Five-feature: SII only"="SII5", "CLARA K-medoids"="CLARA", "Gaussian mixture"="GMM")
dataset_code <- c("NHANES"="NHANES", "MIMIC-IV"="MIMIC", "eICU"="EICU")
st <- read.csv(file.path(root, "output", "cross_database_structural_sensitivity_2026-08-30", "Table83D_outcome_models.csv"))
st <- st[st$comparison == "P1 vs P3",]
for(i in seq_len(nrow(st))) add(paste0("STRUCT_", dataset_code[st$dataset[i]], "_", variant_code[st$variant[i]], "_P1"), st$estimate[i], st$lower_95[i], st$upper_95[i], st$n[i], st$events[i])
nh_counts <- read.csv(file.path(root, "output", "nhanes_albumin_amendment17_2026-08-30", "Table81A_full_phenotype_counts.csv"))
add("NHANES_FULL_N", sum(nh_counts$n), NA_real_, NA_real_, sum(nh_counts$n), sum(nh_counts$deaths))
nh_ipw <- read.csv(file.path(root, "output", "nhanes_albumin_amendment17_2026-08-30", "selection_bias", "Table33D_selection_IPW_Cox.csv"))
nh_ipw <- nh_ipw[grepl("selection IPW", nh_ipw$model, fixed=TRUE), ]
# The IPW table does not export N/events; both fits use the same complete-case
# model cohort (confirmed in the frozen selection-IPW module).
nh_model_counts <- read.csv(file.path(root, "output", "nhanes_albumin_amendment17_2026-08-30", "Table81B_primary_model_counts.csv"))
for(i in seq_len(nrow(nh_ipw))) add(paste0("NHANES_IPW_", sub(" vs P3", "", nh_ipw$comparison[i])), nh_ipw$HR[i], nh_ipw$lower_95[i], nh_ipw$upper_95[i], sum(nh_model_counts$n), sum(nh_model_counts$deaths))
mi_counts <- read.csv(file.path(root, "output", "mimic_albumin_primary_2026-08-30", "Table80A_MIMIC_albumin_counts.csv"))
add("MIMIC_FULL_N", sum(mi_counts$n), NA_real_, NA_real_, sum(mi_counts$n), sum(mi_counts$deaths_365d))
mi_sofa <- read.csv(file.path(root, "output", "mimic_albumin_primary_2026-08-30", "Table80D_MIMIC_albumin_SOFA_models.csv"))
mi_sofa <- mi_sofa[mi_sofa$model == "Selection-IPW official first-day SOFA quartile-stratified Cox", ]
for(i in seq_len(nrow(mi_sofa))) add(paste0("MIMIC_SOFA_IPW_", sub(" vs P3", "", mi_sofa$comparison[i])), mi_sofa$estimate[i], mi_sofa$lower_95[i], mi_sofa$upper_95[i], mi_sofa$n[i], mi_sofa$events[i])
ei_counts <- read.csv(file.path(root, "output", "eicu_albumin_amendment17_2026-08-30", "Table82B_phenotype_counts.csv"))
add("EICU_FEATURE_COMPLETE_N", sum(ei_counts$n), NA_real_, NA_real_, sum(ei_counts$n), NA_real_)
ei_meta <- read.csv(file.path(root, "output", "eicu_albumin_amendment17_2026-08-30", "Table82J_REML_HK_meta_analysis.csv"))
stopifnot(nrow(ei_meta) == 1L)
add("EICU_META_P1", ei_meta$pooled_or, ei_meta$lower_95, ei_meta$upper_95, ei_meta$participants_p1_p3, ei_meta$deaths_p1_p3)
obs <- do.call(rbind, observed)
if(anyDuplicated(expected$result_id) || anyDuplicated(obs$result_id) || !setequal(expected$result_id, obs$result_id)) stop("Missing, extra, or duplicated result identifiers.", call.=FALSE)
check <- merge(expected, obs, by="result_id", suffixes=c("_expected", "_observed"))
num_cols <- c("estimate", "lower_95", "upper_95", "n", "events")
check$passed <- TRUE
for(v in num_cols) check$passed <- check$passed & (is.na(check[[paste0(v,"_expected")]]) | abs(check[[paste0(v,"_expected")]] - check[[paste0(v,"_observed")]]) < 1e-10)
out_dir <- file.path(root, "output", "release_verification"); dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)
write.csv(check, file.path(out_dir, "KEY_RESULT_VERIFICATION.csv"), row.names=FALSE, na="")
if(nrow(check) != nrow(expected) || !isTRUE(all(check$passed))) stop("Key-result verification failed.", call.=FALSE)
message("Key-result verification passed: ", nrow(check), "/", nrow(check))

