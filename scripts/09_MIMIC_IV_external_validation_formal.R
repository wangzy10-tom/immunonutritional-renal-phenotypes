# ==============================================================================
# Script: 09_MIMIC_IV_external_validation_formal.R
# Purpose:
#   External transportability validation of NHANES-derived phenotypes in
#   MIMIC-IV v3.1 older ICU patients.
#
# Design:
#   1) Reconstruct the NHANES scaling parameters and phenotype centroids from
#      the local Freeze V2 projection-reference object built by script 09A.
#   2) Extract first ICU stays aged >=65 years from local MIMIC-IV v3.1.
#   3) Use audited MIMIC lab itemids and a strict 0h/+24h window after ICU admission.
#   4) Run two analyses:
#      - Strict total-protein validation: exact NHANES feature analogue.
#      - Albumin-proxy validation: larger, biologically adjacent sensitivity
#        analysis, explicitly labelled as non-identical marker validation.
#
# Default output: output/mimic_24h_validation under PROJECT_ROOT.
# ==============================================================================

options(stringsAsFactors = FALSE)

# ---- User-configurable paths --------------------------------------------------
project_root <- normalizePath(
  Sys.getenv("PROJECT_ROOT", unset = getwd()), winslash = "/", mustWork = TRUE
)
mimic_root <- Sys.getenv(
  "MIMIC_ROOT", unset = file.path(project_root, "data", "mimic-iv-3.1")
)
nhanes_rdata_path <- Sys.getenv(
  "NHANES_CLUSTERED_RDATA",
  unset = file.path(
    project_root, "output", "nhanes_projection_reference",
    "NHANES_FreezeV2_projection_reference.RData"
  )
)
output_dir <- Sys.getenv(
  "MIMIC_VALIDATION_OUTPUT",
  unset = file.path(project_root, "output", "mimic_24h_validation")
)

lab_window_before_hours <- as.integer(Sys.getenv("LAB_WINDOW_BEFORE_HOURS", unset = "0"))
lab_window_after_hours <- as.integer(Sys.getenv("LAB_WINDOW_AFTER_HOURS", unset = "24"))
bmi_lookback_days <- as.integer(Sys.getenv("BMI_LOOKBACK_DAYS", unset = "365"))

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Package handling ---------------------------------------------------------
required_packages <- c("DBI", "duckdb", "dplyr", "ggplot2", "survival")

install_if_missing <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) == 0) {
    return(invisible(TRUE))
  }

  auto_install <- tolower(Sys.getenv("AUTO_INSTALL_R_PACKAGES", unset = "true")) %in%
    c("true", "1", "yes", "y")
  if (!auto_install) {
    stop(
      "Missing R packages: ", paste(missing, collapse = ", "),
      "\nSet AUTO_INSTALL_R_PACKAGES=true or install them manually in RStudio.",
      call. = FALSE
    )
  }

  cran_mirror <- Sys.getenv("CRAN_MIRROR", unset = "https://cloud.r-project.org")
  options(repos = c(CRAN = cran_mirror))
  pkg_type <- if (.Platform$OS.type == "windows") "binary" else "source"

  for (pkg in missing) {
    message("Installing missing R package: ", pkg)
    tryCatch(
      install.packages(pkg, dependencies = TRUE, type = pkg_type),
      error = function(e) {
        message("Package installation failed for ", pkg, ": ", conditionMessage(e))
      }
    )
  }

  still_missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(still_missing) > 0) {
    stop(
      "These R packages are still missing after installation attempt: ",
      paste(still_missing, collapse = ", "),
      "\nTry manually in RStudio:\n",
      "install.packages(c('", paste(still_missing, collapse = "','"), "'), dependencies = TRUE)\n",
      "If duckdb fails behind a hospital network, switch CRAN_MIRROR or install the Windows binary.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

install_if_missing(required_packages)
suppressPackageStartupMessages({
  library(DBI)
  library(duckdb)
  library(dplyr)
  library(ggplot2)
  library(survival)
})

# ---- Small helpers ------------------------------------------------------------
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

safe_path <- function(path) {
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

# ---- Reconstruct NHANES reference model --------------------------------------
message("Loading NHANES clustered data and reconstructing reference centroids...")
if (!file.exists(nhanes_rdata_path)) {
  stop("NHANES clustered RData not found: ", nhanes_rdata_path, call. = FALSE)
}

loaded_objects <- load(nhanes_rdata_path)
if (!"nhanes_clustered" %in% loaded_objects && !exists("nhanes_clustered")) {
  stop(
    "The RData file does not contain an object named nhanes_clustered. ",
    "A valid phenotype projection requires the final NHANES clustered dataset.",
    call. = FALSE
  )
}

nhanes_features <- c("NLR", "SII", "LBXHGB", "LBXSTP", "BMXBMI", "LBXSCR")
missing_nhanes_cols <- setdiff(c(nhanes_features, "Phenotype"), names(nhanes_clustered))
if (length(missing_nhanes_cols) > 0) {
  stop(
    "NHANES clustered data is missing columns: ",
    paste(missing_nhanes_cols, collapse = ", "),
    call. = FALSE
  )
}

nhanes_ref <- nhanes_clustered |>
  dplyr::select(dplyr::all_of(c("Phenotype", nhanes_features))) |>
  dplyr::filter(stats::complete.cases(dplyr::across(dplyr::all_of(nhanes_features))))

nhanes_center <- vapply(nhanes_ref[nhanes_features], mean, numeric(1), na.rm = TRUE)
nhanes_scale <- vapply(nhanes_ref[nhanes_features], stats::sd, numeric(1), na.rm = TRUE)
if (any(!is.finite(nhanes_scale) | nhanes_scale <= 0)) {
  stop("Invalid NHANES scaling vector; at least one feature has zero or missing SD.", call. = FALSE)
}

nhanes_z <- as.data.frame(scale(nhanes_ref[nhanes_features], center = nhanes_center, scale = nhanes_scale))
nhanes_z$Phenotype <- as.character(nhanes_ref$Phenotype)
nhanes_centers <- nhanes_z |>
  dplyr::group_by(.data$Phenotype) |>
  dplyr::summarise(dplyr::across(dplyr::all_of(nhanes_features), mean), .groups = "drop") |>
  dplyr::arrange(.data$Phenotype)

write_csv_utf8(
  data.frame(feature = nhanes_features, mean = nhanes_center, sd = nhanes_scale),
  "NHANES_reference_scaling_used.csv"
)
write_csv_utf8(nhanes_centers, "NHANES_reference_centers_used.csv")

center_matrix <- as.matrix(nhanes_centers[nhanes_features])
rownames(center_matrix) <- nhanes_centers$Phenotype

# ---- MIMIC extraction through DuckDB -----------------------------------------
message("Connecting to local MIMIC-IV v3.1 files through DuckDB...")
mimic_files <- list(
  patients = file.path(mimic_root, "hosp", "patients.csv.gz"),
  admissions = file.path(mimic_root, "hosp", "admissions.csv.gz"),
  d_labitems = file.path(mimic_root, "hosp", "d_labitems.csv.gz"),
  labevents = file.path(mimic_root, "hosp", "labevents.csv.gz"),
  omr = file.path(mimic_root, "hosp", "omr.csv.gz"),
  icustays = file.path(mimic_root, "icu", "icustays.csv.gz")
)

missing_files <- mimic_files[!vapply(mimic_files, file.exists, logical(1))]
if (length(missing_files) > 0) {
  stop(
    "Missing MIMIC-IV files:\n",
    paste(names(missing_files), unlist(missing_files), sep = ": ", collapse = "\n"),
    call. = FALSE
  )
}

connect_duckdb <- function() {
  candidate_dbdirs <- c(
    ":memory:",
    file.path(tempdir(), paste0("mimic_validation_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".duckdb"))
  )

  last_error <- NULL
  for (dbdir in candidate_dbdirs) {
    con_try <- try(DBI::dbConnect(duckdb::duckdb(), dbdir = dbdir), silent = TRUE)
    if (inherits(con_try, "try-error")) {
      last_error <- as.character(con_try)
      next
    }
    if (!DBI::dbIsValid(con_try)) {
      last_error <- paste("DuckDB connection object is invalid for dbdir:", dbdir)
      try(DBI::dbDisconnect(con_try, shutdown = TRUE), silent = TRUE)
      next
    }

    smoke <- try(DBI::dbGetQuery(con_try, "SELECT 1 AS duckdb_connection_test"), silent = TRUE)
    if (!inherits(smoke, "try-error")) {
      message("DuckDB connection OK: ", dbdir)
      return(con_try)
    }

    last_error <- as.character(smoke)
    try(DBI::dbDisconnect(con_try, shutdown = TRUE), silent = TRUE)
  }

  stop(
    "DuckDB connection failed after memory and temporary-file attempts.\n",
    "Last error: ", last_error, "\n",
    "Recommended quick fix: close RStudio completely, reopen it, and rerun the script. ",
    "If it still fails, use the Python fallback validation script.",
    call. = FALSE
  )
}

con <- connect_duckdb()
on.exit(try(DBI::dbDisconnect(con, shutdown = TRUE), silent = TRUE), add = TRUE)

csv_sql <- function(path, all_varchar = FALSE) {
  paste0(
    "read_csv_auto('", safe_path(path), "', sample_size=200000, all_varchar=",
    ifelse(all_varchar, "true", "false"),
    ")"
  )
}

DBI::dbExecute(con, paste0("CREATE OR REPLACE VIEW patients AS SELECT * FROM ", csv_sql(mimic_files$patients)))
DBI::dbExecute(con, paste0("CREATE OR REPLACE VIEW admissions AS SELECT * FROM ", csv_sql(mimic_files$admissions)))
DBI::dbExecute(con, paste0("CREATE OR REPLACE VIEW d_labitems AS SELECT * FROM ", csv_sql(mimic_files$d_labitems)))
DBI::dbExecute(con, paste0("CREATE OR REPLACE VIEW labevents AS SELECT * FROM ", csv_sql(mimic_files$labevents)))
DBI::dbExecute(con, paste0("CREATE OR REPLACE VIEW omr AS SELECT * FROM ", csv_sql(mimic_files$omr, all_varchar = TRUE)))
DBI::dbExecute(con, paste0("CREATE OR REPLACE VIEW icustays AS SELECT * FROM ", csv_sql(mimic_files$icustays)))

message("Extracting first-ICU cohort and baseline feature matrix...")
lab_sql <- sprintf(
  "
CREATE OR REPLACE TEMP TABLE mimic_validation_matrix AS
WITH first_icu AS (
  SELECT
    *,
    ROW_NUMBER() OVER (PARTITION BY subject_id ORDER BY CAST(intime AS TIMESTAMP), stay_id) AS rn
  FROM icustays
),
cohort AS (
  SELECT
    f.subject_id,
    f.hadm_id,
    f.stay_id,
    CAST(f.intime AS TIMESTAMP) AS intime,
    CAST(f.outtime AS TIMESTAMP) AS outtime,
    f.first_careunit,
    f.last_careunit,
    p.gender,
    p.anchor_age,
    p.anchor_year_group,
    CAST(p.dod AS TIMESTAMP) AS dod,
    a.race,
    CAST(a.admittime AS TIMESTAMP) AS admittime,
    CAST(a.dischtime AS TIMESTAMP) AS dischtime,
    a.hospital_expire_flag,
    CASE
      WHEN p.dod IS NOT NULL
      THEN DATE_DIFF('day', CAST(f.intime AS TIMESTAMP), CAST(p.dod AS TIMESTAMP))
      ELSE NULL
    END AS death_days_from_icu,
    CASE
      WHEN p.dod IS NOT NULL
       AND DATE_DIFF('day', CAST(f.intime AS TIMESTAMP), CAST(p.dod AS TIMESTAMP)) BETWEEN 0 AND 28
      THEN 1 ELSE 0
    END AS mortality_28d,
    CASE
      WHEN p.dod IS NOT NULL
       AND DATE_DIFF('day', CAST(f.intime AS TIMESTAMP), CAST(p.dod AS TIMESTAMP)) BETWEEN 0 AND 90
      THEN 1 ELSE 0
    END AS mortality_90d,
    CASE
      WHEN p.dod IS NOT NULL
       AND DATE_DIFF('day', CAST(f.intime AS TIMESTAMP), CAST(p.dod AS TIMESTAMP)) BETWEEN 0 AND 365
      THEN 1 ELSE 0
    END AS mortality_365d,
    CASE
      WHEN p.dod IS NOT NULL
       AND DATE_DIFF('day', CAST(f.intime AS TIMESTAMP), CAST(p.dod AS TIMESTAMP)) BETWEEN 0 AND 365
      THEN GREATEST(DATE_DIFF('day', CAST(f.intime AS TIMESTAMP), CAST(p.dod AS TIMESTAMP)), 0.5)
      ELSE 365
    END AS survival_days_365
  FROM first_icu f
  INNER JOIN patients p ON f.subject_id = p.subject_id
  INNER JOIN admissions a ON f.hadm_id = a.hadm_id
  WHERE f.rn = 1
    AND p.anchor_age >= 65
),
selected_labs AS (
  SELECT * FROM (
    VALUES
      (51256, 'neutrophils'),
      (51244, 'lymphocytes'),
      (51265, 'platelets'),
      (51222, 'haemoglobin'),
      (50912, 'creatinine'),
      (50976, 'total_protein'),
      (50862, 'albumin')
  ) AS t(itemid, feature)
),
labs_ranked AS (
  SELECT
    c.stay_id,
    s.feature,
    TRY_CAST(l.valuenum AS DOUBLE) AS valuenum,
    CAST(l.charttime AS TIMESTAMP) AS charttime,
    ABS(DATE_DIFF('minute', CAST(c.intime AS TIMESTAMP), CAST(l.charttime AS TIMESTAMP))) AS abs_minutes_from_intime,
    ROW_NUMBER() OVER (
      PARTITION BY c.stay_id, s.feature
      ORDER BY
        ABS(DATE_DIFF('minute', CAST(c.intime AS TIMESTAMP), CAST(l.charttime AS TIMESTAMP))),
        CAST(l.charttime AS TIMESTAMP)
    ) AS rn
  FROM cohort c
  INNER JOIN labevents l
    ON c.subject_id = l.subject_id
   AND c.hadm_id = l.hadm_id
  INNER JOIN selected_labs s
    ON l.itemid = s.itemid
  WHERE TRY_CAST(l.valuenum AS DOUBLE) IS NOT NULL
    AND CAST(l.charttime AS TIMESTAMP)
      BETWEEN c.intime - INTERVAL '%d hours'
          AND c.intime + INTERVAL '%d hours'
),
labs_wide AS (
  SELECT
    stay_id,
    MAX(CASE WHEN feature = 'neutrophils' AND rn = 1 THEN valuenum END) AS neutrophils,
    MAX(CASE WHEN feature = 'lymphocytes' AND rn = 1 THEN valuenum END) AS lymphocytes,
    MAX(CASE WHEN feature = 'platelets' AND rn = 1 THEN valuenum END) AS platelets,
    MAX(CASE WHEN feature = 'haemoglobin' AND rn = 1 THEN valuenum END) AS haemoglobin,
    MAX(CASE WHEN feature = 'creatinine' AND rn = 1 THEN valuenum END) AS creatinine,
    MAX(CASE WHEN feature = 'total_protein' AND rn = 1 THEN valuenum END) AS total_protein,
    MAX(CASE WHEN feature = 'albumin' AND rn = 1 THEN valuenum END) AS albumin
  FROM labs_ranked
  WHERE rn = 1
  GROUP BY stay_id
),
bmi_ranked AS (
  SELECT
    c.stay_id,
    TRY_CAST(o.result_value AS DOUBLE) AS bmi,
    CAST(o.chartdate AS DATE) AS bmi_chartdate,
    ABS(DATE_DIFF('day', CAST(c.intime AS DATE), CAST(o.chartdate AS DATE))) AS bmi_abs_days_from_intime,
    ROW_NUMBER() OVER (
      PARTITION BY c.stay_id
      ORDER BY ABS(DATE_DIFF('day', CAST(c.intime AS DATE), CAST(o.chartdate AS DATE))), CAST(o.chartdate AS DATE)
    ) AS rn
  FROM cohort c
  INNER JOIN omr o
    ON c.subject_id = TRY_CAST(o.subject_id AS BIGINT)
  WHERE regexp_matches(lower(o.result_name), 'bmi|body mass index')
    AND TRY_CAST(o.result_value AS DOUBLE) BETWEEN 10 AND 80
    AND CAST(o.chartdate AS DATE)
      BETWEEN CAST(c.intime AS DATE) - INTERVAL '%d days'
          AND CAST(c.intime AS DATE) + INTERVAL '1 day'
),
bmi AS (
  SELECT stay_id, bmi, 'direct_bmi' AS bmi_source, bmi_chartdate, bmi_abs_days_from_intime
  FROM bmi_ranked
  WHERE rn = 1
)
SELECT
  c.subject_id,
  c.hadm_id,
  c.stay_id,
  c.intime,
  c.outtime,
  c.first_careunit,
  c.last_careunit,
  c.gender,
  c.anchor_age,
  c.anchor_year_group,
  c.race,
  c.admittime,
  c.dischtime,
  c.hospital_expire_flag,
  c.death_days_from_icu,
  c.mortality_28d,
  c.mortality_90d,
  c.mortality_365d,
  c.survival_days_365,
  l.neutrophils,
  l.lymphocytes,
  l.platelets,
  l.haemoglobin,
  l.creatinine,
  l.total_protein,
  l.albumin,
  b.bmi,
  b.bmi_source,
  b.bmi_chartdate,
  b.bmi_abs_days_from_intime
FROM cohort c
LEFT JOIN labs_wide l ON c.stay_id = l.stay_id
LEFT JOIN bmi b ON c.stay_id = b.stay_id
",
  lab_window_before_hours,
  lab_window_after_hours,
  bmi_lookback_days
)

DBI::dbExecute(con, lab_sql)
mimic_raw <- DBI::dbGetQuery(con, "SELECT * FROM mimic_validation_matrix")

message("MIMIC first-ICU cohort rows: ", nrow(mimic_raw))

mimic_raw <- mimic_raw |>
  dplyr::mutate(
    gender = factor(.data$gender),
    nlr = dplyr::if_else(.data$lymphocytes > 0, .data$neutrophils / .data$lymphocytes, NA_real_),
    sii = dplyr::if_else(.data$lymphocytes > 0, .data$neutrophils * .data$platelets / .data$lymphocytes, NA_real_),
    protein_proxy = dplyr::coalesce(.data$total_protein, .data$albumin),
    protein_proxy_source = dplyr::case_when(
      !is.na(.data$total_protein) ~ "total_protein",
      is.na(.data$total_protein) & !is.na(.data$albumin) ~ "albumin_proxy",
      TRUE ~ NA_character_
    )
  )

write_csv_utf8(
  mimic_raw,
  "MIMIC_first_ICU_feature_availability_dataset.csv"
)

write_csv_utf8(
  data.frame(
    n = nrow(mimic_raw),
    hospital_deaths = sum(mimic_raw$hospital_expire_flag == 1, na.rm = TRUE),
    mortality_28d = sum(mimic_raw$mortality_28d == 1, na.rm = TRUE),
    mortality_90d = sum(mimic_raw$mortality_90d == 1, na.rm = TRUE),
    mortality_365d = sum(mimic_raw$mortality_365d == 1, na.rm = TRUE),
    mortality_365d_percent = round(100 * mean(mimic_raw$mortality_365d == 1, na.rm = TRUE), 2)
  ),
  "MIMIC_first_ICU_cohort_summary.csv"
)

write_csv_utf8(
  data.frame(
    feature = c(
      "nlr",
      "sii",
      "haemoglobin",
      "total_protein",
      "albumin",
      "protein_proxy",
      "bmi",
      "creatinine"
    ),
    n_nonmissing = c(
      sum(!is.na(mimic_raw$nlr)),
      sum(!is.na(mimic_raw$sii)),
      sum(!is.na(mimic_raw$haemoglobin)),
      sum(!is.na(mimic_raw$total_protein)),
      sum(!is.na(mimic_raw$albumin)),
      sum(!is.na(mimic_raw$protein_proxy)),
      sum(!is.na(mimic_raw$bmi)),
      sum(!is.na(mimic_raw$creatinine))
    ),
    percent_nonmissing = round(100 * c(
      mean(!is.na(mimic_raw$nlr)),
      mean(!is.na(mimic_raw$sii)),
      mean(!is.na(mimic_raw$haemoglobin)),
      mean(!is.na(mimic_raw$total_protein)),
      mean(!is.na(mimic_raw$albumin)),
      mean(!is.na(mimic_raw$protein_proxy)),
      mean(!is.na(mimic_raw$bmi)),
      mean(!is.na(mimic_raw$creatinine))
    ), 2)
  ),
  "MIMIC_feature_missingness.csv"
)

# ---- Projection and modelling -------------------------------------------------
project_to_nhanes <- function(df, dataset_label, protein_col) {
  feature_data <- data.frame(
    NLR = df$nlr,
    SII = df$sii,
    LBXHGB = df$haemoglobin,
    LBXSTP = df[[protein_col]],
    BMXBMI = df$bmi,
    LBXSCR = df$creatinine
  )

  keep <- stats::complete.cases(feature_data) &
    !is.na(df$survival_days_365) &
    !is.na(df$mortality_365d) &
    !is.na(df$anchor_age) &
    !is.na(df$gender)

  out <- df[keep, , drop = FALSE]
  feature_data <- feature_data[keep, , drop = FALSE]

  if (nrow(out) == 0) {
    stop("No complete cases for dataset: ", dataset_label, call. = FALSE)
  }

  z <- sweep(feature_data[nhanes_features], 2, nhanes_center, FUN = "-")
  z <- sweep(z, 2, nhanes_scale, FUN = "/")

  distances <- sapply(seq_len(nrow(center_matrix)), function(i) {
    rowSums((as.matrix(z) - matrix(center_matrix[i, ], nrow = nrow(z), ncol = ncol(z), byrow = TRUE))^2)
  })
  colnames(distances) <- rownames(center_matrix)

  assigned <- colnames(distances)[max.col(-distances, ties.method = "first")]
  out$validation_dataset <- dataset_label
  out$cloned_phenotype_raw <- assigned
  out$cloned_phenotype <- stats::relevel(factor(assigned, levels = c("1", "2", "3")), ref = "3")
  out$plot_phenotype <- factor(
    assigned,
    levels = c("1", "2", "3"),
    labels = c("Phenotype 1", "Phenotype 2", "Phenotype 3 (Ref)")
  )

  for (feature in nhanes_features) {
    out[[paste0("z_", feature)]] <- z[[feature]]
  }
  out$projection_distance <- apply(distances, 1, min)
  out
}

strict_data <- project_to_nhanes(
  mimic_raw,
  dataset_label = "Strict total-protein validation",
  protein_col = "total_protein"
)

proxy_data <- project_to_nhanes(
  mimic_raw,
  dataset_label = "Albumin-proxy validation",
  protein_col = "protein_proxy"
)

write_csv_utf8(strict_data, "MIMIC_projected_strict_total_protein_dataset.csv")
write_csv_utf8(proxy_data, "MIMIC_projected_albumin_proxy_dataset.csv")

dataset_summary <- dplyr::bind_rows(strict_data, proxy_data) |>
  dplyr::group_by(.data$validation_dataset) |>
  dplyr::summarise(
    n = dplyr::n(),
    deaths_365d = sum(.data$mortality_365d == 1, na.rm = TRUE),
    mortality_365d_percent = round(100 * mean(.data$mortality_365d == 1, na.rm = TRUE), 2),
    age_median_iqr = iqr_text(.data$anchor_age, digits = 1),
    bmi_median_iqr = iqr_text(.data$bmi, digits = 1),
    total_protein_nonmissing = sum(!is.na(.data$total_protein)),
    albumin_nonmissing = sum(!is.na(.data$albumin)),
    .groups = "drop"
  )
write_csv_utf8(dataset_summary, "Table4A_MIMIC_validation_dataset_summary.csv")

phenotype_counts <- dplyr::bind_rows(strict_data, proxy_data) |>
  dplyr::group_by(.data$validation_dataset, .data$cloned_phenotype_raw) |>
  dplyr::summarise(
    n = dplyr::n(),
    deaths_365d = sum(.data$mortality_365d == 1, na.rm = TRUE),
    mortality_365d_percent = round(100 * mean(.data$mortality_365d == 1, na.rm = TRUE), 2),
    .groups = "drop"
  ) |>
  dplyr::group_by(.data$validation_dataset) |>
  dplyr::mutate(percent = round(100 * .data$n / sum(.data$n), 2)) |>
  dplyr::ungroup() |>
  dplyr::select(
    .data$validation_dataset,
    .data$cloned_phenotype_raw,
    .data$n,
    .data$percent,
    .data$deaths_365d,
    .data$mortality_365d_percent
  )
write_csv_utf8(phenotype_counts, "Table4B_MIMIC_validation_phenotype_counts.csv")

feature_profiles <- dplyr::bind_rows(strict_data, proxy_data) |>
  dplyr::group_by(.data$validation_dataset, .data$cloned_phenotype_raw) |>
  dplyr::summarise(
    n = dplyr::n(),
    nlr = iqr_text(.data$nlr),
    sii = iqr_text(.data$sii),
    haemoglobin = iqr_text(.data$haemoglobin),
    total_protein = iqr_text(.data$total_protein),
    albumin = iqr_text(.data$albumin),
    protein_proxy = iqr_text(.data$protein_proxy),
    bmi = iqr_text(.data$bmi),
    creatinine = iqr_text(.data$creatinine),
    .groups = "drop"
  )
write_csv_utf8(feature_profiles, "Table4C_MIMIC_validation_feature_profiles.csv")

fit_cox <- function(df, dataset_label) {
  if (length(unique(df$cloned_phenotype_raw)) < 2) {
    return(data.frame(
      validation_dataset = dataset_label,
      variable = "Model not fitted",
      hazard_ratio_95ci = NA_character_,
      p_value = NA_character_,
      note = "Fewer than two phenotypes were present."
    ))
  }

  fit <- survival::coxph(
    survival::Surv(survival_days_365, mortality_365d) ~ cloned_phenotype + anchor_age + gender,
    data = df
  )
  s <- summary(fit)
  ci <- as.data.frame(s$conf.int)
  co <- as.data.frame(s$coefficients)

  row_names <- rownames(ci)
  pretty <- row_names
  pretty <- sub("^cloned_phenotype1$", "Phenotype 1 vs Phenotype 3", pretty)
  pretty <- sub("^cloned_phenotype2$", "Phenotype 2 vs Phenotype 3", pretty)
  pretty <- sub("^anchor_age$", "Age, per year", pretty)
  pretty <- sub("^genderM$", "Male vs female", pretty)

  out <- data.frame(
    validation_dataset = dataset_label,
    variable = pretty,
    hazard_ratio_95ci = sprintf("%.2f (%.2f-%.2f)", ci[, "exp(coef)"], ci[, "lower .95"], ci[, "upper .95"]),
    p_value = format_p(co[, "Pr(>|z|)"]),
    note = "",
    row.names = NULL
  )

  ph <- survival::cox.zph(fit)
  ph_out <- data.frame(
    validation_dataset = dataset_label,
    variable = rownames(ph$table),
    chisq = ph$table[, "chisq"],
    df = ph$table[, "df"],
    p = ph$table[, "p"],
    row.names = NULL
  )
  write_csv_utf8(
    ph_out,
    if (identical(dataset_label, "Strict total-protein validation")) {
      "Cox_PH_check_strict_total_protein.csv"
    } else {
      "Cox_PH_check_albumin_proxy.csv"
    }
  )

  out
}

cox_table <- dplyr::bind_rows(
  fit_cox(strict_data, "Strict total-protein validation"),
  fit_cox(proxy_data, "Albumin-proxy validation")
)
write_csv_utf8(cox_table, "Table4D_MIMIC_validation_Cox_365d.csv")

plot_km <- function(df, dataset_label, filename) {
  if (length(unique(df$plot_phenotype)) < 2) {
    message("Skipping KM plot for ", dataset_label, ": fewer than two phenotypes present.")
    return(invisible(NULL))
  }

  fit <- survival::survfit(
    survival::Surv(survival_days_365, mortality_365d) ~ plot_phenotype,
    data = df
  )
  surv_sum <- summary(fit)
  km_df <- data.frame(
    time = surv_sum$time,
    survival = surv_sum$surv,
    lower = surv_sum$lower,
    upper = surv_sum$upper,
    plot_phenotype = sub("^plot_phenotype=", "", surv_sum$strata)
  )

  start_rows <- data.frame(
    time = 0,
    survival = 1,
    lower = 1,
    upper = 1,
    plot_phenotype = levels(df$plot_phenotype)
  )
  km_df <- dplyr::bind_rows(start_rows, km_df)
  km_df$plot_phenotype <- factor(km_df$plot_phenotype, levels = levels(df$plot_phenotype))

  logrank <- survival::survdiff(
    survival::Surv(survival_days_365, mortality_365d) ~ plot_phenotype,
    data = df
  )
  p <- stats::pchisq(logrank$chisq, df = length(logrank$n) - 1, lower.tail = FALSE)
  p_label <- paste0("Log-rank P ", ifelse(p < 0.001, "<0.001", paste0("= ", sprintf("%.3f", p))))

  p_obj <- ggplot(km_df, aes(x = .data$time, y = .data$survival, color = .data$plot_phenotype)) +
    geom_step(linewidth = 0.9) +
    scale_color_manual(values = c("Phenotype 1" = "#C73E3A", "Phenotype 2" = "#2F6F9F", "Phenotype 3 (Ref)" = "#3E8B4E")) +
    scale_x_continuous(limits = c(0, 365), breaks = c(0, 90, 180, 270, 365)) +
    scale_y_continuous(limits = c(0, 1), labels = function(x) paste0(round(100 * x), "%")) +
    labs(
      title = dataset_label,
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
    ggplot2::ggsave(file.path(output_dir, filename), plot = p_obj, width = 8, height = 6, dpi = 300)
  }
}

plot_km(strict_data, "MIMIC-IV v3.1 External Validation: Strict Total-Protein Set", "Figure4A_MIMIC_KM_strict_total_protein.png")
plot_km(proxy_data, "MIMIC-IV v3.1 External Validation: Albumin-Proxy Set", "Figure4B_MIMIC_KM_albumin_proxy.png")

# ---- Method note --------------------------------------------------------------
method_note <- c(
  "# MIMIC-IV External Validation Method Note",
  "",
  paste0("- MIMIC root: `", mimic_root, "`"),
  paste0("- NHANES source object: `", nhanes_rdata_path, "`"),
  paste0("- Lab window: ICU intime -", lab_window_before_hours, "h to +", lab_window_after_hours, "h"),
  paste0("- BMI window: ICU intime -", bmi_lookback_days, " days to +1 day; direct OMR BMI only"),
  "- Cohort: first ICU stay per subject, age >=65 years.",
  "- Endpoint: 365-day all-cause mortality after ICU admission; non-events are censored at 365 days.",
  "- Projection: MIMIC variables are standardized using NHANES means and SDs, then assigned to the nearest NHANES phenotype centroid by Euclidean distance.",
  "- Strict total-protein validation uses the exact NHANES nutritional marker analogue, but sample size is expected to be small because total protein is sparse in ICU labs.",
  "- Albumin-proxy validation replaces total protein with albumin only when total protein is unavailable; this is a transportability/sensitivity analysis, not a fully identical external clone.",
  "",
  "## Output Files",
  "",
  "- `NHANES_reference_scaling_used.csv`",
  "- `NHANES_reference_centers_used.csv`",
  "- `MIMIC_feature_missingness.csv`",
  "- `MIMIC_projected_strict_total_protein_dataset.csv`",
  "- `MIMIC_projected_albumin_proxy_dataset.csv`",
  "- `Table4A_MIMIC_validation_dataset_summary.csv`",
  "- `Table4B_MIMIC_validation_phenotype_counts.csv`",
  "- `Table4C_MIMIC_validation_feature_profiles.csv`",
  "- `Table4D_MIMIC_validation_Cox_365d.csv`",
  "- `Figure4A_MIMIC_KM_strict_total_protein.png`",
  "- `Figure4B_MIMIC_KM_albumin_proxy.png`"
)
writeLines(method_note, con = file.path(output_dir, "MIMIC_validation_method_note.md"), useBytes = TRUE)

message("\nDone. Validation outputs saved to: ", output_dir)
message("Strict total-protein n = ", nrow(strict_data))
message("Albumin-proxy n = ", nrow(proxy_data))
message("Main Cox table: ", file.path(output_dir, "Table4D_MIMIC_validation_Cox_365d.csv"))
