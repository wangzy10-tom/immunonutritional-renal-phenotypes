# ==============================================================================
# Add smoking, hypertension, and diabetes covariates to the rebuilt NHANES cohort
# ==============================================================================

required_packages <- c("haven", "dplyr", "readr", "tibble")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop("Missing required package(s): ", paste(missing_packages, collapse = ", "), call. = FALSE)
}

project_root <- normalizePath(
  Sys.getenv("PROJECT_ROOT", unset = getwd()), winslash = "/", mustWork = TRUE
)
input_path <- file.path(
  project_root, "output", "nhanes_2011_2018_rebuild",
  "NHANES_2011_2018_rebuilt_analytical_cohort.rds"
)
output_dir <- file.path(project_root, "output", "nhanes_covariate_upgrade")
cache_dir <- file.path(output_dir, "source_cache")
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(input_path)) {
  stop("Rebuilt NHANES cohort not found. Run script 27 first.", call. = FALSE)
}

cycles <- tibble::tribble(
  ~cycle, ~survey_start_year,
  "G", "2011",
  "H", "2013",
  "I", "2015",
  "J", "2017"
)
components <- c("BPQ", "DIQ", "SMQ")

is_valid_xpt <- function(path) {
  if (!file.exists(path) || file.info(path)$size < 1000) return(FALSE)
  con <- file(path, open = "rb")
  on.exit(close(con), add = TRUE)
  header <- rawToChar(readBin(con, what = "raw", n = 80L))
  grepl("HEADER RECORD", header, fixed = TRUE)
}

download_component <- function(component, cycle, survey_start_year) {
  filename <- paste0(component, "_", cycle, ".xpt")
  destination <- file.path(cache_dir, filename)
  url <- paste0(
    "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/", survey_start_year,
    "/DataFiles/", filename
  )

  if (!is_valid_xpt(destination)) {
    message("Downloading ", filename, "...")
    status <- tryCatch(
      utils::download.file(url, destination, mode = "wb", quiet = TRUE),
      error = function(e) e
    )
    if (inherits(status, "error") || !is_valid_xpt(destination)) {
      stop("Failed to download official NHANES file: ", url, call. = FALSE)
    }
  }

  destination
}

read_component <- function(component, cycle, survey_start_year) {
  path <- download_component(component, cycle, survey_start_year)
  raw <- haven::read_xpt(path)
  keep <- switch(
    component,
    BPQ = c("SEQN", "BPQ020"),
    DIQ = c("SEQN", "DIQ010"),
    SMQ = c("SEQN", "SMQ020", "SMQ040"),
    stop("Unknown component: ", component, call. = FALSE)
  )
  missing <- setdiff(keep, names(raw))
  if (length(missing) > 0L) {
    stop(
      "Missing variable(s) in ", component, "_", cycle, ": ",
      paste(missing, collapse = ", "), call. = FALSE
    )
  }
  raw |>
    dplyr::select(dplyr::all_of(keep)) |>
    dplyr::mutate(Cycle_ID = cycle)
}

cycle_covariates <- lapply(seq_len(nrow(cycles)), function(i) {
  component_data <- lapply(
    components,
    read_component,
    cycle = cycles$cycle[i],
    survey_start_year = cycles$survey_start_year[i]
  )
  Reduce(
    function(x, y) dplyr::full_join(x, y, by = c("SEQN", "Cycle_ID")),
    component_data
  )
})
covariates <- dplyr::bind_rows(cycle_covariates)

if (anyDuplicated(covariates$SEQN) > 0L) {
  stop("Duplicate SEQN detected in augmented covariate data.", call. = FALSE)
}

cohort <- readRDS(input_path)
augmented <- cohort |>
  dplyr::left_join(covariates, by = c("SEQN", "Cycle_ID")) |>
  dplyr::mutate(
    hypertension = dplyr::case_when(
      BPQ020 == 1 ~ "Yes",
      BPQ020 == 2 ~ "No",
      TRUE ~ NA_character_
    ),
    diabetes = dplyr::case_when(
      DIQ010 == 1 ~ "Yes",
      DIQ010 == 3 ~ "Borderline",
      DIQ010 == 2 ~ "No",
      TRUE ~ NA_character_
    ),
    smoking = dplyr::case_when(
      SMQ020 == 2 ~ "Never",
      SMQ020 == 1 & SMQ040 %in% c(1, 2) ~ "Current",
      SMQ020 == 1 & SMQ040 == 3 ~ "Former",
      TRUE ~ NA_character_
    ),
    hypertension = factor(hypertension, levels = c("No", "Yes")),
    diabetes = factor(diabetes, levels = c("No", "Borderline", "Yes")),
    smoking = factor(smoking, levels = c("Never", "Former", "Current"))
  )

if (nrow(augmented) != nrow(cohort)) {
  stop("Row count changed after covariate augmentation.", call. = FALSE)
}

availability <- tibble::tibble(
  variable = c("hypertension", "diabetes", "smoking"),
  nonmissing_n = c(
    sum(!is.na(augmented$hypertension)),
    sum(!is.na(augmented$diabetes)),
    sum(!is.na(augmented$smoking))
  )
) |>
  dplyr::mutate(
    total_n = nrow(augmented),
    nonmissing_percent = 100 * nonmissing_n / total_n
  )

cycle_availability <- augmented |>
  dplyr::group_by(Cycle_ID) |>
  dplyr::summarise(
    n = dplyr::n(),
    hypertension_nonmissing = sum(!is.na(hypertension)),
    diabetes_nonmissing = sum(!is.na(diabetes)),
    smoking_nonmissing = sum(!is.na(smoking)),
    .groups = "drop"
  )

readr::write_csv(availability, file.path(output_dir, "Table27B_A_covariate_availability.csv"))
readr::write_csv(cycle_availability, file.path(output_dir, "Table27B_B_cycle_availability.csv"))
readr::write_csv(augmented, file.path(output_dir, "NHANES_2011_2018_covariate_augmented.csv"))
saveRDS(augmented, file.path(output_dir, "NHANES_2011_2018_covariate_augmented.rds"))

message("NHANES covariate augmentation completed.")
print(availability)
print(cycle_availability)
