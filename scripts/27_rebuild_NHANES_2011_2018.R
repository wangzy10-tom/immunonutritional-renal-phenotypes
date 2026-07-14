# ==============================================================================
# Rebuild NHANES 2011-2018 analytical cohort with validated mortality linkage
#
# This script intentionally writes to a versioned output directory and does not
# overwrite the original project files. It stops if any of the four survey
# cycles or mortality files is missing.
# ==============================================================================

required_packages <- c("haven", "dplyr", "readr", "tibble")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0L) {
  stop(
    "Missing required R package(s): ", paste(missing_packages, collapse = ", "),
    ". Install them before running this script.",
    call. = FALSE
  )
}

project_root <- normalizePath(
  Sys.getenv("PROJECT_ROOT", unset = getwd()), winslash = "/", mustWork = TRUE
)
data_dir <- Sys.getenv("NHANES_DATA_DIR", unset = file.path(project_root, "data", "nhanes"))
output_dir <- file.path(project_root, "output", "nhanes_2011_2018_rebuild")
mortality_cache_dir <- file.path(output_dir, "mortality_cache")
dir.create(mortality_cache_dir, recursive = TRUE, showWarnings = FALSE)

cycles <- tibble::tribble(
  ~cycle, ~survey_years,
  "G", "2011_2012",
  "H", "2013_2014",
  "I", "2015_2016",
  "J", "2017_2018"
)

core_features <- c("NLR", "SII", "LBXHGB", "LBXSTP", "BMXBMI", "LBXSCR")
mcq_variables <- c(
  "MCQ220", "MCQ160A", "MCQ160B", "MCQ160C", "MCQ160D", "MCQ160E",
  "MCQ160F", "MCQ160G", "MCQ160K", "MCQ160L", "MCQ160M", "MCQ160N"
)

required_source_files <- unlist(lapply(cycles$cycle, function(cycle) {
  file.path(data_dir, paste0(c("DEMO_", "CBC_", "BIOPRO_", "BMX_", "MCQ_"), cycle, ".xpt"))
}))
missing_source_files <- required_source_files[!file.exists(required_source_files)]
if (length(missing_source_files) > 0L) {
  stop(
    "Missing NHANES source file(s):\n", paste(missing_source_files, collapse = "\n"),
    call. = FALSE
  )
}

yes_no_flag <- function(x) {
  dplyr::case_when(
    x == 1 ~ 1L,
    x == 2 ~ 0L,
    TRUE ~ NA_integer_
  )
}

read_cycle <- function(cycle) {
  message("Reading NHANES cycle ", cycle, "...")

  demo <- haven::read_xpt(file.path(data_dir, paste0("DEMO_", cycle, ".xpt"))) |>
    dplyr::select(
      SEQN, RIDAGEYR, RIAGENDR, RIDRETH3, DMDEDUC2, DMDMARTL,
      INDFMPIR, WTMEC2YR, SDMVPSU, SDMVSTRA
    )

  cbc <- haven::read_xpt(file.path(data_dir, paste0("CBC_", cycle, ".xpt"))) |>
    dplyr::select(SEQN, LBDNENO, LBDLYMNO, LBXPLTSI, LBXHGB)

  biopro <- haven::read_xpt(file.path(data_dir, paste0("BIOPRO_", cycle, ".xpt"))) |>
    dplyr::select(SEQN, LBXSAL, LBXSCR, LBXSTP)

  bmx <- haven::read_xpt(file.path(data_dir, paste0("BMX_", cycle, ".xpt"))) |>
    dplyr::select(SEQN, BMXWT, BMXHT, BMXBMI)

  mcq <- haven::read_xpt(file.path(data_dir, paste0("MCQ_", cycle, ".xpt"))) |>
    dplyr::select(SEQN, dplyr::any_of(mcq_variables))

  demo |>
    dplyr::left_join(cbc, by = "SEQN") |>
    dplyr::left_join(biopro, by = "SEQN") |>
    dplyr::left_join(bmx, by = "SEQN") |>
    dplyr::left_join(mcq, by = "SEQN") |>
    dplyr::mutate(Cycle_ID = cycle)
}

download_mortality_file <- function(survey_years) {
  filename <- paste0("NHANES_", survey_years, "_MORT_2019_PUBLIC.dat")
  destination <- file.path(mortality_cache_dir, filename)
  url <- paste0(
    "https://ftp.cdc.gov/pub/Health_Statistics/NCHS/datalinkage/linked_mortality/",
    filename
  )

  if (!file.exists(destination) || file.info(destination)$size < 1000) {
    message("Downloading official mortality file: ", filename)
    status <- tryCatch(
      utils::download.file(url, destination, mode = "wb", quiet = TRUE),
      error = function(e) e
    )
    if (inherits(status, "error") || !file.exists(destination) || file.info(destination)$size < 1000) {
      stop("Failed to download mortality file: ", url, call. = FALSE)
    }
  }

  destination
}

read_mortality_file <- function(path, cycle) {
  positions <- readr::fwf_positions(
    start = c(1, 15, 16, 43),
    end = c(6, 15, 16, 45),
    col_names = c("SEQN", "ELIGSTAT", "MORTSTAT", "PERMTH_INT")
  )

  mortality <- readr::read_fwf(
    path,
    col_positions = positions,
    na = ".",
    col_types = readr::cols(
      SEQN = readr::col_integer(),
      ELIGSTAT = readr::col_integer(),
      MORTSTAT = readr::col_integer(),
      PERMTH_INT = readr::col_integer()
    ),
    progress = FALSE
  ) |>
    dplyr::mutate(Mortality_Cycle_ID = cycle)

  if (nrow(mortality) == 0L || anyDuplicated(mortality$SEQN) > 0L) {
    stop("Invalid mortality file for cycle ", cycle, call. = FALSE)
  }

  mortality
}

cycle_data <- lapply(cycles$cycle, read_cycle)
names(cycle_data) <- cycles$cycle
nhanes_merged <- dplyr::bind_rows(cycle_data)

cycle_source_counts <- nhanes_merged |>
  dplyr::count(Cycle_ID, name = "source_n") |>
  dplyr::arrange(Cycle_ID)

if (!identical(sort(unique(nhanes_merged$Cycle_ID)), cycles$cycle)) {
  stop("The merged source data do not contain all G/H/I/J cycles.", call. = FALSE)
}

nhanes_older <- nhanes_merged |>
  dplyr::filter(RIDAGEYR >= 65) |>
  dplyr::mutate(
    NLR = dplyr::if_else(LBDLYMNO > 0, LBDNENO / LBDLYMNO, NA_real_),
    SII = dplyr::if_else(LBDLYMNO > 0, LBDNENO * LBXPLTSI / LBDLYMNO, NA_real_),
    WTMEC8YR = WTMEC2YR / 4
  )

nhanes_complete_features <- nhanes_older |>
  dplyr::filter(dplyr::if_all(dplyr::all_of(core_features), ~ !is.na(.x) & is.finite(.x)))

mortality_paths <- mapply(
  download_mortality_file,
  cycles$survey_years,
  SIMPLIFY = FALSE,
  USE.NAMES = FALSE
)
mortality_data <- Map(read_mortality_file, mortality_paths, cycles$cycle)
names(mortality_data) <- cycles$cycle
all_mortality <- dplyr::bind_rows(mortality_data)

mortality_cycle_counts <- all_mortality |>
  dplyr::count(Mortality_Cycle_ID, name = "mortality_file_n") |>
  dplyr::arrange(Mortality_Cycle_ID)

if (!identical(sort(unique(all_mortality$Mortality_Cycle_ID)), cycles$cycle)) {
  stop("Mortality linkage data do not contain all G/H/I/J cycles.", call. = FALSE)
}

nhanes_older_linked <- nhanes_older |>
  dplyr::left_join(all_mortality, by = "SEQN") |>
  dplyr::filter(
    ELIGSTAT == 1,
    !is.na(MORTSTAT),
    !is.na(PERMTH_INT),
    PERMTH_INT > 0
  ) |>
  dplyr::mutate(
    complete_six_feature = dplyr::if_all(
      dplyr::all_of(core_features),
      ~ !is.na(.x) & is.finite(.x)
    )
  )

nhanes_final_2011_2018 <- nhanes_older_linked |>
  dplyr::filter(complete_six_feature)

if (any(nhanes_final_2011_2018$Cycle_ID != nhanes_final_2011_2018$Mortality_Cycle_ID)) {
  stop("Cycle mismatch detected after mortality linkage.", call. = FALSE)
}
if (!identical(sort(unique(nhanes_final_2011_2018$Cycle_ID)), cycles$cycle)) {
  stop("The final analytical cohort does not contain all G/H/I/J cycles.", call. = FALSE)
}

available_mcq <- intersect(mcq_variables, names(nhanes_final_2011_2018))
nhanes_final_2011_2018 <- nhanes_final_2011_2018 |>
  dplyr::mutate(dplyr::across(dplyr::all_of(available_mcq), yes_no_flag, .names = "flag_{.col}"))

flag_variables <- paste0("flag_", available_mcq)
nhanes_final_2011_2018 <- nhanes_final_2011_2018 |>
  dplyr::mutate(
    Comorbidity_Unknown_Count = rowSums(is.na(dplyr::pick(dplyr::all_of(flag_variables)))),
    Comorbidity_Score_Extended = dplyr::if_else(
      Comorbidity_Unknown_Count == 0,
      rowSums(dplyr::pick(dplyr::all_of(flag_variables))),
      NA_real_
    )
  )

cohort_flow <- tibble::tibble(
  stage = c(
    "All examined participants",
    "Age 65 years or older",
    "Complete six-feature matrix",
    "Mortality-linkage eligible final cohort"
  ),
  n = c(
    nrow(nhanes_merged),
    nrow(nhanes_older),
    nrow(nhanes_complete_features),
    nrow(nhanes_final_2011_2018)
  )
)

cycle_flow <- nhanes_merged |>
  dplyr::count(Cycle_ID, name = "all_examined") |>
  dplyr::left_join(nhanes_older |> dplyr::count(Cycle_ID, name = "age_65_plus"), by = "Cycle_ID") |>
  dplyr::left_join(
    nhanes_complete_features |> dplyr::count(Cycle_ID, name = "complete_features"),
    by = "Cycle_ID"
  ) |>
  dplyr::left_join(
    nhanes_final_2011_2018 |>
      dplyr::group_by(Cycle_ID) |>
      dplyr::summarise(
        final_n = dplyr::n(),
        deaths = sum(MORTSTAT == 1),
        mortality_percent = 100 * mean(MORTSTAT == 1),
        median_followup_months = stats::median(PERMTH_INT),
        .groups = "drop"
      ),
    by = "Cycle_ID"
  ) |>
  dplyr::arrange(Cycle_ID)

validation_summary <- tibble::tibble(
  item = c(
    "Cycles in source data",
    "Cycles in mortality data",
    "Cycles in final cohort",
    "Final cohort size",
    "Deaths",
    "Maximum follow-up months",
    "Duplicate SEQN in final cohort"
  ),
  value = c(
    paste(sort(unique(nhanes_merged$Cycle_ID)), collapse = ","),
    paste(sort(unique(all_mortality$Mortality_Cycle_ID)), collapse = ","),
    paste(sort(unique(nhanes_final_2011_2018$Cycle_ID)), collapse = ","),
    as.character(nrow(nhanes_final_2011_2018)),
    as.character(sum(nhanes_final_2011_2018$MORTSTAT == 1)),
    as.character(max(nhanes_final_2011_2018$PERMTH_INT)),
    as.character(anyDuplicated(nhanes_final_2011_2018$SEQN))
  )
)

readr::write_csv(cycle_source_counts, file.path(output_dir, "Table27A_source_cycle_counts.csv"))
readr::write_csv(mortality_cycle_counts, file.path(output_dir, "Table27B_mortality_file_cycle_counts.csv"))
readr::write_csv(cohort_flow, file.path(output_dir, "Table27C_cohort_flow.csv"))
readr::write_csv(cycle_flow, file.path(output_dir, "Table27D_cycle_specific_flow.csv"))
readr::write_csv(validation_summary, file.path(output_dir, "Table27E_validation_summary.csv"))
readr::write_csv(
  nhanes_final_2011_2018,
  file.path(output_dir, "NHANES_2011_2018_rebuilt_analytical_cohort.csv")
)
saveRDS(
  nhanes_final_2011_2018,
  file.path(output_dir, "NHANES_2011_2018_rebuilt_analytical_cohort.rds")
)
saveRDS(
  nhanes_older_linked,
  file.path(output_dir, "NHANES_2011_2018_age65_mortality_eligible.rds")
)

message("NHANES 2011-2018 rebuild completed successfully.")
print(cohort_flow)
print(cycle_flow)
print(validation_summary)
