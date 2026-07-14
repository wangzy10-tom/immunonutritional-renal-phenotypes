# ==============================================================================
# Build the local NHANES Freeze V2 reference object required by script 09.
#
# This compatibility object contains individual-level NHANES records and is
# written only to the ignored local output directory. It must never be committed
# or distributed. The MIMIC main analysis is de novo; the projection step is an
# extraction/diagnostic bridge and is not claimed as fixed-centroid validation.
# ==============================================================================

required_packages <- c("dplyr", "readr", "tibble")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop("Missing required package(s): ", paste(missing_packages, collapse = ", "), call. = FALSE)
}

project_root <- normalizePath(
  Sys.getenv("PROJECT_ROOT", unset = getwd()), winslash = "/", mustWork = TRUE
)
cohort_path <- file.path(
  project_root, "output", "nhanes_2011_2018_rebuild",
  "NHANES_2011_2018_rebuilt_analytical_cohort.rds"
)
result_path <- file.path(
  project_root, "output", "nhanes_robust_reanalysis",
  "NHANES_robust_reanalysis_results.rds"
)
output_dir <- file.path(project_root, "output", "nhanes_projection_reference")
output_path <- file.path(output_dir, "NHANES_FreezeV2_projection_reference.RData")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (!all(file.exists(c(cohort_path, result_path)))) {
  stop("Run scripts 27 and 28 before building the projection reference.", call. = FALSE)
}

cohort <- readRDS(cohort_path)
assignments <- readRDS(result_path)$assignments
required_cohort <- c("SEQN", "NLR", "SII", "LBXHGB", "LBXSTP", "BMXBMI", "LBXSCR")
required_assignments <- c("SEQN", "phenotype_log6_winsor")
if (!all(required_cohort %in% names(cohort)) ||
    !all(required_assignments %in% names(assignments))) {
  stop("Freeze V2 source objects do not contain the required projection fields.", call. = FALSE)
}

nhanes_clustered <- cohort |>
  dplyr::select(dplyr::all_of(required_cohort)) |>
  dplyr::inner_join(
    assignments |>
      dplyr::select(dplyr::all_of(required_assignments)),
    by = "SEQN"
  ) |>
  dplyr::transmute(
    SEQN,
    NLR,
    SII,
    LBXHGB,
    LBXSTP,
    BMXBMI,
    LBXSCR,
    Phenotype = sub("^P", "", as.character(phenotype_log6_winsor))
  ) |>
  dplyr::filter(
    stats::complete.cases(dplyr::across(dplyr::all_of(required_cohort[-1L]))),
    Phenotype %in% c("1", "2", "3")
  )

if (nrow(nhanes_clustered) != 4636L || anyDuplicated(nhanes_clustered$SEQN)) {
  stop("Unexpected Freeze V2 reference cohort size or duplicate SEQN.", call. = FALSE)
}

save(nhanes_clustered, file = output_path, compress = "xz")
counts <- nhanes_clustered |>
  dplyr::count(Phenotype, name = "n") |>
  dplyr::mutate(percent = 100 * n / sum(n))
readr::write_csv(counts, file.path(output_dir, "NHANES_FreezeV2_projection_reference_counts.csv"))

message("Local projection reference created: ", output_path)
