# Immunonutritional-renal vulnerability phenotypes: Freeze V2 code

This repository contains the path-normalized R code release for the manuscript
"Data-driven immunonutritional-renal vulnerability phenotypes and mortality in
older adults across community and critical-care settings."

## Release status

This is the public Freeze V2 release (`v1.0.0-freeze-v2`) corresponding to the
submission-stage analysis. The exact GitHub release is intended for archival in
Zenodo; see `RELEASE_STATUS.md` for the DOI status.

## Scientific scope

The code implements the Freeze V2 analyses in three settings:

- NHANES 2011-2018 discovery and survey-weighted mortality analyses;
- MIMIC-IV v3.1 first-24-hour conceptual replication with official OASIS and
  official first-day SOFA derived from MIT-LCP `mimic-code` v3.0.1;
- eICU-CRD v2.0 first-24-hour multicenter stress test with APACHE IVa adjustment.

The ICU clusters are refitted independently. They are not presented as a
fixed-centroid external validation. The locally reconstructed OASIS-like score
is retained only as a transparent transition dependency; the official OASIS
analysis supersedes it for all Freeze V2 scientific claims.

## Data are not included

This repository must not contain source data, extracts, row-level intermediate
files, credentials, or small disclosure-prone cells. NHANES source files are
publicly downloadable. MIMIC-IV and eICU require credentialed PhysioNet access
and must remain on the authorized user's local system. See
`documentation/DATA_ACCESS_AND_DISCLOSURE.md`.

## Setup

Use R 4.6.0. Install the packages checked at the beginning of each script. The
principal packages are `survey`, `survival`, `dplyr`, `readr`, `lme4`,
`metafor`, `duckdb`, `DBI`, `cluster`, `mclust`, and `glmnet`.

1. Copy `config/.env.example` outside the repository or set the variables in
   the shell. The scripts do not read a committed `.env` file.
2. Set `PROJECT_ROOT` to the repository root.
3. Set `NHANES_DATA_DIR`, `MIMIC_ROOT`, and `EICU_DIR` as needed.
4. Clone MIT-LCP `mimic-code`, check out tag `v3.0.1`, and set
   `MIMIC_CODE_ROOT` to its `mimic-iv/concepts_duckdb` directory.
5. Run scripts from the repository root. Outputs are written below `output/`,
   which is ignored by version control.

Figure writing is disabled by default. Set `WRITE_FIGURES=true` only when
figures are intentionally required; this is not needed to reproduce the frozen
numeric results.

## Execution order

`scripts/00_run_pipeline.R` provides stage-specific execution. The default
stage is `nhanes`; restricted ICU stages run only when the corresponding local
datasets are available. The expanded dependency order is documented in
`documentation/RELEASE_PIPELINE_MANIFEST.csv`.

Examples:

```text
Rscript scripts/00_run_pipeline.R nhanes
Rscript scripts/00_run_pipeline.R mimic
Rscript scripts/00_run_pipeline.R eicu
Rscript scripts/00_run_pipeline.R cross_database
Rscript scripts/00_run_pipeline.R postfreeze
Rscript scripts/00_run_pipeline.R exploratory_tool
```

The public outputs expected from the frozen analysis are listed in
`aggregate_results/EXPECTED_KEY_RESULTS.csv`. That table contains aggregate
values only and is not a substitute for rerunning the models.

## Deliberate exclusions

- No restricted external-cohort data or derived asset is included.
- The Table S19 model and historical transport scripts use NHANES data only.
  Their neutralized release copies reproduced the locked 6.67% formula and all
  manuscript-facing numeric outputs exactly.
- Figures, PDFs, Word files, and manuscript layout assets are not included.
- Third-party MIT-LCP SQL files are not vendored; users fetch the pinned release.

## Reproducibility boundary

The permitted NHANES, MIMIC-IV, eICU, cross-database, post-freeze, and
exploratory NHANES stages were rerun in a clean local project root using
authorized local data. The rerun identified and corrected input-row-order
dependence in tied OASIS/SOFA quartiles and a two-patient eICU APACHE-join
omission; these changes are locked in Amendment 15. All corrected table and
result-dictionary checks passed. This source-locked GitHub release is ready for
archival; the version-specific Zenodo DOI is recorded after deposit.

## License and citation

Project code is released under the MIT License. Third-party data and code retain
their own access terms and licenses. Citation metadata and the public repository
URL are in `CITATION.cff`; cite the version-specific Zenodo DOI once available.
