# Release status: public GitHub release

Version/tag: `v1.0.0-freeze-v2`  
Freeze basis: Freeze V2 plus source-locked Amendments 14 and 15  
Public repository: `https://github.com/wangzy10-tom/immunonutritional-renal-phenotypes`  
Zenodo DOI: pending archival of the GitHub release

## Completed in this candidate

- Copied the frozen core and sensitivity scripts without changing the original
  project files.
- Replaced personal drive paths with `PROJECT_ROOT` and data-root variables.
- Added explicit support steps for the NHANES projection reference, the
  transitional OASIS-like file, and eICU APACHE extraction.
- Kept MIMIC-IV/eICU source data and all row-level intermediates outside the
  candidate package.
- Excluded all restricted external-cohort data and results, credentials,
  figures, PDFs, and Word files.
- Added disclosure rules, citation metadata, license, pipeline manifest, and
  frozen aggregate expectations.
- Rebuilt the Table S19 locked model and historical transport as an NHANES-only
  chain. Fourteen model outputs and ten historical-transport outputs were
  byte-identical to the frozen numeric sources; coefficients, predictions, and
  the 6.67% threshold were unchanged.
- Completed a clean rerun of the permitted NHANES, MIMIC-IV, eICU,
  cross-database, post-freeze, and exploratory NHANES stages without generating
  figures or PDFs.
- Added deterministic score-quartile tie-breaking and a direct eICU APACHE
  join under Amendment 15. Corrected publication-table QA passed 53/53 and
  final result-dictionary QA passed 15/15.

## Remaining archival step

1. Archive the exact GitHub tag `v1.0.0-freeze-v2` in Zenodo.
2. After the DOI exists, replace the manuscript's repository placeholder with
   the GitHub URL, release tag, commit hash, and version-specific Zenodo DOI.

The GitHub repository may be cited as the public code location. It should not
be described as a DOI-backed permanent archive until the Zenodo record is live
and its files have been checked against this tag.
