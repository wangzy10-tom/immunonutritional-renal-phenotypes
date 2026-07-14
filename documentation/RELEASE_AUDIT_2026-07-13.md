# Public-release audit — updated 14 July 2026

## Overall determination

**Static packaging audit: PASS.**  
**Authorization for public GitHub upload: CONFIRMED BY THE AUTHOR.**

The package completed a clean rerun of the permitted analysis pipeline and is
approved for upload to the author-controlled public repository. A Zenodo DOI
remains a separate archival step after creation of the exact GitHub release tag.

## One-pass audit results

| Check | Result |
|---|---:|
| Files covered by the hash manifest | 54 |
| R scripts parsed | 36/36 |
| Release-manifest scripts present | 36/36 |
| Duplicate script names in release manifest | 0 |
| Personal/local absolute path hits | 0 |
| Credential or secret-assignment hits | 0 |
| Raw or row-level data files included | 0 |
| Excluded restricted-dataset keyword/file hits | 0 |
| Corrected publication-table / result-dictionary QA | 53/53; 15/15 |
| Unresolved repository placeholders | 0 |

The `repository-code` field in `CITATION.cff` points to
`https://github.com/wangzy10-tom/immunonutritional-renal-phenotypes`.

## Scientific and execution checks

- The portable release scripts were synchronized back to the local project
  after the clean rerun so the project and release candidate use the same
  deterministic implementation.
- The release copies use `PROJECT_ROOT`, `NHANES_DATA_DIR`, `MIMIC_ROOT`,
  `EICU_DIR`, and `MIMIC_CODE_ROOT` rather than personal drive paths.
- A local-only NHANES projection-reference builder was added to remove the
  unpublished legacy-object dependency. Its row-level output is written below
  ignored `output/` and is prohibited from distribution.
- Transitional OASIS-like and APACHE support steps are explicitly listed in the
  expanded pipeline manifest. The official MIT-LCP OASIS analysis remains the
  Freeze V2 MIMIC primary severity model.
- The locked Table S19 model was rebuilt using NHANES data only. Fourteen core
  model tables were byte-identical to their frozen sources; coefficients,
  development/validation predictions, and the 6.67% threshold had zero numeric
  difference. The historical NHANES transport was then rerun with 1,000
  bootstrap repetitions, and ten manuscript-facing tables were byte-identical.
- The superseded script that produced a different locked threshold was excluded
  from the release candidate; its outputs are not distributed.
- The full permitted stages were executed in a clean local project root using
  authorized local data. NHANES manuscript-facing outputs were reproduced;
  MIMIC-IV and eICU were rebuilt from database-level sources. No figure or PDF
  was generated.
- The clean rerun exposed input-row-order dependence in tied OASIS/SOFA
  quartiles and a two-patient eICU APACHE-join omission. Amendment 15 locked a
  non-outcome-based `stay_id` tie-breaker and the direct APACHE join. Corrected
  publication-table QA passed 53/53 and result-dictionary QA passed 15/15.

## Remaining release gates

1. Create and inspect the public GitHub tag `v1.0.0-freeze-v2`, then archive
   that exact release in Zenodo and add the DOI to the manuscript.

Until the Zenodo gate passes, the manuscript may cite the working GitHub
repository but must not claim a DOI-backed permanent archive.
