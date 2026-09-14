# Multisystem profiles and mortality: analysis-ready code package

Package version: `v2.0.0-albumin`. Prepared for the manuscript "Multisystem profiles and mortality in older adults across community and critical-care cohorts: a multidatabase observational study". Source and saved-output checks were completed on 4 September 2026; publication metadata were prepared on 6 September 2026. See the repository's Releases page for the actual publication record. No archival DOI is assigned in this package.

## Scope

Features are NLR, SII/SII-like, hemoglobin, serum albumin, BMI, and creatinine. Clustering and label mapping are outcome-free, while overall study development and the explicitly identified later modules were outcome-aware. The primary P1 estimates are NHANES HR 2.08, MIMIC-IV HR 1.81, and eICU OR 1.43. These effect measures must not be pooled. The code is for research reproduction, not patient-level deployment.

The unchanged prior albumin code covers primary, selection-weighted, and structural analyses. This package additionally includes PNI/HALP rescaling, K=2/K=4, K-selection diagnostics, renal-component adjustment, same-cohort no-BMI analyses, and ICU landmarks, with aggregate reference tables and development protocols.

## Data and execution

No source data, patient-level extract, fitted row-level object, credentials, or real `.env` is included. Use only your own authorized NHANES/PhysioNet inputs. See `documentation/ANALYSIS_READY_DATA_CONTRACT.md` and `config/.env.example`. Scripts do not automatically load that example.

From the package root, set the input environment variables, then run an explicitly selected stage:

```text
Rscript scripts/00_run_pipeline.R nhanes
Rscript scripts/00_run_pipeline.R mimic
Rscript scripts/00_run_pipeline.R eicu
Rscript scripts/00_run_pipeline.R structural
Rscript scripts/00_run_pipeline.R verify
Rscript scripts/00_run_pipeline.R pni
Rscript scripts/00_run_pipeline.R resolution
Rscript scripts/00_run_pipeline.R kselection
Rscript scripts/00_run_pipeline.R reviewer
```

The four later stages require the completed primary outputs; the K-resolution stage requires the authorized raw APACHE IVa file. K-selection is computationally intensive. `all` executes every stage and refits models; it is not a read-only audit. Outputs under `output/` must remain private. Source/protocol names retain development identifiers to preserve provenance.

For a read-only integrity check: `Rscript scripts/99_audit_release.R`.

## Verification boundary

The base albumin release had a prior clean rerun with 22 matching aggregate anchors. On 4 September 2026, the analysis scripts were checked against that frozen source, and the entire source of each of the four later modules was reconciled against its original script. Differences in those four modules are limited to explicit input-path portability changes; all 39 top-level functions are unchanged.

The current publication check matched all 30 expected result anchors against saved aggregate outputs and confirmed byte-identical copies of all 34 later-module reference tables. The verifier was extended from 22 to 30 anchors and now rejects missing or duplicate result identifiers. Its numerical tolerance is 1e-10 for non-missing expected estimates, confidence limits, sample sizes, and event counts. The stage launcher now quotes script paths containing spaces. Neither change alters an analysis model.

This check did not rerun models. The authors elected to rely on the prior analyses and these source/output checks rather than repeat a complete computation. No new clean end-to-end rerun of the expanded package is claimed. The package begins at documented analysis-ready denominators, not a validated raw-database extraction workflow. See `documentation/PREPUBLICATION_CHECK_2026-09-04.md` for the exact evidence and boundaries.

## Public repository

Repository: https://github.com/wangzy10-tom/immunonutritional-renal-phenotypes.

This version is intended to supersede the earlier `v1.0.0-freeze-v2` code for the current manuscript. The earlier tag and commit c9c3602d6fcfce59fc572a97dc03d7b2c27708b0 remain historical records and do not represent the current six-feature albumin analysis. Cite the exact commit containing this package, together with its release identifier after publication. Historical audit files describe their status on the stated audit dates.
