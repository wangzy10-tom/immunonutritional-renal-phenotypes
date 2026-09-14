# Prepublication source/output check, 2026-09-04

No models were fitted and no raw or participant-level data were read for this check.

- All 75 original manifest entries matched before editing.
- All 13 base R files other than the expanded launcher were byte-identical to both the prior local base release and its isolated-rerun copy. Subsequently, only the verifier and self-audit received checking improvements; the analysis files remain unchanged.
- All four later modules matched their original whole-file source after explicitly enumerated input-path substitutions. Their 39 top-level functions were identical.
- All 34 later-module CSV reference tables were byte-identical to the corresponding saved aggregate tables.
- All 18 R scripts parsed. The current key-result verifier checks the complete set of 30 expected anchors, not only the prior 22. All 30 saved-output comparisons passed for the non-missing expected estimates, confidence limits, counts, and event counts at tolerance 1e-10.
- NHANES IPW table counts are not exported separately. Its N/events are checked against the primary model cohort counts; the frozen code confirms that primary and IPW fits use the same complete-case rows.
- CSV headers, file types, local-path patterns, credential patterns, manifest coverage, and archive contents were screened. No restricted participant-level files are included.

## Packaging corrections

The old verifier omitted eight expected anchors and did not require equality of the expected/observed identifier sets. This check adds those eight and rejects missing or duplicate identifiers. The launcher quotes script paths containing spaces. Two truncated documentation passages were restored; the stage manifest now lists the four later modules. These changes do not alter cohort construction, features, clustering, seeds, covariates, endpoints, or fitted results.

## Evidence boundary

Reference results were already saved, not recalculated for publication. The prior clean rerun covered the base package. This expanded package has not undergone another complete isolated rerun; the author explicitly chose source/output reconciliation instead. Input reconstruction starts from the documented analysis-ready denominators. A fully validated from-raw extraction pipeline is not included.

## Publication status

The intended GitHub repository was read successfully. Its main branch still pointed to c9c3602d6fcfce59fc572a97dc03d7b2c27708b0, and the current connection reported no push permission. The local command-line client was also unauthenticated. No public write was attempted. A public commit, tag, or DOI must not be asserted until created and verified.
