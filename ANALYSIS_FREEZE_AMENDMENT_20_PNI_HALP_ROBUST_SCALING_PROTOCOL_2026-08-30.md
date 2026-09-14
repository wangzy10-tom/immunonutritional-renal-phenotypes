# Analysis Freeze Amendment 20 — Robust PNI/HALP Benchmark Scaling

Date locked: 30 August 2026  
Scope: NHANES conventional-index benchmark module only

## Reason for amendment

Post-result review identified that the original PNI comparator was standardized directly on its raw scale. A small number of very high absolute lymphocyte counts produced extreme high PNI values and an excessively large standardized range. This made the reported per-1-SD lower-PNI contrast difficult to interpret. HALP uses the same lymphocyte component but had already been protected by winsorization and log transformation; it is nevertheless included in this audit because the underlying concern applies to both lymphocyte-based indices.

This amendment is a transparent post-result data-distribution correction. It is not prospectively registered and will not be used to alter the clustering solution, phenotype labels, primary analysis population, primary covariate model, or interpretation of the primary P1-versus-P3 result.

## Locked population and model

- Source: the frozen Amendment 17 NHANES primary model data.
- Required sample: 3,979 participants and 720 deaths.
- Survey design: NHANES examination weights, strata, and primary sampling units, unchanged from Amendment 17.
- Covariates: age, sex, race/ethnicity, family income-to-poverty ratio, expanded comorbidity score, survey cycle, smoking, hypertension, and diabetes.
- Outcome: all-cause mortality with follow-up in months.

## Locked index definitions

- PNI = 10 × albumin (g/dL) + 5 × absolute lymphocyte count (10^9/L).
- HALP = haemoglobin (g/L) × albumin (g/L) × absolute lymphocyte count (10^9/L) / platelet count (10^9/L).

## Locked robust scaling

- PNI: winsorize the raw PNI distribution at the 1st and 99th percentiles, standardize the winsorized value, and multiply by −1 so that a 1-SD increase represents lower-PNI risk.
- HALP: winsorize the raw HALP distribution at the 1st and 99th percentiles, apply the natural logarithm, standardize, and multiply by −1 so that a 1-SD increase represents lower-HALP risk.
- Quantile limits, counts affected by winsorization, and pre/post-transformation ranges will be reported.

## Locked analyses

1. Separate survey-weighted Cox models for corrected PNI and corrected HALP, each reporting the HR per 1-SD higher risk direction, 95% CI, and P value.
2. Separate models adding the corrected PNI or corrected HALP risk score to the frozen phenotype model, reporting the adjusted P1-versus-P3 HR and the global phenotype Wald P value.
3. A deterministic QA check confirming unchanged N, event count, phenotype assignments, survey variables, and covariate completeness.

## Interpretation boundary

These analyses compare association scales only. They do not establish that one index or the clustering framework has superior discrimination, calibration, clinical utility, or causal relevance. The original uncorrected PNI benchmark will be superseded in the current manuscript and supplement but retained in the historical Amendment 17 files for auditability.
