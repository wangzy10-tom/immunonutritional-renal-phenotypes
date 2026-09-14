# Analysis Freeze Amendment 19 — K-Resolution Sensitivity Protocol

**Date:** 30 August 2026  
**Status:** Locked before execution of the K=2 and K=4 outcome models  
**Purpose:** Evaluate whether the adverse biological and mortality-association direction depends on retaining K=3 as the primary clustering resolution.

## 1. Timing and analysis status

This module was added after the K=3 phenotype profiles and mortality results were known. It is therefore a post-result sensitivity analysis and is not described as prospectively preregistered. K=3 remains the primary descriptive resolution. K=2 and K=4 are supportive analyses and will not replace the primary analysis according to whether their results are favorable.

## 2. Locked cohorts and features

The module will use the same outcome-blind six-feature-complete clustering cohorts as the unified-albumin analysis:

- NHANES: 4,637 participants;
- MIMIC-IV: 1,100 first ICU stays;
- eICU: 15,242 first ICU stays.

The six clustering features remain NLR, SII/SII-like, hemoglobin, albumin, BMI, and creatinine. Within each database, preprocessing will be identical to the primary analysis: 1st/99th percentile winsorization, natural-log transformation of NLR, SII/SII-like, and creatinine, and within-database Z standardization.

## 3. Clustering specifications

K-means will be fitted independently within each database at K=2 and K=4 using Lloyd's algorithm, `nstart=100`, `iter.max=500`, random seed 20260710. No outcome, follow-up time, severity score, covariate, or event count will be used to fit clusters or map biological labels.

## 4. Outcome-blind label mapping

Cluster profiles will be calculated from within-cluster medians of the six original-scale features.

### K=2

1. Calculate the P1 score as the sum of standardized log-median NLR, log-median SII/SII-like, and log-median creatinine across the two raw clusters.
2. Label the cluster with the highest P1 score as P1.
3. Label the other cluster as P3.
4. No P2 label will be created at K=2.

### K=4

1. Label the cluster with the highest P1 score as P1.
2. Among the remaining clusters, label as P2 the cluster with the highest reserve-deficit score, defined as the negative sum of standardized median hemoglobin, albumin, and BMI.
3. Among the two remaining clusters, label as P3 the cluster with the highest favorable-reference score, defined as standardized median hemoglobin + albumin + BMI minus standardized log-median NLR, log-median SII/SII-like, and log-median creatinine.
4. Label the remaining cluster as P4.

If an exact score tie occurs, the smaller deterministic raw cluster number will be selected. The mapping will not be revised after outcome models are examined.

## 5. Locked outcome models

The outcome models and model populations will match the unified-albumin primary analyses:

- NHANES: complex-survey fully adjusted Cox model, n=3,979 with 720 deaths;
- MIMIC-IV: sex-adjusted Cox model stratified by deterministic official OASIS quartiles, n=1,100 with 550 deaths;
- eICU: age-, sex-, and APACHE IVa-adjusted hospital random-intercept logistic model, n=13,234 with 1,990 in-hospital deaths across 166 hospitals.

The common core comparison is P1 versus P3. At K=4, P2 versus P3 and P4 versus P3 will be reported descriptively but will not replace the common core comparison.

## 6. Prespecified outputs

For each database and K, the module will report:

- cluster sizes and median feature profiles;
- the complete outcome model estimates, 95% confidence intervals, and P values;
- P1-versus-P3 estimates as the common sensitivity result;
- adjusted Rand index and P1 recall, precision, and Jaccard overlap relative to the primary K=3 solution;
- K-means objective, iterations, and convergence status;
- eICU mixed-model warnings, convergence messages, and singularity status;
- outcome-blind biological-direction indicators comparing P1 with P3.

No cross-database pooling will be performed because the settings, outcome horizons, and effect measures differ.

## 7. Interpretation rules

The analysis will be considered supportive of resolution-insensitive adverse direction only if P1 retains higher inflammation, lower hemoglobin and albumin, and higher creatinine than P3 in the named resolution/database combination. BMI direction will be reported but is not required.

Mortality estimates will be reported regardless of direction or statistical significance. A confidence interval crossing 1 will not be concealed. Even if all P1 estimates exceed 1, this module cannot establish a unique optimal K, invariant membership, algorithm independence, a diagnostic category, or clinical deployment.

## 8. Reproducibility and privacy

Only aggregate tables, diagnostics, manifests, and a summary-only RDS may be written. No participant-level assignments, restricted extracts, credentials, or fitted model objects may be exported. The final script will be executed twice, and generated aggregate-file hashes must match before the module is considered complete.
