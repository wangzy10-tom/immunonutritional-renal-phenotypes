# Analysis Freeze Amendment 21 — Cross-Database Clusterability and K-Selection Audit

**Date:** 30 August 2026  
**Status:** Locked before execution of the new diagnostics  
**Execution script:** `87_K_selection_diagnostics_2026-08-30.R`  
**Planned result directory:** `output/k_selection_diagnostics_amendment21_2026-08-30`

## 1. Timing, purpose, and boundary

The primary K=3 phenotype results, the K=2/K=4 resolution sensitivity results, and the earlier NHANES gap-statistic results were already known before this amendment. This is therefore a **post-result, outcome-blind methodological audit**, not a prospectively registered analysis.

The audit answers two separate questions:

1. Do the six-feature data in each database support a **unique natural cluster number**?
2. If not, can K=3 still be retained transparently as the study's **descriptive resolution of a continuous vulnerability structure**?

No outcome, follow-up, treatment, severity score, or event count is loaded or used. The audit does not redefine the biological labels or refit outcome models.

## 2. Locked inputs and preprocessing

| Database | Locked input | Expected n |
|---|---|---:|
| NHANES 2011–2018 | `output/nhanes_albumin_amendment17_2026-08-30/INTERNAL_NHANES_albumin_complete_cohort.rds` | 4,637 |
| MIMIC-IV v3.1 | `output/mimic_albumin_primary_2026-08-30/MIMIC_albumin_primary_results.rds` (`$analysis`) | 1,100 |
| eICU-CRD v2.0 | `output/eicu_24h_extraction/eICU_first_ICU_feature_availability_dataset.csv`, filtered exactly as in Amendment 19 | 15,242 |

The six features remain NLR, SII/SII-like, haemoglobin, albumin, BMI, and creatinine. Within each database, preprocessing is locked to 1st/99th percentile winsorisation, natural-log transformation of NLR, SII/SII-like, and creatinine, and within-database Z standardisation. K-means uses Euclidean distance, Lloyd's algorithm, `iter.max=500`, and seed 20260710 unless a diagnostic-specific seed is stated.

## 3. Candidate values and diagnostics

### 3.1 Gap statistic: primary selector that permits K=1

- Candidate range: K=1–6.
- `cluster::clusGap`, B=100 reference datasets, `nstart=20`.
- The primary selection rule is Tibshirani's 1-SE rule (`Tibs2001SEmax`).
- The global maximum and `firstSEmax` selections are recorded as supporting descriptions.
- K=1 is retained as a valid result and means that this diagnostic does not support more than one natural cluster.

### 3.2 Average silhouette width: separation selector conditional on K≥2

- Candidate range: K=2–6.
- Full-cohort K-means solutions use `nstart=100`.
- Silhouette values are evaluated on a fixed outcome-blind sample of at most 5,000 observations per database (all observations when n≤5,000; seed 20260711 otherwise).
- Selection is the K with the highest mean silhouette width; ties within 0.001 are resolved in favour of the smaller K.
- Because silhouette cannot evaluate K=1, it cannot override a gap-statistic selection of K=1.

### 3.3 Consensus clustering: ambiguity selector conditional on K≥2

- Candidate range: K=2–6.
- A fixed outcome-blind audit sample of at most 2,000 observations per database is used (seed 20260712).
- For each K, perform 100 repetitions of 80% subsampling without replacement and K-means with `nstart=20`.
- Record the consensus CDF, area under the CDF, change in area, and the proportion of ambiguous clustering (PAC) for consensus values between 0.1 and 0.9.
- Selection is the K with the lowest PAC; ties within 0.01 are resolved in favour of the smaller K.
- PAC is supporting evidence only and cannot by itself establish a natural K.

### 3.4 Subsampling stability: admissibility gate, not a K selector

- Candidate range: K=2–6.
- For each K, fit the full-cohort reference solution with `nstart=100`, then perform 200 repetitions of 80% subsampling without replacement and K-means with `nstart=20`.
- Compare each subsample solution with the reference labels on the same sampled observations using adjusted Rand index (ARI).
- For each reference cluster, calculate its maximum Jaccard overlap with any subsample cluster in every repetition and report the mean and distribution across repetitions.
- A K is **stable** only when every reference cluster has mean Jaccard ≥0.75. Mean Jaccard 0.60–<0.75 is labelled uncertain and <0.60 unstable. These cut-offs apply only to cluster-wise Jaccard, not to ARI.
- ARI is reported descriptively without imposing the Jaccard 0.75 threshold on it.

### 3.5 Minimum-size gate

A K is structurally admissible only if every full-cohort cluster contains at least 5% of that database's clustering cohort. This gate prevents selection driven by very small partitions.

## 4. Locked decision rules

### 4.1 Database-level conclusion

A database supports a **unique natural K** only when all of the following hold:

1. the primary gap-statistic 1-SE rule selects K≥2;
2. the silhouette and PAC selectors choose the same K as the gap 1-SE rule;
3. that K passes the cluster-wise Jaccard stability gate; and
4. that K passes the minimum-size gate.

If the gap 1-SE rule selects K≥2 and exactly one of silhouette or PAC agrees, with both gates passed, the result is labelled **partially supported K**, not unique K. All other patterns are labelled **no uniquely supported natural K**.

### 4.2 Cross-database conclusion

A common natural K is claimed only if all three databases independently meet the unique-K rule for the same value. Lesser agreement is described without promotion to a common natural K.

### 4.3 Consequence for the manuscript's primary resolution

- If all three databases independently support the same unique K, the manuscript's main resolution will be reconsidered and changed to that K, with transparent disclosure that the decision followed this post-result audit.
- Otherwise, K=3 remains the main **descriptive resolution**, because it displays the inflammatory–renal, low-reserve, and favourable-reference regions used by the study question; it will not be called statistically optimal, naturally occurring, or prospectively prespecified.
- K=2 and K=4 remain resolution-sensitivity analyses. No outcome result may influence this decision.

This rule makes a clear distinction between the number of natural clusters supported by the data and the resolution used to describe a continuous structure.

## 5. QA and reporting

The script must abort if cohort sizes differ from 4,637/1,100/15,242, if transformed values are non-finite, if any required repetition is missing, or if outcome-like fields are included in a clustering matrix. It writes aggregate diagnostics only, records package versions and seeds, and produces SHA-256 manifests. Participant-level assignments are not written.

The report will retain all results, including K=1 selections, discordant selectors, unstable K values, and results unfavourable to K=3.
