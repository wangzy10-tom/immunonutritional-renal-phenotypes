# Analysis Freeze Amendment 14 — Locked protocol

Date locked: 13 July 2026  
Status at lock: no Amendment 14 outcome model had been run  
Scope: bounded NHANES sensitivity analyses requested during pre-submission peer review

## Scientific question

The amendment tests two specific vulnerabilities of the frozen NHANES discovery analysis without changing its primary status:

1. whether the adverse multidomain structure remains visible when the prespecified cluster count is changed from K=3 to K=2 or K=6; and
2. whether the P1 mortality association remains directionally similar when the mathematically dependent NLR/SII pair is represented by NLR alone or SII alone.

The six-feature, robust-transformed K=3 solution remains the sole primary NHANES phenotype analysis regardless of the amendment results.

## Frozen inputs

- Phenotype cohort: `output/nhanes_covariate_upgrade/NHANES_2011_2018_covariate_augmented.rds` (expected n=4,636; expected deaths=831).
- Exact final model cohort and covariates: `output/nhanes_albumin_benchmarks/NHANES_albumin_benchmark_results.rds`, object `model_data` (expected n=3,979; expected deaths=720).
- Frozen primary assignments: `output/nhanes_robust_reanalysis/NHANES_2011_2018_variant_assignments.csv`.
- Frozen numerical reference: `output/final_analysis_freeze_v2/UNIQUE_RESULT_DICTIONARY.csv`.

Historical K=2/K=4 mortality tables and historical NLR-only/SII-only model tables were generated with an earlier 3,992-person adjustment set. They are not admissible evidence for this amendment and will not be imported.

## Locked preprocessing and clustering

- Clustering is fitted in all 4,636 phenotype-complete participants, before outcome modelling.
- Each feature is winsorised at the cohort-specific 1st and 99th percentiles.
- NLR, SII, and creatinine are log transformed; haemoglobin, total protein, and BMI remain on their winsorised original scales.
- Selected features are standardised to Z scores.
- K-means uses Lloyd's algorithm, 100 random starts, and a maximum of 500 iterations.
- Fixed seeds are 20260710 for K=3 and the single-index variants, 20260712 for K=2, and 20260716 for K=6.
- Mortality, follow-up time, survey weights, and covariates are excluded from clustering and from phenotype labelling.

## Locked labels

### K=3 verification and single-index variants

The existing outcome-blind rule is retained. P1 is the cluster with the largest standardised inflammation-plus-creatinine centroid score. Among the other two clusters, P2 has the larger low-haemoglobin/low-total-protein/low-BMI score; the remaining cluster is P3.

### K=2 and K=6

Clusters are ranked without outcomes by the sum of six standardised adverse centroid directions: higher log NLR, higher log SII, lower haemoglobin, lower total protein, lower BMI, and higher log creatinine.

- K=2: `Higher vulnerability` versus `Lower vulnerability`.
- K=6: V1 through V6, ordered from highest to lowest multidomain vulnerability; V6 is the model reference.
- A K=6 `P1-like` cluster is separately designated, before outcome analysis, as the cluster with the largest higher-log-NLR/higher-log-SII/higher-log-creatinine score.

All K=6 contrasts will be reported; no cluster will be selected for reporting according to its mortality estimate.

## Locked outcome model

Assignments are joined by SEQN to the exact 3,979-person final model cohort. Every amendment contrast uses the frozen complex-survey Cox specification:

`Surv(PERMTH_INT, MORTSTAT) ~ phenotype + RIDAGEYR + male + race + INDFMPIR + Comorbidity_Score_Extended + cycle + smoking + hypertension + diabetes`

The survey design uses SDMVPSU, SDMVSTRA, WTMEC8YR, nesting, and `survey.lonely.psu = "adjust"`. Effect estimates are hazard ratios with 95% confidence intervals. P values are descriptive and unadjusted because all amendment analyses are explicitly sensitivity analyses.

## Mandatory reproduction gate

Before any new result is accepted, the new program must:

1. reproduce the frozen K=3 assignments exactly for all 4,636 participants;
2. reproduce the frozen P1/P2/P3 counts of 817/1,894/1,925; and
3. reproduce the frozen final P1 and P2 survey-weighted Cox estimates to an absolute tolerance of 1e-10.

Failure of any gate invalidates the entire amendment.

## Reporting decision rules

- The amendment cannot promote K=2, K=6, NLR-only, or SII-only to a co-primary analysis.
- The harmonised-albumin analysis remains a sensitivity analysis and will not be promoted because of its favourable estimate.
- Evidence across K is described as preservation, subdivision, attenuation, or failure to preserve a biomarker-defined structure; it is not described as validation of a unique K.
- Cross-database agreement continues to mean direction of the P1 association only, not identical membership, effect magnitude, protein measurement, effect measure, or mortality horizon.
- Results will enter the manuscript only after an independent audit passes.

