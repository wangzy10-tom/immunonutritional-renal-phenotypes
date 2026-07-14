# Analysis Freeze Amendment 14 — NHANES cluster-count and inflammation-axis sensitivities

Date: 13 July 2026  
Status: accepted after independent audit  
Scope: source-locked expansion of Table S3 and bounded manuscript clarification; Freeze V2 primary results remain primary

## Why this amendment was performed

Pre-submission review identified two unresolved methodological questions:

1. whether the NHANES high-risk structure remained visible when K was changed from the prespecified K=3 to K=2 or K=6; and
2. whether simultaneous inclusion of the mathematically related NLR and SII measures was necessary for the P1 mortality association.

The analysis was locked before outcome modelling in `ANALYSIS_FREEZE_AMENDMENT_14_PROTOCOL_2026-07-13.md`. Historical K and single-index mortality outputs based on an earlier 3,992-person adjustment set were explicitly excluded.

## Locked population and model

- Clustering population: final NHANES phenotype cohort, n=4,636 with 831 deaths during linked follow-up.
- Outcome-model population: exact final common covariate-complete cohort, n=3,979 with 720 deaths.
- Primary preprocessing: 1st/99th percentile winsorization; natural-log transformation of NLR, SII, and creatinine; within-cohort Z standardization.
- Primary outcome model: the frozen complex-survey Cox model adjusted for age, sex, race/ethnicity, income-to-poverty ratio, extended comorbidity score, survey cycle, smoking, hypertension, and diabetes.
- Mortality and follow-up were excluded from clustering and cluster labelling.

## Mandatory reproduction gate

Before the new analyses were accepted, the program reproduced:

- all 4,636 frozen K=3 assignments exactly;
- the frozen P1/P2/P3 counts of 817/1,894/1,925; and
- both frozen primary phenotype hazard-ratio coefficients with absolute difference 0.000e+00.

All mandatory gates passed.

## Results

### Cluster-count sensitivity

| Solution | Locked comparison | Fully adjusted survey HR (95% CI) | Interpretation |
|---|---|---:|---|
| K=2 | Higher versus lower vulnerability | 1.10 (0.91-1.31) | Direction positive but interval crossed 1. |
| K=6 | V1 versus V6 | 1.14 (0.86-1.50) | V1 was the outcome-blind highest vulnerability rank and P1-like flag; interval crossed 1. |
| K=6 | V2 versus V6 | 1.47 (1.08-2.01) | V2 contained 494/506 (97.6%) primary-P1 members. |
| K=6 | V3 versus V6 | 0.79 (0.58-1.07) | No ordered gradient. |
| K=6 | V4 versus V6 | 0.76 (0.50-1.16) | No ordered gradient. |
| K=6 | V5 versus V6 | 0.68 (0.49-0.96) | No ordered gradient. |

The adjusted Rand index relative to primary K=3 was 0.515 for K=2 and 0.215 for K=6. K=6 subdivided the primary structure: V6 contained only primary P3 members, V2 was predominantly primary P1, and V1 was predominantly primary P2 despite having the largest prespecified centroid vulnerability score. The K=6 risk ordering was therefore nonmonotonic.

### Symmetric NLR-only and SII-only sensitivity

| Variant | P1 versus P3 HR (95% CI) | P2 versus P3 HR (95% CI) | Adjusted Rand index versus primary K=3 |
|---|---:|---:|---:|
| NLR retained; SII excluded | 1.90 (1.56-2.33) | 1.13 (0.87-1.47) | 0.491 |
| SII retained; NLR excluded | 1.71 (1.34-2.19) | 0.99 (0.77-1.28) | 0.354 |

The P1 mortality direction remained positive in both symmetric variants, whereas P2 remained non-core and both P2 intervals crossed 1. Individual membership nevertheless changed materially.

## Interpretation decision

1. The six-feature K=3 solution remains the sole primary phenotype analysis. No sensitivity solution is promoted to co-primary status.
2. The new analyses reduce, but do not eliminate, concern that the P1 association is an artefact of jointly entering NLR and SII: either measure alone preserved the P1 direction.
3. The K=2 and K=6 results do not validate K=3 as the unique natural taxonomy. They show that K=3 is resampling stable and interpretable under the frozen design, while cluster count, membership, and risk ordering remain analysis dependent.
4. The K=6 results must be reported in full. V2 cannot be selected post hoc as the only K=6 result simply because its interval excluded 1.
5. The harmonized-albumin analysis remains a sensitivity analysis rather than a co-primary analysis.
6. eICU is described as a multicenter robustness audit/stress test, not confirmatory external validation.

## Manuscript impact

- The Abstract and Introduction now identify the contribution as a biological stratification framework rather than a deployable prediction tool.
- Methods and Table S3 report the locked K=2, K=6, NLR-only, and SII-only analyses.
- Results disclose the non-significant K=2 contrast, nonmonotonic K=6 findings, and positive P1 directions in both single-index variants.
- Discussion and Conclusion explicitly state cluster-count and membership dependence.
- The eICU role is narrowed to a multicenter robustness audit because its selection-weighted and hospital-summary intervals crossed 1.
- The exploratory 36-month tool is reduced to one main-text sentence and remains detailed only in the Supplementary Material.

No main table, figure, reference, primary sample size, event count, primary effect estimate, or frozen phenotype assignment was changed.

## Audit and provenance

- Internal analysis QA: 13/13 checks passed.
- Independent audit: 12/12 checks passed.
- Independent sensitivity-model reproduction: maximum numerical difference 2.220e-16.
- Analysis source manifest: `output/amendment14_nhanes_cluster_sensitivity/AMENDMENT14_SOURCE_MANIFEST.csv`.
- Independent audit manifest: `output/amendment14_nhanes_cluster_sensitivity/AMENDMENT14_AUDIT_MANIFEST.csv`.
- Complete aggregate results: `output/amendment14_nhanes_cluster_sensitivity/`.

## Decision

Accept Amendment 14 under Freeze V2. Supersede any earlier manuscript wording that implies K=3 is uniquely optimal, that the phenotype is cluster-count invariant, that NLR and SII are independent features, or that eICU provides confirmatory multicenter validation.

