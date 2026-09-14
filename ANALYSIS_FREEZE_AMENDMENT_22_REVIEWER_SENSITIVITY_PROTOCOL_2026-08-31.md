# Analysis Freeze Amendment 22 — Reviewer-Driven Renal, BMI, and ICU Landmark Sensitivities

**Locked:** 31 August 2026, 01:01:50 +08:00  
**Status:** Locked before the formal Amendment 22 reproducibility executions  
**Execution script:** `87_post_result_reviewer_sensitivities_2026-08-30.R`  
**Locked script SHA-256:** `6207467294e8b298b3336003e7a1d75e05eefcffb0c195a072713a1d9f924652`  
**Formal result directory:** `output/amendment22_reviewer_sensitivities_2026-08-31`

## 1. Chronology and analysis identity

The primary unified-albumin analyses, earlier sensitivity analyses, and reviewer critiques were already known before this amendment. In addition, a preliminary execution of the analyses described below was performed and viewed on 30 August 2026 before this Amendment 22 document was created. Those preliminary aggregate outputs are retained in `output/post_result_reviewer_sensitivities_2026-08-30` with an explicit chronology notice.

Consequently, Amendment 22 is a **post-result, reviewer-driven exploratory sensitivity specification**. It is not prospectively registered, is not outcome-blind at the investigator level, and cannot be made prospective by repeating the analysis. The formal runs following this lock test deterministic reproducibility and enforce complete reporting; they are not independent confirmation.

No Amendment 22 result may replace the frozen primary clustering resolution, phenotype definitions, primary cohorts, or primary outcome models because it is more favourable. All listed results will be retained regardless of direction, confidence interval, or P value.

## 2. Questions addressed

Amendment 22 addresses three reviewer questions without redefining the main study:

1. Is the NHANES P1 mortality association completely reducible to the creatinine component used in P1 mapping?
2. Does including BMI materially determine the K=3 partition or P1 mortality direction?
3. Do the ICU associations persist after conditioning on survival and continued observation beyond 1, 3, and 7 days?

The amendment does not claim to eliminate residual confounding, reverse causation, selection bias, measurement error, or hospital-level heterogeneity.

## 3. Locked input populations

| Database/module | Locked source population | Expected n/events |
|---|---|---:|
| NHANES clustering | Unified-albumin six-feature-complete cohort | 4,637 / 831 deaths |
| NHANES outcome models | Frozen fully adjusted survey-Cox cohort | 3,979 / 720 deaths |
| MIMIC-IV clustering and outcome models | Unified-albumin six-feature-complete first-ICU-stay cohort | 1,100 / 550 deaths within 365 days |
| eICU clustering | Unified-albumin six-feature-complete first-ICU-stay cohort | 15,242 stays |
| eICU primary outcome models | APACHE IVa-complete hospital-mortality cohort | 13,234 / 1,990 in-hospital deaths |

No sample expansion is permitted in the no-BMI module. In particular, deleting BMI from the clustering matrix does not allow participants with missing BMI to enter; sample membership must remain identical to the six-feature primary cohort.

## 4. Module A: NHANES creatinine-component adjustment

### 4.1 Locked renal representation

Creatinine will be represented exactly as in the clustering input:

1. winsorise serum creatinine at the 1st and 99th percentiles in the full 4,637-person clustering cohort;
2. apply the natural logarithm;
3. standardise to mean 0 and standard deviation 1 in that clustering cohort;
4. join the fixed standardised value to the 3,979-person outcome-model cohort by participant identifier.

### 4.2 Locked models

Both models will use the NHANES complex-survey design with examination weights, strata, and primary sampling units.

- Reference model: phenotype plus age, sex, race/ethnicity, income-to-poverty ratio, expanded comorbidity score, survey cycle, smoking, hypertension, and diabetes.
- Renal-component model: the identical reference model plus standardised log-creatinine.

Both P1 versus P3 and P2 versus P3 hazard ratios, 95% confidence intervals, P values, sample sizes, and events will be reported.

### 4.3 Interpretation rule

Because creatinine helped define P1, this is a part-whole/non-reducibility analysis rather than ordinary confounder control. Attenuation will be described without claiming that renal disease has been eliminated. A confidence interval crossing 1 will be reported without changing the model.

## 5. Module B: fixed-cohort K=3 clustering after omitting BMI

### 5.1 Locked clustering inputs

Within each unchanged six-feature-complete cohort, the five clustering variables will be:

- NLR;
- SII or SII-like;
- haemoglobin;
- albumin;
- creatinine.

BMI will be deliberately omitted from distance calculation and label mapping. The remaining preprocessing is unchanged: 1st/99th percentile winsorisation, natural-log transformation of NLR, SII/SII-like, and creatinine, and within-database Z standardisation.

### 5.2 Locked clustering and label mapping

- K-means Lloyd algorithm;
- K=3;
- `nstart=100`;
- `iter.max=500`;
- random seed 20260710.

No outcome, follow-up time, event count, severity score, treatment, or model covariate may be used in clustering or label mapping.

P1 is the cluster maximising the sum of standardised log-median NLR, log-median SII/SII-like, and log-median creatinine. Among the remaining clusters, P2 is the cluster maximising reserve depletion defined by lower median haemoglobin and albumin only. The final cluster is P3. BMI cannot be reintroduced to break a tie or revise an unfavourable mapping.

### 5.3 Locked outcome models

- NHANES: fully adjusted complex-survey Cox model used in the primary analysis.
- MIMIC-IV: sex-adjusted 365-day Cox model stratified by deterministic official OASIS quartiles.
- eICU: age-, sex-, and APACHE IVa-adjusted hospital random-intercept logistic model.

For numerical stability only, eICU age and APACHE scores are Z-standardised before fitting. This linear rescaling does not change the covariate space or substantive adjustment.

### 5.4 Required outputs

For every database, report:

- K-means objective and iteration count;
- all three cluster sizes and original-scale median profiles, including BMI as a descriptive variable not used for clustering;
- adjusted Rand index, exact label agreement, and P1 recall, precision, and Jaccard overlap versus the primary six-feature K=3 solution;
- P1 versus P3 and P2 versus P3 effect estimates, 95% confidence intervals, and P values;
- eICU convergence messages, warnings, and singularity status.

The no-BMI solution remains a sensitivity analysis and cannot replace the six-feature primary solution according to its result.

## 6. Module C: ICU landmark analyses

### 6.1 Fixed landmark times

All analyses will be performed at 1, 3, and 7 days. No time point may be added, removed, or promoted after examining estimates.

### 6.2 MIMIC-IV

For each landmark, retain patients with recorded 365-day survival time strictly greater than the landmark and reset analysis time to survival time minus the landmark.

Fit both:

1. phenotype plus sex, stratified by deterministic official OASIS quartile;
2. phenotype plus age and sex, stratified by deterministic official first-day SOFA quartile.

The endpoint remains 365-day all-cause mortality. Report P1 versus P3 and P2 versus P3 results, included sample size and events, early deaths excluded, and any non-events no longer observed before the landmark.

### 6.3 eICU

Two operational endpoints will be retained at every landmark:

1. in-hospital mortality among stays with hospital length of stay strictly greater than the landmark;
2. ICU mortality among stays with ICU length of stay strictly greater than the landmark.

Thus, early deaths and early live discharges are both excluded. Each model will include phenotype, standardised age, sex, standardised APACHE IVa, and a hospital random intercept. Report all P1 versus P3 and P2 versus P3 odds ratios, included sample sizes and events, hospitals represented, early deaths excluded, early live discharges excluded, warnings, convergence messages, and singularity status.

### 6.4 Interpretation rule

Landmark restriction can reduce sensitivity to deaths occurring shortly after baseline but cannot remove reverse causation. In eICU it conditions strongly on remaining in hospital or ICU; later-landmark estimates therefore cannot be interpreted as stronger causal effects or greater biological persistence merely because the point estimate increases.

## 7. Complete-reporting rules

The following must appear in the aggregate output and remain available for manuscript/supplement integration:

- both P1 and P2 contrasts in every fitted model;
- all 1-, 3-, and 7-day time points;
- both MIMIC-IV severity-adjustment approaches;
- both eICU endpoints;
- every confidence interval crossing 1;
- all sample and event losses at each landmark;
- all convergence warnings or model failures;
- all cluster-agreement metrics and profiles.

No cross-database effect pooling will be performed because outcomes, horizons, and effect measures differ.

## 8. QA, reproducibility, and privacy

The formal script must abort if expected source cohort sizes or event counts are not reproduced, any required estimate is missing or non-finite, eICU mixed models are singular or retain convergence warnings, or landmark sample sizes fail to decrease monotonically.

Only aggregate tables, summaries, diagnostics, result dictionaries, and SHA-256 manifests may be written. No participant-level assignments, restricted extracts, fitted model objects, credentials, or local database contents may be exported.

The locked script will be executed twice after this amendment. The first execution's aggregate output manifest will be saved outside the result directory, the second execution will overwrite the same formal directory, and the two manifest contents must match exactly. A mismatch requires investigation and must not be resolved by choosing the more favourable run.

## 9. Manuscript status

No manuscript or supplement wording will be changed until the formal reruns and manifest comparison are complete. Even if the formal outputs reproduce the preliminary results exactly, the Methods and limitations must disclose that this module was specified after the primary results and after a preliminary execution had been viewed.
