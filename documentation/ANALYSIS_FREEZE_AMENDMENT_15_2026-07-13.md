# Analysis Freeze Amendment 15

Date: 13 July 2026

Status: completed; supersedes only the affected exact Freeze V2 values below.

## What changed

### 1. Deterministic OASIS and SOFA quartiles

OASIS and SOFA quartiles are now assigned after sorting by ascending score and then by `stay_id`. Reversing the source-row order produced identical assignments. Both scores yielded group sizes of 287, 286, 286, and 286 in the 1,145-patient MIMIC-IV cohort.

The primary MIMIC-IV P1 result remains 1.67 (95% CI 1.36-2.05). The affected displayed sensitivity estimates are:

| Analysis | Previous display | Amendment 15 display |
|---|---:|---:|
| MIMIC-IV selection-IPW, P1 vs P3 | HR 1.72 (1.30-2.28) | HR 1.71 (1.29-2.25) |
| MIMIC-IV selection-IPW, P2 vs P3 | HR 1.07 (0.77-1.49) | HR 1.06 (0.77-1.47) |
| MIMIC-IV 24-h landmark, P1 vs P3 | HR 1.68 (1.36-2.08) | HR 1.68 (1.36-2.07) |
| SOFA-quartile model, P1 vs P3 | HR 1.57 (1.28-1.93) | HR 1.56 (1.27-1.92) |
| SOFA selection-IPW, P1 vs P3 | HR 1.66 (1.26-2.19) | HR 1.65 (1.25-2.17) |

### 2. eICU APACHE-complete cohort

A direct join of the strict first-24-hour phenotype cohort to the raw APACHE IVa table identified two additional eligible patients who had been omitted because the earlier object inherited APACHE availability from a superseded wider-window cohort. No additional deaths were introduced.

| Quantity | Previous | Amendment 15 |
|---|---:|---:|
| APACHE-complete primary-model n | 12,546 | 12,548 |
| In-hospital deaths | 1,884 | 1,884 |
| Primary P1 vs P3 | OR 1.28 (1.11-1.48) | OR 1.28 (1.11-1.48) |
| Selection-IPW P1 vs P3 | OR 1.13 (0.98-1.31) | OR 1.13 (0.98-1.31) |
| Primary-model inclusion | 18.24% | 18.24% |

## Interpretation

The amendments are reproducibility corrections, not evidence-strengthening choices. NHANES results are unchanged. The MIMIC-IV primary result, eICU displayed primary effect, event counts, direction of all principal comparisons, and the manuscript's conclusion remain unchanged. The eICU analysis remains a multicenter robustness audit rather than confirmatory external validation.

## Verification

- Corrected publication-table QA: 53/53 passed.
- Final Freeze V2 result-dictionary QA: 15/15 passed.
- OASIS row-order invariance: passed.
- SOFA row-order invariance: passed.
- No image or PDF was generated.
