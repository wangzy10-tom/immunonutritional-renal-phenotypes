# Analysis Freeze Amendment 15 Protocol

Date: 13 July 2026

## Purpose

This amendment prospectively locks two reproducibility corrections identified during a clean rerun of the Freeze V2 code package. It does not introduce a new hypothesis, outcome, phenotype definition, database, or outcome-driven model choice.

## Correction 1: deterministic MIMIC-IV severity-score quartiles

Official OASIS, official first-day SOFA, and the historical OASIS-like sensitivity previously used rank-based four-group allocation without an explicit tie-breaker. Because these scores are integer valued, patients with the same score could be assigned differently when the input-row order changed.

The amended rule is:

1. sort by the severity score in ascending order;
2. break score ties by ascending analysis-unit `stay_id`;
3. allocate the ordered cohort into four approximately equal groups;
4. verify that reversing the input-row order produces identical assignments.

The tie-breaker contains no mortality or other outcome information. The analysis cohort, scores, covariates, endpoints, phenotype labels, and model specifications remain unchanged.

## Correction 2: direct eICU APACHE join

The strict first-24-hour eICU phenotype cohort must be joined directly to the raw APACHE IVa table used by the current analysis. It must not inherit APACHE availability from an earlier wider-window cohort object. Eligibility remains determined by the prespecified first-ICU, age, feature-completeness, and APACHE-completeness rules.

## Required checks

- OASIS quartile assignment is invariant to reversed input-row order.
- SOFA quartile assignment is invariant to reversed input-row order.
- Primary NHANES estimates are unchanged.
- MIMIC-IV and eICU outcomes are not used in either correction.
- All corrected publication-table checks and the final result-dictionary checks pass.
- No individual-level data are included in any public or manuscript-facing artifact.

## Reporting rule

All affected exact values in the manuscript, supplementary methods, supplementary tables, submission metadata, code assertions, and result dictionary must be synchronized. Earlier Freeze V2 values remain historical provenance but are superseded where Amendment 15 explicitly supplies a corrected value.
