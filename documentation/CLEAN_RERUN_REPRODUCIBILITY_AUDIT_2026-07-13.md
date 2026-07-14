# Freeze V2 Clean-Rerun Reproducibility Audit

Date: 13 July 2026

## Scope

The Freeze V2 analysis was rerun in a clean temporary project root from database-level source inputs and local public NHANES caches. The clean root is outside the intended public release and contains restricted row-level derivatives that must not be uploaded.

## Results by database

### NHANES

- Final phenotype cohort: 4,636 participants and 831 deaths.
- Covariate-complete primary model: 3,979 participants and 720 deaths.
- Seven manuscript-facing Table 29 files were byte-identical to the frozen outputs.
- Three RDS objects differed only in stored derived columns or removal of a legacy comparison object; all common scientific columns were identical and downstream estimates were unchanged.

### MIMIC-IV

- Eligible denominator: 33,671.
- Strict first-24-hour phenotype cohort: 1,145 with 574 deaths.
- Phenotype counts: P1 515, P2 253, P3 377.
- Official OASIS values were identical to the prior extraction.
- The clean rerun exposed input-row-order dependence in tied-score quartile allocation. Amendment 15 introduced a deterministic `stay_id` tie-breaker and verified order invariance for OASIS and SOFA.

### eICU

- Eligible denominator: 68,798.
- Strict total-protein six-feature cohort: 14,241.
- Direct raw APACHE joining produced 12,548 complete primary-model patients and 1,884 deaths.
- Two additional eligible surviving patients explained the difference from 12,546. Their inclusion had no effect on the displayed primary P1 estimate.

## Cross-database and post-freeze checks

- Corrected publication-table QA: 53/53 passed.
- Final result-dictionary QA: 15/15 passed.
- Alternative clustering QA: 8/8 passed.
- Nonlinearity QA: 12/12 passed.
- Amendment 14 reproduction audit: 12/12 passed.
- Locked 36-month NHANES model: 13/13 passed.
- Historical 2005-2010 transport analysis: 14/14 passed.

## Publication impact

The clean rerun did not change the central scientific interpretation. It required small numerical updates to MIMIC-IV tied-score sensitivity estimates and a two-patient correction to the eICU APACHE-complete denominator. The manuscript, supplementary methods, and text tables were synchronized under Amendment 15. Figures and PDFs were deliberately not generated or modified.
