# Analysis-ready input contract

The current six-feature cohorts are 4,637 NHANES participants, 1,100 MIMIC-IV stays, and 15,242 eICU stays. Primary outcome samples are 3,979/720, 1,100/550, and 13,234/1,990 (n/deaths), respectively. Inputs must reproduce these locked populations; do not silently relax eligibility to make a run finish.

## NHANES

`NHANES_ELIGIBLE_RDS`: the 5,418-person 2011-2018 mortality-eligible age-65-or-older denominator. Required features are NLR, SII, LBXHGB, LBXSAL, BMXBMI, LBXSCR; original LBDNENO/LBDLYMNO/LBXPLTSI are needed for benchmark scores. Include SEQN, Cycle_ID, RIDAGEYR, RIAGENDR, RIDRETH3, INDFMPIR, PERMTH_INT, MORTSTAT, SDMVPSU, SDMVSTRA, WTMEC2YR, and all MCQ components listed in the accompanying field dictionary. `NHANES_COVARIATE_CACHE` contains BPQ, DIQ, and SMQ XPT files for cycles G, H, I, and J. Pipeline assertions are the authoritative required-column checks.

The cluster features come from the MEC examination but current follow-up is PERMTH_INT (interview origin). PERMTH_EXM is not substituted. Eight-year MEC weights equal WTMEC2YR/4. NLR uses absolute counts; SII equals platelets times NLR. Hemoglobin/albumin are g/dL, BMI kg/m^2, creatinine mg/dL, counts 10^3 cells/uL.

## MIMIC-IV

`MIMIC_DENOMINATOR_CSV`: 33,671 candidate first ICU stays aged at least 65 under the original anchor_age rule. Include stay_id, subject_id, age/sex/race fields, first_careunit, anchor_year_group, nlr, sii, haemoglobin, albumin, bmi, creatinine, mortality_365d, survival_days_365, and hospital_expire_flag. OASIS/SOFA inputs contain stay_id plus the corresponding official score. Required ancillary columns are asserted in the scripts.

NLR uses differential percentages (51256/51244), not absolute counts; platelet itemid is 51265. Components are selected separately in 0-24 h and are not specimen-paired. Hemoglobin=51222, albumin=50862, creatinine=50912. BMI is direct OMR BMI nearest admission within -365 to +1 days. The explicit eICU BMI 10-80 filter must not be represented as a common three-database filter. Primary and SOFA sensitivity model adjustment sets differ as specified in the scripts and supplement.

## eICU

`EICU_DENOMINATOR_CSV`: 68,798 candidate first ICU records; include patientunitstayid, hospitalid, age_num, gender_model, ethnicity, unittype, nlr, sii_like, haemoglobin, albumin, bmi, creatinine, hospital_mortality, icu_mortality, hospitaldischargeoffset, unitdischargeoffset, and discharge status/time fields needed by landmark modules. Do not restrict this input to APACHE-complete patients before clustering.

`EICU_DIR` contains apachePatientResult.csv.gz. Only IVa rows are used and negative apachescore is missing. The optional `EICU_APACHE_AUDIT_RDS` has a model_data element with patientunitstayid and apachescore; the K-resolution module currently requires the raw source. Lab medians are separate within-window summaries after plausibility filtering. BMI uses admission kg and cm and is missing outside 10-80 kg/m^2. Full field and range definitions are in `SOURCE_FIELD_DICTIONARY.md`.

## Reproduction boundary

This package starts from documented analysis-ready denominators. It does not yet provide a validated from-raw end-to-end extractor for every database. Independent users must reconstruct these inputs under the stated rules and source-data access terms. Do not claim that the candidate has passed a new isolated end-to-end rerun until such a run is completed.
