### Table S1C. Source fields for the six clustering measurements

| Measurement | NHANES field and unit | MIMIC-IV field and unit | eICU field and unit |
|---|---|---|---|
| NLR | LBDNENO/LBDLYMNO; absolute neutrophil and lymphocyte counts in 10³ cells/µL | labevents itemid 51256/51244; differential percentages | labname -polys/-lymphs; separately summarized differential percentages |
| SII/SII-like | LBXPLTSI × NLR; platelets in 10³ cells/µL | labevents itemid 51265 × NLR; platelets in 10³ cells/µL | labname platelets x 1000 × NLR; platelets in 10³ cells/µL |
| Hemoglobin | LBXHGB; g/dL | labevents itemid 51222; g/dL | labname Hgb; g/dL |
| Albumin | LBXSAL; g/dL | labevents itemid 50862; g/dL | labname albumin; g/dL |
| BMI | BMXBMI; kg/m² | omr result_name matching BMI/body mass index, numeric result_value; kg/m² | admissionweight/[admissionheight/100]²; weight in kg, height in cm |
| Creatinine | LBXSCR; mg/dL | labevents itemid 50912; mg/dL | labname creatinine; mg/dL |

**Aggregation and cleaning:** NHANES required finite values for all six features and positive NLR, SII, albumin, and creatinine. MIMIC-IV selected each nonmissing laboratory component independently within 0–24 h of ICU admission, ordered by time from admission and then charttime; direct BMI was selected from the stated OMR window. Its six-feature cohort likewise required finite features and positive NLR, SII, albumin, and creatinine. The explicit 10–80 kg/m² BMI filter was applied in eICU, not uniformly across all three pipelines. For eICU, values were filtered before within-window medians: neutrophil/lymphocyte percentages >0 and ≤100, platelets >0 and ≤2,000 ×10³ cells/µL, hemoglobin 3–25 g/dL, albumin 0.5–8 g/dL, and creatinine 0.1–25 mg/dL. BMI outside 10–80 kg/m² was set to missing. Final eICU clustering additionally required finite features, 0<NLR≤100, 0<SII-like≤200,000, and positive albumin and creatinine. These extraction filters precede database-specific 1st/99th-percentile winsorization; they are not diagnostic thresholds.

### Table S1D. Primary-model covariate and outcome coding

| Database | Variable | Operational definition |
|---|---|---|
| NHANES | Age and sex | RIDAGEYR, continuous; male=1 for RIAGENDR=1 and 0 for RIAGENDR=2 |
| NHANES | Race/ethnicity | RIDRETH3 as a categorical variable: 1 Mexican American, 2 other Hispanic, 3 non-Hispanic White, 4 non-Hispanic Black, 6 non-Hispanic Asian, 7 other/multiracial |
| NHANES | Income-to-poverty ratio and cycle | INDFMPIR, continuous; Cycle_ID G/H/I/J as categorical indicators for 2011–2012/2013–2014/2015–2016/2017–2018 |
| NHANES | Smoking | Never if SMQ020=2; current if SMQ020=1 and SMQ040=1 or 2; former if SMQ020=1 and SMQ040=3; other combinations missing; never is the reference |
| NHANES | Hypertension | BPQ020=1 yes, 2 no; other responses missing |
| NHANES | Diabetes | DIQ010=1 yes, 2 no, 3 borderline; other responses missing; no is the reference |
| NHANES | Expanded comorbidity score | Complete sum of the 12 indicators listed in Table S1B |
| NHANES | Survey design and mortality | WTMEC8YR=WTMEC2YR/4, strata SDMVSTRA, primary sampling unit SDMVPSU; MORTSTAT=1 for death, follow-up PERMTH_INT in months |
| MIMIC-IV | Age, sex, and severity | patients.anchor_age determined eligibility; male derived from gender; official OASIS and first-day SOFA joined by stay_id. Primary Cox: sex plus OASIS-quartile strata. SOFA sensitivities additionally adjusted for continuous anchor_age and sex |
| MIMIC-IV | Mortality and follow-up | patients.dod linked by subject_id; death within 0–365 days after ICU intime counted as the 365-day event; time was the integer day difference, with same-day deaths assigned 0.5 day and survivors censored at 365 days. In-hospital mortality used hospital_expire_flag |
| eICU | Age, sex, and severity | Numeric age; top-coded strings beginning with > assigned 90 years; male=1 for Male and 0 for Female, other gender values missing. apachescore from apachePatientResult restricted to apacheversion=IVa, negative scores set to missing, one record per patientunitstayid |
| eICU | Mortality | hospitaldischargestatus Expired=1, Alive=0, other/missing status missing; ICU sensitivity used unitdischargestatus analogously. hospitalid defined the random intercept |

**Within-database joins and first-stay rules:** NHANES cycle files and linked mortality were joined by SEQN. MIMIC-IV patients, admissions, and ICU stays were joined using subject_id/hadm_id, with the earliest ICU intime per subject_id selected and stay_id breaking ties; laboratory and OMR records were linked to this stay's admission window. eICU patient/laboratory/APACHE records were joined by patientunitstayid. Candidate eICU records had unitvisitnumber=1 or missing; among records satisfying the age and first-visit rule, one record per uniquepid was selected in ascending patientunitstayid order, with patientunitstayid used when uniquepid was missing. These operational rules do not constitute linkage between the three databases. Native fields and derived indices were not adjudicated as new clinical diagnoses.

**Grouping and labeling:** All primary K-means fits used the Lloyd algorithm, K=3, nstart=100, iter.max=500, and seed 20260710. P1 mapping scores were standardized across all three cluster medians. For P2, hemoglobin, albumin, and BMI medians were standardized across the two remaining clusters before their reversed Z scores were summed. The mortality outcome did not enter either mapping step. The OASIS and SOFA strata used equal-sized rank groups ordered by score and then stay_id; identical boundary scores can therefore occur in adjacent groups (Tables S9B and S11B).
