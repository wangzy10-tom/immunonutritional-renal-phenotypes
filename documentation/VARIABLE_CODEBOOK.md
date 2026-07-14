# Core variable codebook

| Construct | NHANES | MIMIC-IV | eICU |
|---|---|---|---|
| Inflammatory ratio | NLR from absolute neutrophil and lymphocyte counts | NLR from first-24-hour absolute counts | NLR-like ratio from first-24-hour differential percentages |
| Systemic inflammation | SII = platelets × neutrophils / lymphocytes | Same count-based formula | SII-like = platelets × polys percentage / lymphocyte percentage |
| Hematologic reserve | Hemoglobin | Hemoglobin | Hemoglobin |
| Protein marker | Total protein in the primary discovery analysis | Predominantly albumin-based protein proxy; harmonized albumin sensitivity reported separately | Total protein in the primary ICU analysis; harmonized albumin sensitivity reported separately |
| Body reserve | BMI | BMI from locally defined pre-ICU/early-ICU window | BMI |
| Renal marker | Creatinine | Creatinine | Creatinine |
| Main outcome | Long-term all-cause mortality | 365-day all-cause mortality | In-hospital mortality |
| Main severity handling | Multivariable survey-weighted adjustment | Official OASIS quartile-stratified Cox model | APACHE IVa-adjusted logistic mixed model with hospital random intercept |

All clustering labels are assigned without mortality information. Cross-database
agreement refers to recovery of a related biological structure and direction,
not identical membership, effect magnitude, or outcome definition.
