# Historical clean-rerun audit: base albumin package, 2026-08-30

This note restores a truncated documentation file. It summarizes the retained local freeze report and the accompanying `CLEAN_RERUN_KEY_RESULT_VERIFICATION_2026-08-30.csv`; it is not evidence of a new run.

The prior isolated run completed the NHANES base suite (including both 1,000-replicate bootstrap modules), MIMIC-IV analyses, eICU analyses, and cross-database structural sensitivities. The retained aggregate comparison contains 22 matching result rows, covering the verifier's then-implemented primary, IPW, and structural comparisons. The expected-anchor list contained 30 rows; the old verifier did not cover the remaining eight. All 30 are covered by the updated saved-output check dated 2026-09-04.

The prior base run does not establish a clean rerun of the four subsequently added modules or a validated raw-source extraction workflow. No participant-level outputs from the run are distributed.
