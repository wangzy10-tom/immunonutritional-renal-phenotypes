# Albumin-only stage runner. Run from the release root.
root <- normalizePath(Sys.getenv("PROJECT_ROOT", unset=getwd()), winslash="/", mustWork=TRUE)
setwd(root); Sys.setenv(PROJECT_ROOT=root)
rscript <- file.path(R.home("bin"), "Rscript"); if(.Platform$OS.type=="windows" && file.exists(paste0(rscript,".exe"))) rscript <- paste0(rscript,".exe")
stages <- list(nhanes="01_nhanes_albumin_analysis.R", mimic="02_mimic_albumin_analysis.R", eicu="03_eicu_albumin_analysis.R", structural="04_cross_database_structural_sensitivity.R", verify="05_verify_key_results.R", pni="../82_NHANES_PNI_HALP_robust_benchmark.R", resolution="../86_K_resolution_sensitivity_2026-08-30.R", kselection="../87_K_selection_diagnostics_2026-08-30.R", reviewer="../87_post_result_reviewer_sensitivities_2026-08-30.R")
args <- commandArgs(trailingOnly=TRUE); stage <- if(length(args)) args[1] else "verify"
order <- if(stage=="all") names(stages) else stage
if(any(!order %in% names(stages))) stop("Choose nhanes, mimic, eicu, structural, verify, pni, resolution, kselection, reviewer, or all.", call.=FALSE)
for(s in order){ path <- file.path(root,"scripts",stages[[s]]); message("Running ",s,"..."); status <- system2(rscript,shQuote(path)); if(!identical(status,0L)) stop("Stage failed: ",s,call.=FALSE) }
message("Completed stage: ",stage)

