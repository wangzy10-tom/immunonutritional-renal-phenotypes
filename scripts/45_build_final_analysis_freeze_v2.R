# ==============================================================================
# Build final analysis freeze v2 and the unique manuscript result dictionary
# ==============================================================================

required_packages <- c("dplyr", "readr", "tibble")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop("Missing required package(s): ", paste(missing_packages, collapse = ", "), call. = FALSE)
}

root <- normalizePath(
  Sys.getenv("PROJECT_ROOT", unset = getwd()), winslash = "/", mustWork = TRUE
)
table_dir <- file.path(root, "output", "corrected_publication_tables")
output_dir <- file.path(root, "output", "final_analysis_freeze_v2")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

read_table <- function(name) {
  readr::read_csv(file.path(table_dir, name), show_col_types = FALSE)
}

table1 <- read_table("Table1_corrected_weighted_baseline.csv")
table2 <- read_table("Table2_corrected_NHANES_models.csv")
table2b <- read_table("Table2B_internal_validation.csv")
table3b <- read_table("Table3B_cross_database_sensitivity.csv")
table4a <- read_table("Table4A_adjusted_absolute_risk_RMST.csv")
table4b <- read_table("Table4B_adjusted_risk_RMST_contrasts.csv")
table5a <- read_table("Table5A_ICU_selection_flow.csv")
table5b <- read_table("Table5B_ICU_selection_IPW_models.csv")
table6a <- read_table("Table6A_harmonised_albumin_models.csv")
table6b <- read_table("Table6B_harmonised_albumin_agreement.csv")
table7d <- read_table("Table7D_official_SOFA_models.csv")
table8a <- read_table("Table8A_eICU_hospital_meta_analysis.csv")
table8b <- read_table("Table8B_eICU_leave_one_hospital_out.csv")
table8e <- read_table("Table8E_eICU_hospital_effect_directions.csv")
table9a <- read_table("Table9A_age_sex_global_interactions.csv")
table9b <- read_table("Table9B_P1_interaction_ratios.csv")
table9c <- read_table("Table9C_P1_conditional_effects.csv")

cindex <- readr::read_csv(
  file.path(root, "output", "nhanes_bootstrap_cindex", "Table34B_bootstrap_Cindex_summary.csv"),
  show_col_types = FALSE
)
stability <- readr::read_csv(
  file.path(root, "output", "nhanes_cluster_number_validation", "Table28B_C_subsampling_stability_summary.csv"),
  show_col_types = FALSE
)
variant_agreement <- readr::read_csv(
  file.path(root, "output", "nhanes_robust_reanalysis", "Table28F_variant_agreement.csv"),
  show_col_types = FALSE
)
loco_agreement <- readr::read_csv(
  file.path(root, "output", "nhanes_leave_one_cycle_out", "Table35A_cycle_projection_agreement.csv"),
  show_col_types = FALSE
)
loco_interaction <- readr::read_csv(
  file.path(root, "output", "nhanes_leave_one_cycle_out", "Table35E_cycle_interaction.csv"),
  show_col_types = FALSE
)
time_effects <- readr::read_csv(
  file.path(root, "output", "nhanes_corrected_bootstrap", "Table32E_time_varying_score_effects.csv"),
  show_col_types = FALSE
)
time_comparison <- readr::read_csv(
  file.path(root, "output", "nhanes_corrected_bootstrap", "Table32F_time_varying_model_comparison.csv"),
  show_col_types = FALSE
)

dictionary <- tibble::tibble(
  result_id = character(),
  tier = character(),
  dataset = character(),
  analysis = character(),
  n = integer(),
  events = integer(),
  outcome = character(),
  comparison_or_metric = character(),
  effect_measure = character(),
  estimate = double(),
  lower_95 = double(),
  upper_95 = double(),
  p_value = double(),
  multiplicity_adjusted_p = double(),
  display_value = character(),
  destination = character(),
  source_file = character(),
  allowed_statement_cn = character(),
  prohibited_statement_cn = character()
)

add_result <- function(
  result_id, tier, dataset, analysis, n = NA_integer_, events = NA_integer_,
  outcome = "Not applicable", comparison_or_metric, effect_measure,
  estimate = NA_real_, lower_95 = NA_real_, upper_95 = NA_real_,
  p_value = NA_real_, multiplicity_adjusted_p = NA_real_, display_value,
  destination, source_file, allowed_statement_cn, prohibited_statement_cn
) {
  dictionary <<- dplyr::bind_rows(
    dictionary,
    tibble::tibble(
      result_id = result_id,
      tier = tier,
      dataset = dataset,
      analysis = analysis,
      n = as.integer(n),
      events = as.integer(events),
      outcome = outcome,
      comparison_or_metric = comparison_or_metric,
      effect_measure = effect_measure,
      estimate = as.numeric(estimate),
      lower_95 = as.numeric(lower_95),
      upper_95 = as.numeric(upper_95),
      p_value = as.numeric(p_value),
      multiplicity_adjusted_p = as.numeric(multiplicity_adjusted_p),
      display_value = display_value,
      destination = destination,
      source_file = source_file,
      allowed_statement_cn = allowed_statement_cn,
      prohibited_statement_cn = prohibited_statement_cn
    )
  )
}

effect_row <- function(data, ...) {
  row <- dplyr::filter(data, ...)
  if (nrow(row) != 1L) stop("Expected exactly one source row.", call. = FALSE)
  row
}

parse_formatted_effect <- function(value) {
  parts <- strsplit(gsub("[()]", "", value), "[[:space:]-]+")[[1]]
  numbers <- suppressWarnings(as.numeric(parts[nzchar(parts)]))
  if (length(numbers) != 3L || any(!is.finite(numbers))) {
    stop("Could not parse formatted effect: ", value, call. = FALSE)
  }
  names(numbers) <- c("estimate", "lower_95", "upper_95")
  numbers
}

participant_row <- effect_row(table1, Characteristic == "Participants, unweighted n")
death_row <- effect_row(table1, Characteristic == "Deaths during follow-up, n (weighted %)")
add_result(
  "NH-COHORT-001", "Descriptive", "NHANES 2011-2018", "Corrected analytical cohort",
  n = participant_row$Overall, events = 831, outcome = "All-cause mortality",
  comparison_or_metric = "Six-feature complete cohort", effect_measure = "n",
  estimate = participant_row$Overall, display_value = "n=4,636; 831 deaths",
  destination = "Methods cohort flow and Table 1",
  source_file = "output/corrected_publication_tables/Table1_corrected_weighted_baseline.csv",
  allowed_statement_cn = "NHANES 六特征完整且具死亡随访资格的队列为4,636人，共831例死亡。",
  prohibited_statement_cn = "不得再使用3,386人作为NHANES最终队列。"
)
for (phenotype in c("P1", "P2", "P3")) {
  add_result(
    paste0("NH-PHENOTYPE-", phenotype), "Descriptive", "NHANES 2011-2018",
    "Frozen robust K=3 phenotype counts", n = participant_row$Overall, events = 831,
    outcome = "All-cause mortality", comparison_or_metric = paste0(phenotype, " count"),
    effect_measure = "n", estimate = as.numeric(participant_row[[phenotype]]),
    display_value = paste0(phenotype, " n=", participant_row[[phenotype]], "; ", death_row[[phenotype]]),
    destination = "Results phenotype description and Table 1",
    source_file = "output/corrected_publication_tables/Table1_corrected_weighted_baseline.csv",
    allowed_statement_cn = paste0(phenotype, "样本量及复杂抽样加权死亡比例按Table 1报告。"),
    prohibited_statement_cn = "不得把聚类变量的组间P值当作独立验证。"
  )
}

nh_primary_p1 <- effect_row(table2, model_role == "Primary", comparison == "P1 vs P3")
nh_primary_p2 <- effect_row(table2, model_role == "Primary", comparison == "P2 vs P3")
for (row in list(nh_primary_p1, nh_primary_p2)) {
  is_p1 <- row$comparison == "P1 vs P3"
  add_result(
    if (is_p1) "NH-PRIMARY-P1" else "NH-SECONDARY-P2",
    if (is_p1) "Primary" else "Secondary", "NHANES 2011-2018",
    row$model, row$n, row$events, "All-cause mortality", row$comparison, "HR",
    row$HR, row$lower_95, row$upper_95, row$p_value, display_value = row$effect_95ci,
    destination = if (is_p1) "Abstract, main Results, Table 2" else "Main Results, Table 2",
    source_file = "output/corrected_publication_tables/Table2_corrected_NHANES_models.csv",
    allowed_statement_cn = if (is_p1) {
      "复杂抽样完全调整后，P1与P3相比全因死亡风险更高，HR 1.68（95% CI 1.30-2.17）。"
    } else {
      "P2与P3的复杂抽样完全调整关联不明确，HR 0.97（95% CI 0.78-1.20）。"
    },
    prohibited_statement_cn = if (is_p1) {
      "不得使用旧aHR 3.08，也不得声称因果关系。"
    } else {
      "不得把P2描述为稳定的独立高风险表型。"
    }
  )
}

for (model_name in c("Complex survey x selection IPW", "Leave-one-cycle-out complex-survey Cox")) {
  row <- effect_row(table2, model == model_name, comparison == "P1 vs P3")
  add_result(
    if (grepl("selection", model_name, ignore.case = TRUE)) "NH-SENS-IPW-P1" else "NH-SENS-LOCO-P1",
    "Key sensitivity", "NHANES 2011-2018", row$model, row$n, row$events,
    "All-cause mortality", row$comparison, "HR", row$HR, row$lower_95, row$upper_95,
    row$p_value, display_value = row$effect_95ci, destination = "Supplementary Results",
    source_file = "output/corrected_publication_tables/Table2_corrected_NHANES_models.csv",
    allowed_statement_cn = if (grepl("selection", model_name, ignore.case = TRUE)) {
      "选择IPW与NHANES抽样权重合并后，P1关联基本不变。"
    } else {
      "留一周期外推后的复杂抽样模型保持P1风险方向。"
    },
    prohibited_statement_cn = "敏感性结果不得取代主要复杂抽样模型。"
  )
}

bootstrap_result <- table2b |>
  dplyr::filter(validation == "Cluster-aware bootstrap P1 HR")
add_result(
  "NH-SENS-BOOTSTRAP-P1", "Key sensitivity", "NHANES 2011-2018",
  "Outcome-blind cluster-aware bootstrap with reclustering", n = 3979, events = 720,
  outcome = "All-cause mortality", comparison_or_metric = "P1 vs P3",
  effect_measure = "Bootstrap HR", estimate = 1.73, lower_95 = 1.35, upper_95 = 2.15,
  display_value = "Median HR 1.73 (percentile 95% interval 1.35-2.15)",
  destination = "Supplementary Results",
  source_file = "output/corrected_publication_tables/Table2B_internal_validation.csv",
  allowed_statement_cn = "每次重聚类的1,000次bootstrap支持P1效应稳定。",
  prohibited_statement_cn = "不得把bootstrap百分位区间写成独立外部验证。"
)
add_result(
  "NH-DIAG-PH", "Model diagnostic", "NHANES 2011-2018",
  "Fully adjusted phenotype Cox PH test", n = 3979, events = 720,
  outcome = "All-cause mortality", comparison_or_metric = "Global PH test",
  effect_measure = "P value", estimate = 0.480, p_value = 0.480, display_value = "Global P=0.480",
  destination = "Methods and Supplementary Results",
  source_file = "output/corrected_publication_tables/Table2B_internal_validation.csv",
  allowed_statement_cn = "表型主要Cox模型未发现明显比例风险假设违背。",
  prohibited_statement_cn = "不得把P>0.05写成比例风险假设已被证明绝对成立。"
)

for (metric_name in c(
  "60-month risk difference, percentage points",
  "RMST difference through 60 months, months"
)) {
  row <- effect_row(table4b, comparison == "P1 vs P3", metric == metric_name)
  add_result(
    if (grepl("risk difference", metric_name)) "NH-CLINICAL-RISK60" else "NH-CLINICAL-RMST60",
    "Clinical interpretation", "NHANES 2011-2018", "Survey-bootstrap adjusted absolute risk and RMST",
    n = 3979, events = 720, outcome = "All-cause mortality", comparison_or_metric = metric_name,
    effect_measure = if (grepl("risk difference", metric_name)) "Percentage-point difference" else "Month difference",
    estimate = row$estimate, lower_95 = row$lower_95, upper_95 = row$upper_95,
    display_value = row$estimate_95ci, destination = "Main or Supplementary clinical interpretation",
    source_file = "output/corrected_publication_tables/Table4B_adjusted_risk_RMST_contrasts.csv",
    allowed_statement_cn = if (grepl("risk difference", metric_name)) {
      "P1相对P3的调整后60个月死亡绝对风险高7.39个百分点。"
    } else {
      "P1相对P3截至60个月的调整后RMST少1.99个月。"
    },
    prohibited_statement_cn = "不得将调整后绝对风险差解释为治疗效应或因果效应。"
  )
}

k3 <- effect_row(stability, k == 3)
add_result(
  "NH-CLUSTER-K3-STABILITY", "Method robustness", "NHANES 2011-2018",
  "K=3 subsampling stability", n = 4636, comparison_or_metric = "Median adjusted Rand index",
  effect_measure = "ARI", estimate = k3$median_ari, lower_95 = k3$q1_ari, upper_95 = k3$q3_ari,
  display_value = sprintf("Median ARI %.3f (IQR %.3f-%.3f)", k3$median_ari, k3$q1_ari, k3$q3_ari),
  destination = "Methods and Supplementary Results",
  source_file = "output/nhanes_cluster_number_validation/Table28B_C_subsampling_stability_summary.csv",
  allowed_statement_cn = "K=3在重抽样中稳定且具有临床可解释性。",
  prohibited_statement_cn = "不得称K=3为所有诊断一致支持的唯一数学最优解。"
)
raw_robust <- effect_row(variant_agreement, reference == "raw6_original", candidate == "log6_winsor")
add_result(
  "NH-CLUSTER-RAW-ROBUST-AGREEMENT", "Limitation", "NHANES 2011-2018",
  "Raw versus robust transformed K=3 agreement", n = 4636,
  comparison_or_metric = "Adjusted Rand index", effect_measure = "ARI",
  estimate = raw_robust$adjusted_rand_index,
  display_value = sprintf("ARI %.3f", raw_robust$adjusted_rand_index),
  destination = "Discussion limitation and Supplementary Results",
  source_file = "output/nhanes_robust_reanalysis/Table28F_variant_agreement.csv",
  allowed_statement_cn = "原始值与稳健变换后的成员一致性较低，表型并非算法无关。",
  prohibited_statement_cn = "不得声称表型完全算法无关或成员完全一致。"
)

for (predictor_name in c("phenotype_total_protein", "GNRI", "domain_balanced_score")) {
  row <- effect_row(cindex, predictor == predictor_name)
  add_result(
    paste0("NH-OOB-CINDEX-", toupper(predictor_name)), "Exploratory", "NHANES 2011-2018",
    "Cluster-aware out-of-bag C-index bootstrap", n = 3979, events = 720,
    outcome = "All-cause mortality", comparison_or_metric = predictor_name,
    effect_measure = "Delta C-index", estimate = row$median_delta_c_index,
    lower_95 = row$lower_95_delta, upper_95 = row$upper_95_delta,
    display_value = sprintf("Median delta %.4f (95%% bootstrap interval %.4f to %.4f)", row$median_delta_c_index, row$lower_95_delta, row$upper_95_delta),
    destination = "Supplementary Results only",
    source_file = "output/nhanes_bootstrap_cindex/Table34B_bootstrap_Cindex_summary.csv",
    allowed_statement_cn = if (predictor_name == "phenotype_total_protein") {
      "总蛋白表型带来幅度较小但重抽样中稳定为正的区分度增量。"
    } else if (predictor_name == "GNRI") {
      "GNRI的区分度增量不低于总蛋白表型。"
    } else {
      "连续结构评分区分度增量较大，但存在明显时间变化效应，仅作探索。"
    },
    prohibited_statement_cn = "不得宣称表型全面优于所有传统营养指标。"
  )
}

for (index in seq_len(nrow(time_effects))) {
  row <- time_effects[index, ]
  add_result(
    paste0("NH-SCORE-TIME-", row$followup_month, "M"), "Exploratory", "NHANES 2011-2018",
    "Time-varying continuous vulnerability score", n = 3979, events = 720,
    outcome = "All-cause mortality", comparison_or_metric = paste0(row$followup_month, " months; per 1 SD"),
    effect_measure = "HR", estimate = row$HR_per_1_SD, lower_95 = row$lower_95,
    upper_95 = row$upper_95, display_value = row$hazard_ratio_95ci,
    destination = "Supplementary Results only",
    source_file = "output/nhanes_corrected_bootstrap/Table32E_time_varying_score_effects.csv",
    allowed_statement_cn = "连续脆弱评分的关联随随访时间减弱。",
    prohibited_statement_cn = "不得用单一恒定HR概括连续评分，也不得取代表型主分析。"
  )
}

for (dataset_name in c("MIMIC-IV", "eICU")) {
  flow <- effect_row(table5a, dataset == dataset_name)
  add_result(
    paste0(ifelse(dataset_name == "MIMIC-IV", "MIMIC", "EICU"), "-FLOW-PRIMARY"),
    "Limitation", dataset_name, "Complete-case selection flow",
    n = flow$denominator_n, comparison_or_metric = "Primary model inclusion",
    effect_measure = "Percent", estimate = flow$primary_model_percent,
    display_value = sprintf("%s/%s (%.2f%%); effective sample size %.1f", flow$primary_model_n, flow$denominator_n, flow$primary_model_percent, flow$effective_sample_size),
    destination = "Methods, flow diagram, and limitations",
    source_file = "output/corrected_publication_tables/Table5A_ICU_selection_flow.csv",
    allowed_statement_cn = if (dataset_name == "MIMIC-IV") {
      "MIMIC主要模型纳入率为3.40%，选择IPW有效样本量约631，必须报告选择限制。"
    } else {
      "eICU主要模型纳入率为18.24%，选择IPW有效样本量约8,283。"
    },
    prohibited_statement_cn = "不得忽略完整病例选择比例或声称排除了未测量选择偏倚。"
  )
}

add_icu_effect <- function(
  result_id, tier, dataset_name, data, model_name, comparison, outcome,
  destination, allowed, prohibited
) {
  row <- effect_row(
    data,
    model == .env$model_name,
    comparison == .env$comparison
  )
  add_result(
    result_id, tier, dataset_name, row$model, row$n, row$events, outcome,
    row$comparison, row$effect_measure, row$estimate, row$lower_95, row$upper_95,
    row$p_value, display_value = row$effect_95ci, destination = destination,
    source_file = if (identical(data, table5b)) {
      "output/corrected_publication_tables/Table5B_ICU_selection_IPW_models.csv"
    } else {
      "output/corrected_publication_tables/Table7D_official_SOFA_models.csv"
    },
    allowed_statement_cn = allowed, prohibited_statement_cn = prohibited
  )
}

add_icu_effect(
  "MIMIC-PRIMARY-P1", "Primary conceptual replication", "MIMIC-IV v3.1", table5b,
  "Official OASIS quartile-stratified Cox", "P1 vs P3", "365-day all-cause mortality",
  "Abstract, main Results, Table 3",
  "严格首24小时队列中，官方OASIS四分位分层后P1与P3相比365天死亡风险更高，HR 1.67（95% CI 1.36-2.05）。",
  "不得写成72小时结果、OASIS-like、原质心严格外部验证或因果效应。"
)
add_icu_effect(
  "MIMIC-SECONDARY-P2", "Secondary", "MIMIC-IV v3.1", table5b,
  "Official OASIS quartile-stratified Cox", "P2 vs P3", "365-day all-cause mortality",
  "Main Results, Table 3", "MIMIC中P2风险升高，但其选择IPW结果不稳定，因此仅作次要比较。",
  "不得把P2称为跨数据库稳定的核心阳性结果。"
)
add_icu_effect(
  "MIMIC-SENS-IPW-P1", "Key sensitivity", "MIMIC-IV v3.1", table5b,
  "Selection-IPW official OASIS quartile-stratified Cox", "P1 vs P3", "365-day all-cause mortality",
  "Supplementary Results", "确定性OASIS四分位及选择IPW后MIMIC P1关联保持，HR 1.71（95% CI 1.29-2.25）。",
  "不得声称选择IPW消除了所有选择偏倚。"
)
add_icu_effect(
  "MIMIC-SENS-IPW-P2", "Key sensitivity", "MIMIC-IV v3.1", table5b,
  "Selection-IPW official OASIS quartile-stratified Cox", "P2 vs P3", "365-day all-cause mortality",
  "Supplementary Results", "MIMIC P2在选择IPW后区间跨1，提示其稳定性不足。",
  "不得隐去P2的衰减或只报告未加权结果。"
)

for (analysis_name in c("365-day mortality, 24-hour landmark", "In-hospital mortality")) {
  row <- effect_row(table3b, database == "MIMIC-IV", analysis == analysis_name)
  parsed <- parse_formatted_effect(row$P1_vs_P3)
  add_result(
    if (grepl("landmark", analysis_name)) "MIMIC-SENS-LANDMARK-P1" else "MIMIC-SENS-HOSPITAL-P1",
    "Key sensitivity", "MIMIC-IV v3.1", analysis_name, row$n, row$events,
    if (grepl("landmark", analysis_name)) "365-day all-cause mortality" else "In-hospital mortality",
    "P1 vs P3", row$effect_measure, parsed[["estimate"]], parsed[["lower_95"]], parsed[["upper_95"]],
    display_value = row$P1_vs_P3,
    destination = "Supplementary Results",
    source_file = "output/corrected_publication_tables/Table3B_cross_database_sensitivity.csv",
    allowed_statement_cn = if (grepl("landmark", analysis_name)) {
      "24小时landmark分析保持MIMIC P1风险方向。"
    } else {
      "官方OASIS调整后MIMIC P1与院内死亡关联，OR 2.25（95% CI 1.56-3.23）。"
    },
    prohibited_statement_cn = "敏感性结局不得取代365天死亡主结局。"
  )
}

for (model_name in c(
  "Official first-day SOFA quartile-stratified Cox",
  "Selection-IPW official first-day SOFA quartile-stratified Cox"
)) {
  add_icu_effect(
    if (grepl("Selection", model_name)) "MIMIC-SENS-SOFA-IPW-P1" else "MIMIC-SENS-SOFA-P1",
    "Severity sensitivity", "MIMIC-IV v3.1", table7d, model_name, "P1 vs P3",
    "365-day all-cause mortality", "Supplementary Results",
    if (grepl("Selection", model_name)) {
      "选择IPW并按官方首日SOFA四分位分层后，P1关联仍高于1。"
    } else {
      "按官方首日SOFA四分位分层后，P1 HR为1.57（95% CI 1.28-1.93）。"
    },
    "SOFA含肾脏分量，属于可能过度校正的敏感性，不替代官方OASIS主模型。"
  )
}

add_icu_effect(
  "EICU-PRIMARY-P1", "Primary conceptual replication", "eICU v2.0", table5b,
  "Frozen APACHE-adjusted mixed model", "P1 vs P3", "In-hospital mortality",
  "Abstract, main Results, Table 3",
  "严格首24小时、APACHE调整并纳入医院随机截距后，P1与P3相比院内死亡OR为1.28（95% CI 1.11-1.48）。",
  "不得使用旧OR 2.07，不得写成所有医院一致显著。"
)
add_icu_effect(
  "EICU-SECONDARY-P2", "Secondary", "eICU v2.0", table5b,
  "Frozen APACHE-adjusted mixed model", "P2 vs P3", "In-hospital mortality",
  "Main Results, Table 3", "eICU中P2与P3的调整关联不明确。",
  "不得把P2描述为跨数据库稳定阳性。"
)
add_icu_effect(
  "EICU-SENS-IPW-P1", "Key sensitivity", "eICU v2.0", table5b,
  "Selection-IPW APACHE-adjusted mixed model", "P1 vs P3", "In-hospital mortality",
  "Supplementary Results", "eICU选择IPW后P1方向一致但置信区间跨1，OR 1.13（95% CI 0.98-1.31）。",
  "不得隐去选择IPW后的衰减或称完全稳健。"
)

for (analysis_name in c("ICU mortality", "In-hospital mortality, 24-hour landmark")) {
  row <- effect_row(table3b, database == "eICU", analysis == analysis_name)
  parsed <- parse_formatted_effect(row$P1_vs_P3)
  add_result(
    if (analysis_name == "ICU mortality") "EICU-SENS-ICU-P1" else "EICU-SENS-LANDMARK-P1",
    "Key sensitivity", "eICU v2.0", analysis_name, row$n, row$events,
    if (analysis_name == "ICU mortality") "ICU mortality" else "In-hospital mortality",
    "P1 vs P3", row$effect_measure, parsed[["estimate"]], parsed[["lower_95"]], parsed[["upper_95"]],
    display_value = row$P1_vs_P3,
    destination = "Supplementary Results",
    source_file = "output/corrected_publication_tables/Table3B_cross_database_sensitivity.csv",
    allowed_statement_cn = "eICU替代结局或24小时landmark分析保持P1风险方向。",
    prohibited_statement_cn = "敏感性结局不得取代院内死亡主结局。"
  )
}

meta <- table8a[1, ]
add_result(
  "EICU-HOSPITAL-META-P1", "Transportability audit", "eICU v2.0",
  meta$model, n = meta$participants_p1_p3, events = meta$deaths_p1_p3,
  outcome = "In-hospital mortality", comparison_or_metric = "P1 vs P3 across 17 high-information hospitals",
  effect_measure = "Random-effects OR", estimate = meta$pooled_or, lower_95 = meta$lower_95,
  upper_95 = meta$upper_95, p_value = meta$p_value, display_value = paste0(meta$effect_95ci, "; prediction interval ", meta$prediction_interval),
  destination = "Supplementary Results and Discussion limitation",
  source_file = "output/corrected_publication_tables/Table8A_eICU_hospital_meta_analysis.csv",
  allowed_statement_cn = "医院特异汇总方向大于1，但置信区间和预测区间均跨1，提示中心间精度和可迁移性有限。",
  prohibited_statement_cn = "不得声称各医院均一、全部显著或完全可迁移。"
)
direction <- table8e[1, ]
add_result(
  "EICU-HOSPITAL-DIRECTION", "Transportability audit", "eICU v2.0",
  "High-information hospital effect directions", comparison_or_metric = "Hospitals with OR above one",
  effect_measure = "Count", estimate = direction$hospitals_or_above_one,
  display_value = sprintf("%d/%d hospitals above one", direction$hospitals_or_above_one, direction$estimable_hospitals),
  destination = "Supplementary Results",
  source_file = "output/corrected_publication_tables/Table8E_eICU_hospital_effect_directions.csv",
  allowed_statement_cn = "17家高信息量医院中13家P1调整OR大于1，但多数单中心区间跨1。",
  prohibited_statement_cn = "不得按显著性筛选医院或只展示阳性中心。"
)

albumin_models <- list(
  list("ALBUMIN-NH-P1", "NHANES", "Complex-survey albumin phenotype Cox"),
  list("ALBUMIN-MIMIC-IPW-P1", "MIMIC-IV", "Selection-IPW official OASIS albumin phenotype Cox"),
  list("ALBUMIN-EICU-IPW-P1", "eICU", "Selection-IPW APACHE-adjusted albumin phenotype mixed model")
)
for (specification in albumin_models) {
  row <- effect_row(table6a, dataset == specification[[2]], model == specification[[3]], comparison == "P1 vs P3")
  add_result(
    specification[[1]], "Harmonised biomarker sensitivity", specification[[2]], row$model,
    row$n, row$events, if (specification[[2]] == "NHANES") "All-cause mortality" else if (specification[[2]] == "MIMIC-IV") "365-day all-cause mortality" else "In-hospital mortality",
    row$comparison, row$effect_measure, row$estimate, row$lower_95, row$upper_95, row$p_value,
    display_value = row$effect_95ci, destination = "Supplementary Results",
    source_file = "output/corrected_publication_tables/Table6A_harmonised_albumin_models.csv",
    allowed_statement_cn = "统一albumin变量口径的敏感性中，P1方向在三库保持一致。",
    prohibited_statement_cn = "不得把albumin版与总蛋白版称为完全相同表型或替代NHANES主分析。"
  )
}

for (index in seq_len(nrow(table9a))) {
  row <- table9a[index, ]
  prefix <- if (grepl("NHANES", row$dataset)) "NH" else if (grepl("MIMIC", row$dataset)) "MIMIC" else "EICU"
  add_result(
    paste0(prefix, "-INTERACTION-", toupper(row$modifier)), "Exploratory interaction",
    row$dataset, row$model, row$n, row$events,
    if (prefix == "NH") "All-cause mortality" else if (prefix == "MIMIC") "365-day all-cause mortality" else "In-hospital mortality",
    paste0("Global phenotype-by-", tolower(row$modifier), " interaction"), "Global interaction P value",
    estimate = row$wald_chisq, p_value = row$p_value,
    multiplicity_adjusted_p = row$p_value_bh,
    display_value = sprintf("raw P=%.3f; BH P=%.3f", row$p_value, row$p_value_bh),
    destination = "Supplementary Results only",
    source_file = "output/corrected_publication_tables/Table9A_age_sex_global_interactions.csv",
    allowed_statement_cn = "六项年龄/性别全局交互均未通过BH校正，未观察到明确效应修饰证据。",
    prohibited_statement_cn = "不得因某个条件效应显著或MIMIC年龄原始P<0.05而宣称存在确定亚组差异。"
  )
}

for (index in seq_len(nrow(table9c))) {
  row <- table9c[index, ]
  prefix <- if (grepl("NHANES", row$dataset)) "NH" else if (grepl("MIMIC", row$dataset)) "MIMIC" else "EICU"
  level_id <- toupper(gsub("[^A-Za-z0-9]+", "-", row$level))
  add_result(
    paste0(prefix, "-CONDITIONAL-", toupper(row$modifier), "-", level_id),
    "Exploratory conditional estimate", row$dataset, paste0(row$modifier, " interaction model"),
    row$n, row$events,
    if (prefix == "NH") "All-cause mortality" else if (prefix == "MIMIC") "365-day all-cause mortality" else "In-hospital mortality",
    paste0(row$comparison, " at ", row$level), row$effect_measure,
    row$estimate, row$lower_95, row$upper_95, display_value = row$effect_95ci,
    destination = "Supplementary Results only",
    source_file = "output/corrected_publication_tables/Table9C_P1_conditional_effects.csv",
    allowed_statement_cn = "预设年龄点和性别条件估计仅用于解释交互模型，方向总体一致。",
    prohibited_statement_cn = "不得用一个条件估计显著、另一个不显著来证明亚组间差异。"
  )
}

dictionary <- dictionary |>
  dplyr::arrange(
    factor(tier, levels = c(
      "Primary", "Primary conceptual replication", "Secondary", "Descriptive",
      "Key sensitivity", "Severity sensitivity", "Clinical interpretation",
      "Harmonised biomarker sensitivity", "Method robustness", "Model diagnostic",
      "Transportability audit", "Limitation", "Exploratory", "Exploratory interaction",
      "Exploratory conditional estimate"
    )),
    dataset, result_id
  )

legacy_blacklist <- tibble::tribble(
  ~legacy_item, ~required_replacement, ~reason,
  "NHANES final n=3,386", "Corrected six-feature cohort n=4,636; primary model n=3,979", "Old mortality-linkage restriction was incorrect",
  "NHANES P1 aHR=3.08", "Primary survey-weighted HR=1.68 (95% CI 1.30-2.17)", "Old cohort and preprocessing",
  "MIMIC 72-hour P1 aHR=3.21", "Strict 0-24-hour official OASIS HR=1.67 (95% CI 1.36-2.05)", "Old window and model",
  "MIMIC OASIS-like primary HR=1.64", "MIT-LCP v3.0.1 official OASIS primary HR=1.67 (95% CI 1.36-2.05)", "Official severity score supersedes local approximation",
  "eICU 72-hour P1 OR=2.07", "Strict 0-24-hour APACHE mixed-model OR=1.28 (95% CI 1.11-1.48)", "Old window and insufficient centre adjustment",
  "Pure nutritional-inflammatory phenotype", "Immunonutritional-renal vulnerability phenotype", "Creatinine axis materially contributes",
  "P1 is the highest-inflammation phenotype", "P1 has elevated inflammation with greatest renal-haematological vulnerability", "NHANES P2 has higher NLR/SII than P1",
  "Perfect or strict external validation", "Conceptual replication or cross-setting reproducibility", "ICU phenotypes are de novo and settings differ",
  "Identical phenotype membership across databases", "Related biological structure with non-identical membership", "Different biomarkers and distributions",
  "K=3 is the unique mathematical optimum", "K=3 is stable, clinically interpretable, and pre-specified, but not uniquely optimal", "Gap statistic does not uniquely select K=3",
  "Phenotype is algorithm-independent", "High-risk structure is recoverable under selected alternatives", "Raw-versus-robust ARI is low",
  "Phenotype outperforms all conventional indices", "Phenotype adds modest information; GNRI performs at least comparably", "OOB C-index results",
  "Obesity paradox is proven", "No obesity-paradox claim", "BMI pattern alone cannot establish the paradox",
  "All eICU hospitals show the same significant effect", "Direction is not driven by one hospital, but centre-specific intervals are imprecise", "Meta-analysis CI and prediction interval cross one",
  "MIMIC age interaction is significant", "No age/sex interaction survives six-test BH correction", "Raw age interaction P=0.040 but BH P=0.243",
  "Pool NHANES/MIMIC HRs with eICU OR", "Report database-specific HRs and ORs without a common-effect pool", "Effect measures and settings differ"
)

write_markdown_table <- function(data, path) {
  display <- data
  display[] <- lapply(display, function(value) {
    value <- as.character(value)
    value[is.na(value)] <- ""
    gsub("\\|", "\\\\|", value)
  })
  header <- paste0("| ", paste(names(display), collapse = " | "), " |")
  separator <- paste0("| ", paste(rep("---", ncol(display)), collapse = " | "), " |")
  rows <- apply(display, 1, function(row) paste0("| ", paste(row, collapse = " | "), " |"))
  writeLines(c(header, separator, rows), path, useBytes = TRUE)
}

readr::write_csv(dictionary, file.path(output_dir, "UNIQUE_RESULT_DICTIONARY.csv"))
dictionary_compact <- dictionary |>
  dplyr::select(
    result_id, tier, dataset, analysis, n, events, comparison_or_metric,
    effect_measure, display_value, destination, source_file, allowed_statement_cn
  )
write_markdown_table(dictionary_compact, file.path(output_dir, "UNIQUE_RESULT_DICTIONARY.md"))
readr::write_csv(legacy_blacklist, file.path(output_dir, "FORBIDDEN_LEGACY_RESULTS.csv"))
write_markdown_table(legacy_blacklist, file.path(output_dir, "FORBIDDEN_LEGACY_RESULTS.md"))

source_files <- sort(unique(dictionary$source_file))
source_paths <- file.path(root, source_files)
source_manifest <- tibble::tibble(
  source_file = source_files,
  exists = file.exists(source_paths),
  bytes = ifelse(file.exists(source_paths), file.info(source_paths)$size, NA_real_),
  md5 = ifelse(file.exists(source_paths), unname(tools::md5sum(source_paths)), NA_character_)
)
readr::write_csv(source_manifest, file.path(output_dir, "RESULT_SOURCE_MANIFEST.csv"))
write_markdown_table(source_manifest, file.path(output_dir, "RESULT_SOURCE_MANIFEST.md"))

nh_p1 <- dictionary |> dplyr::filter(result_id == "NH-PRIMARY-P1")
mimic_p1 <- dictionary |> dplyr::filter(result_id == "MIMIC-PRIMARY-P1")
eicu_p1 <- dictionary |> dplyr::filter(result_id == "EICU-PRIMARY-P1")
mimic_ipw <- dictionary |> dplyr::filter(result_id == "MIMIC-SENS-IPW-P1")
eicu_ipw <- dictionary |> dplyr::filter(result_id == "EICU-SENS-IPW-P1")

freeze_lines <- c(
  "# Final Statistical Analysis Freeze v2.0",
  "",
  "冻结日期：2026-07-10",
  "",
  "## 一、版本地位",
  "",
  "本文件完全取代 `ANALYSIS_FREEZE_V1_2026-07-10.md` 及其后所有零散口头数字。V1 中 MIMIC OASIS-like 主模型已由 MIT-LCP `mimic-code v3.0.1` 官方 OASIS 模型取代。除发现程序错误、数据错误或审稿人提出必要要求外，不再新增结局驱动模型、亚组或阈值搜索。",
  "",
  "所有稿件数字必须从 `output/final_analysis_freeze_v2/UNIQUE_RESULT_DICTIONARY.csv` 提取；不得凭记忆、旧Word稿或旧截图录入。",
  "",
  "## 二、冻结的研究定位",
  "",
  "研究定位为：在 NHANES 社区老年人中发现免疫营养-肾脏脆弱表型，并在 MIMIC-IV 和 eICU 重症老年人中进行跨场景概念性复现。",
  "",
  "- 允许：`conceptual replication`、`cross-setting reproducibility`、`biological structure re-emerged`。",
  "- 不允许：原质心严格外部验证、完美复现、完全相同表型、算法无关。",
  "- P1 对 P3 是核心比较；P2 对 P3 是次要比较。",
  "- NHANES/MIMIC 报告 HR，eICU 报告 OR，不进行共同效应合并。",
  "",
  "## 三、冻结的表型与命名",
  "",
  "NHANES 主表型使用 NLR、SII、血红蛋白、血清总蛋白、BMI 和肌酐，经1%/99% winsorisation、NLR/SII/肌酐自然对数转换、Z标准化后进行 K-means（K=3，nstart=100，iter.max=500）。",
  "",
  "- P1：`renal-haematological vulnerability with elevated inflammation`，中文为“炎症升高的肾脏-造血脆弱组”。",
  "- P2：`inflammation-dominant lower-reserve phenotype`，中文为“炎症占优的低储备组”。",
  "- P3：`low-inflammation reference phenotype`，中文为“低炎症参考组”。",
  "",
  "不得把P1写成炎症绝对最高组；NHANES中P2的NLR/SII高于P1。K=3是稳定、可解释且预先冻结的方案，不是唯一数学最优解。",
  "",
  "## 四、主要推断结果",
  "",
  "| 数据库 | 主样本 | 结局 | 冻结模型 | P1 vs P3 | 解释层级 |",
  "|---|---:|---|---|---|---|",
  sprintf("| NHANES 2011-2018 | %s（%s例死亡） | 全因死亡 | 复杂抽样完全调整 Cox | HR %s | 主要推断 |", nh_p1$n, nh_p1$events, nh_p1$display_value),
  sprintf("| MIMIC-IV v3.1 | %s（%s例死亡） | 365天全因死亡 | 官方OASIS四分位分层 Cox | HR %s | 主要概念复现 |", mimic_p1$n, mimic_p1$events, mimic_p1$display_value),
  sprintf("| eICU v2.0 | %s（%s例死亡） | 院内死亡 | APACHE调整、医院随机截距 logistic | OR %s | 多中心概念复现 |", eicu_p1$n, eicu_p1$events, eicu_p1$display_value),
  "",
  "NHANES六特征完整队列为4,636人、831例死亡；完全调整共同样本为3,979人、720例死亡。MIMIC和eICU均使用ICU入科后0-24小时窗口。",
  "",
  "## 五、必须同时报告的稳健性边界",
  "",
  sprintf("- NHANES选择IPW与留一周期外推均保持P1方向；1,000次每次重聚类bootstrap中位HR为1.73（1.35-2.15）。"),
  sprintf("- MIMIC选择IPW后P1为HR %s；完整病例纳入率仅3.40%%，必须报告。", mimic_ipw$display_value),
  sprintf("- eICU选择IPW后P1为OR %s，区间跨1；不能写成完全稳健。", eicu_ipw$display_value),
  sprintf("- eICU 17家高信息量医院随机效应汇总OR为%s，预测区间%s；不能写成医院间均一。", meta$effect_95ci, meta$prediction_interval),
  "- 官方SOFA、24小时landmark、统一albumin变量口径属于敏感性分析，不替代三库各自主模型。",
  "- 连续结构评分存在明显时间变化效应，只能作为探索性补充。",
  "",
  "## 六、临床解释结果",
  "",
  "NHANES中，P1相对P3的调整后60个月死亡绝对风险差为7.39个百分点（95% CI 4.74-10.34），截至60个月RMST少1.99个月（95% CI 1.26-2.83个月）。这些是模型化临床解释，不是因果治疗效应。",
  "",
  "## 七、效应修饰",
  "",
  "仅检验年龄和性别。三数据库共六项2自由度全局交互统一BH校正后均不显著。MIMIC年龄交互原始P=0.040、BH P=0.243，只能作为提示，不得宣称确定年龄亚组差异。",
  "",
  "## 八、报告禁区",
  "",
  paste0("- ", legacy_blacklist$legacy_item, " -> ", legacy_blacklist$required_replacement),
  "",
  "## 九、结果与来源管理",
  "",
  "- `UNIQUE_RESULT_DICTIONARY.csv`：稿件唯一取数来源。",
  "- `FORBIDDEN_LEGACY_RESULTS.csv`：旧结果和过强措辞黑名单。",
  "- `RESULT_SOURCE_MANIFEST.csv`：每个来源文件的大小和MD5。",
  "- `FREEZE_V2_QA_REPORT.csv`：冻结一致性自动校验。",
  "",
  "后续可以重写稿件和调整呈现顺序，但不得更换冻结的主要样本、时间窗、结局、模型、表型命名边界或主要效应值。"
)
writeLines(freeze_lines, file.path(root, "FINAL_ANALYSIS_FREEZE_V2_2026-07-10.md"), useBytes = TRUE)

checks <- list()
add_check <- function(check, passed, detail) {
  checks[[length(checks) + 1L]] <<- tibble::tibble(
    check = check, passed = isTRUE(passed), detail = as.character(detail)
  )
}
add_check("Result IDs are unique", !anyDuplicated(dictionary$result_id), paste0(nrow(dictionary), " rows"))
add_check("All result source files exist", all(source_manifest$exists), paste(source_manifest$source_file[!source_manifest$exists], collapse = "; "))
add_check("All source files have MD5 hashes", all(nchar(source_manifest$md5) == 32L), paste0(nrow(source_manifest), " sources"))
add_check("Required dictionary fields are complete", all(nzchar(dictionary$result_id)) && all(nzchar(dictionary$tier)) && all(nzchar(dictionary$source_file)) && all(nzchar(dictionary$allowed_statement_cn)), "Required text fields")
effect_rows <- dictionary |> dplyr::filter(!is.na(lower_95), !is.na(estimate), !is.na(upper_95))
add_check("All numeric intervals are ordered", all(effect_rows$lower_95 <= effect_rows$estimate & effect_rows$estimate <= effect_rows$upper_95), paste0(nrow(effect_rows), " intervals"))
add_check("NHANES primary result is frozen", nrow(nh_p1) == 1L && abs(nh_p1$estimate - 1.6830536847349045) < 1e-12 && nh_p1$n == 3979, nh_p1$display_value)
add_check("MIMIC primary result uses deterministic official OASIS quartiles", nrow(mimic_p1) == 1L && abs(mimic_p1$estimate - 1.6670612438483148) < 1e-12 && !grepl("OASIS-like", mimic_p1$analysis), mimic_p1$display_value)
add_check("eICU primary result preserves OR", nrow(eicu_p1) == 1L && eicu_p1$effect_measure == "OR" && abs(eicu_p1$estimate - 1.281767915094283) < 1e-12, eicu_p1$display_value)
add_check("Primary P1 results use expected model cohorts", identical(c(nh_p1$n, mimic_p1$n, eicu_p1$n), c(3979L, 1145L, 12548L)), paste(c(nh_p1$n, mimic_p1$n, eicu_p1$n), collapse = ", "))
add_check("All six interaction tests remain BH nonsignificant", nrow(table9a) == 6L && all(table9a$p_value_bh > 0.05), paste0("minimum BH P=", round(min(table9a$p_value_bh), 3)))
add_check("Legacy blacklist covers core obsolete estimates", all(c("NHANES P1 aHR=3.08", "MIMIC 72-hour P1 aHR=3.21", "eICU 72-hour P1 OR=2.07") %in% legacy_blacklist$legacy_item), paste0(nrow(legacy_blacklist), " banned items"))
legacy_point_estimate <- grepl("^(3\\.08|3\\.21|2\\.07) \\(", dictionary$display_value)
legacy_sample_size <- grepl("n=3,386", dictionary$display_value, fixed = TRUE)
add_check(
  "Dictionary contains no obsolete core point estimates",
  !any(legacy_point_estimate | legacy_sample_size),
  "Legacy point estimates and sample size absent; confidence-limit values are not misclassified"
)
add_check("P2 is not assigned a primary tier", !any(grepl("P2", dictionary$result_id) & dictionary$tier %in% c("Primary", "Primary conceptual replication")), "P2 retained as secondary")
add_check("Freeze v2 file exists", file.exists(file.path(root, "FINAL_ANALYSIS_FREEZE_V2_2026-07-10.md")), "Root freeze document")
add_check("Freeze v2 explicitly supersedes v1", any(grepl("完全取代", freeze_lines)) && any(grepl("官方 OASIS|官方OASIS", freeze_lines)), "Supersession and official OASIS recorded")

qa <- dplyr::bind_rows(checks)
readr::write_csv(qa, file.path(output_dir, "FREEZE_V2_QA_REPORT.csv"))
qa_lines <- c(
  "Final analysis freeze v2 QA",
  "",
  paste0("Checks passed: ", sum(qa$passed), "/", nrow(qa)),
  "",
  paste(capture.output(print(qa)), collapse = "\n")
)
writeLines(qa_lines, file.path(output_dir, "FREEZE_V2_QA_REPORT.txt"), useBytes = TRUE)
cat(paste(qa_lines, collapse = "\n"), "\n")

if (!all(qa$passed)) {
  stop("Final analysis freeze v2 QA failed.", call. = FALSE)
}
