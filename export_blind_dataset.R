## =============================================================================
## export_blind_dataset.R
##
## Writes the packet for the independent Figure 1 coding exercise:
##   blind_coding_dataset.csv   one row per participant, goal columns REMOVED
##   blind_coding_dictionary.csv  what each column is, and what NHANES lacks
##
## Giving the coder the derived clinical variables rather than raw NHANES files
## is deliberate. It isolates the step we want tested -- the reading of Figure 1
## -- from cohort construction and variable derivation, which are not in
## question. A disagreement in the returned file is then attributable to the
## reading and to nothing else.
##
## RUN AFTER nhanes_apob_analysis.R, in the same session.
## =============================================================================

.need <- c("primary_df", "secondary_df")
.missing <- .need[!vapply(.need, exists, logical(1), envir = globalenv())]
if (length(.missing)) {
  stop("\n\n export_blind_dataset.R cannot run. Missing: ",
       paste(.missing, collapse = ", "),
       "\n Run source(\"nhanes_apob_analysis.R\") first, in this session.\n\n",
       call. = FALSE)
}

suppressPackageStartupMessages(library(dplyr))

## Columns the coder needs to evaluate Figure 1's conditions, and nothing that
## reveals how we resolved them.
keep <- c(
  "SEQN", "cycle",
  "RIDAGEYR", "RIAGENDR",
  "ascvd_hx", "MCQ160C", "MCQ160D", "MCQ160E", "MCQ160F",
  "diabetes", "ckd", "egfr", "uacr",
  "avg_sbp", "on_antihtn", "current_smoker", "on_statin",
  "LBXTC", "LBDHDD", "LBXTR", "ldl_primary", "LBXAPB",
  "BMXBMI", "prevent_ascvd_10yr"
)

blind <- bind_rows(
  primary_df   %>% mutate(cohort = "primary prevention"),
  secondary_df %>% mutate(cohort = "secondary prevention")
) %>%
  select(cohort, any_of(keep)) %>%
  arrange(cohort, SEQN)

## Guard: nothing that encodes our reading may leave this script.
forbidden <- c("apob_goal", "ldl_goal", "has_apob_goal", "above_apob_goal",
               "at_ldl_goal", "discordant", "goal_cell", "goal_row",
               "risk_tier", "g_ascvd", "g_dm", "g_shc", "g_htg", "g_pp",
               "firstmatch", "lenient", "severe_hc", "has_ascvd_rf",
               "ascvd_rf_count", "ascvd_very_high", "tg_in_goal_band",
               "dm_risk_modifier", "htn_flag", "n_rows_matched")
leak <- names(blind)[vapply(names(blind), function(n)
  any(vapply(forbidden, function(f) grepl(f, n, fixed = TRUE), logical(1))),
  logical(1))]
if (length(leak)) {
  stop("ABORTED: these columns would reveal our operationalisation: ",
       paste(leak, collapse = ", "), call. = FALSE)
}

write.csv(blind, "blind_coding_dataset.csv", row.names = FALSE)

dict <- tibble::tribble(
  ~column,               ~meaning,
  "cohort",              "Primary prevention (no established ASCVD) or secondary prevention",
  "SEQN",                "NHANES participant identifier",
  "cycle",               "Survey cycle: D=2005-06, E=2007-08, F=2009-10, G=2011-12, H=2013-14, I=2015-16",
  "RIDAGEYR",            "Age in years",
  "RIAGENDR",            "1 = male, 2 = female",
  "ascvd_hx",            "1 if any of MCQ160C, D, E, F is yes",
  "MCQ160C",             "Ever told you had coronary heart disease (1 = yes)",
  "MCQ160D",             "Ever told you had angina (1 = yes)",
  "MCQ160E",             "Ever told you had a heart attack (1 = yes)",
  "MCQ160F",             "Ever told you had a stroke (1 = yes)",
  "diabetes",            "1 if self-reported diagnosis or HbA1c >= 6.5%",
  "ckd",                 "1 if eGFR < 60 or UACR >= 30",
  "egfr",                "mL/min/1.73 m2, CKD-EPI 2021 without race coefficient",
  "uacr",                "Urine albumin-to-creatinine ratio, mg/g",
  "avg_sbp",             "Mean systolic blood pressure, mm Hg, across available readings",
  "on_antihtn",          "1 if a prescription antihypertensive was recorded",
  "current_smoker",      "1 if currently smoking every day or some days",
  "on_statin",           "1 if a prescription statin was recorded",
  "LBXTC",               "Total cholesterol, mg/dL",
  "LBDHDD",              "HDL cholesterol, mg/dL",
  "LBXTR",               "Triglycerides, mg/dL",
  "ldl_primary",         "LDL-C, mg/dL, Sampson/NIH equation 2",
  "LBXAPB",              "Apolipoprotein B, mg/dL, calibrated to a common assay scale",
  "BMXBMI",              "Body mass index, kg/m2",
  "prevent_ascvd_10yr",  "10-year total ASCVD risk, AHA PREVENT equations, as a proportion (0-1)"
)

notavail <- tibble::tribble(
  ~item, ~status,
  "Coronary artery calcium score",              "Not collected by NHANES in any cycle",
  "Acute coronary syndrome within 12 months",   "Timing of events is not collected",
  "Symptomatic peripheral arterial disease",    "Not collected",
  "Coronary revascularisation",                 "Not collected",
  "History of heart failure",                   "Not collected",
  "LDL-C above 100 mg/dL on maximally tolerated statin plus ezetimibe", "Drug dose and tolerance are not recorded",
  "Genetic confirmation of familial hypercholesterolaemia", "Not available",
  "Diabetes duration, retinopathy, neuropathy", "Not collected"
)

write.csv(dict,     "blind_coding_dictionary.csv",   row.names = FALSE)
write.csv(notavail, "blind_coding_unavailable.csv",  row.names = FALSE)

message("\n=== BLIND CODING PACKET ===")
message("  blind_coding_dataset.csv     ", nrow(blind), " rows x ", ncol(blind), " columns")
message("  blind_coding_dictionary.csv  ", nrow(dict), " variables described")
message("  blind_coding_unavailable.csv ", nrow(notavail), " items NHANES does not collect")
message("\n  No column encodes our reading of Figure 1 (checked against ",
        length(forbidden), " forbidden patterns).")
message("  Send these three files with BLIND_CODING_BRIEF.md and the guideline PDF.")
message("  Do NOT send the supplement, the manuscript, or any analysis script.\n")
