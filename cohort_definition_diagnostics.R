## =============================================================================
## cohort_definition_diagnostics.R
##
## What the LBDLDL eligibility filter excludes and why, and how the
## hypertension definition behaves under first-reading versus averaged blood
## pressure. Reads only; writes one CSV. Changes nothing about the analysis.
##
## RUN AFTER nhanes_apob_analysis.R, in the same session.
## =============================================================================

suppressPackageStartupMessages({ library(dplyr); library(survey) })

.need <- c("primary_raw", "primary_df", "primary_design")
.missing <- .need[!vapply(.need, exists, logical(1), envir = globalenv())]
if (length(.missing)) {
  stop("\n\n cohort_definition_diagnostics.R cannot run. Missing: ",
       paste(.missing, collapse = ", "),
       "\n Run source(\"nhanes_apob_analysis.R\") first, in this session.\n\n",
       call. = FALSE)
}

out <- list()
add <- function(finding, metric, value) {
  out[[length(out) + 1]] <<- data.frame(finding = finding, metric = metric,
                                        value = as.character(value))
}

message("\n=== REVIEW DIAGNOSTICS ===")

## -----------------------------------------------------------------------------
## FINDING 1. What does requiring non-missing LBDLDL cost?
##
## NHANES blanks LBDLDL when TG > 400 mg/dL, because Friedewald is invalid
## there. The cohort filter requires !is.na(LBDLDL), so those adults are gone --
## even though Sampson (the primary equation) stays valid to roughly TG 800 and
## the guideline's primary-prevention band runs to 499.
## -----------------------------------------------------------------------------
message("\n-- F1. Cost of the !is.na(LBDLDL) cohort filter --")

base <- derive_vars(primary_raw) %>%
  filter(RIDAGEYR >= CONFIG$age_min, RIDAGEYR <= CONFIG$age_max,
         pregnant == 0, ascvd_hx == 0,
         !is.na(LBXAPB), !is.na(LBXTR), WTSAF2YR > 0)

n_base      <- nrow(base)
n_kept      <- sum(!is.na(base$LBDLDL))
lost        <- base %>% filter(is.na(LBDLDL))
n_lost      <- nrow(lost)
lost_tg400  <- sum(lost$LBXTR > 400, na.rm = TRUE)
lost_sampok <- sum(!is.na(lost$ldl_sampson))
# how many of the lost sit inside the guideline's stated 150-499 band
lost_inband <- sum(lost$LBXTR >= 150 & lost$LBXTR <= 499, na.rm = TRUE)

message("  Eligible before the LBDLDL requirement : ", n_base)
message("  Retained                               : ", n_kept)
message("  DROPPED by !is.na(LBDLDL)              : ", n_lost,
        "  (", sprintf("%.1f%%", 100*n_lost/n_base), ")")
message("    of those, TG > 400 mg/dL             : ", lost_tg400)
message("    of those, Sampson LDL-C computable   : ", lost_sampok)
message("    of those, TG inside the 150-499 band : ", lost_inband)
if (n_lost > 0) {
  message("  Median TG among the dropped            : ",
          sprintf("%.0f mg/dL", median(lost$LBXTR, na.rm = TRUE)))
  message("  Median apoB among the dropped          : ",
          sprintf("%.1f mg/dL", median(lost$LBXAPB, na.rm = TRUE)))
  message("  Median apoB among the retained         : ",
          sprintf("%.1f mg/dL", median(base$LBXAPB[!is.na(base$LBDLDL)], na.rm = TRUE)))
}
message("  Max TG in the analytic cohort          : ",
        sprintf("%.0f mg/dL", max(primary_df$LBXTR, na.rm = TRUE)))
message("  --> The manuscript states a 150-499 band throughout. Compare that",
        "\n      against the maximum above before deciding what to disclose.")

add("F1 LBDLDL filter", "eligible_before_filter", n_base)
add("F1 LBDLDL filter", "dropped_by_filter", n_lost)
add("F1 LBDLDL filter", "dropped_pct", sprintf("%.2f", 100*n_lost/n_base))
add("F1 LBDLDL filter", "dropped_with_TG_over_400", lost_tg400)
add("F1 LBDLDL filter", "dropped_recoverable_by_sampson", lost_sampok)
add("F1 LBDLDL filter", "dropped_inside_150_499_band", lost_inband)
add("F1 LBDLDL filter", "max_TG_in_analytic_cohort", round(max(primary_df$LBXTR, na.rm = TRUE)))

## -----------------------------------------------------------------------------
## FINDING 2. Treated hypertension is invisible to htn_flag.
##
## htn_flag is BPXSY1 >= 140 only. on_antihtn is loaded and is already passed to
## PREVENT as bp_tx, so the risk model treats these people as hypertensive while
## the guideline's risk-factor count does not. htn_flag feeds ascvd_rf_count
## (which sets goal STRINGENCY) and ascvd_highrisk_cond_count (which sets
## very-high-risk ASCVD status).
## -----------------------------------------------------------------------------
message("\n-- F2. Treated hypertension not counted as a risk factor --")

pd <- primary_df
htn_now <- pd$htn_flag
htn_alt <- as.integer(pd$htn_flag == 1 | pd$on_antihtn == 1)

message("  Flagged hypertensive now (SBP1 >=140)        : ", sum(htn_now),
        "  (", sprintf("%.1f%%", 100*mean(htn_now)), ")")
message("  Flagged if treated hypertension also counted : ", sum(htn_alt),
        "  (", sprintf("%.1f%%", 100*mean(htn_alt)), ")")
message("  Additional people captured                   : ", sum(htn_alt) - sum(htn_now))

rf_now <- pd$ascvd_rf_count
rf_alt <- rf_now - htn_now + htn_alt
message("  has_ascvd_rf now                             : ", sum(rf_now >= 1))
message("  has_ascvd_rf under the wider definition      : ", sum(rf_alt >= 1))

## Does it change who has a goal, or only how strict the goal is?
goal_now <- pd$has_apob_goal
message("  Adults with an apoB goal now                 : ", sum(goal_now))
message("  --> has_ascvd_rf gates the 70-vs-90 value in the diabetes and",
        "\n      hypertriglyceridaemia rows, so expect the goal-assigned COUNT to be",
        "\n      stable and the DISTRIBUTION of goal values to shift stricter.")
message("  apoB goal value distribution now:")
print(table(pd$apob_goal, useNA = "ifany"))

## sensitivity: also check the first BP reading against the 3-reading average
if ("avg_sbp" %in% names(pd)) {
  d_sbp <- pd$BPXSY1 - pd$avg_sbp
  message("  Mean (SBP reading 1 - average of readings)    : ",
          sprintf("%+.2f mmHg", mean(d_sbp, na.rm = TRUE)))
  message("  Hypertensive by reading 1 but not by average  : ",
          sum(pd$BPXSY1 >= 140 & pd$avg_sbp < 140, na.rm = TRUE))
  message("  Hypertensive by average but not by reading 1  : ",
          sum(pd$avg_sbp >= 140 & pd$BPXSY1 < 140, na.rm = TRUE))
  add("F2 hypertension", "mean_reading1_minus_average_mmHg",
      sprintf("%.2f", mean(d_sbp, na.rm = TRUE)))
  add("F2 hypertension", "htn_by_reading1_only",
      sum(pd$BPXSY1 >= 140 & pd$avg_sbp < 140, na.rm = TRUE))
  add("F2 hypertension", "htn_by_average_only",
      sum(pd$avg_sbp >= 140 & pd$BPXSY1 < 140, na.rm = TRUE))
}

add("F2 hypertension", "flagged_untreated_definition", sum(htn_now))
add("F2 hypertension", "flagged_incl_treated", sum(htn_alt))
add("F2 hypertension", "has_ascvd_rf_now", sum(rf_now >= 1))
add("F2 hypertension", "has_ascvd_rf_wider", sum(rf_alt >= 1))

## -----------------------------------------------------------------------------
## FINDING 3. Korn-Graubard vs logit intervals.
##
## The manuscript states Korn-Graubard intervals. svyciprop(method = "beta") is
## the Korn-Graubard-style exact interval; method = "logit" is not. Most of the
## reported estimates use logit. The two diverge most where the proportion is
## small, which is exactly where the headline discordance estimate sits.
## -----------------------------------------------------------------------------
message("\n-- F3. Korn-Graubard (beta) vs logit intervals on the headline numbers --")

apob_d <- subset(primary_design, has_apob_goal == 1)

cmp <- function(label, formula, design) {
  b <- svyciprop(formula, design, method = "beta")
  l <- svyciprop(formula, design, method = "logit")
  data.frame(
    estimate = label,
    pct      = sprintf("%.3f", 100*as.numeric(b)),
    beta_KG  = sprintf("%.3f-%.3f", 100*confint(b)[1], 100*confint(b)[2]),
    logit    = sprintf("%.3f-%.3f", 100*confint(l)[1], 100*confint(l)[2]),
    width_beta  = sprintf("%.3f", 100*(confint(b)[2] - confint(b)[1])),
    width_logit = sprintf("%.3f", 100*(confint(l)[2] - confint(l)[1]))
  )
}

f3 <- bind_rows(
  cmp("goal-assigned share",        ~I(has_apob_goal == 1), primary_design),
  cmp("above apoB goal",            ~I(above_apob_goal == 1), apob_d),
  cmp("discordant (headline)",      ~I(discordant == 1),      apob_d)
)
at_goal_d <- subset(apob_d, at_ldl_goal == 1)
f3 <- bind_rows(f3, cmp("reclassified at LDL-C goal", ~I(discordant == 1), at_goal_d))
print(f3)

for (i in seq_len(nrow(f3))) {
  add("F3 interval method", paste0(f3$estimate[i], " | beta(KG)"), f3$beta_KG[i])
  add("F3 interval method", paste0(f3$estimate[i], " | logit"),    f3$logit[i])
}

message("\n  --> Wherever these two columns differ materially, the manuscript's",
        "\n      'Korn-Graubard' sentence does not describe what produced the number.")

## -----------------------------------------------------------------------------
res <- bind_rows(out)
write.csv(res, "cohort_definition_diagnostics.csv", row.names = FALSE)
message("\n=== COMPLETE === wrote cohort_definition_diagnostics.csv")
message("Nothing in the analysis was changed. Decide on each finding first.")
