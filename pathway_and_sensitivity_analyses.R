## =============================================================================
## pathway_and_sensitivity_analyses.R
## Which Figure 1 row assigns which goal, how the rows overlap, and the
## definitional sensitivities: hypertension threshold, stroke-only ASCVD,
## the PREVENT age restriction, and the ASCVD risk-factor age rule.
##
## RUN AFTER nhanes_apob_analysis.R, in the same session. Reads objects the
## main analysis leaves in the global environment; writes CSVs only, and
## changes nothing about any existing estimate.
##
## WHAT IT PRODUCES
##   r4a_goal_by_pathway.csv        Comment 7  -- which pathway assigns which goal
##   r4a_pathway_overlap.csv        Comment 7  -- how often pathways co-occur
##   r4a_secondary_goal_by_route.csv Comment 7 -- the <55 target and the ASCVD
##                                  route, which primary prevention cannot show
##   r4b_phenotype_at_ldl_goal.csv  Comment 8  -- the 349, against everyone else
##   r4c_htn_threshold_sensitivity.csv  Comment 13a -- SBP >=140 versus >=130
##   r4d_stroke_only_ascvd.csv      Comment 13b -- ASCVD entered via stroke alone
##   r4e_secondary_no_age_limit.csv Comment 3   -- secondary prevention, all ages
##   r4f_ascvd_rf_age_rule.csv      age 55 for all vs men 55 / women 65 as an
##   r4f_ascvd_rf_age_moves.csv     ASCVD risk factor, and who moves between them
## =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(survey); library(purrr); library(tidyr)
})

.need <- c("CONFIG", "primary_df", "secondary_df", "primary_design",
           "secondary_design", "primary_raw", "derive_vars", "assign_goals")
.missing <- .need[!vapply(.need, exists, logical(1), envir = globalenv())]
if (length(.missing))
  stop("\n\n pathway_and_sensitivity_analyses.R cannot run. Missing: ",
       paste(.missing, collapse = ", "),
       "\n Run source(\"nhanes_apob_analysis.R\") first, in this session.\n\n",
       call. = FALSE)

## =============================================================================
## GUARD: is CONFIG still what built the cohort in memory?
## =============================================================================
## R4c below temporarily changes CONFIG$htn_sbp_threshold. If an earlier run of
## THIS script died between changing it and restoring it, the value is left
## dirty, and a later run silently re-derives cohorts at the wrong threshold --
## producing numbers that look plausible and are wrong. That happened once: a
## crashed run left the threshold at 130, and the next run reported the
## secondary-prevention cohort as 856 goal-assigned instead of 850.
##
## Rather than trusting the restore, this recomputes htn_flag for the cohort
## already in memory using the CURRENT CONFIG value and compares it with the
## flag the cohort actually carries. They must match.
local({
  live <- sum(primary_df$htn_flag == 1, na.rm = TRUE)
  recomputed <- if (isTRUE(CONFIG$htn_counts_treatment)) {
    sum(dplyr::coalesce(primary_df$avg_sbp >= CONFIG$htn_sbp_threshold, FALSE) |
        dplyr::coalesce(primary_df$on_antihtn == 1, FALSE), na.rm = TRUE)
  } else {
    sum(dplyr::coalesce(primary_df$BPXSY1 >= CONFIG$htn_sbp_threshold, FALSE),
        na.rm = TRUE)
  }
  if (live != recomputed)
    stop("\n\nCONFIG does not match the cohort in memory.\n",
         "  htn_flag in primary_df counts ", live, " hypertensive adults;\n",
         "  CONFIG$htn_sbp_threshold = ", CONFIG$htn_sbp_threshold,
         " would give ", recomputed, ".\n\n",
         "  A previous run left CONFIG modified. Restart R and re-run the ",
         "pipeline from the top;\n  do not read any output produced in this ",
         "session.\n", call. = FALSE)
  message("  CONFIG check: htn_sbp_threshold = ", CONFIG$htn_sbp_threshold,
          ", consistent with the cohort in memory (", live, " hypertensive).")
})

.written <- character(0)
OUT <- function(x, f) {
  write.csv(x, f, row.names = FALSE); .written <<- c(.written, f)
  message("  wrote ", f)
}

message("\n=== ROUND 4 ANALYSES ===")

## =============================================================================
## R4a. WHICH PATHWAY ASSIGNS WHICH GOAL  (comment 7)
##
## Dani's point: Table 2 says how many adults get a goal but not who they are.
## A reader cannot tell whether the goal-assigned population is mostly people
## with diabetes, mostly people with high triglycerides, or a mixture, nor
## which of them get the stricter <70 target.
##
## The per-row goal columns already exist -- g_dm_apob, g_htg_apob, g_pp_apob,
## g_shc_apob, g_ascvd_apob -- each NA where that row states no apoB goal. So
## "which pathways apply to this person" is readable directly, without
## re-deriving anything, and the union rule means an adult can appear under
## more than one. That is the point: the overlap is the finding.
## =============================================================================

message("\n-- R4a. Goal value by generating pathway --")

pw_design <- subset(primary_design, has_apob_goal == 1)

## Pathway indicators. A pathway "applies" when its own row states a goal.
pw_design <- update(
  pw_design,
  pw_diabetes  = as.integer(!is.na(g_dm_apob)),
  pw_hypertrig = as.integer(!is.na(g_htg_apob)),
  pw_primary   = as.integer(!is.na(g_pp_apob)),
  pw_severe_hc = as.integer(!is.na(g_shc_apob)),
  goal_value   = apob_goal
)

.pw <- c(pw_diabetes  = "Diabetes",
         pw_hypertrig = "Hypertriglyceridaemia (TG 150-499)",
         pw_primary   = "Primary prevention by PREVENT risk")

## For each pathway: how many goal-assigned adults it applies to, and of those
## how many end up at <70 rather than <90. The goal an adult actually carries
## is the most stringent across matched rows, so a pathway that states <90 can
## still sit behind a person whose final goal is <70.
pathway_rows <- bind_rows(lapply(names(.pw), function(v) {
  f  <- as.formula(paste0("~I(", v, " == 1)"))
  pr <- svyciprop(f, pw_design, method = "beta", na.rm = TRUE)
  n  <- sum(pw_design$variables[[v]] == 1, na.rm = TRUE)
  sub <- subset(pw_design, get(v) == 1)
  g70 <- tryCatch(svyciprop(~I(goal_value == 70), sub, method = "beta", na.rm = TRUE),
                  error = function(e) NULL)
  tibble::tibble(
    pathway        = unname(.pw[v]),
    n_unweighted   = n,
    pct_of_goal_assigned = 100 * as.numeric(pr),
    ci_low  = 100 * confint(pr)[1],
    ci_high = 100 * confint(pr)[2],
    pct_with_goal_70 = if (is.null(g70)) NA_real_ else 100 * as.numeric(g70),
    g70_lo  = if (is.null(g70)) NA_real_ else 100 * confint(g70)[1],
    g70_hi  = if (is.null(g70)) NA_real_ else 100 * confint(g70)[2]
  )
}))

## The goal actually carried, which is what a clinician acts on.
goal_value_rows <- bind_rows(lapply(c(70, 90), function(g) {
  pr <- svyciprop(as.formula(paste0("~I(goal_value == ", g, ")")),
                  pw_design, method = "beta", na.rm = TRUE)
  tibble::tibble(
    pathway = paste0("[GOAL CARRIED] apoB <", g, " mg/dL"),
    n_unweighted = sum(pw_design$variables$goal_value == g, na.rm = TRUE),
    pct_of_goal_assigned = 100 * as.numeric(pr),
    ci_low = 100 * confint(pr)[1], ci_high = 100 * confint(pr)[2],
    pct_with_goal_70 = NA_real_, g70_lo = NA_real_, g70_hi = NA_real_)
}))

goal_by_pathway <- bind_rows(goal_value_rows, pathway_rows)
OUT(goal_by_pathway, "r4a_goal_by_pathway.csv")
print(as.data.frame(goal_by_pathway))

## Overlap. Percentages above sum past 100 because the rows are not mutually
## exclusive; this says by how much, and which combinations carry the weight.
ov <- pw_design$variables %>%
  mutate(combo = paste0(
    ifelse(pw_diabetes  == 1, "DM+", ""),
    ifelse(pw_hypertrig == 1, "TG+", ""),
    ifelse(pw_primary   == 1, "PREVENT+", ""))) %>%
  mutate(combo = ifelse(combo == "", "(none -- check)", sub("\\+$", "", combo)),
         n_pathways = pw_diabetes + pw_hypertrig + pw_primary) %>%
  count(combo, n_pathways, name = "n_unweighted") %>%
  arrange(desc(n_unweighted)) %>%
  mutate(pct_unweighted = round(100 * n_unweighted / sum(n_unweighted), 2))
OUT(ov, "r4a_pathway_overlap.csv")
print(as.data.frame(ov))
message("  Goal-assigned adults matching more than one pathway: ",
        sum(ov$n_unweighted[ov$n_pathways > 1]), " of ", sum(ov$n_unweighted))

## ---------------------------------------------------------------------------
## A3 (round 6): the one-row / more-than-one-row split, SURVEY WEIGHTED.
##
## Table 3 printed 47.1% and 52.9% here. Those are unweighted -- 1,567 and
## 1,757 divided by 3,324 -- sitting in a column of otherwise survey-weighted
## percentages, which is exactly the mixing round 6 asked us to stop. Recomputed
## as weighted proportions with Korn-Graubard intervals, matching every other
## percentage in that table.
pw_design <- update(pw_design,
                    multi_row = as.integer((pw_diabetes + pw_hypertrig +
                                            pw_primary) > 1))
.multi <- svyciprop(~I(multi_row == 1), pw_design, method = "beta", na.rm = TRUE)
.one   <- svyciprop(~I(multi_row == 0), pw_design, method = "beta", na.rm = TRUE)
rows_matched <- tibble(
  rows_matched   = c("One row only", "More than one row"),
  n_unweighted   = c(sum(pw_design$variables$multi_row == 0, na.rm = TRUE),
                     sum(pw_design$variables$multi_row == 1, na.rm = TRUE)),
  pct_weighted   = 100 * c(as.numeric(.one),  as.numeric(.multi)),
  ci_low         = 100 * c(attr(.one, "ci")[1], attr(.multi, "ci")[1]),
  ci_high        = 100 * c(attr(.one, "ci")[2], attr(.multi, "ci")[2]),
  ci_method      = "korn-graubard")
OUT(rows_matched, "r4a_rows_matched_weighted.csv")
print(as.data.frame(rows_matched))
message("  Rows matched per adult, survey weighted: one row ",
        sprintf("%.1f%%", rows_matched$pct_weighted[1]),
        ", more than one ", sprintf("%.1f%%", rows_matched$pct_weighted[2]),
        "  (unweighted were 47.1% and 52.9%)")

## -----------------------------------------------------------------------------
## R4a, secondary-prevention arm.
##
## The reviewer asked for a display showing which patients receive each apoB
## target and which pathway generated it, naming three targets -- <90, <70 and
## <55 -- and four pathways including ASCVD. Everything above is primary
## prevention, where <55 cannot occur: it is stated only under clinical ASCVD
## at very high risk. So the <55 column and the ASCVD pathway have no source
## unless secondary prevention is broken down the same way, which is what this
## does. Without it Table 3 answers two thirds of the question asked.
##
## The two routes are mutually exclusive by construction -- very high risk is
## evaluated first and carries <55; everything else with a stated goal carries
## <70 from the ASCVD row with TG >=150 -- so unlike the primary-prevention
## panel above these rows partition rather than overlap.
## -----------------------------------------------------------------------------
message("\n-- R4a (secondary). Goal value by route, established ASCVD --")

sec_ap <- subset(secondary_design, has_apob_goal == 1)
sec_n  <- sum(secondary_df$has_apob_goal == 1, na.rm = TRUE)

sec_goal_rows <- bind_rows(lapply(c(55, 70, 90), function(g) {
  n <- sum(sec_ap$variables$apob_goal == g, na.rm = TRUE)
  if (n == 0) return(NULL)
  pr <- svyciprop(as.formula(paste0("~I(apob_goal == ", g, ")")),
                  sec_ap, method = "beta", na.rm = TRUE)
  tibble::tibble(
    route = paste0("[GOAL CARRIED] apoB <", g, " mg/dL"),
    n_unweighted = n,
    pct_of_goal_assigned = 100 * as.numeric(pr),
    ci_low = 100 * confint(pr)[1], ci_high = 100 * confint(pr)[2])
}))

## Which route put them there. ascvd_very_high already encodes the guideline's
## very-high-risk definition including the CKD branch adjudicated in eMethods 1.
sec_route_rows <- bind_rows(lapply(
  list(c("ascvd_very_high == 1", "Clinical ASCVD at very high risk"),
       c("ascvd_very_high != 1", "Clinical ASCVD, not very high risk")),
  function(x) {
    pr <- svyciprop(as.formula(paste0("~I(", x[1], ")")), sec_ap,
                    method = "beta", na.rm = TRUE)
    n <- sum(with(sec_ap$variables, eval(parse(text = x[1]))), na.rm = TRUE)
    tibble::tibble(route = x[2], n_unweighted = n,
                   pct_of_goal_assigned = 100 * as.numeric(pr),
                   ci_low = 100 * confint(pr)[1],
                   ci_high = 100 * confint(pr)[2])
  }))

sec_goal_by_route <- bind_rows(sec_goal_rows, sec_route_rows)
OUT(sec_goal_by_route, "r4a_secondary_goal_by_route.csv")
print(as.data.frame(sec_goal_by_route))
stopifnot(sum(sec_goal_rows$n_unweighted) == sec_n,
          sum(sec_route_rows$n_unweighted) == sec_n)
message("  Both panels partition the ", sec_n,
        " goal-assigned secondary-prevention adults.")

## =============================================================================
## R4b. WHO THE 349 ARE  (comment 8)
##
## The 22.8% (now recomputed) is a conditional probability within a selected
## subgroup, and Dani is right that a reader can mistake it for a population
## prevalence. The defence is to show who that subgroup is: if the adults at
## their LDL-C goal look systematically different from the rest of the
## goal-assigned cohort, the conditional reading is unavoidable and the table
## makes it visible.
## =============================================================================

message("\n-- R4b. Phenotype of adults at their assigned LDL-C goal --")

ph <- update(subset(primary_design, has_apob_goal == 1),
             at_ldl = factor(ifelse(at_ldl_goal == 1, "At LDL-C goal",
                                    "Above LDL-C goal")))

.cont <- c(RIDAGEYR = "Age, years", BMXBMI = "Body-mass index, kg/m2",
           LBXTR = "Triglycerides, mg/dL", LBXAPB = "Apolipoprotein B, mg/dL",
           ldl_primary = "LDL-C (Sampson), mg/dL", LBDHDD = "HDL-C, mg/dL",
           prevent_ascvd_10yr = "PREVENT 10-year risk, proportion")
.bin  <- c(diabetes = "Diabetes, %", hypertrig = "Triglycerides >=150 mg/dL, %",
           on_statin = "Taking a statin, %", ckd = "Chronic kidney disease, %",
           htn_flag = "Hypertension, %")

pheno <- bind_rows(
  bind_rows(lapply(names(.cont), function(v) {
    m <- svyby(as.formula(paste0("~", v)), ~at_ldl, ph, svymean, na.rm = TRUE)
    tibble::tibble(characteristic = unname(.cont[v]),
                   at_ldl_goal    = sprintf("%.1f (%.2f)", m[[2]][m$at_ldl == "At LDL-C goal"],
                                            m$se[m$at_ldl == "At LDL-C goal"]),
                   above_ldl_goal = sprintf("%.1f (%.2f)", m[[2]][m$at_ldl == "Above LDL-C goal"],
                                            m$se[m$at_ldl == "Above LDL-C goal"]),
                   scale = "mean (SE)")
  })),
  bind_rows(lapply(names(.bin), function(v) {
    m <- svyby(as.formula(paste0("~I(", v, " == 1)")), ~at_ldl, ph,
               svyciprop, vartype = "ci", na.rm = TRUE)
    fmt <- function(i) sprintf("%.1f (%.1f-%.1f)", 100 * m[[2]][i],
                               100 * m[[3]][i], 100 * m[[4]][i])
    ia <- which(m$at_ldl == "At LDL-C goal"); ib <- which(m$at_ldl == "Above LDL-C goal")
    tibble::tibble(characteristic = unname(.bin[v]),
                   at_ldl_goal = fmt(ia), above_ldl_goal = fmt(ib),
                   scale = "% (95% CI)")
  }))
)
n_at  <- sum(primary_df$has_apob_goal == 1 & primary_df$at_ldl_goal == 1, na.rm = TRUE)
n_abv <- sum(primary_df$has_apob_goal == 1 & primary_df$at_ldl_goal == 0, na.rm = TRUE)
pheno <- bind_rows(
  tibble::tibble(characteristic = "Unweighted n",
                 at_ldl_goal = format(n_at, big.mark = ","),
                 above_ldl_goal = format(n_abv, big.mark = ","), scale = "n"),
  pheno)
OUT(pheno, "r4b_phenotype_at_ldl_goal.csv")
print(as.data.frame(pheno))

## =============================================================================
## R4c. HYPERTENSION DEFINITION  (comment 13a)
##
## Re-derives the cohort with the systolic threshold at 130 instead of 140 and
## reports what moves. htn_flag feeds ascvd_rf_count, which sets the 70-vs-90
## goal VALUE in the diabetes and hypertriglyceridaemia rows, so membership of
## the goal-assigned population cannot change and attainment can. The run
## asserts the first half rather than claiming it.
## =============================================================================

message("\n-- R4c. Hypertension threshold: >=140 versus >=130 mmHg --")

## Setting CONFIG has to be done explicitly rather than with `<<-`.
## `CONFIG$x <<- v` works inside a function -- the search starts at the
## enclosing environment and finds CONFIG in the global one -- but at the TOP
## LEVEL of a sourced script the search starts at the global environment's
## PARENT, i.e. the package search path, where CONFIG does not exist. It fails
## with "object 'CONFIG' not found", and it fails at the restore line, after
## the analysis has run and before the result is written. assign() into
## globalenv() behaves the same way in both places.
.set_htn <- function(thresh) {
  cfg <- get("CONFIG", envir = globalenv())
  cfg$htn_sbp_threshold <- thresh
  assign("CONFIG", cfg, envir = globalenv())
}

.orig_thresh <- CONFIG$htn_sbp_threshold
htn_variant <- function(thresh) {
  .set_htn(thresh)
  df <- derive_vars(primary_raw) %>%
    filter(RIDAGEYR >= CONFIG$age_min, RIDAGEYR <= CONFIG$age_max,
           pregnant == 0, ascvd_hx == 0,
           !is.na(LBXAPB), !is.na(ldl_primary), !is.na(LBXTR), WTSAF2YR > 0) %>%
    compute_prevent_risk() %>%
    filter(!is.na(prevent_ascvd_10yr)) %>%
    assign_goals() %>%
    mutate(wt_pooled = WTSAF2YR / length(primary_suffixes_used))
  d  <- svydesign(id = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~wt_pooled,
                  nest = TRUE, data = df)
  ap <- subset(d, has_apob_goal == 1)
  sh <- svyciprop(~I(has_apob_goal == 1), d,  method = "beta", na.rm = TRUE)
  ab <- svyciprop(~I(above_apob_goal == 1), ap, method = "beta", na.rm = TRUE)
  di <- svyciprop(~I(discordant == 1),      ap, method = "beta", na.rm = TRUE)
  tibble::tibble(
    sbp_threshold = thresh,
    n_hypertensive = sum(df$htn_flag == 1, na.rm = TRUE),
    n_goal_assigned = sum(df$has_apob_goal == 1, na.rm = TRUE),
    n_goal_70 = sum(df$apob_goal == 70, na.rm = TRUE),
    n_goal_90 = sum(df$apob_goal == 90, na.rm = TRUE),
    goal_assigned_pct = 100 * as.numeric(sh),
    above_goal_pct = 100 * as.numeric(ab),
    above_lo = 100 * confint(ab)[1], above_hi = 100 * confint(ab)[2],
    discordance_pct = 100 * as.numeric(di),
    disc_lo = 100 * confint(di)[1], disc_hi = 100 * confint(di)[2])
}

## finally = runs whether or not the variants succeed, so a failure inside
## one of them cannot leave CONFIG holding a threshold the rest of the session
## did not ask for.
htn_sens <- tryCatch(
  bind_rows(htn_variant(140), htn_variant(130)),
  finally = .set_htn(.orig_thresh)
)
stopifnot(identical(CONFIG$htn_sbp_threshold, .orig_thresh))

htn_sens$note <- c("As reported", "Contemporary BP definition")
OUT(htn_sens, "r4c_htn_threshold_sensitivity.csv")
print(as.data.frame(htn_sens))

if (htn_sens$n_goal_assigned[1] != htn_sens$n_goal_assigned[2])
  warning("Goal-assigned count changed with the BP threshold. It should not: ",
          "htn_flag sets goal VALUES, not membership. Investigate before ",
          "reporting.", call. = FALSE, immediate. = TRUE)
message("  Adults whose goal VALUE moves from 90 to 70 at the lower threshold: ",
        htn_sens$n_goal_70[2] - htn_sens$n_goal_70[1])

## =============================================================================
## R4d. STROKE-ONLY ASCVD  (comment 13b)
##
## NHANES asks whether a stroke was ever diagnosed (MCQ160F) and does not
## distinguish ischaemic from haemorrhagic. Haemorrhagic stroke is not
## atherosclerotic, so anyone entering the secondary-prevention cohort on a
## stroke answer ALONE is potentially misclassified. This counts them, which
## turns a limitation sentence into a bounded one.
## =============================================================================

message("\n-- R4d. ASCVD entered through stroke alone --")

sd <- secondary_df
stroke_only <- sd %>%
  mutate(only_stroke = as.integer(
    coalesce(MCQ160F == 1, FALSE) &
    !coalesce(MCQ160C == 1, FALSE) &
    !coalesce(MCQ160D == 1, FALSE) &
    !coalesce(MCQ160E == 1, FALSE)))
n_only <- sum(stroke_only$only_stroke == 1, na.rm = TRUE)

sdd <- update(secondary_design, only_stroke = stroke_only$only_stroke)
pr  <- svyciprop(~I(only_stroke == 1), sdd, method = "beta", na.rm = TRUE)

r4d <- tibble::tibble(
  metric = "Secondary-prevention adults whose only qualifying event is stroke",
  n_unweighted = n_only,
  n_cohort = nrow(sd),
  pct_weighted = 100 * as.numeric(pr),
  ci_low = 100 * confint(pr)[1], ci_high = 100 * confint(pr)[2],
  note = paste("NHANES MCQ160F does not distinguish ischaemic from",
               "haemorrhagic stroke; this is the maximum possible",
               "misclassification of the secondary-prevention cohort."))
OUT(r4d, "r4d_stroke_only_ascvd.csv")
print(as.data.frame(r4d))

## =============================================================================
## R4e. SECONDARY PREVENTION WITHOUT THE PREVENT AGE RESTRICTION  (comment 3)
##
## 30-79 is the validated range of the PREVENT equations. PREVENT is not used
## in established ASCVD -- goals there come from the clinical-ASCVD row -- so
## the age restriction is inherited rather than required. This rebuilds the
## secondary cohort across all adult ages and reports the difference.
## =============================================================================

message("\n-- R4e. Secondary prevention, PREVENT age restriction removed --")

sec_variant <- function(lo, hi, label) {
  df <- derive_vars(primary_raw) %>%
    filter(RIDAGEYR >= lo, RIDAGEYR <= hi, pregnant == 0, ascvd_hx == 1,
           !is.na(LBXAPB), !is.na(ldl_primary), !is.na(LBXTR), WTSAF2YR > 0) %>%
    mutate(prevent_ascvd_10yr = NA_real_) %>%
    assign_goals() %>%
    mutate(wt_pooled = WTSAF2YR / length(primary_suffixes_used))
  d  <- svydesign(id = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~wt_pooled,
                  nest = TRUE, data = df)
  ap <- subset(d, has_apob_goal == 1)
  sh <- svyciprop(~I(has_apob_goal == 1),  d,  method = "beta", na.rm = TRUE)
  ab <- svyciprop(~I(above_apob_goal == 1), ap, method = "beta", na.rm = TRUE)
  tibble::tibble(
    cohort = label, age_range = paste0(lo, "-", hi),
    n = nrow(df), n_goal_assigned = sum(df$has_apob_goal == 1, na.rm = TRUE),
    goal_assigned_pct = 100 * as.numeric(sh),
    ga_lo = 100 * confint(sh)[1], ga_hi = 100 * confint(sh)[2],
    above_goal_pct = 100 * as.numeric(ab),
    ag_lo = 100 * confint(ab)[1], ag_hi = 100 * confint(ab)[2])
}

sec_age <- bind_rows(
  sec_variant(CONFIG$age_min, CONFIG$age_max, "As reported (PREVENT range)"),
  sec_variant(20, 150, "All adults 20+"),
  sec_variant(CONFIG$age_min, 150, "30 and older, no upper limit"))
OUT(sec_age, "r4e_secondary_no_age_limit.csv")
print(as.data.frame(sec_age))
message("  PREVENT plays no part in secondary-prevention goal assignment, so ",
        "any difference here is population composition, not method.")


## =============================================================================
## R4f. AGE AS AN ASCVD RISK FACTOR: ONE CUT AT 55, OR THE SEX-SPECIFIC PAIR
##
## Figure 1 gives an apoB goal to several rows only when the patient also has
## "at least 1 ASCVD risk factor", and never defines the term. The document's
## only enumeration is a footnote to the icosapent-ethyl recommendation in
## diabetes, taken from REDUCE-IT, which reads "men >=55 years or women >=65
## years". The primary analysis uses a single cut at 55 for both sexes.
##
## The difference is not cosmetic. A woman of 58 with triglycerides of 200 and
## nothing else has a risk factor under the single cut and none under the
## footnote, and the hypertriglyceridaemia row states an apoB goal only when a
## risk factor is present -- so she moves in or out of the goal-assigned
## population, which is the denominator for every headline estimate.
##
## Reported, not resolved: whether a footnote attached to one recommendation
## governs a figure elsewhere in the document is exactly the kind of question
## this paper says the guideline leaves open.
## =============================================================================

message("\n-- R4f. Age as an ASCVD risk factor: 55 for all vs men 55 / women 65 --")

.set_rf <- function(v) {
  cfg <- get("CONFIG", envir = globalenv())
  cfg$ascvd_rf_age_sex_specific <- v
  assign("CONFIG", cfg, envir = globalenv())
}
.orig_rf <- CONFIG$ascvd_rf_age_sex_specific

rf_variant <- function(sex_specific, label) {
  .set_rf(sex_specific)
  df <- derive_vars(primary_raw) %>%
    filter(RIDAGEYR >= CONFIG$age_min, RIDAGEYR <= CONFIG$age_max,
           pregnant == 0, ascvd_hx == 0,
           !is.na(LBXAPB), !is.na(ldl_primary), !is.na(LBXTR), WTSAF2YR > 0) %>%
    compute_prevent_risk() %>%
    filter(!is.na(prevent_ascvd_10yr)) %>%
    assign_goals() %>%
    mutate(wt_pooled = WTSAF2YR / length(primary_suffixes_used))
  d  <- svydesign(id = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~wt_pooled,
                  nest = TRUE, data = df)
  ap <- subset(d, has_apob_goal == 1)
  sh <- svyciprop(~I(has_apob_goal == 1),   d,  method = "beta", na.rm = TRUE)
  ab <- svyciprop(~I(above_apob_goal == 1), ap, method = "beta", na.rm = TRUE)
  di <- svyciprop(~I(discordant == 1),      ap, method = "beta", na.rm = TRUE)
  list(
    row = tibble::tibble(
      rule = label,
      n_cohort = nrow(df),
      n_with_any_rf   = sum(df$has_ascvd_rf == 1, na.rm = TRUE),
      n_goal_assigned = sum(df$has_apob_goal == 1, na.rm = TRUE),
      n_goal_70 = sum(df$apob_goal == 70, na.rm = TRUE),
      n_goal_90 = sum(df$apob_goal == 90, na.rm = TRUE),
      goal_assigned_pct = 100 * as.numeric(sh),
      ga_lo = 100 * confint(sh)[1], ga_hi = 100 * confint(sh)[2],
      above_goal_pct = 100 * as.numeric(ab),
      ag_lo = 100 * confint(ab)[1], ag_hi = 100 * confint(ab)[2],
      discordance_pct = 100 * as.numeric(di),
      disc_lo = 100 * confint(di)[1], disc_hi = 100 * confint(di)[2]),
    ## keyed by SEQN so the two arms can be compared person by person
    key = df %>% select(SEQN, RIAGENDR, RIDAGEYR, has_ascvd_rf, has_apob_goal,
                        apob_goal))
}

rf_out <- tryCatch({
  a <- rf_variant(FALSE, "Age >=55, both sexes (as reported)")
  b <- rf_variant(TRUE,  "Men >=55 / women >=65 (guideline footnote)")
  list(a = a, b = b)
}, finally = .set_rf(.orig_rf))
stopifnot(identical(CONFIG$ascvd_rf_age_sex_specific, .orig_rf))

rf_sens <- bind_rows(rf_out$a$row, rf_out$b$row)
OUT(rf_sens, "r4f_ascvd_rf_age_rule.csv")
print(as.data.frame(rf_sens))

## Who actually moves, and whether they are the women the footnote predicts.
mv <- rf_out$a$key %>%
  rename(rf_a = has_ascvd_rf, goal_a = has_apob_goal, val_a = apob_goal) %>%
  inner_join(rf_out$b$key %>%
               select(SEQN, rf_b = has_ascvd_rf, goal_b = has_apob_goal,
                      val_b = apob_goal),
             by = "SEQN")
stopifnot(nrow(mv) == rf_out$a$row$n_cohort)

rf_moves <- tibble::tibble(
  quantity = c("Lost a risk factor under the footnote",
               "  of whom women aged 55-64",
               "Lost their apoB goal entirely",
               "Goal value moved from <70 to <90",
               "Any change to the assigned apoB goal"),
  n = c(sum(mv$rf_a == 1 & mv$rf_b == 0),
        sum(mv$rf_a == 1 & mv$rf_b == 0 &
            mv$RIAGENDR == 2 & mv$RIDAGEYR >= 55 & mv$RIDAGEYR < 65),
        sum(mv$goal_a == 1 & mv$goal_b == 0),
        sum(coalesce(mv$val_a == 70 & mv$val_b == 90, FALSE)),
        ## val_a != val_b is NA whenever either arm assigns no goal, so the
        ## coalesce falls back to the has-a-goal comparison: both NA means no
        ## goal either way and no change; one NA means the goal appeared or
        ## disappeared, which is a change.
        sum(coalesce(mv$val_a != mv$val_b, mv$goal_a != mv$goal_b))))
OUT(rf_moves, "r4f_ascvd_rf_age_moves.csv")
print(as.data.frame(rf_moves))

## Nobody should GAIN a risk factor by narrowing the definition.
stopifnot(sum(mv$rf_a == 0 & mv$rf_b == 1) == 0)

message("  Goal-assigned share: ",
        sprintf("%.2f%% -> %.2f%%", rf_sens$goal_assigned_pct[1],
                rf_sens$goal_assigned_pct[2]),
        "   (", rf_sens$n_goal_assigned[1], " -> ",
        rf_sens$n_goal_assigned[2], " adults)")
message("  Report this as it comes out. A shift here belongs in the ",
        "interpretive range\n  alongside the 15.7%-33.4% already reported.")

message("\n=== ROUND 4 ANALYSES COMPLETE ===")
message(length(.written), " file(s) written:")
for (f in .written) message("  - ", f)

## ---------------------------------------------------------------------------
## B2 (round 6): SEVERE HYPERTRIGLYCERIDAEMIA AND THE DISCORDANT GROUP
##
## Eight of the 77 discordant adults have triglycerides above 800 mg/dL. At
## those levels the lipid panel behaves differently -- apoB rises with particle
## number in a way that is not the atherogenic picture the goal framework is
## aimed at, and calculated LDL-C is least reliable. If a tenth of the
## discordant group sits there, a reviewer will want to know the finding is not
## carried by them.
##
## Reported on the SAME three quantities as the main analysis, so the
## comparison is like for like: apoB non-attainment among goal-assigned adults,
## overall discordance, and discordance among those at their LDL-C goal.
local({
  TG_CUT <- 800
  base <- apob_design
  keep <- subset(base, coalesce(LBXTR <= TG_CUT, TRUE))

  n_excl <- sum(coalesce(base$variables$LBXTR > TG_CUT, FALSE), na.rm = TRUE)
  n_excl_disc <- sum(coalesce(base$variables$LBXTR > TG_CUT, FALSE) &
                     coalesce(base$variables$discordant == 1, FALSE),
                     na.rm = TRUE)

  q <- function(d, f) {
    ci <- svyciprop(f, d, method = "beta", na.rm = TRUE)
    c(100 * as.numeric(ci), 100 * attr(ci, "ci")[1], 100 * attr(ci, "ci")[2])
  }
  atldl      <- subset(base, at_ldl_goal == 1)
  atldl_keep <- subset(keep, at_ldl_goal == 1)

  rows <- list(
    c("ApoB non-attainment among goal-assigned",
      q(base, ~I(above_apob_goal == 1)), q(keep, ~I(above_apob_goal == 1))),
    c("Discordance among goal-assigned",
      q(base, ~I(discordant == 1)),      q(keep, ~I(discordant == 1))),
    c("Discordance among those at LDL-C goal",
      q(atldl, ~I(above_apob_goal == 1)),
      q(atldl_keep, ~I(above_apob_goal == 1))))

  out <- tibble(
    quantity     = vapply(rows, function(r) r[1], character(1)),
    all_pct      = as.numeric(vapply(rows, function(r) r[2], character(1))),
    all_lo       = as.numeric(vapply(rows, function(r) r[3], character(1))),
    all_hi       = as.numeric(vapply(rows, function(r) r[4], character(1))),
    excl_pct     = as.numeric(vapply(rows, function(r) r[5], character(1))),
    excl_lo      = as.numeric(vapply(rows, function(r) r[6], character(1))),
    excl_hi      = as.numeric(vapply(rows, function(r) r[7], character(1))),
    tg_cut       = TG_CUT,
    n_excluded   = n_excl,
    n_excluded_discordant = n_excl_disc,
    ci_method    = "korn-graubard")
  OUT(out, "r6_tg800_sensitivity.csv")
  print(as.data.frame(out))
  message("\nB2 excluding triglycerides > ", TG_CUT, " mg/dL: ", n_excl,
          " goal-assigned adults removed, of whom ", n_excl_disc,
          " were discordant.")
  for (i in seq_len(nrow(out)))
    message("  ", format(out$quantity[i], width = 40),
            sprintf("%6.2f%% -> %6.2f%%", out$all_pct[i], out$excl_pct[i]))
})
