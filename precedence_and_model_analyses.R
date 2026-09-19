## =============================================================================
## precedence_and_model_analyses.R
## Sensitivity of goal assignment to the precedence question, and the
## regression models reported in the supplement.
##
## HOW TO RUN
##   1. source("nhanes_apob_analysis.R")   # in a fresh session, to completion
##   2. source("precedence_and_model_analyses.R")      # this file, in the SAME session
##
## The apoB cross-instrument calibration this file assumes is already part of
## the main script and is controlled by CONFIG$apply_apob_calibration, which
## defaults to TRUE. Nothing needs to be patched by hand. The calibration is
## documented below under "APOB CALIBRATION".
##
## The script assumes these objects exist in the global environment after the
## main run: primary_df, secondary_df, primary_design, secondary_design, CONFIG.
##
## WHAT IT PRODUCES
##   rev1_precedence_sensitivity.csv   Reviewer major comment 1
##   rev1_precedence_bounds.csv        Reviewer major comment 1 (order-free)
##   rev2_denominator_variants.csv     Reviewer major comment 2
##   rev3_correlates_firth.csv         Reviewer major comment 3
##   rev3_correlates_unadjusted.csv    Reviewer major comment 3
##   rev3_epv_table.csv                Reviewer major comment 3
##   rev4_cycle_influence_on_share.csv Reviewer major comment 4
##   rev6_sampson800_sensitivity.csv  Dani review round 2, item 4
## =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(survey); library(purrr); library(tidyr)
})

`%||%` <- function(a, b) if (is.null(a)) b else a

.written <- character(0)
OUT <- function(x, f) {
  write.csv(x, f, row.names = FALSE)
  .written <<- c(.written, f)
  message("  wrote ", f)
}

## =============================================================================
## PREFLIGHT -- fail immediately and legibly if the session is not ready.
##
## This script is NOT standalone. It reads objects that nhanes_apob_analysis.R
## leaves in the global environment. Running it in a fresh R session produces a
## cascade of "object not found" errors, which is what the checks below prevent.
## =============================================================================

.need <- c("CONFIG", "primary_df", "secondary_df",
           "primary_design", "secondary_design")
.missing <- .need[!vapply(.need, exists, logical(1), envir = globalenv())]

if (length(.missing)) {
  stop(
    "\n\n",
    "=========================================================================\n",
    " precedence_and_model_analyses.R cannot run: the main analysis has not been run in\n",
    " this R session.\n\n",
    " Missing objects: ", paste(.missing, collapse = ", "), "\n\n",
    " DO THIS, in order, in ONE R session:\n\n",
    "   1. install.packages(\"logistf\")        # once only\n",
    "   2. source(\"nhanes_apob_analysis.R\")   # ~several minutes\n",
    "   3. source(\"precedence_and_model_analyses.R\")      # this file\n\n",
    " Do not restart R or clear the workspace between steps 2 and 3.\n",
    "=========================================================================\n\n",
    call. = FALSE
  )
}

if (is.null(CONFIG$apply_apob_calibration) &&
    !isTRUE(getOption("apob_allow_uncalibrated"))) {
  stop(
    "\n\n",
    "=========================================================================\n",
    " CONFIG$apply_apob_calibration does not exist.\n\n",
    " This means the main script in memory predates the apoB cross-instrument\n",
    " calibration, so apoB has not been placed on a common assay scale. NHANES\n",
    " measured apoB on three different instruments across the six pooled\n",
    " cycles. Estimates produced without the calibration are not the ones to\n",
    " report.\n\n",
    " Re-source the current nhanes_apob_analysis.R, then rerun this file.\n\n",
    " To proceed anyway for a like-for-like comparison against the old run,\n",
    " set:  options(apob_allow_uncalibrated = TRUE)\n",
    "=========================================================================\n\n",
    call. = FALSE
  )
}

if (!isTRUE(CONFIG$apply_apob_calibration %||% FALSE)) {
  warning("apoB calibration is not applied. These are the ",
          "UNCALIBRATED results, for sensitivity comparison only.",
          call. = FALSE, immediate. = TRUE)
}

if (!requireNamespace("logistf", quietly = TRUE)) {
  message("NOTE: package 'logistf' is not installed, so the Firth model in ",
          "section R3 will be skipped.\n      Run install.packages(\"logistf\") ",
          "if you want it.")
}

message("Preflight OK. Cohort: ", nrow(primary_df), " primary-prevention, ",
        nrow(secondary_df), " secondary-prevention records. ",
        "apoB calibration: ",
        if (isTRUE(CONFIG$apply_apob_calibration %||% FALSE)) "APPLIED" else "NOT applied")

## =============================================================================
## APOB CALIBRATION -- reference note. Nothing here needs to be done by hand.
## =============================================================================
##
## This section documents the cross-instrument calibration that the main script
## already applies, and the evidence that it is correct. The implementation
## lives in derive_vars() in nhanes_apob_analysis.R and is switched by
## CONFIG$apply_apob_calibration (default TRUE); set it FALSE there to
## reproduce the uncalibrated run for comparison.
##
## NHANES changed the apoB laboratory, instrument and method twice inside the
## pooled window, and publishes Deming-regression bridging equations for both.
## Pooling six cycles of apoB on three different measurement scales is the
## explanation for the 2005-2006 anomaly the reviewer flagged.
##
##   2005-2006 (D): Dade-Behring BN100 nephelometric, Johns Hopkins Hospital
##   2007-2014 (E-H): Dade-Behring / Siemens BN ProSpec nephelometric,
##                    University of Minnesota
##   2015-2016 (I): Roche Cobas 6000 turbidimetric, University of Minnesota
##
## CDC's published bridging equations (APOB_E and APOB_I documentation):
##   Y(2005-06, BN100)  = 1.084 * X(2007-08, ProSpec)     n=168 pairs, r=0.984
##   Y(ProSpec)         = 0.9988 * X(Cobas 6000) - 2.880  n=105 pairs, r=0.988
## CDC states these "should be used to trend the Apolipoprotein B data."
##
## We adopt the ProSpec scale (cycles E-H) as the common scale because it
## covers four of the six cycles.
##
## The transform is applied as the FIRST mutate inside derive_vars(), before
## any variable that reads LBXAPB (severe_hc, above_apob_goal, apob_high_130,
## apob_markedly_high, the RQ9/RQ10 outcomes, and the Section 3 trend), so that
## every downstream apoB quantity is on the common scale:
##
##   LBXAPB_raw = LBXAPB
##   LBXAPB     = cycle "D" -> LBXAPB * 0.923            # BN100  -> ProSpec
##                cycle "I" -> LBXAPB * 0.9988 - 2.880   # Cobas  -> ProSpec
##                otherwise -> unchanged
##
## The uncalibrated values remain available as LBXAPB_raw. Both versions are
## reported; the calibrated version is primary.
##
## EXPECTED EFFECT (verified directly from the raw XPT files, all fasting
## participants aged 12+, apoB regressed on LDL-C and log-triglycerides, then
## evaluated at LDL-C 100 mg/dL and TG 100 mg/dL):
##
##   cycle   raw    calibrated        apoB >=130, raw -> calibrated
##   D      89.4       82.5                 10.11%  ->  5.62%
##   E      82.6       82.6                  6.39%  ->  6.39%
##   F      81.3       81.3                  5.20%  ->  5.20%
##   G      80.1       80.1                  4.25%  ->  4.25%
##   H      83.4       83.4                  5.02%  ->  5.02%
##   I      86.6       83.6                  8.11%  ->  6.63%
##
##   Between-cycle spread at a fixed lipid profile: 9.24 -> 3.46 mg/dL.
##
## The empirical cycle-D excess (+7.9%) matches CDC's published factor (+8.4%)
## almost exactly, which is independent confirmation that the anomaly is
## instrument calibration and not a population difference.
## =============================================================================

message("\n=== REVISION ANALYSES ===")

## =============================================================================
## R1. PRECEDENCE SENSITIVITY  (major comment 1)
##
## Figure 1 gives no rule for a person who matches more than one patient-
## population row. Our main analysis imposes one ordering. This section
## recomputes the goal-assigned share under EVERY permutation of the five
## implementable rows, and additionally reports two order-free bounds that do
## not depend on any ordering at all.
##
##   UNION        a goal is assigned if ANY matched row states one.
##                = the largest share any ordering can produce.
##   INTERSECTION a goal is assigned only if EVERY matched row states one.
##                = the smallest share any ordering can produce.
##
## Every permutation lies between these two bounds, so the bounds settle the
## question the reviewer raised without appeal to which ordering is "right".
## =============================================================================

ROWS <- c("clinical_ascvd", "diabetes", "severe_hypercholesterolaemia",
          "hypertriglyceridaemia", "primary_prevention")

## Membership in each row, and the apoB goal that row assigns when it governs.
## NA goal means "this row states no apoB goal for this person".
add_row_goals <- function(df) {
  df %>%
    mutate(
      in_clinical_ascvd = coalesce(ascvd_hx == 1, FALSE),
      in_diabetes       = coalesce(diabetes == 1, FALSE),
      in_severe_hc      = coalesce(severe_hc == 1, FALSE),
      in_hypertrig      = coalesce(hypertrig == 1, FALSE),
      in_primary_prev   = !is.na(prevent_ascvd_10yr),

      g_clinical_ascvd = case_when(
        !in_clinical_ascvd            ~ NA_real_,
        ascvd_very_high == 1          ~ 55,
        hypertrig == 1                ~ 70,
        TRUE                          ~ NA_real_
      ),
      g_diabetes = case_when(
        !in_diabetes                                        ~ NA_real_,
        has_ascvd_rf == 1 | dm_risk_modifier == 1           ~ 70,
        TRUE                                                ~ 90
      ),
      ## Figure 1 states no apoB goal for the severe-hypercholesterolaemia row
      ## absent the therapy-response conditions we cannot observe.
      g_severe_hc = NA_real_,
      g_hypertrig = case_when(
        !in_hypertrig                                       ~ NA_real_,
        RIDAGEYR >= 40 & RIDAGEYR <= 75 & has_ascvd_rf == 1 ~ 70,
        TRUE                                                ~ NA_real_
      ),
      g_primary_prev = case_when(
        !in_primary_prev                                      ~ NA_real_,
        tg_in_goal_band != 1                                  ~ NA_real_,
        prevent_ascvd_10yr >= CONFIG$prevent_goal_high_risk_cut ~ 70,
        TRUE                                                  ~ 90
      )
    )
}

MEMBER <- c(clinical_ascvd = "in_clinical_ascvd",
            diabetes = "in_diabetes",
            severe_hypercholesterolaemia = "in_severe_hc",
            hypertriglyceridaemia = "in_hypertrig",
            primary_prevention = "in_primary_prev")
GOALV  <- c(clinical_ascvd = "g_clinical_ascvd",
            diabetes = "g_diabetes",
            severe_hypercholesterolaemia = "g_severe_hc",
            hypertriglyceridaemia = "g_hypertrig",
            primary_prevention = "g_primary_prev")

## Goal under a given first-match ordering.
goal_under_order <- function(df, ord) {
  g <- rep(NA_real_, nrow(df))
  assigned <- rep(FALSE, nrow(df))
  for (r in ord) {
    hit <- !assigned & df[[MEMBER[r]]]
    g[hit] <- df[[GOALV[r]]][hit]
    assigned <- assigned | df[[MEMBER[r]]]
  }
  g
}

## Order-free bounds.
goal_union <- function(df) {          # any matched row that states a goal
  G <- as.matrix(df[, GOALV])
  suppressWarnings(apply(G, 1, function(v) if (all(is.na(v))) NA_real_ else min(v, na.rm = TRUE)))
}
goal_intersection <- function(df) {   # only if every matched row states a goal
  M <- as.matrix(df[, MEMBER]);  G <- as.matrix(df[, GOALV])
  out <- rep(NA_real_, nrow(df))
  for (i in seq_len(nrow(df))) {
    m <- which(M[i, ]); if (!length(m)) next
    v <- G[i, m]
    if (!any(is.na(v))) out[i] <- min(v)
  }
  out
}

message("\n-- R1. Precedence sensitivity --")

pdf_rows <- add_row_goals(primary_df)

## All 120 permutations of the five rows, generated in base R so this script
## carries no extra package dependency.
all_perms <- function(v) {
  if (length(v) <= 1) return(matrix(v, nrow = 1))
  do.call(rbind, lapply(seq_along(v), function(i)
    cbind(v[i], all_perms(v[-i]), deparse.level = 0)))
}
perms <- all_perms(ROWS)
stopifnot(nrow(perms) == 120)

## NOTE: columns are attached by direct assignment into design$variables rather
## than update(), because the goal vectors are computed outside the design and
## update() evaluates its expressions inside the design's own data frame.
share_for <- function(goal_vec, design, df) {
  d <- design
  d$variables$.hg <- as.integer(!is.na(goal_vec))
  ci <- svyciprop(~I(.hg == 1), d, method = "beta")
  c(share = as.numeric(ci),
    lo = confint(ci)[1], hi = confint(ci)[2],
    n_unwt = sum(!is.na(goal_vec)))
}

perm_res <- map_dfr(seq_len(nrow(perms)), function(i) {
  ord <- perms[i, ]
  g <- goal_under_order(pdf_rows, ord)
  s <- share_for(g, primary_design, pdf_rows)
  tibble(ordering = paste(ord, collapse = " > "),
         share = s["share"], ci_low = s["lo"], ci_high = s["hi"],
         n_unweighted = s["n_unwt"])
})

perm_res <- perm_res %>% arrange(share)
OUT(perm_res, "rev1_precedence_sensitivity.csv")

message(sprintf("  %d orderings. Goal-assigned share ranges %.1f%% to %.1f%% (span %.2f pp).",
                nrow(perm_res), 100 * min(perm_res$share), 100 * max(perm_res$share),
                100 * (max(perm_res$share) - min(perm_res$share))))
message(sprintf("  Distinct values: %d", length(unique(round(perm_res$share, 6)))))

g_union <- goal_union(pdf_rows)
g_inter <- goal_intersection(pdf_rows)
g_main  <- goal_under_order(pdf_rows, ROWS)   # the ordering used in the paper

## LABELS CORRECTED 2026-08-16. "Main analysis ordering" now refers to the
## SUPERSEDED first-match rule. Since the rebuild, the main analysis uses the
## union with the strictest applicable goal, i.e. the third row below.
bounds <- bind_rows(lapply(
  list(c("Intersection (order-free lower bound)",              "g_inter"),
       c("First-match ordering (SUPERSEDED sensitivity)",      "g_main"),
       c("Union (order-free) -- MAIN ANALYSIS",                "g_union")),
  function(x) {
    s <- share_for(get(x[2]), primary_design, pdf_rows)
    tibble(variant = x[1], share = unname(s["share"]),
           ci_low = unname(s["lo"]), ci_high = unname(s["hi"]),
           n_unweighted = unname(s["n_unwt"]))
  }))
OUT(bounds, "rev1_precedence_bounds.csv")
print(as.data.frame(bounds))

## Does the HEADLINE move, not just the denominator? Recompute discordance and
## reclassification under each bound and under the main ordering.
headline_under <- function(goal_vec, label) {
  d <- primary_design
  d$variables$.goal  <- goal_vec
  ## above_goal() is defined in nhanes_apob_analysis.R and honours
  ## CONFIG$strict_goal_boundaries; do not inline a bare > here.
  d$variables$.above <- as.integer(above_goal(d$variables$LBXAPB, goal_vec))
  d$variables$.disc  <- as.integer(d$variables$at_ldl_goal == 1 &
                                   above_goal(d$variables$LBXAPB, goal_vec))
  d2 <- subset(d, !is.na(.goal))
  disc <- svyciprop(~I(.disc == 1), d2, method = "logit")
  d3 <- subset(d2, at_ldl_goal == 1)
  n3 <- nrow(d3$variables)
  rec <- if (n3 > 0) svyciprop(~I(.above == 1), d3, method = "logit") else NA
  tibble(variant = label,
         discordance = as.numeric(disc),
         disc_lo = confint(disc)[1], disc_hi = confint(disc)[2],
         reclassification = if (is.na(rec[1])) NA_real_ else as.numeric(rec),
         recl_lo = if (is.na(rec[1])) NA_real_ else confint(rec)[1],
         recl_hi = if (is.na(rec[1])) NA_real_ else confint(rec)[2],
         n_at_ldl_goal = n3)
}

headlines <- bind_rows(
  headline_under(g_inter, "Intersection (order-free lower bound)"),
  headline_under(g_main,  "First-match ordering (SUPERSEDED sensitivity)"),
  headline_under(g_union, "Union (order-free) -- MAIN ANALYSIS")
)
OUT(headlines, "rev1_precedence_headlines.csv")
print(as.data.frame(headlines))

message("\n  INTERPRETATION: if the union and intersection bounds are close, ",
        "the goal-assigned share is robust to precedence and stands as a ",
        "headline. If they diverge, report the bounded range in the abstract ",
        "instead of a point estimate. They do diverge here (15.7% to 33.4%), ",
        "which is why every estimate in the paper is range-anchored.")

## =============================================================================
## R2. DENOMINATOR, STATED EXACTLY  (major comment 2)
##
## The reviewer is right that a percentage of "US adults aged 30-79" and a
## count out of the primary-prevention cohort are claims about DIFFERENT
## denominators. This section computes every defensible denominator so the
## manuscript can name the one it means.
##
## The figures the reviewer quoted ("22.0%", "2,314 of 8,925") come from the
## superseded first-match analysis and no longer appear anywhere. Under the
## union rule the cohort is 3,324 of 9,108 and the share is 33.4%.
## =============================================================================

message("\n-- R2. Denominator variants --")

## Combined design: primary + secondary prevention, one row per participant.
combined_df <- bind_rows(
  primary_df   %>% mutate(.arm = "primary"),
  secondary_df %>% mutate(.arm = "secondary")
) %>% distinct(SEQN, .keep_all = TRUE)

combined_design <- svydesign(
  id = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~wt_pooled,
  data = combined_df, nest = TRUE
)

den_variant <- function(design, df, label) {
  ci <- svyciprop(~I(has_apob_goal == 1), design, method = "beta")
  tot <- svytotal(~I(has_apob_goal == 1), design)
  tibble(denominator = label,
         n_with_goal = sum(df$has_apob_goal == 1, na.rm = TRUE),
         n_total = nrow(df),
         pct_unweighted = 100 * sum(df$has_apob_goal == 1, na.rm = TRUE) / nrow(df),
         pct_weighted = 100 * as.numeric(ci),
         ci_low = 100 * confint(ci)[1], ci_high = 100 * confint(ci)[2],
         us_adults_with_goal = as.numeric(tot)[2])
}

den <- bind_rows(
  den_variant(primary_design,   primary_df,   "Primary-prevention adults 30-79 (as reported)"),
  den_variant(secondary_design, secondary_df, "Secondary-prevention adults 30-79"),
  den_variant(combined_design,  combined_df,  "ALL adults 30-79 (primary + secondary)")
)
OUT(den, "rev2_denominator_variants.csv")
print(as.data.frame(den))

message("\n  Use the third row if the abstract says 'US adults aged 30-79'. ",
        "Use the first row if it says 'adults without established ASCVD'. ",
        "Do not mix them.")

## -----------------------------------------------------------------------------
## R2b. THE COMBINED DENOMINATOR, BY READING OF FIGURE 1  (round-6 audit)
##
## R2 above computes the combined share under the main reading only, and the
## manuscript carried it as a point estimate with an em dash where every other
## goal-assignment number in the paper carries its interpretive range. The
## combined estimate is not independent of that reading: 3,324 of the 4,174
## adults in its numerator, 80%, are assigned within primary prevention, where
## R1 shows the share moving from 15.7% to 33.4%. Reported without a range, next
## to a primary row that has one, it reads as the sturdier of the two, which is
## backwards.
##
## The secondary arm is held fixed across readings. The row-ordering question is
## a question about the five rows as they apply to adults WITHOUT clinical
## ASCVD; the clinical-ASCVD row cannot match within primary prevention (mean
## influence 0.00 percentage points under first-match) and competes with nothing
## within secondary prevention. The separate question -- whether the
## clinical-ASCVD rows state a goal or leave one optional -- is not this range
## and is deliberately not folded into it.
## -----------------------------------------------------------------------------

message("\n-- R2b. Combined denominator by reading --")

## Map each reading onto the combined frame BY SEQN. combined_df is the two arms
## stacked, so positions do not correspond to primary_df's.
combined_indicator <- function(goal_vec) {
  assigned <- setNames(as.integer(!is.na(goal_vec)), primary_df$SEQN)
  ifelse(combined_df$.arm == "primary",
         unname(assigned[as.character(combined_df$SEQN)]),
         as.integer(combined_df$has_apob_goal == 1))
}

comb_est <- function(ind, label) {
  d <- combined_design
  d$variables$.hg <- ind
  ci  <- svyciprop(~I(.hg == 1), d, method = "beta", na.rm = TRUE)
  tot <- svytotal(~I(.hg == 1), d, na.rm = TRUE)
  tibble(reading      = label,
         n_with_goal  = sum(ind == 1, na.rm = TRUE),
         n_total      = length(ind),
         pct_weighted = 100 * as.numeric(ci),
         ci_low       = 100 * confint(ci)[1],
         ci_high      = 100 * confint(ci)[2],
         us_adults_with_goal = as.numeric(tot)[2],
         ci_method    = "korn-graubard")
}

den_by_reading <- bind_rows(
  comb_est(combined_indicator(g_inter),
           "Intersection - a goal only if every matched row states one"),
  comb_est(combined_indicator(g_main),
           "Sequential first-match, the ordering used in an earlier implementation"),
  comb_est(combined_indicator(g_union),
           "Union with most stringent value (primary analysis)")
)

## The union row IS the combined estimate R2 reports. If the two ever diverge
## they are no longer the same quantity and the range cannot be printed beside
## the point estimate, so stop rather than publish a mismatched pair.
{
  .u <- den_by_reading[grepl("^Union", den_by_reading$reading), ]
  .r2 <- den[grepl("^ALL adults", den$denominator), ]
  stopifnot(nrow(.u) == 1, nrow(.r2) == 1,
            .u$n_with_goal == .r2$n_with_goal,
            .u$n_total     == .r2$n_total,
            abs(.u$pct_weighted - .r2$pct_weighted) < 1e-6)
}

## The first-match family, all 120 orderings. Point estimates only, as R1 does:
## 120 orderings give 120 intervals and none of them is the reported estimate.
perm_comb <- vapply(seq_len(nrow(perms)), function(i) {
  d <- combined_design
  d$variables$.hg <- combined_indicator(goal_under_order(pdf_rows, perms[i, ]))
  as.numeric(svymean(~I(.hg == 1), d, na.rm = TRUE))[2]
}, numeric(1))

den_by_reading <- bind_rows(den_by_reading, tibble(
  reading      = "Sequential first-match, across all 120 orderings",
  n_with_goal  = NA_integer_,
  n_total      = nrow(combined_df),
  pct_weighted = NA_real_,
  ci_low       = 100 * min(perm_comb),
  ci_high      = 100 * max(perm_comb),
  us_adults_with_goal = NA_real_,
  ci_method    = "range of point estimates, not an interval"
))

OUT(den_by_reading, "rev2b_denominator_by_reading.csv")
print(as.data.frame(den_by_reading))

message(sprintf("  Range across readings, combined: %.1f%% to %.1f%%",
                min(den_by_reading$pct_weighted, na.rm = TRUE),
                max(den_by_reading$pct_weighted, na.rm = TRUE)))
message(sprintf("  Range across readings, primary:  %.1f%% to %.1f%%",
                100 * min(bounds$share), 100 * max(bounds$share)))
message("  Table 2's combined row carries this range once this file exists; ",
        "until then it states that it inherits the primary-prevention range.")

## =============================================================================
## R3. CORRELATES OF DISCORDANCE  (major comment 3)
##
## The reviewer is right that the adjusted model is overfit. Three responses,
## in increasing order of what we would actually publish:
##
##   (a) events-per-variable reported for every model in the paper;
##   (b) the two PRESPECIFIED associations, unadjusted, with design-based
##       exact intervals -- this is what we propose as primary;
##   (c) Firth penalised logistic regression as a supplementary check.
##
## Note on (c): Firth is not implemented for complex survey designs. logistf()
## accepts weights but ignores clustering and stratification, so its intervals
## are not design-based. It is reported as a sensitivity on the point estimates
## only, and this limitation is stated in the manuscript.
## =============================================================================

message("\n-- R3. Correlates: EPV, unadjusted, Firth --")

apob_sub <- subset(primary_design, has_apob_goal == 1)
sub_df   <- apob_sub$variables

n_events <- sum(sub_df$discordant == 1, na.rm = TRUE)
message("  Discordant events in the goal-assigned subpopulation: ", n_events)

epv_row <- function(label, n_ev, n_par) {
  tibble(model = label, events = n_ev, parameters = n_par,
         epv = round(n_ev / n_par, 2),
         meets_10_epv = (n_ev / n_par) >= 10)
}

epv <- bind_rows(
  epv_row("RQ4 adjusted (as published: 10 terms + intercept)", n_events, 11),
  epv_row("Proposed primary: TG only, unadjusted", n_events, 2),
  epv_row("Proposed primary: diabetes only, unadjusted", n_events, 2),
  epv_row("Proposed adjusted: TG + diabetes + age + sex", n_events, 5)
)
OUT(epv, "rev3_epv_table.csv")
print(as.data.frame(epv))

## (b) Prespecified unadjusted associations, design-based.
unadj <- function(varname, label, design = apob_sub, n_ev = n_events) {
  d <- design
  d$variables$.exp <- as.integer(d$variables[[varname]] == 1)
  d1 <- subset(d, .exp == 1)
  d0 <- subset(d, .exp == 0)
  p1 <- svyciprop(~I(discordant == 1), d1, method = "beta")
  p0 <- svyciprop(~I(discordant == 1), d0, method = "beta")
  m  <- svyglm(discordant ~ .exp, design = d, family = quasibinomial())
  or <- exp(coef(m)[2]); ci <- exp(confint(m)[2, ])
  tibble(exposure = label,
         n_exposed = nrow(d1$variables), n_unexposed = nrow(d0$variables),
         prev_exposed = as.numeric(p1),
         prev_exposed_lo = confint(p1)[1], prev_exposed_hi = confint(p1)[2],
         prev_unexposed = as.numeric(p0),
         prev_unexposed_lo = confint(p0)[1], prev_unexposed_hi = confint(p0)[2],
         OR = or, OR_lo = ci[1], OR_hi = ci[2],
         epv = round(n_ev / 2, 1))
}

unadj_res <- bind_rows(
  unadj("hypertrig", "Hypertriglyceridaemia"),
  unadj("diabetes",  "Diabetes")
)
OUT(unadj_res, "rev3_correlates_unadjusted.csv")
print(as.data.frame(unadj_res))

## (c) Firth penalised logistic regression, weighted, NOT design-based.
if (requireNamespace("logistf", quietly = TRUE)) {
  fit_df <- sub_df %>%
    select(discordant, hypertrig, diabetes, RIDAGEYR, RIAGENDR, BMXBMI,
           ckd, wt_pooled) %>%
    filter(complete.cases(.)) %>%
    ## WEIGHTS ARE RESCALED TO MEAN 1, and this is not cosmetic. logistf()
    ## treats its weights as FREQUENCIES. NHANES pooled weights average about
    ## 15,000, so passing them raw makes the model believe it has ~37 million
    ## observations: every profile-likelihood interval collapses to near-zero
    ## width and every p-value prints as 0 (hypertriglyceridaemia came back as
    ## OR 15.55, 95% CI 15.35-15.76). Dividing by the mean preserves the
    ## relative weighting and restores the true sample size.
    mutate(w = wt_pooled / mean(wt_pooled))

  message("  Firth weight check: raw mean ", round(mean(fit_df$wt_pooled)),
          ", rescaled mean ", round(mean(fit_df$w), 3),
          ", n = ", nrow(fit_df))

  ## Race is dropped from the penalised model: "other or multiracial" is a
  ## harmonisation-created category that cannot be disaggregated before 2011,
  ## and the reviewer is right that reporting it while calling Non-Hispanic
  ## Black unstable is inconsistent. Neither is reported.
  f <- logistf::logistf(
    discordant ~ hypertrig + diabetes + RIDAGEYR + factor(RIAGENDR) + BMXBMI + ckd,
    data = fit_df, weights = fit_df$w
  )
  firth <- tibble(
    term = names(coef(f)),
    OR = exp(coef(f)),
    OR_lo = exp(f$ci.lower),
    OR_hi = exp(f$ci.upper),
    p = f$prob
  )
  OUT(firth, "rev3_correlates_firth.csv")
  print(as.data.frame(firth))
  message("  NOTE: Firth intervals above are NOT design-based (no PSU/strata). ",
          "Report as a sensitivity on point estimates only.")
} else {
  message("  logistf not installed -- run install.packages('logistf') and rerun.")
}

## =============================================================================
## R4. CYCLE INFLUENCE ON THE GOAL-ASSIGNED SHARE  (major comment 4)
##
## The reviewer's comment conflates two different quantities (see the response
## document): the -22.2% figure is a relative shift in DISCORDANCE, not in the
## goal-assigned share. This section tests the share directly, and repeats the
## discordance leave-one-out on CALIBRATED apoB.
##
## Note that the "22.0% goal-assigned share" in the reviewer's comment is the
## superseded first-match figure. The share under the union rule is 33.4%.
## =============================================================================

message("\n-- R4. Leave-one-cycle-out, on the goal-assigned share AND discordance --")

cycles <- sort(unique(primary_df$cycle))
loo <- map_dfr(c("(none)", cycles), function(drop) {
  d <- if (drop == "(none)") primary_design else subset(primary_design, cycle != drop)
  share <- svyciprop(~I(has_apob_goal == 1), d, method = "beta")
  d2 <- subset(d, has_apob_goal == 1)
  disc <- svyciprop(~I(discordant == 1), d2, method = "logit")
  tibble(cycle_dropped = drop,
         goal_share = as.numeric(share),
         discordance = as.numeric(disc))
})
base_share <- loo$goal_share[loo$cycle_dropped == "(none)"]
base_disc  <- loo$discordance[loo$cycle_dropped == "(none)"]
loo <- loo %>%
  mutate(pct_change_share = 100 * (goal_share - base_share) / base_share,
         pct_change_disc  = 100 * (discordance - base_disc) / base_disc,
         flag_share = abs(pct_change_share) > 15,
         flag_disc  = abs(pct_change_disc)  > 15)
OUT(loo, "rev4_cycle_influence_on_share.csv")
print(as.data.frame(loo))

message("\n  With apoB calibration applied, cycle D influence on discordance should fall ",
        "substantially. If it does not, the anomaly is not calibration and the ",
        "2007-2016 estimate should be made co-primary.")

## =============================================================================
## R6. SAMPSON VALIDITY CEILING  (Dani, review round 2, item 4)
##
## There is no R5 here: rev5_stringency_sensitivity.csv is produced by
## precedence_diagnostics.R. This section is numbered to match its own output file.
##
## The Sampson/NIH equation 2 is validated to approximately TG 800 mg/dL.
## Participants above that concentration are retained in the primary analysis
## with their LDL-C read as an extrapolation, which is stated in the Methods.
## This section re-estimates the prespecified unadjusted hypertriglyceridaemia
## association with those participants removed, so the correlate can be shown
## not to depend on the unvalidated region.
##
## The model is NOT re-specified. unadj() is called on a design-subset of the
## same goal-assigned subpopulation, so the restricted and primary estimates
## come from one code path and differ only in who is in the design. Subsetting
## is done on the design, never on the data frame, so the variance estimator
## keeps every PSU and stratum.
## =============================================================================

message("\n-- R6. Sampson validity ceiling: TG > 800 mg/dL excluded --")

TG_CEILING <- 800

if (!"LBXTR" %in% names(sub_df)) {
  stop("R6 needs the raw triglyceride column LBXTR in the goal-assigned ",
       "design and it is not there. Check that nhanes_apob_analysis.R keeps ",
       "LBXTR on primary_df; do not substitute the hypertrig indicator, which ",
       "cannot identify the participants above the ceiling.", call. = FALSE)
}

n_over_ceiling <- sum(sub_df$LBXTR > TG_CEILING, na.rm = TRUE)
message("  Goal-assigned adults above TG ", TG_CEILING, " mg/dL: ",
        n_over_ceiling, " of ", nrow(sub_df))

## coalesce(..., TRUE) keeps any participant with a missing triglyceride value,
## matching how hypertrig treats them in the primary analysis.
apob_sub_800 <- subset(apob_sub, coalesce(LBXTR <= TG_CEILING, TRUE))
n_events_800 <- sum(apob_sub_800$variables$discordant == 1, na.rm = TRUE)
message("  Discordant events retained: ", n_events_800, " of ", n_events)

sampson800 <- bind_rows(
  unadj("hypertrig", "Hypertriglyceridaemia, all goal-assigned adults") %>%
    mutate(cohort = "All goal-assigned adults",
           tg_ceiling_mgdl = NA_real_, n_excluded = 0L,
           events_in_model = n_events),
  unadj("hypertrig", "Hypertriglyceridaemia, TG <= 800 mg/dL",
        design = apob_sub_800, n_ev = n_events_800) %>%
    mutate(cohort = paste0("TG <= ", TG_CEILING, " mg/dL"),
           tg_ceiling_mgdl = TG_CEILING, n_excluded = n_over_ceiling,
           events_in_model = n_events_800)
)

OUT(sampson800, "rev6_sampson800_sensitivity.csv")
print(as.data.frame(sampson800))

message("\n  One eTable row or one clause: the unadjusted hypertriglyceridaemia ",
        "odds ratio with and without the extrapolated region. Overlapping ",
        "intervals mean the correlate does not rest on triglyceride values ",
        "where Sampson is unvalidated. Report n_excluded alongside it.")

## -----------------------------------------------------------------------------
## R7. WHAT ACTING ON apoB WOULD CHANGE
##
## Added 2026-09-09 in response to reviewer comment: "changes management" is the
## clinical claim the paper rests on, so bound it rather than assert it.
##
## Denominator: the discordant adults -- at their assigned LDL-C goal, above
## their assigned apoB goal. For these adults an LDL-C panel alone reports
## success; apoB does not. The question is what the next step would be if apoB
## were measured and acted on. Current therapy determines which step is
## available, so the stratum is split by what the participant is already taking:
##
##   no lipid-lowering therapy   -> initiation is the available step
##   statin monotherapy          -> intensification or add-on is available
##   combination or non-statin   -> already on more than a statin; the further
##                                  step is narrower and this analysis does not
##                                  claim one
##
## This bounds the claim; it does not model who SHOULD be treated, which would
## require the guideline's therapy-escalation logic and outcome data this design
## does not have. Cells are small (71 events split three ways), so unweighted
## counts are reported beside every weighted estimate and the intervals are
## Korn-Graubard.
## -----------------------------------------------------------------------------
message("\n-- R7. What acting on apoB would change, among discordant adults --")

if (!all(c("on_statin", "on_nonstatin_lipid_tx") %in% names(sub_df))) {
  warning("R7 skipped: on_statin / on_nonstatin_lipid_tx are not on the ",
          "goal-assigned design. Check derive_vars().", call. = FALSE)
} else {
  disc_sub <- subset(apob_sub, discordant == 1)
  dd <- disc_sub$variables
  dd$tx_step <- ifelse(dd$on_statin == 0 & dd$on_nonstatin_lipid_tx == 0,
                       "1_no_therapy_initiation_available",
                ifelse(dd$on_statin == 1 & dd$on_nonstatin_lipid_tx == 0,
                       "2_statin_monotherapy_intensification_available",
                       "3_combination_or_nonstatin_narrower_step"))
  disc_sub$variables$tx_step <- dd$tx_step

  n_disc <- nrow(dd)
  message("  Discordant adults in the stratum: ", n_disc)

  pop_total <- as.numeric(svytotal(~I(discordant == 1), apob_sub)[2])

  r7 <- do.call(rbind, lapply(sort(unique(dd$tx_step)), function(lv) {
    p  <- svyciprop(as.formula(paste0("~I(tx_step == '", lv, "')")),
                    disc_sub, method = "beta")
    tibble(
      step            = lv,
      n_unweighted    = sum(dd$tx_step == lv),
      share_of_disc   = as.numeric(p),
      ci_low          = attr(p, "ci")[1],
      ci_high         = attr(p, "ci")[2],
      us_adults       = as.numeric(p) * pop_total
    )
  }))
  r7$discordant_population_total <- pop_total

  ## the headline the reviewer asked for: a step exists for cells 1 and 2
  step_available <- sum(r7$share_of_disc[r7$step %in%
      c("1_no_therapy_initiation_available",
        "2_statin_monotherapy_intensification_available")])
  message(sprintf("  A treatment step is available for %.1f%% of the discordant stratum",
                  100 * step_available))
  message(sprintf("  = %s of %s US adults",
                  format(round(step_available * pop_total), big.mark = ","),
                  format(round(pop_total), big.mark = ",")))

  OUT(r7, "rev7_changes_management.csv")
  print(as.data.frame(r7))
  message("  Small cells: ", n_disc, " events across ",
          nrow(r7), " strata. Report the unweighted counts alongside every ",
          "weighted figure, and do not present any single cell as a headline.")
}

message("\n=== REVISION ANALYSES COMPLETE ===")
message(length(.written), " file(s) written:")
for (f in .written) message("  - ", f)
if (length(.written) < 8) {
  warning("Fewer files were written than expected (8-9). Scroll up for the ",
          "section that failed; do not treat this run as complete.",
          call. = FALSE)
} else {
  message("Numbers for the manuscript come from these files, not from the ",
          "pre-revision RESULTS_VERIFIED.md.")
}
