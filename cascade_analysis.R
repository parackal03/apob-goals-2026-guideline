## =============================================================================
## cascade_analysis.R -- RQ13: the measurement cascade, overall / primary /
## secondary prevention.
## =============================================================================
## ADDED October 2026 for the reviewer's comment, which asks for one flow through the
## whole cohort and then the same flow repeated for each arm:
##
##   adults aged 30-79 years
##     -> in an apoB measurement group (a row of Figure 1 states an apoB value)
##       -> above their stated apoB value
##         -> at BOTH their LDL-C and paired non-HDL-C goals: above apoB value
##         -> at their LDL-C goal (whatever the non-HDL-C): above apoB value
##         -> at their LDL-C goal but NOT their non-HDL-C goal: above apoB value
##         -> at EITHER goal: above apoB value  (the literal "and/or" reading)
##
## for three cohorts: OVERALL (primary + secondary), PRIMARY, SECONDARY.
##
## WHAT IS NEW. The overall column, and the goal-pair rows for secondary
## prevention, have not been computed before. Everything for primary
## prevention, and the measurement-group and above-value rows for secondary
## prevention and overall, already exist in other outputs; this script
## recomputes them from the same objects and STOPS if any point estimate
## differs from the published one, so the new table cannot drift from the old
## numbers.
##
## METHODS, the same as everywhere else in the pipeline:
##   * designs built on the cohort data frames, as primary_design and
##     secondary_design are; the overall design is the one
##     precedence_and_model_analyses.R uses for the combined denominator;
##   * each cell's design is built from the SEQNs of the very mask that
##     produced its unweighted n (as in RQ11), so n and estimate describe the
##     same people;
##   * Korn-Graubard intervals (svyciprop, method = "beta");
##   * every cell gated by nchs_present(), defined in nhanes_apob_analysis.R.
##     Small secondary-prevention cells are expected to fail it; that is
##     recorded, never acted on here. Suppression is a document decision.
##
## NEEDS, in the global environment: primary_df, secondary_df and
## nchs_present(), all from nhanes_apob_analysis.R. Run by run_all.R after
## precedence_and_model_analyses.R, which writes the combined-denominator file
## this script validates against.
##
## WRITES rq13_cascade_by_cohort.csv.
## =============================================================================

suppressPackageStartupMessages({ library(dplyr); library(survey) })

.need13 <- c("primary_df", "secondary_df", "nchs_present")
.miss13 <- .need13[!vapply(.need13, exists, logical(1), envir = globalenv())]
if (length(.miss13))
  stop("cascade_analysis.R needs ", paste(.miss13, collapse = ", "),
       " from nhanes_apob_analysis.R. Run it through run_all.R.", call. = FALSE)

## Lonely-PSU handling is inherited from nhanes_apob_analysis.R ("adjust"),
## which sets it for the whole run; this script does not change it.

## ---- cohorts -------------------------------------------------------------------
.cols13 <- c("SEQN", "SDMVPSU", "SDMVSTRA", "wt_pooled", "has_apob_goal",
             "above_apob_goal", "at_ldl_goal", "goal_pair_cell")
.cohort13 <- function(df, arm) {
  miss <- setdiff(.cols13, names(df))
  if (length(miss)) stop("cascade_analysis.R: ", arm, " cohort lacks ",
                         paste(miss, collapse = ", "), call. = FALSE)
  df[, .cols13] %>% mutate(.arm = arm)
}
.prim13 <- .cohort13(primary_df,   "primary")
.sec13  <- .cohort13(secondary_df, "secondary")
if (length(intersect(.prim13$SEQN, .sec13$SEQN)))
  stop("cascade_analysis.R: a participant is in both cohorts.", call. = FALSE)
.all13 <- bind_rows(.prim13, .sec13)

.design13 <- function(d) svydesign(id = ~SDMVPSU, strata = ~SDMVSTRA,
                                   weights = ~wt_pooled, nest = TRUE, data = d)
.cohorts13 <- list(overall = .all13, primary = .prim13, secondary = .sec13)

## ---- one cell ------------------------------------------------------------------
## Proportion with `num` among those with `den`, both logical vectors over the
## cohort frame `d`, estimated on `des` (the cohort's own design).
.cell13 <- function(cohort, step, d, des, den, num) {
  den <- coalesce(den, FALSE); num <- coalesce(num, FALSE) & den
  n_den <- sum(den); n_num <- sum(num)
  est <- lo <- hi <- dff <- eff <- dfree <- tot <- NA_real_
  if (n_den > 0) {
    dsub <- subset(des, SEQN %in% d$SEQN[den])
    dsub <- update(dsub, .num13 = as.integer(dsub$variables$SEQN %in% d$SEQN[num]))
    ci <- try(svyciprop(~I(.num13 == 1), dsub, method = "beta"), silent = TRUE)
    if (!inherits(ci, "try-error")) {
      est <- as.numeric(ci); lo <- confint(ci)[1]; hi <- confint(ci)[2]
      m <- try(suppressWarnings(svymean(~I(.num13 == 1), dsub, deff = TRUE)),
               silent = TRUE)
      if (!inherits(m, "try-error")) {
        dff <- as.numeric(deff(m))[2]
        if (is.finite(dff) && dff > 0) eff <- n_den / dff
      }
      dfree <- suppressWarnings(as.numeric(degf(dsub)))
      tt <- try(svytotal(~I(.num13 == 1), dsub), silent = TRUE)
      if (!inherits(tt, "try-error")) tot <- as.numeric(coef(tt))[2]
    }
  }
  tibble(cohort = cohort, step = step,
         n_denominator = n_den, n_numerator = n_num,
         estimate = est, ci_low = lo, ci_high = hi,
         ci_width_abs = hi - lo, design_effect = dff, effective_n = eff,
         df = dfree, us_adults_numerator = tot, ci_method = "korn-graubard",
         nchs_reportable = nchs_present(est, lo, hi, n_den, eff,
                                        df = ifelse(is.na(dfree), 94, dfree),
                                        n_events = n_num))
}

## ---- the cascade -----------------------------------------------------------------
rq13_cascade <- bind_rows(lapply(names(.cohorts13), function(cn) {
  d   <- .cohorts13[[cn]]
  des <- .design13(d)
  G   <- d$has_apob_goal == 1
  A   <- d$above_apob_goal == 1
  gp  <- d$goal_pair_cell
  bind_rows(
    .cell13(cn, "1_in_measurement_group",       d, des, rep(TRUE, nrow(d)), G),
    .cell13(cn, "2_above_apob_value",           d, des, G, A),
    .cell13(cn, "3_both_goals_above_apob",      d, des, G & gp %in% "both", A),
    .cell13(cn, "4_ldl_goal_any_above_apob",    d, des, G & d$at_ldl_goal == 1, A),
    .cell13(cn, "5_ldl_goal_only_above_apob",   d, des, G & gp %in% "ldl_only", A),
    .cell13(cn, "6_either_goal_above_apob",     d, des,
            G & gp %in% c("both", "ldl_only", "nonhdl_only"), A)
  )
}))

## ---- internal consistency ------------------------------------------------------------
## both + ldl_only is exactly "at LDL-C goal", in every cohort; and the overall
## counts are the sums of the arms. Either failing means a definition drifted.
local({
  for (cn in names(.cohorts13)) {
    r <- rq13_cascade[rq13_cascade$cohort == cn, ]
    s <- function(st, col) r[[col]][r$step == st]
    if (s("3_both_goals_above_apob", "n_denominator") +
        s("5_ldl_goal_only_above_apob", "n_denominator") !=
        s("4_ldl_goal_any_above_apob", "n_denominator"))
      stop("cascade_analysis.R: both + ldl_only != at LDL-C goal in ", cn,
           call. = FALSE)
  }
  for (st in unique(rq13_cascade$step)) for (col in c("n_denominator", "n_numerator")) {
    v <- setNames(rq13_cascade[[col]][rq13_cascade$step == st],
                  rq13_cascade$cohort[rq13_cascade$step == st])
    if (v[["overall"]] != v[["primary"]] + v[["secondary"]])
      stop("cascade_analysis.R: overall ", col, " for ", st,
           " is not primary + secondary.", call. = FALSE)
  }
})

## ---- agreement with the published estimates ------------------------------------------
## Point estimates only: they do not depend on the interval method. Each file
## is from the same run. A missing file is reported, not fatal, so this
## script can be developed against a partial run; a MISMATCH is fatal.
local({
  chk <- function(file, pick, cohort, step, label) {
    if (!file.exists(file)) {
      message("  rq13 check skipped (", file, " not found): ", label); return()
    }
    want <- pick(read.csv(file, stringsAsFactors = FALSE))
    got  <- rq13_cascade$estimate[rq13_cascade$cohort == cohort &
                                  rq13_cascade$step == step]
    if (length(want) != 1 || is.na(got) || abs(as.numeric(want) - got) > 1e-6)
      stop("cascade_analysis.R: ", label, " = ", signif(got, 8),
           " here but ", signif(as.numeric(want), 8), " in ", file, ".",
           call. = FALSE)
    message("  rq13 agrees with ", file, ": ", label)
  }
  chk("rev2_denominator_variants.csv",
      function(x) x$pct_weighted[grepl("^ALL adults", x$denominator)] / 100,
      "overall", "1_in_measurement_group", "overall measurement-group share")
  chk("rev2_denominator_variants.csv",
      function(x) x$pct_weighted[grepl("^Primary", x$denominator)] / 100,
      "primary", "1_in_measurement_group", "primary measurement-group share")
  chk("rev2_denominator_variants.csv",
      function(x) x$pct_weighted[grepl("^Secondary", x$denominator)] / 100,
      "secondary", "1_in_measurement_group", "secondary measurement-group share")
  chk("rq1b_overall_apob_nonattainment.csv", function(x) x$estimate[1],
      "primary", "2_above_apob_value", "primary above apoB value")
  chk("rq6b_secondary_prevention_apob_nonattainment.csv", function(x) x$estimate[1],
      "secondary", "2_above_apob_value", "secondary above apoB value")
  chk("rq1c_reclassification_at_ldl_goal.csv", function(x) x$estimate[1],
      "primary", "4_ldl_goal_any_above_apob", "primary at LDL-C goal, above apoB")
  chk("rq11_goal_pair_partition.csv",
      function(x) x$estimate[x$population == "all" & x$cell == "both"],
      "primary", "3_both_goals_above_apob", "primary at both goals, above apoB")
  chk("rq11_goal_pair_partition.csv",
      function(x) x$estimate[x$population == "all" & x$cell == "ldl_only"],
      "primary", "5_ldl_goal_only_above_apob", "primary LDL-C goal only, above apoB")
  chk("rq11_goal_pair_partition.csv",
      function(x) x$estimate[x$population == "all" & x$cell == "either"],
      "primary", "6_either_goal_above_apob", "primary either goal, above apoB")
})

write.csv(rq13_cascade, "rq13_cascade_by_cohort.csv", row.names = FALSE)
message("RQ13: measurement cascade by cohort written to rq13_cascade_by_cohort.csv")
print(as.data.frame(rq13_cascade[, c("cohort", "step", "n_denominator",
                                     "n_numerator", "estimate", "ci_low",
                                     "ci_high", "nchs_reportable")]),
      digits = 3)
if (exists("DIAG", mode = "function"))
  for (i in seq_len(nrow(rq13_cascade)))
    DIAG(paste0("rq13_", rq13_cascade$cohort[i], "_", rq13_cascade$step[i], "_n"),
         rq13_cascade$n_denominator[i], "RQ13 cascade denominator")
