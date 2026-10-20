## =============================================================================
## reviewer_strata_analysis.R -- RQ14: residual apoB elevation by triglyceride
## category, by number of Figure 1 rows matched, and by lipid-lowering therapy.
## =============================================================================
## ADDED 16 October 2026 in response to reviewer comments, which asked whether
## residual apoB elevation is concentrated in particular groups:
##   * by triglyceride band (<150, 150-199, 200-299, 300-499, >=500 mg/dL)
##   * by whether an adult matched one guideline row or more than one
##   * by lipid-lowering therapy (taking vs not taking)
##
## Primary prevention, adults in an apoB measurement group. For each stratum:
##   above_apob      above the stated apoB value
##   both_goals      above it, among those at BOTH LDL-C and non-HDL-C goals
##   ldl_goal_any    above it, among those at the LDL-C goal (any non-HDL-C)
##   discordance     at LDL-C goal AND above the apoB value, among everyone in
##                   the stratum (the 2.54% outcome, by stratum)
##
## METHODS as in cascade_analysis.R: each cell's design is the primary design
## subset to the SEQNs of the mask that produced its n; Korn-Graubard
## intervals; every cell gated by nchs_present(). Small cells are EXPECTED to
## fail the gate; that is recorded, not acted on.
##
## Rows matched uses the definition behind the 56.3% in eFigure 3
## (pathway_and_sensitivity_analyses.R): the diabetes, hypertriglyceridemia
## and PREVENT rows, each counted when that row states a value.
##
## CHECKS: the strata of each variable partition the measurement group, so
## their denominators must sum to the whole (3,324 / 259 / 349). A mismatch
## stops the run.
##
## NEEDS primary_df and nchs_present() from nhanes_apob_analysis.R. Run by
## run_all.R after cascade_analysis.R. Also writes analysis_frames_cache.rds
## (primary_df, secondary_df, nchs_present) so that a later follow-up
## analysis can be run without the 60-minute rebuild:
##   x <- readRDS("analysis_frames_cache.rds"); list2env(x, globalenv())
##   source("reviewer_strata_analysis.R")
## The cache is an output, not code; it is not part of the repository.
##
## WRITES rq14_reviewer_strata.csv
## =============================================================================

suppressPackageStartupMessages({ library(dplyr); library(survey) })

.need14 <- c("primary_df", "nchs_present")
.miss14 <- .need14[!vapply(.need14, exists, logical(1), envir = globalenv())]
if (length(.miss14) && file.exists("analysis_frames_cache.rds")) {
  list2env(readRDS("analysis_frames_cache.rds"), globalenv())
  .miss14 <- .need14[!vapply(.need14, exists, logical(1), envir = globalenv())]
}
if (length(.miss14))
  stop("reviewer_strata_analysis.R needs ", paste(.miss14, collapse = ", "),
       " from nhanes_apob_analysis.R. Run it through run_all.R.", call. = FALSE)

if (exists("secondary_df", envir = globalenv()))
  saveRDS(list(primary_df = primary_df, secondary_df = secondary_df,
               nchs_present = nchs_present), "analysis_frames_cache.rds")

.cols14 <- c("SEQN", "SDMVPSU", "SDMVSTRA", "wt_pooled", "has_apob_goal",
             "above_apob_goal", "at_ldl_goal", "goal_pair_cell", "discordant",
             "LBXTR", "g_dm_apob", "g_htg_apob", "g_pp_apob", "on_any_lipid_tx")
.m14 <- setdiff(.cols14, names(primary_df))
if (length(.m14)) stop("reviewer_strata_analysis.R: primary_df lacks ",
                       paste(.m14, collapse = ", "), call. = FALSE)

d14 <- primary_df[, .cols14] %>%
  mutate(
    tg_band = cut(LBXTR, c(-Inf, 150, 200, 300, 500, Inf), right = FALSE,
                  labels = c("<150", "150-199", "200-299", "300-499", ">=500")),
    rows_matched = ifelse((!is.na(g_dm_apob)) + (!is.na(g_htg_apob)) +
                            (!is.na(g_pp_apob)) > 1, "more_than_one", "one"),
    lipid_tx = ifelse(coalesce(on_any_lipid_tx == 1, FALSE),
                      "on_lipid_tx", "not_on_lipid_tx"))
des14 <- svydesign(id = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~wt_pooled,
                   nest = TRUE, data = d14)

.cell14 <- function(variable, level, outcome, den, num) {
  den <- coalesce(den, FALSE); num <- coalesce(num, FALSE) & den
  n_den <- sum(den); n_num <- sum(num)
  est <- lo <- hi <- eff <- dfree <- NA_real_
  if (n_den > 0) {
    dsub <- subset(des14, SEQN %in% d14$SEQN[den])
    dsub <- update(dsub, .num14 = as.integer(dsub$variables$SEQN %in% d14$SEQN[num]))
    ci <- try(svyciprop(~I(.num14 == 1), dsub, method = "beta"), silent = TRUE)
    if (!inherits(ci, "try-error")) {
      est <- as.numeric(ci); lo <- confint(ci)[1]; hi <- confint(ci)[2]
      m <- try(suppressWarnings(svymean(~I(.num14 == 1), dsub, deff = TRUE)),
               silent = TRUE)
      if (!inherits(m, "try-error")) {
        dff <- as.numeric(deff(m))[2]
        if (is.finite(dff) && dff > 0) eff <- n_den / dff
      }
      dfree <- suppressWarnings(as.numeric(degf(dsub)))
    }
  }
  tibble(variable = variable, level = level, outcome = outcome,
         n_denominator = n_den, n_numerator = n_num,
         estimate = est, ci_low = lo, ci_high = hi, effective_n = eff,
         df = dfree, ci_method = "korn-graubard",
         nchs_reportable = nchs_present(est, lo, hi, n_den, eff,
                                        df = ifelse(is.na(dfree), 94, dfree),
                                        n_events = n_num))
}

G <- d14$has_apob_goal == 1
A <- d14$above_apob_goal == 1
.strata14 <- list(
  tg_band      = levels(d14$tg_band),
  rows_matched = c("one", "more_than_one"),
  lipid_tx     = c("on_lipid_tx", "not_on_lipid_tx"))

rq14 <- bind_rows(lapply(names(.strata14), function(v) {
  bind_rows(lapply(.strata14[[v]], function(lv) {
    S <- G & as.character(d14[[v]]) == lv
    bind_rows(
      .cell14(v, lv, "above_apob",   S, A),
      .cell14(v, lv, "both_goals",   S & d14$goal_pair_cell %in% "both", A),
      .cell14(v, lv, "ldl_goal_any", S & d14$at_ldl_goal == 1, A),
      .cell14(v, lv, "discordance",  S, d14$discordant == 1))
  }))
}))

## ---- partition checks --------------------------------------------------------
local({
  tot <- c(above_apob   = sum(G, na.rm = TRUE),
           both_goals   = sum(G & d14$goal_pair_cell %in% "both", na.rm = TRUE),
           ldl_goal_any = sum(G & d14$at_ldl_goal == 1, na.rm = TRUE))
  for (v in names(.strata14)) for (o in names(tot)) {
    s <- sum(rq14$n_denominator[rq14$variable == v & rq14$outcome == o])
    if (s != tot[[o]])
      stop("reviewer_strata_analysis.R: ", v, " strata sum to ", s, " for ",
           o, ", not ", tot[[o]], ".", call. = FALSE)
  }
  message("  rq14 strata partition the measurement group: ",
          paste(names(tot), tot, sep = " = ", collapse = ", "))
})

write.csv(rq14, "rq14_reviewer_strata.csv", row.names = FALSE)
message("RQ14: reviewer strata written to rq14_reviewer_strata.csv")
print(as.data.frame(rq14[, c("variable", "level", "outcome", "n_denominator",
                             "n_numerator", "estimate", "ci_low", "ci_high",
                             "nchs_reportable")]), digits = 3)
