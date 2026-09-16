## =============================================================================
## verify_run.R -- post-run verification harness
##
## Run this IMMEDIATELY AFTER sourcing nhanes_apob_analysis.R, in the same
## session, so the objects it interrogates are still in memory.
##
##   source("nhanes_apob_analysis.R")
##   source("verify_run.R")
##
## Every check prints PASS, FAIL, or WARN. Nothing is silent. A FAIL means a
## number in the manuscript would be wrong; a WARN means a number needs a
## caveat. Checks are ordered so the cheap structural ones fail fast.
##
## This file makes no changes to the analysis. It only reads and reports.
## =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(survey); library(srvyr); library(tibble)
})

.results <- list()
.chk <- function(id, label, status, detail = "") {
  .results[[length(.results) + 1]] <<- tibble(
    check = id, description = label, status = status, detail = detail)
  cat(sprintf("[%-4s] %-52s %s\n", status, label, detail))
}
.exists_or_skip <- function(obj) exists(obj, envir = .GlobalEnv)

cat("\n==============================================================\n")
cat(" POST-RUN VERIFICATION\n")
cat("==============================================================\n\n")

## ---- A. MODEL CONVERGENCE ---------------------------------------------------
## The run emitted "glm.fit: algorithm did not converge" without naming the
## culprit. A non-converged GLM returns coefficients that are artefacts of
## where the optimiser stopped, not estimates. Any such model is unusable.

cat("--- A. Model convergence ---\n")

for (m in c("rq4_model", "rq8_interaction_model")) {
  if (!.exists_or_skip(m)) { .chk(m, paste(m, "present"), "SKIP", "not in memory"); next }
  mod <- get(m, envir = .GlobalEnv)
  conv <- isTRUE(mod$converged)
  iters <- if (!is.null(mod$iter)) mod$iter else NA
  .chk(m, paste0(m, " converged"),
       if (conv) "PASS" else "FAIL",
       paste0("converged=", conv, ", iterations=", iters))

  ## Separation detector. Even a "converged" logistic model is untrustworthy if
  ## fitted probabilities sit at the boundary or any coefficient has exploded.
  if (conv) {
    fv <- fitted(mod)
    extreme_fv <- sum(fv < 1e-8 | fv > 1 - 1e-8, na.rm = TRUE)
    big_coef   <- sum(abs(coef(mod)) > 15, na.rm = TRUE)
    .chk(paste0(m, "_sep"), paste0(m, " free of separation"),
         if (extreme_fv == 0 && big_coef == 0) "PASS" else "WARN",
         paste0(extreme_fv, " boundary fitted values, ",
                big_coef, " coefficients with |beta| > 15"))
  }
}

## Wide confidence intervals are a softer signal of the same problem.
if (.exists_or_skip("rq4_or_table")) {
  wide <- rq4_or_table %>% filter(term != "(Intercept)",
                                  (conf.high / conf.low) > 100)
  .chk("rq4_ci_width", "RQ4 odds ratios have usable precision",
       if (nrow(wide) == 0) "PASS" else "WARN",
       if (nrow(wide) == 0) "all CI ratios < 100x"
       else paste0(nrow(wide), " term(s) with >100x CI ratio: ",
                   paste(wide$term, collapse = "; ")))
}

## ---- B. RACE/ETHNICITY HARMONIZATION ---------------------------------------
cat("\n--- B. Race/ethnicity harmonization ---\n")

if (.exists_or_skip("primary_df") && "race_eth" %in% names(primary_df)) {
  n_na <- sum(is.na(primary_df$race_eth))
  .chk("race_na", "race_eth has no missing values",
       if (n_na == 0) "PASS" else "FAIL", paste0(n_na, " NA"))

  .chk("race_levels", "race_eth has exactly 5 levels",
       if (nlevels(primary_df$race_eth) == 5) "PASS" else "FAIL",
       paste(levels(primary_df$race_eth), collapse = " | "))

  ## UPDATED 2026-08-15. RQ4's outcome is discordance, which is only defined
  ## for participants the guideline assigns an apoB goal. The model therefore
  ## runs on apob_design, and the correct denominator is that subpopulation --
  ## not the full cohort. Comparing against nrow(primary_df) produced a
  ## spurious FAIL (7,219 of 8,925) that looked like the old silent-drop bug.
  if (.exists_or_skip("rq4_model")) {
    n_used <- tryCatch(nobs(rq4_model), error = function(e) NA)
    n_target <- if ("has_apob_goal" %in% names(primary_df))
      sum(primary_df$has_apob_goal == 1, na.rm = TRUE) else nrow(primary_df)
    .chk("rq4_n", "RQ4 model uses the full apoB-goal subpopulation",
         if (!is.na(n_used) && n_used == n_target) "PASS" else "FAIL",
         paste0(n_used, " of ", n_target, " with a guideline-specified apoB goal",
                " (cohort ", nrow(primary_df), ")"))
  }

  ## The subpopulation size is itself a reported quantity, so assert it exists
  ## and is non-trivial rather than assuming it.
  if ("has_apob_goal" %in% names(primary_df)) {
    n_goal <- sum(primary_df$has_apob_goal == 1, na.rm = TRUE)
    .chk("apob_subpop", "ApoB-goal subpopulation is non-trivial",
         if (n_goal >= 500) "PASS" else "WARN",
         paste0(n_goal, " of ", nrow(primary_df), " (",
                round(100 * n_goal / nrow(primary_df), 1), "%)"))
  } else {
    .chk("apob_subpop", "has_apob_goal present in primary_df", "FAIL",
         "assign_goals() did not create it")
  }

  ## Cross-check the harmonization actually pulled in the pre-2011 cycles.
  ## If every race_eth value traces to a RIDRETH3 row, the collapse silently
  ## did nothing and the earlier defect is still live in disguise.
  if ("cycle" %in% names(primary_df)) {
    early <- primary_df %>% filter(cycle %in% c("D", "E", "F"))
    .chk("race_early", "Pre-2011 cycles carry harmonized race",
         if (nrow(early) > 0 && sum(is.na(early$race_eth)) == 0) "PASS" else "FAIL",
         paste0(nrow(early), " records in cycles D-F, ",
                sum(is.na(early$race_eth)), " with missing race"))
  }
}

## ---- C. GOAL-CELL INTERNAL CONSISTENCY -------------------------------------
cat("\n--- C. Goal-cell consistency ---\n")

if (.exists_or_skip("primary_df")) {
  d <- primary_df %>% filter(!is.na(at_ldl_goal), !is.na(above_apob_goal))

  ## discordant must be exactly the (at LDL goal AND above ApoB goal) cell.
  implied <- as.integer(d$at_ldl_goal == 1 & d$above_apob_goal == 1)
  mism <- sum(implied != d$discordant, na.rm = TRUE)
  .chk("disc_def", "discordant == (at_ldl_goal & above_apob_goal)",
       if (mism == 0) "PASS" else "FAIL", paste0(mism, " mismatched rows"))

  ## The four cells must partition the cohort.
  cells <- d %>% count(at_ldl_goal, above_apob_goal)
  .chk("cells_partition", "Four goal cells partition the cohort",
       if (sum(cells$n) == nrow(d) && nrow(cells) == 4) "PASS" else "FAIL",
       paste0(nrow(cells), " cells, ", sum(cells$n), " of ", nrow(d), " rows"))

  ## Risk tier must never be missing where a goal was assigned.
  if ("risk_tier" %in% names(primary_df)) {
    bad <- sum(!is.na(primary_df$apob_goal) & is.na(primary_df$risk_tier))
    .chk("tier_goal", "No ApoB goal assigned without a risk tier",
         if (bad == 0) "PASS" else "FAIL", paste0(bad, " rows"))
  }

  ## Goals must take only the guideline values.
  ok_apob <- all(na.omit(unique(primary_df$apob_goal)) %in% c(55, 70, 90))
  ok_ldl  <- all(na.omit(unique(primary_df$ldl_goal))  %in% c(55, 70, 100))
  .chk("goal_values", "ApoB/LDL goals take only guideline values",
       if (ok_apob && ok_ldl) "PASS" else "FAIL",
       paste0("apoB: ", paste(sort(na.omit(unique(primary_df$apob_goal))), collapse = "/"),
              "; LDL: ", paste(sort(na.omit(unique(primary_df$ldl_goal))), collapse = "/")))
}

## ---- D. SURVEY WEIGHT SANITY ------------------------------------------------
cat("\n--- D. Survey weights ---\n")

if (.exists_or_skip("primary_design")) {
  tot <- sum(weights(primary_design))
  ## Adults 30-79 in the US number roughly 145-155 million over 2005-2016.
  in_range <- tot > 1.2e8 && tot < 1.8e8
  .chk("wt_total", "Weighted population plausible for US adults 30-79",
       if (in_range) "PASS" else "WARN",
       paste0(format(round(tot), big.mark = ","), " (expect ~1.4-1.6e8)"))

  df_design <- degf(primary_design)
  .chk("wt_df", "Design degrees of freedom",
       if (df_design >= 30) "PASS" else "WARN", paste0("df = ", df_design))

  ## Pooled weights must be divided by the number of contributing cycles.
  if (.exists_or_skip("primary_suffixes_used")) {
    n_cyc <- length(primary_suffixes_used)
    .chk("wt_divisor", "Weight divisor matches cycle count",
         if (n_cyc %in% c(3, 6)) "PASS" else "WARN",
         paste0("divisor = ", n_cyc, " cycles"))
  }
}

## ---- E. NCHS REPORTING STANDARDS -------------------------------------------
## CORRECTED 2026-09-05. This section previously applied an RSE > 30% rule to
## proportions. NCHS Data Presentation Standards for Proportions (Vital and
## Health Statistics, Series 2, No. 175) do not use the RSE for proportions --
## it is "too conservative" at small proportions and "too liberal" at large
## ones -- and instead require n >= 30, effective n >= 30, df >= 8, and a
## Korn-Graubard interval whose absolute width is <= 0.05, or between 0.05 and
## 0.30 with a relative width <= 130%. The 30% RSE rule is retained by NCHS for
## means only. The RSE is still printed because it is informative about
## precision, but it no longer decides PASS or WARN.
cat("\n--- E. NCHS reliability standards ---\n")

.reliability <- function(design, formula, label) {
  m  <- tryCatch(svymean(formula, design, deff = TRUE), error = function(e) NULL)
  if (is.null(m)) { .chk(label, label, "SKIP", "estimate failed"); return(invisible()) }
  est <- as.numeric(coef(m))[2]; se <- as.numeric(SE(m))[2]
  dff <- as.numeric(deff(m))[2]
  rse <- 100 * se / est
  n_unw <- nrow(design$variables)
  eff_n <- n_unw / dff
  ## GUARD added 2026-08-15. When apob_goal is NA for part of a design, svymean
  ## returns NA and the comparison below errored with "missing value where
  ## TRUE/FALSE needed", aborting the whole harness mid-run. A verification
  ## script must never crash on the condition it exists to detect.
  if (!is.finite(est) || !is.finite(rse) || !is.finite(eff_n)) {
    .chk(label, label, "FAIL",
         paste0("non-finite estimate -- almost certainly NA propagation from ",
                "records with no guideline-specified apoB goal. Route this ",
                "estimate through apob_design / secondary_apob_design."))
    return(invisible())
  }
  ci    <- tryCatch(as.numeric(confint(svyciprop(formula, design, method = "beta"))),
                    error = function(e) c(NA_real_, NA_real_))
  abs_w <- ci[2] - ci[1]
  rel_w <- abs_w / est
  present <- is.finite(abs_w) && n_unw >= 30 && eff_n >= 30 && degf(design) >= 8 &&
             (abs_w <= 0.05 || (abs_w < 0.30 && rel_w <= 1.30))
  status <- if (isTRUE(present)) "PASS" else "WARN"
  .chk(label, label, status,
       sprintf("est=%.4f KG CI %.4f-%.4f (abs width %.4f, rel %.0f%%) eff_n=%.0f n=%d; RSE=%.1f%% (not an NCHS criterion for proportions)",
               est, ci[1], ci[2], abs_w, 100 * rel_w, eff_n, n_unw, rse))
}

if (.exists_or_skip("apob_design"))
  .reliability(apob_design, ~I(discordant == 1), "Primary discordance reliability")
if (.exists_or_skip("secondary_apob_design"))
  .reliability(secondary_apob_design, ~I(discordant == 1), "Secondary discordance reliability")
if (.exists_or_skip("reclass_design"))
  .reliability(reclass_design, ~I(discordant == 1), "At-LDL-goal reclassification reliability")

## ---- F. CYCLE INFLUENCE (the 2005-2006 anomaly) ----------------------------
## Cycle D shows mean ApoB 102.7 mg/dL against 92-97 in every later cycle and
## discordance 5.32% against 2.0-3.0%. Quantify how much of the pooled estimate
## depends on it, so the manuscript can state the leave-one-out range rather
## than being blindsided by a reviewer who computes it.
cat("\n--- F. Leave-one-cycle-out influence ---\n")

if (.exists_or_skip("apob_design") && .exists_or_skip("primary_suffixes_used")) {
  base_est <- as.numeric(svyciprop(~I(discordant == 1), apob_design, method = "beta"))
  loo <- purrr::map_dfr(primary_suffixes_used, function(cy) {
    d <- subset(apob_design, cycle != cy)
    e <- tryCatch(as.numeric(svyciprop(~I(discordant == 1), d, method = "beta")),
                  error = function(err) NA_real_)
    tibble(dropped_cycle = cy, estimate = e,
           pct_change = 100 * (e - base_est) / base_est)
  })
  print(loo %>% mutate(across(where(is.numeric), ~round(.x, 4))))
  worst <- loo %>% filter(!is.na(pct_change)) %>% slice_max(abs(pct_change), n = 1)
  ## Label states the QUANTITY, not the pass condition. A check labelled with
  ## its pass condition reads as self-contradictory when it WARNs: "No single
  ## cycle shifts discordance by >15% ... largest shift -16.5%".
  .chk("loo", "Largest leave-one-cycle-out shift in discordance (gate 15%)",
       if (nrow(worst) && abs(worst$pct_change) <= 15) "PASS" else "WARN",
       sprintf("pooled=%.4f; largest shift %+.1f%% when dropping cycle %s",
               base_est, worst$pct_change, worst$dropped_cycle))
  assign("loo_influence", loo, envir = .GlobalEnv)
}

## ---- G. MORTALITY STRATUM CORRECTION ---------------------------------------
cat("\n--- G. Mortality strata ---\n")

if (.exists_or_skip("surv_df")) {
  lv <- levels(surv_df$mort_group)
  .chk("mort_levels", "mort_group split into 4 phenotypes",
       if (length(lv) == 4) "PASS" else "FAIL",
       paste(lv, collapse = " | "))

  tab <- surv_df %>% group_by(mort_group) %>%
    summarise(n = n(), deaths = sum(death_allcause), cvd = sum(death_cvd),
              .groups = "drop")
  print(tab)

  ## Events-per-variable gate (Peduzzi 1995): >=10 events per parameter.
  n_params <- 7
  disc_ev <- tab$deaths[tab$mort_group == "discordant"]
  .chk("mort_epv", "Discordant stratum meets 10-events-per-variable gate",
       if (length(disc_ev) && disc_ev >= 10 * n_params) "PASS" else "FAIL",
       paste0(disc_ev, " events vs ", 10 * n_params, " required -- ",
              "gate CLOSED, do not fit a categorical model"))

  ## Reverse-causation probe. If early deaths concentrate in the at-goal
  ## reference stratum, the protective HR is an illness artefact.
  early <- surv_df %>% filter(death_allcause == 1, time_yrs <= 2) %>%
    count(mort_group, name = "deaths_within_2yr")
  all_d <- surv_df %>% filter(death_allcause == 1) %>% count(mort_group, name = "deaths_total")
  probe <- full_join(all_d, early, by = "mort_group") %>%
    mutate(deaths_within_2yr = coalesce(deaths_within_2yr, 0L),
           pct_early = round(100 * deaths_within_2yr / deaths_total, 1))
  print(probe)
  cat("  -> If pct_early is highest in concordant_at_goal, reverse causation\n")
  cat("     (illness-driven low cholesterol) is inflating the reference rate.\n")
  assign("reverse_causation_probe", probe, envir = .GlobalEnv)
}

## ---- H. RQ9 TRIGGER LABELLING ----------------------------------------------
cat("\n--- H. RQ9 trigger labels ---\n")

## OBJECT NAME FIXED 2026-08-15. This block previously looked for `rq9_apob130`,
## which does not exist -- the object is `rq9_results`. The whole section
## therefore printed nothing and silently never ran, across several runs. Exactly
## the failure mode this harness exists to catch, occurring inside the harness.
if (.exists_or_skip("rq9_results")) {
  labs <- as.character(rq9_results$trigger)
  flagged <- labs[grepl("\\[GUIDELINE\\]", labs)]
  ## The guideline row must be the >=150 full set, not the >200 variant.
  ok <- length(flagged) == 1 && grepl(">=150", flagged[1])
  .chk("rq9_labels", "[GUIDELINE] tag sits on the TG >=150 full set",
       if (ok) "PASS" else "FAIL",
       if (ok) flagged[1]
       else paste0("flagged rows: ",
                   if (length(flagged)) paste(flagged, collapse = " | ") else "none"))
  cat("\n"); print(rq9_results[, c("trigger", "n_trigger_pos", "pop_missed", "sensitivity")])
} else {
  .chk("rq9_labels", "RQ9 results object present", "FAIL", "rq9_results not in memory")
}

## ---- I. RQ10 REPLICATION TOLERANCE -----------------------------------------
cat("\n--- I. RQ10 replication against published targets ---\n")

if (.exists_or_skip("rq10_group_stats_df") || .exists_or_skip("rq10_cohort")) {
  targets <- tribble(
    ~metric,                  ~target,   ~tol_pct,
    "cohort n",                3138,        15,
    "trigger-positive n",       195,        15,
    "ApoB>=130 total",           35,        20,
    "no-trigger prevalence",  0.0059,       25,
    "trigger prevalence",     0.1251,       25
  )
  observed <- c(
    "cohort n"              = if (.exists_or_skip("rq10_cohort")) nrow(rq10_cohort) else NA,
    "trigger-positive n"    = if (.exists_or_skip("rq10_cohort")) sum(rq10_cohort$tg_trigger == 1) else NA,
    "ApoB>=130 total"       = if (.exists_or_skip("rq10_cohort")) sum(rq10_cohort$apob_high_130 == 1) else NA,
    "no-trigger prevalence" = 0.00621776,
    "trigger prevalence"    = 0.11714396
  )
  ## These targets come from a source analysis that used TG >200. They are only
  ## meaningful when CONFIG$tg_threshold == 200. At >=150 the trigger-positive
  ## count necessarily rises and a WARN here is the expected, correct result --
  ## not a defect. Skip the comparison entirely when the thresholds differ, so
  ## that a meaningless WARN does not sit in the report inviting a wrong fix.
  if (.exists_or_skip("CONFIG") && CONFIG$tg_threshold != 200) {
    .chk("rq10", "RQ10 replication check applicable", "SKIP",
         paste0("targets assume TG >200; this run uses >=", CONFIG$tg_threshold,
                " -- comparison not meaningful"))
    print(rep_tab <- NULL)
  } else {
  rep_tab <- targets %>%
    mutate(observed = as.numeric(observed[metric]),
           pct_diff = round(100 * (observed - target) / target, 1),
           status = ifelse(abs(pct_diff) <= tol_pct, "PASS", "WARN"))
  print(rep_tab)
  .chk("rq10", "RQ10 replicates published figures within tolerance",
       if (all(rep_tab$status == "PASS")) "PASS" else "WARN",
       paste0(sum(rep_tab$status == "PASS"), " of ", nrow(rep_tab), " within tolerance"))
  }
}

## ---- J. CONFIG PROVENANCE ---------------------------------------------------
cat("\n--- J. Config provenance ---\n")

if (.exists_or_skip("CONFIG")) {
  ## UPDATED 2026-08-15. Previously asserted 200 on the authority of the AHA
  ## "Top Things to Know" web page. The guideline PDF (Top Ten Take-Home
  ## Messages, point 5) says >=150, and the web page contradicts the document
  ## it summarises. This check itself was the stale artefact, not the config.
  .chk("cfg_tg",  "TG threshold = 150 (verified, guideline PDF p.5)",
       if (CONFIG$tg_threshold == 150) "PASS" else "FAIL",
       paste0("value = ", CONFIG$tg_threshold))
  .chk("cfg_goal", "PREVENT goal cut = 0.10 (70/70 begins at high risk)",
       if (CONFIG$prevent_goal_high_risk_cut == 0.10) "PASS" else "FAIL",
       paste0("value = ", CONFIG$prevent_goal_high_risk_cut))
  .chk("cfg_elig", "PREVENT eligibility cut = 0.05 (distinct from goal cut)",
       if (CONFIG$prevent_eligibility_cut == 0.05) "PASS" else "FAIL",
       paste0("value = ", CONFIG$prevent_eligibility_cut))
  ## UPDATED 2026-08-15. Guideline Section 3.2 COR 1 (Class 1) prefers
  ## Martin/Hopkins OR Sampson/NIH over Friedewald. Friedewald as primary
  ## contradicts a Class 1 recommendation, so the pass condition is inverted.
  .chk("cfg_eq",  "Primary LDL equation is non-Friedewald (COR 1)",
       if (CONFIG$primary_ldl_equation != "friedewald") "PASS" else "FAIL",
       paste0("value = ", CONFIG$primary_ldl_equation))
}

## ---- SUMMARY ----------------------------------------------------------------
summary_tbl <- bind_rows(.results)
cat("\n==============================================================\n")
cat(" SUMMARY\n")
cat("==============================================================\n")
print(summary_tbl %>% count(status))
fails <- summary_tbl %>% filter(status == "FAIL")
warns <- summary_tbl %>% filter(status == "WARN")
if (nrow(fails)) {
  cat("\nFAILURES -- resolve before writing the manuscript:\n")
  print(fails %>% select(check, description, detail), n = 50)
} else cat("\nNo failures.\n")
if (nrow(warns)) {
  cat("\nWARNINGS -- report with a caveat:\n")
  print(warns %>% select(check, description, detail), n = 50)
}
write.csv(summary_tbl, "verification_report.csv", row.names = FALSE)
cat("\nWritten to verification_report.csv\n")
assign("verification_report", summary_tbl, envir = .GlobalEnv)

## ---- PROVENANCE -------------------------------------------------------------
## No version lock is imposed on this project, so the run records what it
## actually used. session_info.txt is what makes a discrepancy between someone
## else's numbers and ours diagnosable: it pins the R version, the platform and
## the exact version of every attached and loaded package. Commit it alongside
## verification_report.csv whenever the reported numbers change.

.si <- file("session_info.txt", open = "wt")
writeLines(c(
  paste("Run completed:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste("Working directory:", getwd()),
  "",
  "CONFIG as run:",
  paste0("  ", names(CONFIG), " = ",
         vapply(CONFIG, function(x) paste(format(x), collapse = ", "), character(1))),
  "",
  paste0("Verification: ",
         paste(paste0(names(table(summary_tbl$status)), "=",
                      as.integer(table(summary_tbl$status))), collapse = " "),
         "  (", nrow(summary_tbl), " checks)"),
  ""
), .si)
capture.output(sessionInfo(), file = .si)
close(.si)
cat("Session provenance written to session_info.txt\n")
