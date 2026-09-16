## =============================================================================
## precedence_diagnostics.R
##
## Two diagnostics on the precedence question. Run AFTER
## nhanes_apob_analysis.R and precedence_and_model_analyses.R, in the same session.
##
##   D1. Which Figure 1 row drives the precedence span, and by what mechanism.
##       Writes rev1_precedence_driver.csv.
##
##   D2. Cycle-level influence, and the stringency sensitivity.
##       Writes rev5_stringency_sensitivity.csv.
##
## HISTORY. This file previously also carried two corrections -- a malformed
## rev1_precedence_bounds.csv and a Firth model fitted on raw survey weights.
## Both are now fixed at source in precedence_and_model_analyses.R, so both were removed
## from here. Nothing in this file overwrites another script's output any
## more; if you skip it you lose only the two diagnostics named above.
## =============================================================================
suppressPackageStartupMessages({library(dplyr); library(survey); library(purrr)})

stopifnot(exists("primary_design"), exists("pdf_rows"),
          exists("g_union"), exists("g_inter"), exists("g_main"))

## ---- D1: what drives the precedence span -------------------------------------
## The run gave a 14.7-32.3% span across 120 orderings. This asks which row is
## responsible, because the answer determines whether the span is a real
## ambiguity in the guideline or an artefact of first-match semantics.
message("\n-- D1. Source of the precedence span --")

perm_csv <- read.csv("rev1_precedence_sensitivity.csv", stringsAsFactors = FALSE)
pos_of <- function(ord, nm) match(nm, strsplit(ord, " > ")[[1]])
d1 <- expand.grid(row = ROWS, position = 1:5, stringsAsFactors = FALSE) %>%
  rowwise() %>%
  mutate(mean_share = mean(perm_csv$share[
    vapply(perm_csv$ordering, pos_of, integer(1), nm = row) == position])) %>%
  ungroup() %>%
  tidyr::pivot_wider(names_from = position, values_from = mean_share,
                     names_prefix = "pos") %>%
  mutate(span_pp = 100 * (pmax(pos1,pos2,pos3,pos4,pos5) -
                          pmin(pos1,pos2,pos3,pos4,pos5)))
write.csv(d1, "rev1_precedence_driver.csv", row.names = FALSE)
print(as.data.frame(d1))

## How many people are absorbed by a row that states no apoB goal, when another
## row they also match does state one? This is the entire mechanism.
absorbed <- pdf_rows %>%
  mutate(
    any_row_states_goal  = !is.na(g_union),
    main_states_goal     = !is.na(g_main),
    absorbed_by_silent_row = any_row_states_goal & !main_states_goal
  )
message("  Records matching >=1 row that states an apoB goal: ",
        sum(absorbed$any_row_states_goal))
message("  Of these, denied a goal under the main ordering because a SILENT ",
        "row was reached first: ", sum(absorbed$absorbed_by_silent_row))

## Which silent row absorbed them?
sil <- absorbed %>% filter(absorbed_by_silent_row) %>%
  summarise(by_severe_hc = sum(in_severe_hc & !in_diabetes),
            by_hypertrig = sum(in_hypertrig & !in_diabetes & !in_severe_hc),
            both         = sum(in_severe_hc & in_hypertrig))
print(as.data.frame(sil))

## Is clinical_ascvd inert here? (primary_df excludes ASCVD by construction)
message("  Records in the clinical_ascvd row within primary prevention: ",
        sum(pdf_rows$in_clinical_ascvd),
        " -- if 0, that row is a no-op and the real permutation set is 4! = 24.")

## ---- D2: cycle G ------------------------------------------------------------
## Calibration fixed cycle D (-22.2% -> -6.0%) but cycle G now shows -17.3%.
## Check whether that is an apoB level effect or a composition effect.
message("\n-- D2. Cycle G influence --")
cyc <- svyby(~LBXAPB, ~cycle, primary_design, svymean, vartype = "ci")
print(cyc)
cyc_raw <- svyby(~LBXAPB_raw, ~cycle, primary_design, svymean, vartype = "ci")
print(cyc_raw)
cyc_tg <- svyby(~I(hypertrig == 1), ~cycle, primary_design, svymean, vartype = "ci")
print(cyc_tg)
cyc_n  <- primary_df %>% count(cycle, discordant) %>%
  tidyr::pivot_wider(names_from = discordant, values_from = n, names_prefix = "disc")
print(as.data.frame(cyc_n))
message("  Cycle G has the lowest calibrated mean apoB. If its discordant cell ",
        "count is in single figures, the -17.3% is small-cell noise rather ",
        "than a measurement problem, and should be reported as such.")

## ---- D3: does the STRINGENCY rule matter? -----------------------------------
## Membership is settled without any ordering (union). Where two matched rows
## state different goals we take the stricter. That is a single stated rule
## rather than an ordering, but it is still a choice, and a reviewer will ask
## how much it moves the answer. It cannot move the DENOMINATOR -- the same
## people have goals either way -- only the goal values, and therefore
## above-goal status and discordance.
message("\n-- D3. Stringency rule: strictest vs most lenient applicable goal --")

stopifnot("apob_goal_lenient" %in% names(primary_df))

d_str <- primary_design
## above_goal() comes from nhanes_apob_analysis.R and honours
## CONFIG$strict_goal_boundaries. Both arms must use it, or the stringency
## comparison mixes two boundary conventions.
d_str$variables$.above_strict  <- as.integer(
  above_goal(d_str$variables$LBXAPB, d_str$variables$apob_goal))
d_str$variables$.above_lenient <- as.integer(
  above_goal(d_str$variables$LBXAPB, d_str$variables$apob_goal_lenient))
d_str$variables$.disc_strict   <- as.integer(
  d_str$variables$at_ldl_goal == 1 & d_str$variables$.above_strict  == 1)
d_str$variables$.disc_lenient  <- as.integer(
  d_str$variables$at_ldl_goal == 1 & d_str$variables$.above_lenient == 1)

d_sub <- subset(d_str, has_apob_goal == 1)

grab <- function(design, var) {
  ci <- svyciprop(as.formula(paste0("~I(", var, " == 1)")), design, method = "logit")
  c(est = as.numeric(ci), lo = confint(ci)[1], hi = confint(ci)[2])
}
d_ldl <- subset(d_sub, at_ldl_goal == 1)

strin <- bind_rows(
  tibble(rule = "Strictest applicable (adopted)",
         above_goal = grab(d_sub, ".above_strict")["est"],
         above_lo   = grab(d_sub, ".above_strict")["lo"],
         above_hi   = grab(d_sub, ".above_strict")["hi"],
         discordance = grab(d_sub, ".disc_strict")["est"],
         reclassification = grab(d_ldl, ".above_strict")["est"]),
  tibble(rule = "Most lenient applicable",
         above_goal = grab(d_sub, ".above_lenient")["est"],
         above_lo   = grab(d_sub, ".above_lenient")["lo"],
         above_hi   = grab(d_sub, ".above_lenient")["hi"],
         discordance = grab(d_sub, ".disc_lenient")["est"],
         reclassification = grab(d_ldl, ".above_lenient")["est"])
)
write.csv(strin, "rev5_stringency_sensitivity.csv", row.names = FALSE)
print(as.data.frame(strin))

n_diff <- sum(primary_df$apob_goal != primary_df$apob_goal_lenient, na.rm = TRUE)
message("  Participants whose goal VALUE differs between the two rules: ",
        n_diff, " of ", sum(primary_df$has_apob_goal == 1),
        " with a goal (", round(100 * n_diff / sum(primary_df$has_apob_goal == 1), 1), "%).")
message("  Denominator is identical under both rules by construction -- ",
        "the union settles membership, stringency only sets the value.")

message("\n=== FIXES AND DIAGNOSTICS COMPLETE ===")
