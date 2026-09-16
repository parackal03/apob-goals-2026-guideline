## =============================================================================
## compare_blind_coding.R
##
## Compares the second coder's independent reading of Figure 1 against ours.
##
## WRITTEN BEFORE THE CODER'S FILE ARRIVED. Every quantity reported here was
## chosen without sight of the result. Do not add a statistic after reading the
## output; if something further is needed, say in the paper that it was post hoc.
##
## RUN AFTER nhanes_apob_analysis.R, in the same session.
##
## Expects the returned file to carry, at minimum:
##   SEQN, apob_goal_blind, ldl_goal_blind
## Goals in mg/dL; blank / NA means "their reading assigns none".
## =============================================================================

## Looks in the packet folder first, then the project root.
RETURNED <- local({
  cand <- c(file.path("BLIND_CODING_PACKET", "blind_coding_returned.csv"),
            "blind_coding_returned.csv")
  hit <- cand[file.exists(cand)]
  if (!length(hit)) cand[1] else hit[1]
})

.need <- c("primary_df", "secondary_df", "primary_design")
.missing <- .need[!vapply(.need, exists, logical(1), envir = globalenv())]
if (length(.missing)) {
  stop("\n\n compare_blind_coding.R cannot run. Missing: ",
       paste(.missing, collapse = ", "),
       "\n Run source(\"nhanes_apob_analysis.R\") first, in this session.\n\n",
       call. = FALSE)
}
if (!file.exists(RETURNED)) {
  stop("\n\n Cannot find '", RETURNED, "'. Rename the coder's file to that ",
       "and put it in this folder.\n\n", call. = FALSE)
}

suppressPackageStartupMessages({ library(dplyr); library(survey) })

.out <- character(0)
OUT <- function(x, f) { write.csv(x, f, row.names = FALSE)
                        .out <<- c(.out, f); message("  wrote ", f) }

## ---------------------------------------------------------------- ingest ---
blind <- read.csv(RETURNED, stringsAsFactors = FALSE)
names(blind)[names(blind) == "seqn"] <- "SEQN"

req <- c("SEQN", "apob_goal_blind", "ldl_goal_blind")
miss <- setdiff(req, names(blind))
if (length(miss)) stop("Returned file is missing: ", paste(miss, collapse = ", "),
                       call. = FALSE)

ours <- bind_rows(
  primary_df   %>% mutate(cohort = "primary prevention"),
  secondary_df %>% mutate(cohort = "secondary prevention")
) %>%
  select(SEQN, cohort, goal_row, apob_goal, ldl_goal, ldl_primary, LBXAPB)

## ------------------------------------------------------ integrity checks ---
message("\n=== INTEGRITY OF THE RETURNED FILE ===")
dup <- sum(duplicated(blind$SEQN))
only_theirs <- setdiff(blind$SEQN, ours$SEQN)
only_ours   <- setdiff(ours$SEQN, blind$SEQN)
message("  rows returned          : ", nrow(blind))
message("  duplicate SEQN         : ", dup)
message("  SEQN only in their file: ", length(only_theirs))
message("  SEQN missing from it   : ", length(only_ours))

## Values outside the guideline's set are REPORTED, never silently corrected --
## a different reading may legitimately produce a different value.
odd_a <- sort(setdiff(unique(na.omit(blind$apob_goal_blind)), c(55, 70, 90)))
odd_l <- sort(setdiff(unique(na.omit(blind$ldl_goal_blind)),  c(55, 70, 100)))
if (length(odd_a)) message("  NOTE apoB values we did not use: ", paste(odd_a, collapse = ", "))
if (length(odd_l)) message("  NOTE LDL-C values we did not use: ", paste(odd_l, collapse = ", "))
if (dup || length(only_theirs) || length(only_ours))
  warning("The returned file does not align one-to-one with ours. Resolve this ",
          "before reading any agreement statistic below.", call. = FALSE)

cmp <- ours %>% inner_join(blind[, req], by = "SEQN")
message("  matched for comparison : ", nrow(cmp))

## ------------------------------------------------------------ statistics ---
## Cohen's kappa, computed directly so no extra package is needed.
kappa_of <- function(a, b) {
  lv <- union(sort(unique(a)), sort(unique(b)))
  tb <- table(factor(a, lv), factor(b, lv))
  n  <- sum(tb); po <- sum(diag(tb)) / n
  pe <- sum(rowSums(tb) * colSums(tb)) / n^2
  k  <- (po - pe) / (1 - pe)
  se <- sqrt(po * (1 - po) / (n * (1 - pe)^2))     # Fleiss approximation
  list(table = tb, n = n, po = po, pe = pe, kappa = k,
       lo = k - 1.96 * se, hi = k + 1.96 * se)
}
## Linear-weighted kappa, for the ordered goal VALUES.
wkappa_of <- function(a, b) {
  lv <- sort(union(unique(a), unique(b))); k <- length(lv)
  if (k < 2) return(NULL)
  tb <- table(factor(a, lv), factor(b, lv)); n <- sum(tb)
  w  <- outer(seq_len(k), seq_len(k), function(i, j) 1 - abs(i - j) / (k - 1))
  p  <- tb / n; pr <- rowSums(p); pc <- colSums(p)
  po <- sum(w * p); pe <- sum(w * outer(pr, pc))
  list(table = tb, n = n, kappa_w = (po - pe) / (1 - pe))
}

report <- list()
for (coh in c("primary prevention", "secondary prevention", "ALL")) {
  d <- if (coh == "ALL") cmp else filter(cmp, cohort == coh)

  for (target in c("apob", "ldl")) {
    ourg   <- if (target == "apob") d$apob_goal       else d$ldl_goal
    theirg <- if (target == "apob") d$apob_goal_blind else d$ldl_goal_blind

    ## (A) does a goal apply at all
    a <- kappa_of(as.integer(!is.na(ourg)), as.integer(!is.na(theirg)))
    ## (B) among those where BOTH assign one, is the value the same
    both <- !is.na(ourg) & !is.na(theirg)
    val_agree <- if (any(both)) mean(ourg[both] == theirg[both]) else NA_real_
    wk <- if (sum(both) > 1) wkappa_of(ourg[both], theirg[both]) else NULL

    report[[length(report) + 1]] <- tibble(
      cohort = coh, goal = target, n = a$n,
      both_assign_n       = sum(both),
      agree_applies_pct   = 100 * a$po,
      kappa_applies       = a$kappa,
      kappa_lo            = a$lo,
      kappa_hi            = a$hi,
      expected_agree_pct  = 100 * a$pe,
      ours_assign_pct     = 100 * mean(!is.na(ourg)),
      theirs_assign_pct   = 100 * mean(!is.na(theirg)),
      agree_value_pct     = 100 * val_agree,
      kappa_value_weighted = if (is.null(wk)) NA_real_ else wk$kappa_w
    )
  }
}
agreement <- bind_rows(report)
OUT(agreement, "blind_agreement_summary.csv")
print(as.data.frame(agreement))

message("\n  Read kappa alongside ours_assign_pct and theirs_assign_pct. When both ",
        "\n  margins are extreme, kappa can be low despite high raw agreement; ",
        "\n  report the raw percentage and the marginals with it, not kappa alone.")

## -------------------------------------------- where the disagreements are ---
disagree <- cmp %>%
  mutate(
    apob_applies_differs = is.na(apob_goal) != is.na(apob_goal_blind),
    apob_value_differs   = !is.na(apob_goal) & !is.na(apob_goal_blind) &
                            apob_goal != apob_goal_blind,
    ldl_applies_differs  = is.na(ldl_goal)  != is.na(ldl_goal_blind),
    ldl_value_differs    = !is.na(ldl_goal) & !is.na(ldl_goal_blind) &
                            ldl_goal != ldl_goal_blind
  ) %>%
  filter(apob_applies_differs | apob_value_differs |
         ldl_applies_differs  | ldl_value_differs)

by_row <- disagree %>%
  count(cohort, goal_row,
        apob_applies_differs, apob_value_differs,
        ldl_applies_differs,  ldl_value_differs, name = "n") %>%
  arrange(desc(n))
OUT(by_row, "blind_disagreements_by_row.csv")
OUT(disagree %>% select(SEQN, cohort, goal_row, apob_goal, apob_goal_blind,
                        ldl_goal, ldl_goal_blind),
    "blind_disagreements_participants.csv")
message("\n  ", nrow(disagree), " participants differ on at least one goal.")
print(as.data.frame(head(by_row, 15)))

## ------------------------- what the headline would be under their reading ---
## Mirrors the primary definitions exactly: at_ldl_goal uses <=, above_apob
## uses >, discordant is the conjunction.
pb <- primary_df %>%
  select(SEQN) %>%
  left_join(blind[, req], by = "SEQN")
stopifnot(identical(pb$SEQN, primary_df$SEQN))

pd <- update(primary_design,
             has_apob_goal_blind = as.integer(!is.na(pb$apob_goal_blind)),
             ## Same boundary convention as the primary analysis, or the
             ## agreement statistic measures the convention rather than the
             ## coding. at_goal()/above_goal() come from the main script.
             discordant_blind = as.integer(
               at_goal(primary_df$ldl_primary, pb$ldl_goal_blind) &
               above_goal(primary_df$LBXAPB,   pb$apob_goal_blind)))

share_blind <- svyciprop(~I(has_apob_goal_blind == 1), pd, method = "beta")
disc_blind  <- svyciprop(~I(discordant_blind == 1),
                         subset(pd, has_apob_goal_blind == 1), method = "beta")

## OUR side must be COMPUTED, not typed. It was typed -- c(33.4, 2.38) -- and
## that survived the strict-goal-boundary correction, because a literal cannot
## go stale in any way R can see. The blind side was recomputed and moved from
## 2.25% to 2.41%; ours stayed at a 2.38% that no longer existed anywhere else
## in the analysis, and the supplement printed the pair as though both had been
## produced under the same rule. Both sides now come from the same designs, in
## the same session, through the same comparators.
share_ours <- svyciprop(~I(has_apob_goal == 1), primary_design, method = "beta")
disc_ours  <- svyciprop(~I(discordant == 1),
                        subset(primary_design, has_apob_goal == 1),
                        method = "beta")

impact <- tibble(
  quantity = c("Goal-assigned share, primary prevention",
               "Discordance among goal-assigned"),
  ours_pct = 100 * c(as.numeric(share_ours), as.numeric(disc_ours)),
  ours_lo  = 100 * c(confint(share_ours)[1], confint(disc_ours)[1]),
  ours_hi  = 100 * c(confint(share_ours)[2], confint(disc_ours)[2]),
  theirs_pct = 100 * c(as.numeric(share_blind), as.numeric(disc_blind)),
  theirs_lo  = 100 * c(confint(share_blind)[1], confint(disc_blind)[1]),
  theirs_hi  = 100 * c(confint(share_blind)[2], confint(disc_blind)[2])
)

## Each estimate must sit inside the other reading's interval. That is the
## claim eTable 2 makes; if it ever stops being true, the build should not be
## the place it is discovered.
stopifnot(all(impact$ours_pct >= impact$theirs_lo &
              impact$ours_pct <= impact$theirs_hi))
OUT(impact, "blind_headline_impact.csv")
print(as.data.frame(impact))

message("\n=== BLIND CODING COMPARISON COMPLETE ===")
message(length(.out), " file(s) written:")
for (f in .out) message("  - ", f)
message("\n  Report the result as it came out. A disagreement is a finding about ",
        "\n  the legibility of Figure 1 and belongs in the paper.")
