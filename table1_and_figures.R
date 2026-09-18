## =============================================================================
## table1_and_figures.R
##
## Produces the two things the manuscript is missing:
##   TABLE 1  Baseline characteristics -- never computed in any previous run
##   FIGURE 1 Goal-assigned share under alternative readings of Figure 1
##   FIGURE 2 Participant flow (data for the diagram; draw in PowerPoint)
##   FIGURE S1 Apolipoprotein B by cycle, uncalibrated vs calibrated
##
## RUN AFTER nhanes_apob_analysis.R, in the same session.
## =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(survey); library(ggplot2); library(tidyr); library(purrr)
})

## ---- PREFLIGHT ---------------------------------------------------------------
## This script is NOT standalone. It reads objects that nhanes_apob_analysis.R
## leaves in the global environment, and two CSVs that precedence_and_model_analyses.R
## writes. Running it in a fresh session produces a cascade of "object not
## found" errors, which the checks below prevent.

.need <- c("primary_design", "secondary_design", "primary_df", "secondary_df")
.missing <- .need[!vapply(.need, exists, logical(1), envir = globalenv())]
.csv  <- c("rev1_precedence_bounds.csv", "rev1_precedence_sensitivity.csv")
.nocsv <- .csv[!file.exists(.csv)]

if (length(.missing) || length(.nocsv)) {
  stop(
    "\n\n=========================================================================\n",
    " table1_and_figures.R cannot run.\n\n",
    if (length(.missing))
      paste0(" Missing objects: ", paste(.missing, collapse = ", "),
             "\n   -> the main analysis has not been run in this session.\n") else "",
    if (length(.nocsv))
      paste0(" Missing files:   ", paste(.nocsv, collapse = ", "),
             "\n   -> either precedence_and_model_analyses.R has not been run, or the working\n",
             "      directory is wrong. Current: ", getwd(), "\n") else "",
    "\n DO THIS, in order, in ONE R session:\n\n",
    "   setwd(\"<the folder containing these scripts>\")\n",
    "   source(\"nhanes_apob_analysis.R\")\n",
    "   source(\"precedence_and_model_analyses.R\")\n",
    "   source(\"table1_and_figures.R\")\n\n",
    " Do not restart R or clear the workspace between them.\n",
    "=========================================================================\n\n",
    call. = FALSE
  )
}

.written <- character(0)
.note <- function(f) { .written <<- c(.written, f); message("  wrote ", f) }

message("\n=== TABLE 1 AND FIGURES ===")
message("Preflight OK. ", nrow(primary_df), " primary-prevention and ",
        nrow(secondary_df), " secondary-prevention records in memory.")

## -----------------------------------------------------------------------------
## TABLE 1. Baseline characteristics
##
## Four columns, because the comparison that matters is between adults the
## guideline assigns a goal and adults it does not. That contrast is the
## paper's thesis made concrete: the goal-assigned group should be visibly
## enriched for hypertriglyceridaemia, diabetes and elevated risk.
## -----------------------------------------------------------------------------

d_all  <- primary_design
d_goal <- subset(primary_design, has_apob_goal == 1)
d_none <- subset(primary_design, has_apob_goal == 0)
d_sec  <- secondary_design

## continuous: weighted mean (SE)
cont <- function(design, var) {
  f <- as.formula(paste0("~", var))
  m <- svymean(f, design, na.rm = TRUE)
  sprintf("%.1f (%.1f)", coef(m)[1], SE(m)[1])
}
## binary: weighted percent (95% CI)
##
## Structural-zero guard, added 2026-09-07. An adult with diabetes is assigned
## an apoB goal by definition, so the "Without an apoB goal" column contains no
## diabetic participants at all. svyciprop(method = "logit") fits a
## quasibinomial GLM, which cannot converge on a cell with no events, and R
## reported "glm.fit: algorithm did not converge" at the end of every run. The
## estimate is 0 with a zero-width interval by construction, so it is returned
## directly. This changes no output -- the cell reads 0.0 (0.0-0.0) either way,
## which build_jama.py renders as an em dash -- it only stops a benign warning
## from looking like a modelling failure.
binr <- function(design, expr) {
  f <- as.formula(paste0("~I(", expr, ")"))
  v <- eval(parse(text = expr), design$variables)
  if (sum(v %in% TRUE, na.rm = TRUE) == 0) return("0.0 (0.0-0.0)")
  p <- svyciprop(f, design, method = "logit", na.rm = TRUE)
  sprintf("%.1f (%.1f-%.1f)", 100*as.numeric(p), 100*confint(p)[1], 100*confint(p)[2])
}
## categorical level percent
catl <- function(design, var, lvl) {
  f <- as.formula(paste0("~I(", var, " == '", lvl, "')"))
  p <- svyciprop(f, design, method = "logit", na.rm = TRUE)
  sprintf("%.1f", 100*as.numeric(p))
}

designs <- list(`All primary prevention` = d_all,
                `With an apoB goal`      = d_goal,
                `Without an apoB goal`   = d_none,
                `Secondary prevention`   = d_sec)

row_cont <- function(label, var)
  c(Characteristic = label, vapply(designs, cont, character(1), var = var))
row_binr <- function(label, expr)
  c(Characteristic = label, vapply(designs, binr, character(1), expr = expr))
row_catl <- function(label, var, lvl)
  c(Characteristic = label, vapply(designs, catl, character(1), var = var, lvl = lvl))

n_row <- c(Characteristic = "Unweighted n",
           vapply(designs, function(d) format(nrow(d$variables), big.mark = ","), character(1)))

tab1 <- bind_rows(
  as.list(n_row),
  as.list(row_cont("Age, years, mean (SE)", "RIDAGEYR")),
  as.list(row_binr("Female, % (95% CI)", "RIAGENDR == 2")),
  as.list(row_catl("Non-Hispanic White, %", "race_eth", "Non-Hispanic White")),
  as.list(row_catl("Non-Hispanic Black, %", "race_eth", "Non-Hispanic Black")),
  as.list(row_catl("Mexican American, %",   "race_eth", "Mexican American")),
  as.list(row_catl("Other Hispanic, %",     "race_eth", "Other Hispanic")),
  as.list(row_catl("Other or multiracial, %","race_eth", "Other or multiracial")),
  as.list(row_cont("Body-mass index, kg/m2, mean (SE)", "BMXBMI")),
  as.list(row_binr("Current smoking, % (95% CI)", "smoker_flag == 1")),
  as.list(row_cont("Systolic blood pressure, mmHg, mean (SE)", "BPXSY1")),
  as.list(row_binr("Diabetes, % (95% CI)", "diabetes == 1")),
  as.list(row_binr("Chronic kidney disease, % (95% CI)", "ckd == 1")),
  as.list(row_cont("Total cholesterol, mg/dL, mean (SE)", "LBXTC")),
  as.list(row_cont("HDL-C, mg/dL, mean (SE)", "LBDHDD")),
  as.list(row_cont("LDL-C (Sampson), mg/dL, mean (SE)", "ldl_primary")),
  as.list(row_cont("Triglycerides, mg/dL, mean (SE)", "LBXTR")),
  as.list(row_cont("Apolipoprotein B, mg/dL, mean (SE)", "LBXAPB")),
  as.list(row_binr("Triglycerides >=150 mg/dL, % (95% CI)", "hypertrig == 1")),
  as.list(row_binr("Taking a statin, % (95% CI)", "on_statin == 1"))
)

## PREVENT risk only where calculable (primary prevention)
## PREVENT is stored as a proportion; it is reported as a percentage, so it is
## rescaled here rather than left for hand-conversion at the pasting stage.
pr <- vapply(designs, function(d) {
  v <- d$variables
  if (!"prevent_ascvd_10yr" %in% names(v) || all(is.na(v$prevent_ascvd_10yr))) return("-")
  m <- svymean(~prevent_ascvd_10yr, d, na.rm = TRUE)
  sprintf("%.1f (%.1f)", 100*coef(m)[1], 100*SE(m)[1])
}, character(1))
tab1 <- bind_rows(tab1, as.list(c(Characteristic = "PREVENT 10-year ASCVD risk, %, mean (SE)", pr)))

write.csv(tab1, "table1_baseline_characteristics.csv", row.names = FALSE)
.note("table1_baseline_characteristics.csv")
print(as.data.frame(tab1))

message("\n  NOTE: the 0.0% diabetes cell in the 'Without an apoB goal' column is ",
        "structural, not empty: an adult with diabetes is assigned a goal by ",
        "definition, so the cell is rendered as an em dash. The triglyceride cell ",
        "in that column is 0.2%, NOT zero -- the primary-prevention band stops at ",
        "499 mg/dL, so an adult above that can match the row and still be assigned ",
        "no goal. Do not describe the two cells as if they were the same thing.")

## -----------------------------------------------------------------------------
## FIGURE 1. Goal-assigned share under alternative readings of Figure 1
##
## This is the highest-value figure in the paper: it makes the range-anchoring
## visual and answers the reviewer's first question before it is asked.
## -----------------------------------------------------------------------------

fig1_dat <- read.csv("rev1_precedence_bounds.csv", stringsAsFactors = FALSE)
perm     <- read.csv("rev1_precedence_sensitivity.csv", stringsAsFactors = FALSE)

## x is kept NUMERIC, not a factor. The shaded band spanning the ordering range
## is drawn with numeric xmin/xmax, which a discrete scale rejects; the labels
## are restored through scale_x_continuous(breaks = , labels = ) instead.
fig1_lab <- c("Intersection\n(strictest reading)",
              "Sequential first-match\n(one ordering)",
              "Union, most stringent goal\n(primary analysis)")
fig1_dat$variant <- fig1_lab
fig1_dat$xpos    <- seq_len(nrow(fig1_dat))

p1 <- ggplot(fig1_dat, aes(x = xpos, y = 100*share)) +
  annotate("rect", xmin = 0.4, xmax = 3.6,
           ymin = 100*min(perm$share), ymax = 100*max(perm$share),
           fill = "grey85", alpha = 0.55) +
  annotate("text", x = 0.55, y = 100*max(perm$share) - 1.2,
           label = "range across all 120 sequential orderings",
           hjust = 0, size = 3.1, colour = "grey30") +
  geom_errorbar(aes(ymin = 100*ci_low, ymax = 100*ci_high), width = 0.12, linewidth = 0.6) +
  geom_point(size = 3.4) +
  ## label is anchored above the upper CI cap, not above the point, otherwise it
  ## is drawn through the whisker
  geom_text(aes(y = 100*ci_high + 0.7, label = sprintf("%.1f%%", 100*share)),
            vjust = 0, size = 3.6) +
  scale_x_continuous(breaks = fig1_dat$xpos, labels = fig1_lab,
                     limits = c(0.4, 3.6), expand = c(0, 0)) +
  scale_y_continuous(limits = c(10, 38), breaks = seq(10, 35, 5),
                     labels = function(x) paste0(x, "%")) +
  labs(x = NULL,
       y = "US adults aged 30-79 without established ASCVD\nassigned an apolipoprotein B goal",
       caption = "Points are survey-weighted estimates with 95% Korn-Graubard intervals.") +
  theme_classic(base_size = 11) +
  theme(axis.text.x = element_text(size = 9.5),
        plot.caption = element_text(hjust = 0, size = 8.5, colour = "grey35"))

ggsave("figure1_interpretation_range.png", p1, width = 6.5, height = 4.6, dpi = 400)
ggsave("figure1_interpretation_range.tiff", p1, width = 6.5, height = 4.6, dpi = 400,
       compression = "lzw")
.note("figure1_interpretation_range.png"); .note("figure1_interpretation_range.tiff")

## -----------------------------------------------------------------------------
## FIGURE 2 (manuscript). Cascade: three conditional proportions on shrinking
## denominators.
##
## The paper's whole argument is that the denominator changes at each step, so
## the figure makes the denominator the visible thing: each bar spans the group
## named at its left, and the funnel shows that group being carved out of the
## one above. Every value is read from the verified CSVs, not typed in.
##
## Colours are one hue, light -> dark (ordinal). Encoding is lightness, not hue,
## so it survives greyscale printing and colour-vision deficiency.
## -----------------------------------------------------------------------------

FILL <- "#1c5cab"; TRACK <- "#86b6ef"; FUNNEL <- "#dce9fc"
INK  <- "#0b0b0b"; MUTED <- "#52514e"

.pct <- function(x) 100 * as.numeric(x)

## --- values, all sourced ------------------------------------------------------
.bounds <- read.csv("rev1_precedence_bounds.csv", stringsAsFactors = FALSE)
.union  <- .bounds[grepl("^Union", .bounds$variant), ]
.inter  <- .bounds[grepl("^Intersection", .bounds$variant), ]
.above  <- read.csv("rq1b_overall_apob_nonattainment.csv", stringsAsFactors = FALSE)
.recl   <- read.csv("rq1c_reclassification_at_ldl_goal.csv", stringsAsFactors = FALSE)

n_pp        <- nrow(primary_df)
n_goal      <- sum(primary_df$has_apob_goal == 1, na.rm = TRUE)
n_atgoal    <- as.integer(.recl$n_unweighted[1])
pct_atgoal  <- 100 * n_atgoal / n_goal

cascade <- data.frame(
  row   = 1:3,
  pct   = c(.pct(.union$share),   .pct(.above$estimate), .pct(.recl$estimate)),
  lo    = c(.pct(.union$ci_low),  .pct(.above$ci_low),   .pct(.recl$ci_low)),
  hi    = c(.pct(.union$ci_high), .pct(.above$ci_high),  .pct(.recl$ci_high)),
  title = c("US adults 30-79 y,\nno established ASCVD",
            "Adults assigned\nan apoB goal",
            "Of those, adults at\ntheir LDL-C goal"),
  nlab  = c(sprintf("n = %s", format(n_pp,     big.mark = ",")),
            sprintf("n = %s", format(n_goal,   big.mark = ",")),
            sprintf("n = %s", format(n_atgoal, big.mark = ","))),
  sub   = c("assigned an apoB goal", "above their apoB goal",
            "above apoB goal with LDL-C at goal"),
  stringsAsFactors = FALSE
)

BAR_L <- 29; BAR_R <- 67; BAR_W <- BAR_R - BAR_L; BAR_H <- 12
ROWY  <- c(70, 39, 8)
cascade$y  <- ROWY
cascade$cy <- cascade$y + BAR_H / 2
xat <- function(p) BAR_L + BAR_W * p / 100

## funnel polygons: the shaded slice of one bar becomes the whole of the next
w1 <- BAR_W * cascade$pct[1] / 100
funnels <- rbind(
  data.frame(g = 1, x = c(BAR_L, BAR_L + w1, BAR_R, BAR_L),
             y = c(ROWY[1], ROWY[1], ROWY[2] + BAR_H, ROWY[2] + BAR_H)),
  data.frame(g = 2, x = c(BAR_L, xat(pct_atgoal), BAR_R, BAR_L),
             y = c(ROWY[2], ROWY[2], ROWY[3] + BAR_H, ROWY[3] + BAR_H))
)

p2 <- ggplot() +
  geom_polygon(data = funnels, aes(x, y, group = g), fill = FUNNEL) +
  geom_rect(data = cascade, aes(xmin = BAR_L, xmax = BAR_R, ymin = y, ymax = y + BAR_H),
            fill = TRACK) +
  geom_rect(data = cascade, aes(xmin = BAR_L, xmax = xat(pct), ymin = y, ymax = y + BAR_H),
            fill = FILL) +
  ## CI whisker, drawn white over black so it reads on fill and on track alike
  geom_segment(data = cascade, aes(x = xat(lo), xend = xat(hi), y = cy, yend = cy),
               colour = "#0b0b0b", linewidth = 1.5) +
  geom_segment(data = cascade, aes(x = xat(lo), xend = xat(lo), y = cy - 2.1, yend = cy + 2.1),
               colour = "#0b0b0b", linewidth = 1.5) +
  geom_segment(data = cascade, aes(x = xat(hi), xend = xat(hi), y = cy - 2.1, yend = cy + 2.1),
               colour = "#0b0b0b", linewidth = 1.5) +
  geom_segment(data = cascade, aes(x = xat(lo), xend = xat(hi), y = cy, yend = cy),
               colour = "white", linewidth = 0.8) +
  geom_segment(data = cascade, aes(x = xat(lo), xend = xat(lo), y = cy - 2.1, yend = cy + 2.1),
               colour = "white", linewidth = 0.8) +
  geom_segment(data = cascade, aes(x = xat(hi), xend = xat(hi), y = cy - 2.1, yend = cy + 2.1),
               colour = "white", linewidth = 0.8) +
  ## left column
  geom_text(data = cascade, aes(x = BAR_L - 2.2, y = cy + 2.4, label = title),
            hjust = 1, vjust = 0.5, size = 3.0, fontface = "bold", colour = INK, lineheight = 1.05) +
  geom_text(data = cascade, aes(x = BAR_L - 2.2, y = cy - 4.2, label = nlab),
            hjust = 1, vjust = 0.5, size = 2.7, colour = MUTED) +
  ## right column
  geom_text(data = cascade, aes(x = BAR_R + 2.2, y = cy + 3.4, label = sprintf("%.1f%%", pct)),
            hjust = 0, vjust = 0.5, size = 5.6, fontface = "bold", colour = INK) +
  geom_text(data = cascade, aes(x = BAR_R + 2.2, y = cy - 1.4,
                                label = sprintf("95%% CI, %.1f-%.1f", lo, hi)),
            hjust = 0, vjust = 0.5, size = 2.6, colour = MUTED) +
  geom_text(data = cascade, aes(x = BAR_R + 2.2, y = cy - 6.6, label = sub),
            hjust = 0, vjust = 0.5, size = 2.7, colour = MUTED) +
  ## notes between the bars
  annotate("text", x = BAR_L, y = ROWY[1] - 2.2, hjust = 0, vjust = 1, size = 2.5,
           colour = MUTED, fontface = "italic",
           label = sprintf("%.1f%%-%.1f%% across readings of Figure 1",
                           .pct(.inter$share), .pct(.union$share))) +
  annotate("text", x = BAR_L, y = ROWY[2] - 2.2, hjust = 0, vjust = 1, size = 2.5,
           colour = MUTED, fontface = "italic",
           label = sprintf("%s (%.1f%%) of them reached their LDL-C goal",
                           format(n_atgoal, big.mark = ","), pct_atgoal)) +
  scale_x_continuous(limits = c(0, 100), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 84), expand = c(0, 0)) +
  labs(caption = paste0("Survey-weighted estimates with 95% Korn-Graubard intervals, ",
                        "NHANES 2005-2016.\n",
                        "Each bar spans the group named at its left; shading is that ",
                        "percentage of the group.")) +
  theme_void(base_size = 11) +
  theme(plot.caption = element_text(hjust = 0, size = 7.4, colour = MUTED,
                                    margin = margin(t = 6)),
        plot.margin = margin(6, 6, 4, 6))

ggsave("figure2_cascade.png", p2, width = 7.6, height = 3.9, dpi = 400)
ggsave("figure2_cascade.tiff", p2, width = 7.6, height = 3.9, dpi = 400, compression = "lzw")
.note("figure2_cascade.png"); .note("figure2_cascade.tiff")

message("\n  Figure 2 cascade values, for checking against RESULTS_VERIFIED.md:")
print(cascade[, c("pct", "lo", "hi", "nlab")])

## -----------------------------------------------------------------------------
## FIGURE 2. Participant flow -- numbers for the diagram
##
## Draw the box-and-arrow diagram in PowerPoint or Illustrator; this prints the
## counts so nothing is transcribed by hand.
## -----------------------------------------------------------------------------

flow <- tibble::tribble(
  ~step, ~n,
  "NHANES 2005-2016, fasting subsample, aged 30-79", nrow(primary_df) + nrow(secondary_df),
  "Primary prevention (no established ASCVD)",       nrow(primary_df),
  "  - assigned an apoB goal",                       sum(primary_df$has_apob_goal == 1),
  "  - no apoB goal specified",                      sum(primary_df$has_apob_goal == 0),
  "  - with a goal, at their LDL-C goal",            sum(primary_df$has_apob_goal == 1 &
                                                        primary_df$at_ldl_goal == 1, na.rm = TRUE),
  "  - with a goal, at LDL-C goal, above apoB goal", sum(primary_df$discordant == 1, na.rm = TRUE),
  "Secondary prevention (established ASCVD)",        nrow(secondary_df),
  "  - assigned an apoB goal",                       sum(secondary_df$has_apob_goal == 1)
)
write.csv(flow, "figure2_participant_flow_counts.csv", row.names = FALSE)
.note("figure2_participant_flow_counts.csv")
print(as.data.frame(flow))

## -----------------------------------------------------------------------------
## FIGURE S1. Apolipoprotein B by cycle, uncalibrated vs calibrated
##
## Supplementary. Makes the assay-calibration argument visible rather than
## asking the reader to accept two coefficients on trust.
## -----------------------------------------------------------------------------

if ("LBXAPB_raw" %in% names(primary_df)) {
  cyc_raw <- svyby(~LBXAPB_raw, ~cycle, primary_design, svymean, vartype = "ci", na.rm = TRUE)
  cyc_cal <- svyby(~LBXAPB,     ~cycle, primary_design, svymean, vartype = "ci", na.rm = TRUE)
  names(cyc_raw) <- c("cycle","mean","lo","hi"); cyc_raw$scale <- "Uncalibrated"
  names(cyc_cal) <- c("cycle","mean","lo","hi"); cyc_cal$scale <- "Calibrated"
  cy <- bind_rows(cyc_raw, cyc_cal)
  cy$years <- factor(cy$cycle, levels = c("D","E","F","G","H","I"),
                     labels = c("2005-06","2007-08","2009-10","2011-12","2013-14","2015-16"))
  ## the two series coincide exactly in cycles E-H; ordering the factor this way
  ## draws the calibrated series last so it is the one visible where they overlap
  cy$scale <- factor(cy$scale, levels = c("Uncalibrated", "Calibrated"))
  cy <- cy[order(cy$scale), ]

  ## Round-6 audit: the eFigure 2 note quotes the endpoint cycle means and
  ## the range across cycles, and those numbers existed only in the console
  ## log. One of them ("spread falls from 9.24 to 3.46") was stale and
  ## nothing could have caught it, because there was nothing to check it
  ## against. Written out so the supplement can be built from a file and
  ## guarded like every other quoted value.
  cyc_out <- cy[, c("scale", "cycle", "years", "mean", "lo", "hi")]
  cyc_out$mean <- as.numeric(cyc_out$mean)
  write.csv(cyc_out, "rq3_apob_by_cycle_calibration.csv", row.names = FALSE)
  .rng <- tapply(cyc_out$mean, cyc_out$scale, function(v) max(v) - min(v))
  write.csv(data.frame(scale = names(.rng),
                       range_mg_dl = as.numeric(.rng),
                       min_mg_dl = as.numeric(tapply(cyc_out$mean, cyc_out$scale, min)),
                       max_mg_dl = as.numeric(tapply(cyc_out$mean, cyc_out$scale, max))),
            "rq3_apob_cycle_range.csv", row.names = FALSE)
  message("apoB cycle means written; range across cycles: ",
          paste(sprintf("%s %.2f", names(.rng), .rng), collapse = " | "))

  pS1 <- ggplot(cy, aes(years, mean, group = scale, colour = scale, shape = scale)) +
    geom_line(aes(linetype = scale), linewidth = 0.6) +
    geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.1, linewidth = 0.5) +
    geom_point(size = 2.6) +
    scale_colour_manual(values = c(Uncalibrated = "grey55", Calibrated = "black")) +
    scale_linetype_manual(values = c(Uncalibrated = "22", Calibrated = "solid")) +
    labs(linetype = NULL) +
    labs(x = NULL, y = "Mean apolipoprotein B, mg/dL", colour = NULL, shape = NULL,
         caption = "Survey-weighted means with 95% confidence intervals. Calibration places all cycles on the BN ProSpec scale.") +
    theme_classic(base_size = 11) +
    theme(legend.position = "top",
          plot.caption = element_text(hjust = 0, size = 8.5, colour = "grey35"))

  ggsave("figureS1_apob_by_cycle.png", pS1, width = 6.5, height = 4.2, dpi = 400)
  .note("figureS1_apob_by_cycle.png")
} else {
  message("  LBXAPB_raw not found - rerun the main script so the uncalibrated ",
          "values are retained, then rerun this file for Figure S1.")
}

message("\n=== COMPLETE ===")
message(length(.written), " file(s) written:")
for (f in .written) message("  - ", f)
if (length(.written) < 4) warning("Fewer files than expected. Scroll up for the section that failed.", call. = FALSE)
message("Table 1 needs a caption and a footnote defining the abbreviations. ",
        "Figure 1 goes in the main text; Figure S1 in the supplement.")
