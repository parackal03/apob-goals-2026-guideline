## make_cascade_figure_dualgoal.R -- ALTERNATIVE main-text Figure (option B).
##
## Written October 2026 after the reviewer's comments, which make the "apparently
## controlled" estimate the one for adults at BOTH their LDL-C and non-HDL-C
## goals (6.6%), with the LDL-C-only figure (24.3%) alongside it.
##
## The top two bars are the current Figure unchanged. The 349 adults at their
## LDL-C goal are then shown as a bracketed pair of bars, two definitions of
## control for the same group, rather than as two further steps in the
## cascade: the 259 at both goals are a subset of the 349, but the 6.6% is not
## drawn from the 24.3%, and a chained bar would suggest that it was.
##
## Writes figure2_cascade_dualgoal.png/.tiff and
## figure2_cascade_dualgoal_values.csv. It does NOT touch figure2_cascade.png,
## so the current Figure is unchanged until the reviewers approve the switch
## (FIGURE_VARIANT in build_jama.py).
##
##   Rscript make_cascade_figure_dualgoal.R      (from the project folder)

suppressPackageStartupMessages(library(ggplot2))
.note <- function(f) message("  wrote ", f)

FILL  <- "#1c5cab"; TRACK  <- "#86b6ef"; FUNNEL <- "#dce9fc"
FILL2 <- "#5b88c4"; TRACK2 <- "#c4dbf6"            # lighter pair: the secondary estimate
INK   <- "#0b0b0b"; MUTED  <- "#52514e"; RULE <- "#8a8985"
.pct <- function(x) 100 * as.numeric(x)

## --- values, all sourced ------------------------------------------------------
.bounds <- read.csv("rev1_precedence_bounds.csv", stringsAsFactors = FALSE)
.union  <- .bounds[grepl("^Union", .bounds$variant), ]
.inter  <- .bounds[grepl("^Intersection", .bounds$variant), ]
.above  <- read.csv("rq1b_overall_apob_nonattainment.csv", stringsAsFactors = FALSE)
.recl   <- read.csv("rq1c_reclassification_at_ldl_goal.csv", stringsAsFactors = FALSE)
.part   <- read.csv("rq11_goal_pair_partition.csv", stringsAsFactors = FALSE)
.both   <- .part[.part$population == "all" & .part$cell == "both", ]
.ldlo   <- .part[.part$population == "all" & .part$cell == "ldl_only", ]
.t1     <- read.csv("table1_baseline_characteristics.csv", stringsAsFactors = FALSE,
                    check.names = FALSE)

n_pp   <- as.integer(gsub(",", "", .t1[.t1[[1]] == "Unweighted n", 2]))
n_goal <- as.integer(.union$n_unweighted)
n_at   <- as.integer(.recl$n_unweighted[1])
n_both <- as.integer(.both$n_unweighted)
stopifnot(n_both + as.integer(.ldlo$n_unweighted) == n_at)      # 259 + 90 = 349
pct_at <- 100 * n_at / n_goal

r1 <- c(.pct(.union$share),    .pct(.union$ci_low),   .pct(.union$ci_high))
r2 <- c(.pct(.above$estimate), .pct(.above$ci_low),   .pct(.above$ci_high))
rB <- c(.pct(.both$estimate),  .pct(.both$ci_low),    .pct(.both$ci_high))
rL <- c(.pct(.recl$estimate),  .pct(.recl$ci_low),    .pct(.recl$ci_high))

## --- drawing helpers ------------------------------------------------------------
BL <- 29; BR <- 67; BW <- BR - BL; H <- 12
xa <- function(p) BL + BW * p / 100

bar <- function(y, v, fill = FILL, track = TRACK) {
  cy <- y + H / 2
  whisk <- function(lw, col) list(
    annotate("segment", x = xa(v[2]), xend = xa(v[3]), y = cy, yend = cy, colour = col, linewidth = lw),
    annotate("segment", x = xa(v[2]), xend = xa(v[2]), y = cy - 2.1, yend = cy + 2.1, colour = col, linewidth = lw),
    annotate("segment", x = xa(v[3]), xend = xa(v[3]), y = cy - 2.1, yend = cy + 2.1, colour = col, linewidth = lw))
  list(annotate("rect", xmin = BL, xmax = BR,       ymin = y, ymax = y + H, fill = track),
       annotate("rect", xmin = BL, xmax = xa(v[1]), ymin = y, ymax = y + H, fill = fill),
       whisk(1.5, INK), whisk(0.8, "white"))
}
left <- function(y, title, n) { cy <- y + H / 2; list(
  annotate("text", x = BL - 2.2, y = cy + 2.4, label = title, hjust = 1, size = 3.0,
           fontface = "bold", colour = INK, lineheight = 1.05),
  annotate("text", x = BL - 2.2, y = cy - 4.2, label = n, hjust = 1, size = 2.7, colour = MUTED)) }
right <- function(y, v, sub, col = INK) { cy <- y + H / 2; list(
  annotate("text", x = BR + 2.2, y = cy + 3.4, label = sprintf("%.1f%%", v[1]), hjust = 0,
           size = 5.6, fontface = "bold", colour = col),
  annotate("text", x = BR + 2.2, y = cy - 1.4, label = sprintf("95%% CI, %.1f-%.1f", v[2], v[3]),
           hjust = 0, size = 2.6, colour = MUTED),
  annotate("text", x = BR + 2.2, y = cy - 6.6, label = sub, hjust = 0, size = 2.7, colour = MUTED)) }
fmt <- function(n) format(n, big.mark = ",")

## --- layout -----------------------------------------------------------------------
Y1 <- 113; Y2 <- 82; HDR <- 63; Y3 <- 36; Y4 <- 8

p <- ggplot() +
  ## funnels first, so the notes sit on top of them
  annotate("polygon", x = c(BL, xa(r1[1]), BR, BL), y = c(Y1, Y1, Y2 + H, Y2 + H), fill = FUNNEL) +
  annotate("polygon", x = c(BL, xa(pct_at), BR, BL), y = c(Y2, Y2, HDR, HDR), fill = FUNNEL) +
  bar(Y1, r1) + bar(Y2, r2) +
  left(Y1, "US adults 30-79 y,\nno established ASCVD", sprintf("n = %s", fmt(n_pp))) +
  left(Y2, "Adults in an apoB\nmeasurement group",   sprintf("n = %s", fmt(n_goal))) +
  right(Y1, r1, "in an apoB measurement group") +
  right(Y2, r2, "above their apoB value") +
  ## the pair: two definitions of control for the 349 at LDL-C goal
  annotate("segment", x = BL, xend = BR, y = HDR, yend = HDR, colour = RULE, linewidth = 0.4) +
  annotate("text", x = BL, y = HDR - 3.5, hjust = 0, size = 2.8, fontface = "bold", colour = INK,
           label = sprintf("The %s adults at their LDL-C goal, by definition of control", fmt(n_at))) +
  annotate("text", x = BL, y = HDR - 7.4, hjust = 0, size = 2.5, colour = MUTED, fontface = "italic",
           label = sprintf("The %s are a subset of the %s: two definitions of control, not successive steps",
                           fmt(n_both), fmt(n_at))) +
  bar(Y3, rB) + bar(Y4, rL, FILL2, TRACK2) +
  left(Y3, "LDL-C and non-HDL-C\ngoals both met", sprintf("n = %s (primary)", fmt(n_both))) +
  left(Y4, "LDL-C goal met,\nany non-HDL-C",      sprintf("n = %s", fmt(n_at))) +
  right(Y3, rB, "above apoB value") +
  right(Y4, rL, "above apoB value", col = "#3b3a37") +
  annotate("segment", x = 5, xend = 5, y = Y4, yend = Y3 + H, colour = RULE, linewidth = 0.4) +
  annotate("segment", x = 5, xend = 6.2, y = c(Y4, Y3 + H), yend = c(Y4, Y3 + H),
           colour = RULE, linewidth = 0.4) +
  ## notes under the top two bars
  annotate("text", x = BL, y = Y1 - 2.2, hjust = 0, vjust = 1, size = 2.5, colour = MUTED,
           fontface = "italic",
           label = sprintf("%.1f%%-%.1f%% across readings of Figure 1", .pct(.inter$share), r1[1])) +
  annotate("text", x = BL, y = Y2 - 2.2, hjust = 0, vjust = 1, size = 2.5, colour = MUTED,
           fontface = "italic",
           label = sprintf("%s (%.1f%%) of them reached their LDL-C goal", fmt(n_at), pct_at)) +
  scale_x_continuous(limits = c(0, 100), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 127), expand = c(0, 0)) +
  labs(caption = paste0("Survey-weighted estimates with 95% Korn-Graubard intervals, ",
                        "NHANES 2005-2016.\n",
                        "Each bar spans the group named at its left; shading is that ",
                        "percentage of the group. Counts (n) are unweighted.")) +
  theme_void(base_size = 11) +
  theme(plot.caption = element_text(hjust = 0, size = 7.4, colour = MUTED, margin = margin(t = 6)),
        plot.margin = margin(6, 6, 4, 6))

ggsave("figure2_cascade_dualgoal.png",  p, width = 7.6, height = 5.6, dpi = 400)
ggsave("figure2_cascade_dualgoal.tiff", p, width = 7.6, height = 5.6, dpi = 400, compression = "lzw")
.note("figure2_cascade_dualgoal.png"); .note("figure2_cascade_dualgoal.tiff")

## Every number printed on the figure, written out so check_all_numbers.py can
## trace it, including the unweighted 10.5% share in the note and the legend.
vals <- data.frame(
  quantity = c("in measurement group, share", "above apoB value",
               "both LDL-C and non-HDL-C goals met, above apoB value",
               "LDL-C goal met (any non-HDL-C), above apoB value",
               "unweighted share at LDL-C goal (349/3324)"),
  n        = c(n_pp, n_goal, n_both, n_at, n_at),
  pct      = round(c(r1[1], r2[1], rB[1], rL[1], pct_at), 1),
  ci_low   = round(c(r1[2], r2[2], rB[2], rL[2], NA), 1),
  ci_high  = round(c(r1[3], r2[3], rB[3], rL[3], NA), 1))
write.csv(vals, "figure2_cascade_dualgoal_values.csv", row.names = FALSE)
.note("figure2_cascade_dualgoal_values.csv")
print(vals)
