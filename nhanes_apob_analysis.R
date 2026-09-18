## =============================================================================
## NHANES Apolipoprotein B (ApoB) Discordance & Reclassification Analysis
## Companion script to: NHANES_ApoB_Study_Protocol.docx
## Modeled on: Anderson TS, Wilson LM, Sussman JB. Implications of the 2026
##   Dyslipidemia Guideline for Primary Prevention Statin Therapy. JAMA. 2026.
##   doi:10.1001/jama.2026.11246
##
## -----------------------------------------------------------------------------
## HOW TO REPLICATE (read first)
## -----------------------------------------------------------------------------
##   1. Set the working directory to the folder containing this script.
##   2. source("nhanes_apob_analysis.R")     # full analysis, writes CSVs + PNGs
##   3. source("verify_run.R")               # 29-check verification harness
##
## Data are downloaded automatically on first run: NHANES .xpt files into
## nhanes_raw/, and the NCHS Linked Mortality Files (2019 public-use release)
## into the same folder. Nothing needs to be placed by hand. Parsed cohorts are
## cached as .rds in nhanes_cache/; DELETE THAT FOLDER to force a clean re-parse
## from source, otherwise the cache is reused and edits to the loading code will
## appear to have no effect.
##
## Required packages: dplyr, tidyr, purrr, tibble, broom, survey, srvyr,
## haven, preventr, ggplot2, readr. No version lock is in place -- results were
## produced under the versions current in August 2026.
##
## Everything that governs the analysis lives in the CONFIG block below. Change
## it there, never inline.
##
## -----------------------------------------------------------------------------
## DATA AVAILABILITY
## -----------------------------------------------------------------------------
## NHANES measured ApoB (LBXAPB) in SIX consecutive cycles -- 2005-2006,
## 2007-2008, 2009-2010, 2011-2012, 2013-2014, 2015-2016 -- and has NOT measured
## it since. There is no ApoB in the 2017-March 2020 pre-pandemic file
## (P_TRIGLY) or in August 2021-2023. Confirmed against CDC variable lists.
## The most recent nationally representative US ApoB data is therefore 2015-2016,
## and results describe that period rather than a current snapshot.
##
## The PRIMARY analysis pools ALL SIX cycles (CONFIG$primary_cycle_set = "all6").
## An earlier version restricted to 2011-2016 on the stated grounds that those
## cycles were "assay-consistent". That justification was checked against the
## APOB_G documentation and found to be unsupported -- it records "no changes to
## the lab site, lab instrument or methodology" -- so the restriction was
## discarding half the sample for no reason. Set primary_cycle_set = "recent3"
## to reproduce the old three-cycle behaviour.
##
## -----------------------------------------------------------------------------
## LDL-C EQUATION: Sampson/NIH is PRIMARY (Class 1 recommendation)
## -----------------------------------------------------------------------------
## CONFIG$primary_ldl_equation = "sampson". Friedewald (NHANES's supplied
## LBDLDL) is retained as the comparison arm and reported in Section 11b.
##
## Guideline Section 3.2, COR 1, LOE B-NR, verbatim from the PDF: "In adults and
## children who have undergone a standard lipid profile, use of either the
## Martin/Hopkins equation or the Sampson/National Institutes of Health (NIH)
## equation is preferred over calculation by the Friedewald equation to estimate
## LDL-C." Class 1 is the strongest recommendation tier, and it names
## Martin/Hopkins and Sampson as equivalent acceptable alternatives.
##
## HISTORY, because this reversed twice. Sampson was made primary, then reverted
## to Friedewald on the argument that Friedewald is what laboratories actually
## report and therefore what determines real-world goal classification. That
## argument was wrong: the guideline being operationalized explicitly deprecates
## Friedewald at Class 1. Anchoring on it contradicts the document the paper
## claims to test. The equation effect remains a reportable finding in its own
## right -- Section 11b quantifies it -- but Sampson is primary.
##
## Martin/Hopkins is not implemented, and this is NOT a limitation. COR 1 names
## Martin/Hopkins OR Sampson as acceptable alternatives to Friedewald;
## implementing either satisfies it. Do not describe the absence of a second
## arm as a shortcoming in the manuscript, and do not assert anything about the
## availability of the Martin/Hopkins factor table -- it is published (Martin
## et al. JAMA 2013;310:2061-2068) and has been applied to these same NHANES
## cycles by others. The correct and complete statement in the methods is that
## Sampson is one of the two methods the recommendation names.
##
## The Sampson implementation below is independently confirmed correct: CDC's
## published LBDLDLN formula in the TRIGLY_J codebook is character-for-character
## identical to sampson_ldl(). Cite that in the methods -- it means this analysis
## uses the same equation NHANES itself supplies from 2017-2018 onward.
##
## Any figure in older comments from the Friedewald-primary period is obsolete.
## Trust RESULTS_VERIFIED.md, not stray numbers in comments.
##
## -----------------------------------------------------------------------------
## GUIDELINE THRESHOLDS -- verified against the guideline PDF itself
## -----------------------------------------------------------------------------
## Verified 2026-08-15 against Blumenthal RS, et al. 2026 ACC/AHA guideline on
## the management of dyslipidemia. JACC 2026;87(19):2624-2757.
##   * TG trigger is >=150 mg/dL (Top Ten Take-Home Messages, point 5).
##     WARNING: the AHA "Top Things to Know" WEB PAGE states >200 and is
##     inconsistent with the guideline it summarises. An earlier pass of this
##     script "corrected" 150 to 200 on the authority of that web page. That was
##     an error in the wrong direction. Use the PDF.
##   * The recommendation body (Section 3.3, COR 2a) says only "elevated TG"
##     with no number, and scopes to adults ON lipid-lowering therapy with
##     ASCVD, CKM syndrome, type 2 diabetes, and/or elevated TG. Adults NOT on
##     LLT fall under the weaker COR 2b. The cohorts in RQ9/RQ10 are untreated
##     and therefore sit under 2b, not 2a -- state this in any write-up.
##   * The 70/70 goal pair begins at 10% PREVENT 10-year risk, NOT 5%.
## The guideline specifies two DISTINCT PREVENT thresholds that an earlier
## version conflated: eligibility (therapy "may be considered" >=3%, "should be
## considered" >=5%) versus goal assignment (stricter goals from >=10%). These
## are held separately in CONFIG as prevent_eligibility_cut and
## prevent_goal_high_risk_cut.
## =============================================================================

## =============================================================================
## CONFIG -- every tunable analytic choice lives here, in one place.
## Change a value, re-source the script, and every downstream result updates.
## Nothing else in the script hard-codes these.
## =============================================================================
CONFIG <- list(

  ## APOB CROSS-INSTRUMENT CALIBRATION -- added 2026-08-16 (reviewer comment 4).
  ##
  ## NHANES measured apoB on THREE different instrument/laboratory combinations
  ## across the six pooled cycles, and publishes Deming bridging equations for
  ## both transitions. This analysis had been pooling all six cycles without
  ## applying them, which is the cause of the 2005-2006 anomaly that survived
  ## three previous rounds of correction to the goal-assignment logic.
  ##
  ##   2005-2006 (D)   Dade-Behring BN100,     nephelometric, Johns Hopkins
  ##   2007-2014 (E-H) Siemens BN ProSpec,     nephelometric, U. Minnesota
  ##   2015-2016 (I)   Roche Cobas 6000,       TURBIDIMETRIC, U. Minnesota
  ##
  ## CDC (APOB_E and APOB_I documentation, analytic notes):
  ##   Y(2005-06, BN100) = 1.084 * X(2007-08, ProSpec)     n=168, r=0.984
  ##   Y(ProSpec) = 0.9988 * X(Cobas 6000) - 2.880         n=105, r=0.988
  ## CDC states these "should be used" when trending or combining apoB.
  ##
  ## We adopt the ProSpec scale (four of six cycles) as the common scale.
  ## Set FALSE to reproduce the uncalibrated run for the sensitivity comparison.
  apply_apob_calibration = TRUE,

  ## HYPERTENSION DEFINITION -- added 2026-08-27 (script review, finding 2).
  ##
  ## TRUE  = hypertensive if averaged systolic BP >=140 mmHg OR on
  ##         antihypertensive therapy. This is the clinically correct reading
  ##         of the guideline's unqualified "hypertension", and it is what
  ##         PREVENT already assumes -- on_antihtn is passed to the equations
  ##         as bp_tx regardless of this setting.
  ## FALSE = the superseded definition, first systolic reading >=140 only,
  ##         which counted treated-and-controlled adults as normotensive.
  ##
  ## htn_flag feeds ascvd_rf_count (which sets the 70-vs-90 goal VALUE in the
  ## diabetes and hypertriglyceridaemia rows) and ascvd_highrisk_cond_count
  ## (which sets very-high-risk ASCVD status). It does not affect whether an
  ## adult is assigned a goal at all, so the goal-assigned DENOMINATOR is
  ## unchanged by this switch; goal-attainment estimates are not.
  htn_counts_treatment = TRUE,

  ## CLINICAL ASCVD WITH CKD -- added 2026-09-02 (reviewer comment, Dani r2).
  ##
  ## Figure 1 names chronic kidney disease exactly once. It is NOT a patient-
  ## population row. It appears as a bullet in the Clinical ASCVD row, in the
  ## LDL-C <55 / non-HDL-C <85 column, printed as:
  ##
  ##      * At very high risk
  ##          - apoB goal: <55 mg/dL
  ##      * With CKD
  ##
  ## So clinical ASCVD plus CKD is a second route into that column. The earlier
  ## code had no such route and gave these adults LDL-C <70.
  ##
  ## SCOPE, stated precisely because it is narrow:
  ##   - Primary prevention is untouched. ascvd_hx == 1 excludes those adults
  ##     from that cohort before goals are assigned, so every headline estimate
  ##     (goal-assigned share, attainment, discordance) is unchanged by design,
  ##     not by luck. verify_run.R asserts this.
  ##   - ApoB is untouched. The "<55" apoB bullet is nested under "At very high
  ##     risk", not under "With CKD", so CKD carries the column's LDL-C and
  ##     non-HDL-C goals only. g_ascvd_apob is deliberately NOT modified.
  ##     Reading it the other way would require inferring that CKD confers
  ##     very-high-risk status, which the figure does not say.
  ##   - Only secondary-prevention LDL-C goal attainment (eTable 7) moves.
  ##
  ## FALSE reproduces the superseded mapping for the sensitivity comparison.
  ##
  ## ---------------------------------------------------------------------------
  ## SUPERSEDED IN PART, 2026-09-10. READ THIS BEFORE THE PARAGRAPH ABOVE.
  ##
  ## The "ApoB is untouched" bullet above is no longer what this script does.
  ## It was written on 2 September from Figure 1 alone. On 10 September the
  ## guideline PDF was re-read in full (round 4, comment 9) and Figure 11's
  ## entry box states verbatim: "Clinical ASCVD at very high risk (includes
  ## ASCVD with CKD)", with an optional apoB goal of <55 mg/dL. That IS the
  ## document saying what the bullet above says it does not say.
  ##
  ## CONFIG$ascvd_ckd_very_high (below) now carries that reading and is TRUE by
  ## default, so ASCVD + CKD does reach the apoB goal, via ascvd_very_high
  ## rather than via this switch. It assigns an apoB goal to 54 adults who
  ## would otherwise have none.
  ##
  ## The paragraph above is kept, not rewritten, because it is the record of
  ## what was believed when this switch was added and of why the earlier
  ## estimates differ. Its reasoning stands on Figure 1; it was outvoted by
  ## Figure 11. This switch still does only what it says -- it governs the
  ## LDL-C goal and nothing else.
  ## ---------------------------------------------------------------------------
  ascvd_ckd_ldl55 = TRUE,

  ## STRICT GOAL BOUNDARIES -- added 2026-09-10 (round 4, comment 5).
  ##
  ## The guideline states goals as strict inequalities: "apoB goal: <70 mg/dL".
  ## An adult at exactly 70 mg/dL is therefore NOT at that goal. This script
  ## previously used <= for LDL-C attainment and > for apoB non-attainment,
  ## which counts a value sitting exactly on the threshold as at goal on both
  ## sides. NHANES reports apoB and the lipid panel to the nearest whole number
  ## in most cycles, so boundary ties are not hypothetical.
  ##
  ## TRUE  = literal reading.  at goal: value <  goal.  above goal: value >= goal.
  ## FALSE = the superseded inclusive reading, retained for the comparison.
  ##
  ## Direction of the effect is known in advance and is the same on both sides:
  ## fewer adults at their LDL-C goal (the 349 denominator falls) and more above
  ## their apoB goal (the 85.6% rises). The size is an empirical question.
  ## HYPERTENSION SYSTOLIC THRESHOLD -- added 2026-09-10 (round 4, comment 13).
  ##
  ## 140 mmHg is what the 2026 guideline's own very-high-risk criteria box uses
  ## ("Hypertension" as a high-risk condition, with no number given, read here
  ## against the historical >=140 definition). Contemporary US practice defines
  ## hypertension at >=130/80, which is less restrictive and would classify more
  ## adults as carrying an ASCVD risk factor.
  ##
  ## htn_flag feeds ascvd_rf_count, which sets the 70-versus-90 goal VALUE in
  ## the diabetes and hypertriglyceridaemia rows. It does not decide whether an
  ## adult is assigned a goal at all, so the goal-assigned denominator is
  ## invariant to this setting; attainment is not. pathway_and_sensitivity_analyses.R quantifies
  ## the difference rather than leaving it as an assertion.
  htn_sbp_threshold = 140,

  strict_goal_boundaries = TRUE,

  ## DOMAIN ANALYSIS -- added 2026-09-10 (round 4, comment 6).
  ##
  ## TRUE  = build the survey design on every fasting-subsample record with a
  ##         positive weight, then restrict to the analytic cohort with
  ##         subset(), which is the NCHS-recommended domain approach: records
  ##         outside the domain still inform the variance through their strata
  ##         and PSUs.
  ## FALSE = the superseded behaviour, svydesign() on the already-filtered
  ##         data frame.
  ##
  ## Point estimates are expected to be identical. Standard errors may differ
  ## slightly where filtering emptied a PSU. Section 0b prints both so the
  ## comparison is reportable rather than asserted.
  use_domain_design = TRUE,

  ## ASCVD + CKD AS VERY HIGH RISK -- added 2026-09-10 (round 4, comment 9).
  ##
  ## Adjudicated against the guideline PDF on 10 September 2026, and the
  ## document is not self-consistent. Both readings are recorded here.
  ##
  ##   Figure 1, Clinical ASCVD row, LDL-C <55 column, prints two bullets at
  ##   the same indent level:
  ##        * At very high risk
  ##            - apoB goal: <55 mg/dL
  ##        * With CKD
  ##   Read alone, the apoB value is nested under the first bullet only, and
  ##   CKD is a parallel route to the column's LDL-C and non-HDL-C goals. That
  ##   is what this script did.
  ##
  ##   Figure 11, "Secondary ASCVD Prevention for Adults at Very High Risk",
  ##   entry box, verbatim: "Clinical ASCVD at very high risk (includes ASCVD
  ##   with CKD)". The box that follows states the LDL-C <55 / non-HDL-C <85
  ##   goals AND "Optional apoB goal <55 mg/dL". Under this reading the two
  ##   Figure 1 bullets are not parallel: the second is a member of the first,
  ##   so CKD carries the apoB goal too.
  ##
  ##   Figure 10, "Criteria for Defining At Very High Risk", lists the
  ##   qualifying arithmetic -- >=2 major events, or 1 major event plus >=2 of
  ##   (age >=65, bypass or PCI, current smoker, diabetes, heart failure,
  ##   hypertension, LDL-C >=100 despite maximally tolerated statin plus
  ##   ezetimibe). CKD appears nowhere in it.
  ##
  ## Figure 11's parenthetical is the only explicit statement in the document
  ## about how CKD relates to very-high-risk status, and it is unambiguous, so
  ## TRUE is adopted as primary. Figure 10 is read as enumerating the
  ## events-and-conditions route, not as an exhaustive definition.
  ##
  ## SCOPE: secondary prevention only. Primary prevention excludes ascvd_hx == 1
  ## before goals are assigned, so no primary-prevention estimate can move.
  ## verify_run.R asserts this.
  ascvd_ckd_very_high = TRUE,

  ## AGE AS AN ASCVD RISK FACTOR: one threshold for everyone, or the
  ## sex-specific pair the guideline's own footnote gives.
  ##
  ## Several rows of Figure 1 turn on "at least 1 ASCVD risk factor" but the
  ## figure never defines the term. The only enumeration in the document sits
  ## in a footnote to the icosapent-ethyl recommendation in diabetes, taken
  ## from the REDUCE-IT entry criteria, and it reads "men >=55 years or women
  ## >=65 years" -- not a single cut at 55.
  ##
  ## Whether that footnote governs Figure 1 is genuinely unclear: it is
  ## attached to one recommendation and is not offered as a general
  ## definition. FALSE keeps the single 55 cut used throughout, which is the
  ## primary analysis; TRUE applies the footnote's sex-specific pair.
  ## pathway_and_sensitivity_analyses.R (R4f) runs both and reports the difference, because the
  ## choice moves who has a risk factor and therefore, through the
  ## hypertriglyceridaemia row, who is assigned an apoB goal at all.
  ##
  ## Changing this changes the cohort's goals, so it must be set before
  ## assign_goals() runs, not afterwards.
  ascvd_rf_age_sex_specific = FALSE,

  # PRIMARY LDL-C EQUATION: "sampson" (primary) or "friedewald" (comparison).
  #
  #   sampson    = Sampson/NIH equation 2, recomputed here. THE PRIMARY
  #                ANALYSIS USES THIS.
  #   friedewald = NHANES's supplied LBDLDL. Retained as the comparison arm,
  #                reported in Section 11b, because it is what US laboratories
  #                reported during 2005-2016.
  #
  # Why Sampson is primary. Guideline Section 3.2, COR 1 (LOE B-NR), verbatim:
  # "In adults and children who have undergone a standard lipid profile, use of
  # either the Martin/Hopkins equation or the Sampson/National Institutes of
  # Health (NIH) equation is preferred over calculation by the Friedewald
  # equation to estimate LDL-C." That is a Class 1 recommendation naming
  # Martin/Hopkins and Sampson as equivalent acceptable alternatives, so
  # anchoring the primary analysis on Friedewald would contradict the document
  # this study operationalizes. Implementing either named method satisfies
  # COR 1; the absence of a Martin/Hopkins arm is not a limitation.
  #
  # The implementation is verified: CDC's published LBDLDLN formula in the
  # TRIGLY_J codebook is character-for-character identical to sampson_ldl()
  # below, so this analysis uses the same equation NHANES itself supplies from
  # 2017-2018 onward.
  primary_ldl_equation = "sampson",

  # Triglyceride threshold defining hypertriglyceridemia AND the apoB testing
  # trigger. THE VALUE IS 150 mg/dL, taken from the guideline PDF itself:
  # Top Ten Take-Home Messages, point 5, verbatim: "...particularly in those
  # with elevated triglycerides (TG) (>=150 mg/dL), diabetes, or low achieved
  # LDL-C (<70 mg/dL)."
  #
  # DO NOT change this to 200 without quoting guideline text. The AHA "Top
  # Things to Know" web page summarising this guideline states >200 mg/dL and
  # contradicts the document it summarises. An earlier pass of this script
  # "corrected" 150 to 200 on that page's authority; the PDF is authoritative
  # and that change was reverted.
  #
  # Note that the recommendation body itself (Section 3.3, COR 2a) says only
  # "elevated TG" with NO numeric threshold, and lists ASCVD, CKM syndrome and
  # type 2 diabetes -- not low achieved LDL-C. The >=150 figure comes from the
  # Take-Home summary. Both framings are reported, and RQ9 evaluates both the
  # 150 and 200 cuts, so its conclusions do not depend on this setting.
  tg_threshold = 150,

  # Markedly-elevated ApoB threshold for the RQ9 testing-yield analysis.
  # 130 mg/dL is the 2018 ACC/AHA risk-enhancing-factor threshold, stated to
  # correspond to LDL-C >=160. NOT a treatment goal. Retained at 130 for
  # comparability with the published Circulation 2025 analysis.
  apob_high_threshold = 130,

  # Analytic age range. 30-79 matches Anderson et al. (JAMA 2026), aiding
  # comparability with the concurrent NHANES guideline analyses, and matches
  # the range over which the guideline recommends PREVENT-based risk
  # estimation.
  age_min = 30,
  age_max = 79,

  # CYCLE SET FOR THE PRIMARY AND SECONDARY-PREVENTION COHORTS.
  #   "all6"    = 2005-2016, all six ApoB-measuring cycles (DEFAULT)
  #   "recent3" = 2011-2016 only
  #
  # CHANGED 2026-08-02 from "recent3" to "all6". The original restriction was
  # justified in this script as pooling "the three most recent, assay-
  # consistent cycles" -- but that assay-consistency claim was never verified.
  # The CDC documentation for APOB_G (2011-2012) states explicitly: "There no
  # changes to the lab site, lab instrument or methodology." There is therefore
  # no documented methods basis for the 2011 boundary, and the restriction was
  # discarding half the available sample for no reason.
  #
  # UNRESOLVED CAVEAT: 2005-2006 remains empirically anomalous -- mean ApoB
  # 102.7 mg/dL versus 92-97 in every later cycle, and discordance roughly
  # double. Whether that reflects a real population shift or something about
  # that cycle is not established; the 2007-2008 documentation has not been
  # checked. Set to "recent3" to reproduce the earlier restricted analysis, or
  # exclude cycle "D" specifically to test the 2005-2006 anomaly.
  primary_cycle_set = "all6",

  # TWO SEPARATE PREVENT THRESHOLDS -- they do different jobs and the guideline
  # sets them at different values. Conflating them (as this script previously
  # did, using 0.05 for both) mis-assigns treatment goals AND under-counts
  # statin eligibility.
  #
  # GOAL ASSIGNMENT. Verified from the AHA "Top Things to Know" (13 Mar 2026)
  # and the guideline-at-a-glance: risk categories are borderline 3% to <5%
  # and intermediate 5% to <10%, and "for those at borderline or intermediate
  # risk, the LDL-C goal should be less than 100 mg/dL, and for those at high
  # risk, the goal is less than 70 mg/dL." So the stricter 70/70 goal begins at
  # 10%, NOT 5%. Everyone below 10% takes the 100/90 pair.
  prevent_goal_high_risk_cut = 0.10,
  #
  # TREATMENT ELIGIBILITY. Separate question, lower bar: LLT "can be
  # considered" at 3% to <5% and "should be considered" at 5% to <10% after a
  # clinician-patient discussion. Used only by RQ9 to define who the guideline
  # does not recommend for therapy. Set to 0.05 ("should be considered") as
  # the conservative reading; 0.03 is the permissive alternative and is worth
  # running as a sensitivity check, since it materially changes the RQ9
  # denominator.
  prevent_eligibility_cut = 0.05,

  # Caching. Reading and merging the .xpt files across six cycles is the
  # slowest step in the script. When TRUE, the merged raw cohorts are written
  # to disk on first run and read back thereafter. Config changes do NOT
  # invalidate this cache, because the cached objects are pre-derivation --
  # they contain only raw merged NHANES variables. Delete the cache directory
  # to force a clean reload.
  use_cache = TRUE,
  cache_dir = "nhanes_cache"
)

stopifnot(CONFIG$primary_ldl_equation %in% c("friedewald", "sampson"))
dir.create(CONFIG$cache_dir, showWarnings = FALSE)

## ---- 0. Setup --------------------------------------------------------------

## Packages that must be ATTACHED, because the code below calls them without a
## namespace prefix. tibble and scales are used only as tibble:: and scales::,
## so they need to be INSTALLED but not attached; they are listed separately so
## the installation check covers them.
required_pkgs <- c("haven", "dplyr", "tidyr", "purrr", "survey", "srvyr",
                   "ggplot2", "broom", "stringr", "readr", "survival")
namespaced_pkgs <- c("tibble", "scales")

## Nothing is installed without consent. Sourcing a script should not silently
## write to a user's package library, so a missing package stops the run with
## the exact command to fix it. Set options(apob_auto_install = TRUE) before
## sourcing to opt in to automatic installation instead.
new_pkgs <- setdiff(c(required_pkgs, namespaced_pkgs),
                    rownames(installed.packages()))
if (length(new_pkgs)) {
  if (isTRUE(getOption("apob_auto_install"))) {
    install.packages(new_pkgs, repos = "https://cloud.r-project.org")
  } else {
    stop(
      "\n\n",
      "=========================================================================\n",
      " Missing R packages: ", paste(new_pkgs, collapse = ", "), "\n\n",
      " Install them, then re-source this script:\n\n",
      "   install.packages(c(",
      paste0('"', new_pkgs, '"', collapse = ", "), "))\n\n",
      " Or opt in to automatic installation for this session:\n\n",
      "   options(apob_auto_install = TRUE)\n",
      "=========================================================================\n\n",
      call. = FALSE
    )
  }
}
invisible(lapply(required_pkgs, library, character.only = TRUE))

## NHANES's default Taylor-series variance estimator requires >=2 PSUs per
## stratum. That holds for the full primary/trend cohorts (confirmed -- RQ0
## through RQ5 all ran without incident), but breaks for small subgroups: the
## secondary-prevention extension (Section 12b) subsets to just ~500 known-
## ASCVD adults, small enough that at least one stratum ends up represented
## by only a single PSU once every OTHER respondent in that stratum happens
## to lack ASCVD history. Without this option, svyciprop()/svymean()/
## svytotal() throw "Stratum (x) has only one PSU at stage 1" and abort.
## "adjust" is the standard NCHS/NHANES analytic-guidance workaround for
## exactly this situation: it centers the lonely PSU's contribution at the
## sample grand mean instead of erroring, which is a documented, accepted
## (mildly conservative) approximation for public-use survey data. This is a
## global session option, so setting it once here is sufficient for every
## survey call below -- and since no stratum in the full primary_design/
## trend_design/surv_design ever hit this condition, it has zero effect on
## any of the already-validated RQ0-RQ5 results; it only activates for the
## small secondary-prevention subgroup where it's actually needed.
options(survey.lonely.psu = "adjust")

## Round-4 reanalysis switches, echoed so any saved console log states which
## reading produced it. All three change reported numbers.
message("\nRound-4 settings: strict goal boundaries = ",
        isTRUE(CONFIG$strict_goal_boundaries),
        " | domain design = ", isTRUE(CONFIG$use_domain_design),
        " | ASCVD+CKD very high risk = ", isTRUE(CONFIG$ascvd_ckd_very_high),
        " | sex-specific age RF = ", isTRUE(CONFIG$ascvd_rf_age_sex_specific))

## preventr implements the AHA's actual published PREVENT equations (Khan et
## al, Circulation 2023;148:1982-2004 and Circulation 2024;149:430-449) --
## used below (Section 5b) to replace the simplified diabetes/CKD risk-tier
## proxy with the guideline's real 10-year PREVENT-ASCVD risk score. Not
## reliably on the default CRAN mirror; install from its r-universe build,
## which also mirrors CRAN, so this one repos= call covers both cases.
if (!requireNamespace("preventr", quietly = TRUE)) {
  if (isTRUE(getOption("apob_auto_install"))) {
    install.packages("preventr",
      repos = c("https://martingmayer.r-universe.dev", "https://cloud.r-project.org"))
  } else {
    stop(
      "\n\n",
      "=========================================================================\n",
      " Package 'preventr' is required and is not installed. It is not on the\n",
      " default CRAN mirror; install it from its r-universe build:\n\n",
      "   install.packages(\"preventr\",\n",
      "     repos = c(\"https://martingmayer.r-universe.dev\",\n",
      "               \"https://cloud.r-project.org\"))\n\n",
      " Or opt in to automatic installation for this session:\n\n",
      "   options(apob_auto_install = TRUE)\n",
      "=========================================================================\n\n",
      call. = FALSE
    )
  }
}
library(preventr)

data_dir <- "nhanes_raw"
dir.create(data_dir, showWarnings = FALSE)

## ---- 1. Cycle configuration -------------------------------------------------
## ApoB file naming: 2005-2006 ApoB is bundled inside TRIGLY_D; 2007-2016 ApoB
## is a standalone APOB_<suffix> file. WTSAF2YR (fasting subsample weight) is
## present in both the TRIGLY_<suffix> and APOB_<suffix> files for a cycle and
## should be identical -- we pull it from TRIGLY_<suffix> throughout for
## consistency.

cycles <- tibble::tribble(
  ~cycle,        ~suffix, ~apob_has_own_file,
  "2005-2006",   "D",     FALSE,   # ApoB lives inside TRIGLY_D (LBXAPB)
  "2007-2008",   "E",     TRUE,
  "2009-2010",   "F",     TRUE,
  "2011-2012",   "G",     TRUE,
  "2013-2014",   "H",     TRUE,
  "2015-2016",   "I",     TRUE
)

primary_suffixes <- c("G", "H", "I")          # 2011-2016: assay-consistent, most recent
trend_suffixes    <- cycles$suffix             # all 6, for secular trend (2005-2016)

cdc_base <- "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public"

## Map suffix -> 4-digit year folder CDC uses in URLs
year_folder <- c(D = "2005", E = "2007", F = "2009", G = "2011", H = "2013", I = "2015")

## ---- 2. Download helper -----------------------------------------------------

get_xpt <- function(component, suffix) {
  # component e.g. "DEMO", "BMX", "DIQ", "GHB", "BIOPRO", "ALB_CR", "MCQ",
  # "TRIGLY", "APOB", "HDL", "TCHOL", "RXQ_RX"
  fname <- paste0(component, "_", suffix, ".xpt")
  local_path <- file.path(data_dir, fname)
  if (!file.exists(local_path)) {
    url <- sprintf("%s/%s/DataFiles/%s_%s.xpt", cdc_base, year_folder[suffix], component, suffix)
    message("Downloading ", url)
    ok <- tryCatch({ download.file(url, local_path, mode = "wb", quiet = TRUE); TRUE },
                    error = function(e) { message("  -> FAILED: ", e$message); FALSE })
    if (!ok || !file.exists(local_path) || file.size(local_path) == 0) {
      warning("Could not download ", fname, " -- check CDC file naming for this cycle/component.")
      return(NULL)
    }
  }
  haven::read_xpt(local_path)
}

## ---- 3. Load and merge one cycle --------------------------------------------

load_cycle <- function(suffix) {
  message("== Loading cycle suffix: ", suffix, " ==")

  demo    <- get_xpt("DEMO",    suffix)
  bmx     <- get_xpt("BMX",     suffix)
  diq     <- get_xpt("DIQ",     suffix)
  ghb     <- get_xpt("GHB",     suffix)
  biopro  <- get_xpt("BIOPRO",  suffix)   # LBXSCR (serum creatinine)
  albcr   <- get_xpt("ALB_CR",  suffix)   # URXUMA, URXUCR
  mcq     <- get_xpt("MCQ",     suffix)
  trigly  <- get_xpt("TRIGLY",  suffix)   # LBXTR, LBDLDL, WTSAF2YR
  hdl     <- get_xpt("HDL",     suffix)   # LBDHDD
  tchol   <- get_xpt("TCHOL",   suffix)   # LBXTC
  rxq     <- get_xpt("RXQ_RX",  suffix)   # RXDDRUG (drug name text)
  bpx     <- get_xpt("BPX",     suffix)   # BPXSY1-3 (systolic BP readings) -- needed for PREVENT
  smq     <- get_xpt("SMQ",     suffix)   # SMQ020, SMQ040 (smoking history/current) -- needed for PREVENT
  hiq     <- get_xpt("HIQ",     suffix)   # HIQ011 (covered by health insurance) -- needed for disparities analysis

  apob <- if (suffix == "D") {
    trigly %>% select(SEQN, any_of("LBXAPB"))   # ApoB is bundled in TRIGLY_D for 2005-2006
  } else {
    get_xpt("APOB", suffix) %>% select(SEQN, LBXAPB)
  }

  if (is.null(demo) || is.null(trigly)) {
    warning("Missing core files for cycle ", suffix, "; skipping.")
    return(NULL)
  }

  # Statin use: collapse RXQ_RX (long, one row per medication) to person-level flag
  statin_names <- c("ATORVASTATIN", "SIMVASTATIN", "ROSUVASTATIN", "PRAVASTATIN",
                     "LOVASTATIN", "FLUVASTATIN", "PITAVASTATIN")
  statin_flag <- if (!is.null(rxq) && "RXDDRUG" %in% names(rxq)) {
    rxq %>%
      mutate(is_statin = str_detect(toupper(RXDDRUG), paste(statin_names, collapse = "|"))) %>%
      group_by(SEQN) %>%
      summarise(on_statin = as.integer(any(is_statin, na.rm = TRUE)), .groups = "drop")
  } else {
    tibble(SEQN = numeric(0), on_statin = integer(0))
  }

  # Antihypertensive use: same name-matching approach as statin_flag above --
  # a pragmatic proxy covering the most commonly prescribed generic names
  # across the five major antihypertensive classes (ACE inhibitors, ARBs,
  # beta blockers, calcium channel blockers, thiazide/thiazide-like
  # diuretics), not a full lookup against NHANES's Multum Lexicon
  # therapeutic-class crosswalk. Same caveat as the statin list: will miss
  # less common agents.
  antihtn_names <- c(
    "LISINOPRIL","ENALAPRIL","RAMIPRIL","BENAZEPRIL","CAPTOPRIL","QUINAPRIL","FOSINOPRIL","MOEXIPRIL","TRANDOLAPRIL",
    "LOSARTAN","VALSARTAN","IRBESARTAN","OLMESARTAN","CANDESARTAN","TELMISARTAN","AZILSARTAN",
    "AMLODIPINE","NIFEDIPINE","DILTIAZEM","VERAPAMIL","FELODIPINE","NICARDIPINE","ISRADIPINE",
    "METOPROLOL","ATENOLOL","PROPRANOLOL","CARVEDILOL","BISOPROLOL","NADOLOL","LABETALOL","NEBIVOLOL",
    "HYDROCHLOROTHIAZIDE","CHLORTHALIDONE","INDAPAMIDE","FUROSEMIDE","BUMETANIDE","TORSEMIDE",
    "SPIRONOLACTONE","EPLERENONE","TRIAMTERENE",
    "DOXAZOSIN","TERAZOSIN","PRAZOSIN","CLONIDINE","HYDRALAZINE","METHYLDOPA"
  )
  antihtn_flag <- if (!is.null(rxq) && "RXDDRUG" %in% names(rxq)) {
    rxq %>%
      mutate(is_antihtn = str_detect(toupper(RXDDRUG), paste(antihtn_names, collapse = "|"))) %>%
      group_by(SEQN) %>%
      summarise(on_antihtn = as.integer(any(is_antihtn, na.rm = TRUE)), .groups = "drop")
  } else {
    tibble(SEQN = numeric(0), on_antihtn = integer(0))
  }

  # Non-statin lipid-lowering therapy: ezetimibe (cholesterol absorption
  # inhibitor, approved 2002 -- present throughout the study window),
  # PCSK9 inhibitors evolocumab/alirocumab (approved mid-2015, so these can
  # only appear in the 2015-2016 cycle I of this pooled sample), and
  # bempedoic acid / inclisiran included for completeness even though both
  # were approved (2020, 2021) after NHANES's ApoB measurement window ended
  # -- those two flags will read zero in every cycle here, which is
  # expected, not a bug. Captured as a SEPARATE flag from statin_flag above
  # so treatment-adequacy analysis (Section 12b) can distinguish "on no
  # lipid therapy at all" from "on therapy but statin-free."
  nonstatin_lipid_names <- c("EZETIMIBE", "EVOLOCUMAB", "ALIROCUMAB", "INCLISIRAN", "BEMPEDOIC")
  nonstatin_lipid_flag <- if (!is.null(rxq) && "RXDDRUG" %in% names(rxq)) {
    rxq %>%
      mutate(is_nonstatin_lipid = str_detect(toupper(RXDDRUG), paste(nonstatin_lipid_names, collapse = "|"))) %>%
      group_by(SEQN) %>%
      summarise(on_nonstatin_lipid_tx = as.integer(any(is_nonstatin_lipid, na.rm = TRUE)), .groups = "drop")
  } else {
    tibble(SEQN = numeric(0), on_nonstatin_lipid_tx = integer(0))
  }

  # RIDRETH3 was only introduced in the 2011-2012 (G) cycle; 2005-2010
  # (D/E/F) demo files carry RIDRETH1 instead. any_of() keeps this from
  # erroring on the earlier cycles -- downstream code only relies on
  # RIDRETH3 for the primary (G/H/I) regression, where it is always present.
  df <- demo %>%
    select(SEQN, RIDAGEYR, RIAGENDR, any_of(c("RIDRETH3", "RIDRETH1")),
           RIDEXPRG, SDMVPSU, SDMVSTRA, WTINT2YR, WTMEC2YR, INDFMPIR) %>%
    left_join(trigly %>% select(SEQN, LBXTR, LBDLDL, WTSAF2YR), by = "SEQN") %>%
    left_join(apob, by = "SEQN") %>%
    left_join(hdl   %>% select(SEQN, LBDHDD), by = "SEQN") %>%
    left_join(tchol %>% select(SEQN, LBXTC),  by = "SEQN") %>%
    { if (!is.null(bmx))    left_join(., bmx    %>% select(SEQN, BMXBMI), by = "SEQN") else mutate(., BMXBMI = NA_real_) } %>%
    { if (!is.null(diq))    left_join(., diq    %>% select(SEQN, DIQ010), by = "SEQN") else mutate(., DIQ010 = NA_real_) } %>%
    { if (!is.null(ghb))    left_join(., ghb    %>% select(SEQN, LBXGH),  by = "SEQN") else mutate(., LBXGH = NA_real_) } %>%
    { if (!is.null(biopro)) left_join(., biopro %>% select(SEQN, LBXSCR), by = "SEQN") else mutate(., LBXSCR = NA_real_) } %>%
    { if (!is.null(albcr))  left_join(., albcr  %>% select(SEQN, URXUMA, URXUCR), by = "SEQN") else mutate(., URXUMA = NA_real_, URXUCR = NA_real_) } %>%
    { if (!is.null(mcq))    left_join(., mcq %>% select(SEQN, any_of(c("MCQ160C","MCQ160D","MCQ160E","MCQ160F"))), by = "SEQN") else . } %>%
    { if (!is.null(bpx))    left_join(., bpx %>% select(SEQN, any_of(c("BPXSY1","BPXSY2","BPXSY3"))), by = "SEQN") else mutate(., BPXSY1 = NA_real_, BPXSY2 = NA_real_, BPXSY3 = NA_real_) } %>%
    { if (!is.null(smq))    left_join(., smq %>% select(SEQN, any_of(c("SMQ020","SMQ040"))), by = "SEQN") else mutate(., SMQ020 = NA_real_, SMQ040 = NA_real_) } %>%
    { if (!is.null(hiq))    left_join(., hiq %>% select(SEQN, any_of("HIQ011")), by = "SEQN") else mutate(., HIQ011 = NA_real_) } %>%
    left_join(statin_flag, by = "SEQN") %>%
    left_join(antihtn_flag, by = "SEQN") %>%
    left_join(nonstatin_lipid_flag, by = "SEQN") %>%
    mutate(cycle = suffix, on_statin = coalesce(on_statin, 0L), on_antihtn = coalesce(on_antihtn, 0L),
           on_nonstatin_lipid_tx = coalesce(on_nonstatin_lipid_tx, 0L))

  df
}

## ---- 4. Load all cycles ------------------------------------------------------

## Cached loader. The cached objects are PRE-DERIVATION -- they contain only
## raw merged NHANES variables, so they remain valid across CONFIG changes.
## Delete the cache directory to force a clean reload from the .xpt files.
load_cohort_cached <- function(suffixes, cache_name, label) {
  cf <- file.path(CONFIG$cache_dir, paste0(cache_name, ".rds"))
  if (CONFIG$use_cache && file.exists(cf)) {
    message("Loading ", label, " from cache: ", cf)
    return(readRDS(cf))
  }
  message("Loading ", label, " from .xpt files...")
  out <- bind_rows(purrr::map(suffixes, load_cycle))
  if (CONFIG$use_cache) {
    saveRDS(out, cf)
    message("  cached to ", cf)
  }
  out
}

primary_raw_3cyc <- load_cohort_cached(primary_suffixes, "primary_raw", "2011-2016 cycles")
trend_raw        <- load_cohort_cached(trend_suffixes,   "trend_raw",   "2005-2016 cycles")

## Select the cohort feeding primary_df and secondary_df, per CONFIG.
## primary_suffixes_used is carried forward so the pooled-weight divisor
## matches the number of cycles actually pooled -- getting this wrong would
## silently scale every population estimate.
if (CONFIG$primary_cycle_set == "all6") {
  primary_raw <- trend_raw
  primary_suffixes_used <- trend_suffixes
} else {
  primary_raw <- primary_raw_3cyc
  primary_suffixes_used <- primary_suffixes
}
message("Primary/secondary cohorts drawn from ", length(primary_suffixes_used),
        " cycles (", CONFIG$primary_cycle_set, "); pooled weight divisor = ",
        length(primary_suffixes_used), ".")

## ---- 5. Derived variables -----------------------------------------------------

## Sampson/NIH equation 2 (Sampson et al., JAMA Cardiology 2020). Coefficients
## verified against published sources rather than reproduced from memory:
##   LDL-C = TC/0.948 - HDL-C/0.971
##           - (TG/8.56 + TG*nonHDL-C/2140 - TG^2/16100) - 9.44
## All inputs and the result in mg/dL.
##
## WHICH EQUATION IS PRIMARY IS SET BY CONFIG$primary_ldl_equation at the top
## of this file. Whichever is selected becomes ldl_primary and drives every
## downstream result; the other becomes ldl_alt and appears only in the
## Section 11b comparison. Both are always computed, so switching the CONFIG
## value and re-sourcing is sufficient to flip the entire analysis.
##
## Default is "friedewald": it is the value NHANES publishes for every cycle in
## the analytic window, it is what most US laboratories reported during
## 2005-2016 and therefore what clinicians acted on, and it preserves
## comparability with prior NHANES lipid literature. "sampson" is the
## pre-specified sensitivity arm.
##
## WHY THE COMPARISON MATTERS -- what the sensitivity run found:
##   - Friedewald underestimated LDL-C in 88% of records (mean bias ~2.2
##     mg/dL). Small in absolute terms, but it lands right at the goal
##     thresholds where it flips classifications.
##   - 104 people (2.3%) moved OUT of at-LDL-goal status under Sampson;
##     only 2 moved in. Overwhelmingly one-directional.
##   - The misclassification was DIFFERENTIAL by treatment: 3.7% of
##     statin-treated adults changed status vs 1.9% of statin-naive.
##   - Consequence: overall discordance fell from 1.86% to 1.00%, the
##     at-goal reclassification share fell from 7.16% to 4.25%, and the
##     statin discordance gap (RQ8b) collapsed from 2.7x with separated
##     confidence intervals to 1.6x with heavily overlapping ones.
## Roughly 40-46% of what earlier versions of this script called "ApoB
## discordance" was the LDL-C equation rather than lipoprotein biology.
##
## GUIDELINE CONCORDANCE: the 2026 dyslipidemia guideline whose goals this
## script applies preferentially recommends Martin-Hopkins over Friedewald.
## Sampson is far closer to Martin-Hopkins than Friedewald is, and unlike
## Martin-Hopkins it is closed-form -- no 180-cell lookup table whose values
## would have to be transcribed correctly. Using Friedewald with this
## guideline's thresholds was internally inconsistent; using Sampson is not.
sampson_ldl <- function(tc, hdl_c, tg) {
  nonhdl <- tc - hdl_c
  tc / 0.948 - hdl_c / 0.971 -
    (tg / 8.56 + (tg * nonhdl) / 2140 - (tg^2) / 16100) - 9.44
}

## =============================================================================
## GOAL-ATTAINMENT COMPARATORS -- one definition, used everywhere
## =============================================================================
## The guideline states every lipid goal as a strict inequality ("apoB goal:
## <70 mg/dL"), so a value sitting exactly on the threshold is NOT at goal.
## These two functions are the only place that convention is written down.
## Every attainment comparison in this project calls them -- here, in
## precedence_and_model_analyses.R, precedence_diagnostics.R and compare_blind_coding.R -- so that
## CONFIG$strict_goal_boundaries switches all of them together and none can
## drift out of step with the primary analysis.
##
## Do not write a bare `<=` or `>` against a goal anywhere. That is how the
## inclusive convention survived four rounds of review: it was correct in the
## variable that was checked and wrong in four others that were not.
at_goal <- function(value, goal) {
  if (isTRUE(CONFIG$strict_goal_boundaries)) value <  goal else value <= goal
}
above_goal <- function(value, goal) {
  if (isTRUE(CONFIG$strict_goal_boundaries)) value >= goal else value >  goal
}

derive_vars <- function(df) {
  # RIDRETH1/RIDRETH3 harmonization guard. NHANES introduced RIDRETH3 in the
  # 2011-2012 cycle; 2005-2010 demo files carry RIDRETH1 only. Depending on
  # which cycles are pooled, one or both columns may be absent, so both are
  # created as NA if missing before the harmonization below references them.
  if (!"RIDRETH1" %in% names(df)) df$RIDRETH1 <- NA_real_
  if (!"RIDRETH3" %in% names(df)) df$RIDRETH3 <- NA_real_

  ## ---------------------------------------------------------------------------
  ## APOB CROSS-INSTRUMENT CALIBRATION -- must run BEFORE anything reads LBXAPB.
  ##
  ## Everything downstream depends on this: severe_hc (apoB >=140),
  ## above_apob_goal, apob_high_130, apob_markedly_high, the RQ9/RQ10 outcomes
  ## and the Section 3 cycle trend. Placing it here means the cached raw RDS
  ## stays valid -- the correction is applied on every run, not baked into the
  ## cache. LBXAPB_raw preserves the uncalibrated value for diagnostics.
  ##
  ## Verified effect at a fixed lipid profile (LDL-C 100, TG 100 mg/dL):
  ##   cycle D  89.4 -> 82.5 mg/dL,  apoB >=130  10.11% -> 5.62%
  ##   cycle I  86.6 -> 83.6 mg/dL,  apoB >=130   8.11% -> 6.63%
  ##   between-cycle spread 9.24 -> 3.46 mg/dL
  ## The empirical cycle-D excess (+7.9%) matches CDC's published factor
  ## (+8.4%), which is independent confirmation that this is instrument
  ## calibration and not a population difference.
  ## ---------------------------------------------------------------------------
  if (!"cycle" %in% names(df)) {
    stop("derive_vars(): 'cycle' column absent, so the apoB calibration cannot ",
         "be applied. Check load_cycle().")
  }
  df$LBXAPB_raw <- df$LBXAPB
  if (isTRUE(CONFIG$apply_apob_calibration)) {
    df$LBXAPB <- dplyr::case_when(
      df$cycle == "D" ~ df$LBXAPB * 0.923,           # BN100  -> ProSpec
      df$cycle == "I" ~ df$LBXAPB * 0.9988 - 2.880,  # Cobas  -> ProSpec
      TRUE            ~ df$LBXAPB
    )
    message("apoB calibrated to the BN ProSpec scale: ",
            sum(df$cycle == "D", na.rm = TRUE), " records rescaled (cycle D), ",
            sum(df$cycle == "I", na.rm = TRUE), " records rescaled (cycle I).")
  } else {
    message("apoB calibration DISABLED -- uncalibrated values in use ",
            "(sensitivity run only).")
  }

  df %>%
    mutate(
      # eGFR, CKD-EPI 2021 creatinine (race-free), Inker et al NEJM 2021
      egfr = case_when(
        RIAGENDR == 2 & LBXSCR <= 0.7 ~ 142 * (LBXSCR / 0.7)^-0.241 * 0.9938^RIDAGEYR * 1.012,
        RIAGENDR == 2 & LBXSCR >  0.7 ~ 142 * (LBXSCR / 0.7)^-1.200 * 0.9938^RIDAGEYR * 1.012,
        RIAGENDR == 1 & LBXSCR <= 0.9 ~ 142 * (LBXSCR / 0.9)^-0.302 * 0.9938^RIDAGEYR,
        RIAGENDR == 1 & LBXSCR >  0.9 ~ 142 * (LBXSCR / 0.9)^-1.200 * 0.9938^RIDAGEYR,
        TRUE ~ NA_real_
      ),
      uacr = ifelse(!is.na(URXUMA) & !is.na(URXUCR) & URXUCR > 0, (URXUMA / URXUCR) * 100, NA_real_),
      ckd = as.integer(coalesce(egfr < 60, FALSE) | coalesce(uacr >= 30, FALSE)),
      diabetes = as.integer(coalesce(DIQ010 == 1, FALSE) | coalesce(LBXGH >= 6.5, FALSE)),
      obesity = as.integer(coalesce(BMXBMI >= 30, FALSE)),
      hypertrig = as.integer(coalesce(LBXTR >= CONFIG$tg_threshold, FALSE)),
      ascvd_hx = as.integer(
        coalesce(MCQ160C == 1, FALSE) | coalesce(MCQ160D == 1, FALSE) |
        coalesce(MCQ160E == 1, FALSE) | coalesce(MCQ160F == 1, FALSE)
      ),
      nonhdl = LBXTC - LBDHDD,
      pregnant = as.integer(coalesce(RIDEXPRG == 1, FALSE)),

      # --- HARMONIZED RACE/ETHNICITY ------------------------------------
      # ADDED 2026-08-02. Extending the primary cohort to 2005-2016 exposed a
      # silent defect: RQ4's regression included factor(RIDRETH3), which does
      # not exist before 2011-2012. svyglm drops rows with missing predictors,
      # so the model quietly reverted to the 2011-2016 subset and returned
      # results identical to the 3-cycle run while every other section used
      # all six cycles. Nothing errored; the output simply described a
      # different sample than the surrounding text claimed.
      #
      # Harmonization collapses both coding schemes onto the common
      # denominator, which is RIDRETH1's five categories:
      #   RIDRETH3: 1 Mexican American, 2 Other Hispanic, 3 NH White,
      #             4 NH Black, 6 NH Asian, 7 Other/multiracial
      #   RIDRETH1: 1 Mexican American, 2 Other Hispanic, 3 NH White,
      #             4 NH Black, 5 Other incl. multiracial
      # RIDRETH3 codes 6 and 7 both map to 5.
      #
      # COST OF THIS CHOICE, stated plainly: Non-Hispanic Asian adults cannot
      # be separately identified in the pooled cohort, because NHANES did not
      # distinguish them before 2011. Analyses needing that category must be
      # restricted to 2011-2016. This is a genuine loss, accepted in exchange
      # for doubling the sample; it is not a workaround that recovers the
      # missing information.
      race_eth_code = coalesce(
        case_when(RIDRETH3 %in% c(6, 7) ~ 5,
                  !is.na(RIDRETH3)      ~ as.numeric(RIDRETH3),
                  TRUE                  ~ NA_real_),
        as.numeric(RIDRETH1)
      ),
      race_eth = factor(
        case_when(
          race_eth_code == 1 ~ "Mexican American",
          race_eth_code == 2 ~ "Other Hispanic",
          race_eth_code == 3 ~ "Non-Hispanic White",
          race_eth_code == 4 ~ "Non-Hispanic Black",
          race_eth_code == 5 ~ "Other or multiracial",
          TRUE ~ NA_character_
        ),
        levels = c("Non-Hispanic White", "Non-Hispanic Black",
                   "Mexican American", "Other Hispanic", "Other or multiracial")
      ),

      # Sampson-recomputed LDL-C. Whether this or NHANES's Friedewald LBDLDL
      # becomes the primary analytic value is decided immediately below by
      # CONFIG$primary_ldl_equation; both are always computed.
      ldl_sampson = sampson_ldl(LBXTC, LBDHDD, LBXTR),
      ldl_diff_sampson_minus_fried = ldl_sampson - LBDLDL,  # positive = Friedewald underestimated
      # CONFIG-selected primary LDL-C, and the other equation carried in
      # parallel for the Section 11b sensitivity comparison. Defined here
      # rather than in assign_goals() so both are available to the proxy
      # tiering below as well.
      ldl_primary = if (CONFIG$primary_ldl_equation == "sampson") ldl_sampson else LBDLDL,
      ldl_alt     = if (CONFIG$primary_ldl_equation == "sampson") LBDLDL else ldl_sampson,

      # --- Socioeconomic variables for the health-equity extension -------
      # income_tier from INDFMPIR (ratio of family income to the federal
      # poverty line for that survey year): <1 = below poverty, 1-<2 = low
      # income, >=2 = higher income -- a standard three-tier cut used
      # throughout the health-disparities literature. Left as NA (not
      # coalesced to a default tier) when INDFMPIR is missing, since NHANES
      # does not compute INDFMPIR for a meaningful share of respondents
      # (income refused, or reported only as a bracket) and silently
      # assigning those people a tier would misrepresent them.
      # insured from HIQ011 ("covered by health insurance or some other
      # kind of health care plan?"): 1=Yes -> insured=1, 2=No -> insured=0.
      # Refused/don't-know (7/9) and true missing are all left as NA rather
      # than assumed either way -- same reasoning as income_tier above.
      income_tier = case_when(
        is.na(INDFMPIR) ~ NA_character_,
        INDFMPIR < 1     ~ "below_poverty",
        INDFMPIR < 2     ~ "low_income",
        TRUE             ~ "higher_income"
      ),
      insured = case_when(
        HIQ011 == 1 ~ 1L,
        HIQ011 == 2 ~ 0L,
        TRUE ~ NA_integer_
      ),

      # --- New PREVENT-model inputs --------------------------------------
      # Average of available systolic BP readings (BPXSY1-3, stable across
      # all six cycles used here). current_smoker from SMQ040, which NHANES
      # only asks of people who report ever smoking >=100 cigarettes
      # (SMQ020==1); for never-smokers SMQ040 is legitimately skip-coded
      # missing, so coalescing to FALSE (not NA) below is correct, not a
      # data-quality workaround.
      avg_sbp = rowMeans(cbind(BPXSY1, BPXSY2, BPXSY3), na.rm = TRUE),
      current_smoker = as.integer(coalesce(SMQ040 %in% c(1, 2), FALSE)),

      # --- Simplified risk tier (PROXY -- retained only for comparison) --
      # This is the ORIGINAL diabetes/CKD-based proxy from the first version
      # of this script, kept under an explicit _proxy suffix (rather than
      # deleted) so it can be compared head-to-head against the real
      # PREVENT-based tiering computed below in Section 5b/5c -- see the RQ0
      # comparison after primary_design is built.
      risk_tier_proxy = case_when(
        ascvd_hx == 1 ~ "secondary_prevention",
        diabetes == 1 | ckd == 1 ~ "high_risk",
        TRUE ~ "borderline_intermediate"
      ),
      apob_goal_proxy = case_when(
        risk_tier_proxy == "secondary_prevention"    ~ 55,
        risk_tier_proxy == "high_risk"               ~ 70,
        risk_tier_proxy == "borderline_intermediate" ~ 90
      ),
      ldl_goal_proxy = case_when(
        risk_tier_proxy == "secondary_prevention"    ~ 55,
        risk_tier_proxy == "high_risk"               ~ 70,
        risk_tier_proxy == "borderline_intermediate" ~ 100
      ),
      # Uses ldl_primary, matching the main definition. RQ0 compares the proxy
      # TIERING method against PREVENT TIERING, so the LDL-C equation must be
      # held constant across both arms or the comparison would confound two
      # changes at once.
      at_ldl_goal_proxy  = as.integer(at_goal(ldl_primary, ldl_goal_proxy)),
      above_apob_goal_proxy = as.integer(above_goal(LBXAPB, apob_goal_proxy)),
      discordant_proxy = as.integer(at_ldl_goal_proxy == 1 & above_apob_goal_proxy == 1)
    )
}

## ---- 5b. PREVENT-based risk scoring (replaces the proxy going forward) ----
## Runs the AHA's actual published PREVENT-ASCVD equations (via the preventr
## package) batch-scored across the analytic sample. risk_tier/apob_goal/
## ldl_goal/at_ldl_goal/above_apob_goal/discordant are (re)assigned afterward
## in assign_goals() using real 10-year risk instead of the diabetes/CKD
## stand-in -- every RQ1 onward below then automatically runs on the
## PREVENT-based numbers without further code changes, since they all
## reference these same column names.

compute_prevent_risk <- function(df) {
  # preventr's validated input ranges are narrower than NHANES's raw values in
  # a small number of cases (BMI up to 39.9, SBP up to 180, eGFR up to 140).
  # Per the AHA's own PREVENT FAQ, out-of-range values should be clipped to
  # the nearest valid bound rather than dropped; n_clipped below makes that
  # visible rather than silent.
  clip <- function(x, lo, hi) pmin(pmax(x, lo), hi)

  df <- df %>%
    mutate(
      prevent_sbp     = clip(avg_sbp, 90, 180),
      prevent_total_c = clip(LBXTC,   130, 320),
      prevent_hdl_c   = clip(LBDHDD,  20, 100),
      prevent_egfr    = clip(egfr,    15, 140),
      prevent_bmi     = clip(BMXBMI,  18.5, 39.9),
      prevent_sex     = ifelse(RIAGENDR == 2, "female", "male")
    )

  n_clipped <- sum(
    (df$avg_sbp != df$prevent_sbp | df$LBXTC != df$prevent_total_c |
     df$LBDHDD != df$prevent_hdl_c | df$egfr != df$prevent_egfr |
     df$BMXBMI != df$prevent_bmi),
    na.rm = TRUE
  )
  message("PREVENT input clipping: ", n_clipped, " of ", nrow(df),
          " records had >=1 input value clipped to the equation's validated range.")

  # Score only rows with every required PREVENT input present. estimate_risk()'s
  # handling of missing values in batch mode was not independently verified
  # against a raw run in this environment, so completeness is enforced
  # explicitly here instead of relying on it. This is a genuinely NEW source
  # of sample attrition versus the old proxy -- BP and smoking data were
  # never required before, and NHANES's smoking/BP questionnaires have their
  # own missingness -- so primary_df/trend_df will end up somewhat smaller
  # than in the proxy version of this script. That's expected, not a bug.
  scorable <- df %>%
    filter(
      !is.na(RIDAGEYR), !is.na(RIAGENDR), !is.na(prevent_sbp),
      !is.na(on_antihtn), !is.na(prevent_total_c), !is.na(prevent_hdl_c),
      !is.na(on_statin), !is.na(diabetes), !is.na(current_smoker),
      !is.na(prevent_egfr), !is.na(prevent_bmi)
    )

  message("PREVENT-scorable records: ", nrow(scorable), " of ", nrow(df),
          " (", round(100 * nrow(scorable) / nrow(df), 1), "%) -- ",
          "the remainder are missing blood pressure, smoking status, or another required input.")

  scored <- estimate_risk(
    use_dat = scorable %>%
      transmute(
        SEQN, cycle,
        age = RIDAGEYR, sex = prevent_sex, sbp = prevent_sbp,
        bp_tx = as.logical(on_antihtn), total_c = prevent_total_c,
        hdl_c = prevent_hdl_c, statin = as.logical(on_statin),
        dm = as.logical(diabetes), smoking = as.logical(current_smoker),
        egfr = prevent_egfr, bmi = prevent_bmi
      ),
    time = "10yr",
    quiet = TRUE
  )

  out <- df %>%
    left_join(scored %>% select(SEQN, cycle, prevent_ascvd_10yr = ascvd), by = c("SEQN", "cycle"))

  # 30-YEAR PREVENT-ASCVD RISK. Needed only by the RQ10 replication, whose
  # cohort definition uses a 30-year threshold alongside the 10-year one.
  # Wrapped in tryCatch because the 30-year equations are validated over a
  # narrower age band than the 10-year ones (the AHA estimates 30-year risk for
  # younger adults), and preventr's handling of out-of-range ages was not
  # verified here. If the call fails or returns nothing usable, the column is
  # left NA and RQ10 falls back to the 10-year criterion alone -- which is a
  # DEVIATION from the replicated cohort definition and is reported as such
  # rather than passed over silently.
  scored30 <- tryCatch(
    estimate_risk(
      use_dat = scorable %>%
        transmute(
          SEQN, cycle,
          age = RIDAGEYR, sex = prevent_sex, sbp = prevent_sbp,
          bp_tx = as.logical(on_antihtn), total_c = prevent_total_c,
          hdl_c = prevent_hdl_c, statin = as.logical(on_statin),
          dm = as.logical(diabetes), smoking = as.logical(current_smoker),
          egfr = prevent_egfr, bmi = prevent_bmi
        ),
      time = "30yr", quiet = TRUE
    ),
    error = function(e) { message("30-year PREVENT unavailable: ", conditionMessage(e)); NULL }
  )

  if (!is.null(scored30) && "ascvd" %in% names(scored30)) {
    out <- out %>%
      left_join(scored30 %>% select(SEQN, cycle, prevent_ascvd_30yr = ascvd),
                by = c("SEQN", "cycle"))
    message("30-year PREVENT scored for ", sum(!is.na(out$prevent_ascvd_30yr)), " records.")
  } else {
    out$prevent_ascvd_30yr <- NA_real_
    message("30-year PREVENT NOT available -- RQ10 will use the 10-year criterion only.")
  }
  out
}

assign_goals <- function(df) {
  df %>%
    mutate(
      # Guideline-based risk tier. VERIFIED against the AHA "Top Things to
      # Know" summary (updated 13 March 2026):
      #   - Risk categories: borderline 3% to <5%; intermediate 5% to <10%.
      #   - "For those at borderline or intermediate risk, the LDL-C goal
      #     should be less than 100 mg/dL, and for those at high risk, the
      #     goal is less than 70 mg/dL."
      # Therefore the 70/70 goal pair begins at 10% 10-year PREVENT-ASCVD
      # risk. Borderline and intermediate collapse into one goal tier because
      # the guideline assigns them the SAME goal -- this is not a
      # simplification on our part, it is how the guideline is written.
      #
      # CORRECTED 2026-08-01: this previously used a 5% cut, which assigned the
      # stricter 70/70 pair to the entire 5-10% intermediate band. That was
      # wrong and inflated every above-goal estimate in the script.
      #
      # Known ASCVD keeps its own tier regardless of PREVENT score (in
      # practice unreachable in primary_df/trend_df, since ascvd_hx==1 is
      # filtered out upstream).
      ## ===================================================================
      ## REBUILT 2026-08-15 from Figure 1 of the guideline PDF
      ## ("Lipoprotein Goals for ASCVD Risk Reduction", JACC 2026;87:2643).
      ##
      ## WHAT WAS WRONG BEFORE: the previous version routed every participant
      ## through a single PREVENT-based pathway. Figure 1 is organised by
      ## PATIENT POPULATION, not by risk level, and PREVENT appears in only
      ## ONE of its six rows. Diabetes, severe hypercholesterolaemia,
      ## hypertriglyceridaemia, subclinical atherosclerosis, and clinical
      ## ASCVD each have their own column-assignment criteria. Applying the
      ## primary-prevention pathway to a person with diabetes assigns them a
      ## goal by a rule the guideline does not use for them.
      ##
      ## SECOND DEFECT: apoB goals were assigned to everyone. In primary
      ## prevention Figure 1 states the apoB goal only "If TG >=150 mg/dL to
      ## 499 mg/dL". For a participant below that band the guideline specifies
      ## NO apoB goal, so "above apoB goal" is UNDEFINED for them -- not FALSE.
      ## apob_goal is left NA in that case and every apoB analysis must
      ## restrict to !is.na(apob_goal).
      ##
      ## COLUMNS (goal sets) from Figure 1:
      ##   col 1: LDL-C <100, non-HDL-C <130
      ##   col 2: LDL-C <70,  non-HDL-C <100
      ##   col 3: LDL-C <55,  non-HDL-C <85
      ##
      ## PRECEDENCE. Rows are not mutually exclusive (a person may have
      ## diabetes AND hypertriglyceridaemia AND ASCVD). Figure 1 gives no
      ## precedence rule, so one is imposed here, most-specific first, and
      ## must be stated in the manuscript as an analytic decision:
      ##   clinical ASCVD > diabetes > severe hypercholesterolaemia >
      ##   hypertriglyceridaemia > primary prevention (PREVENT)
      ## Section 11c reports a sensitivity using the most stringent applicable
      ## goal instead, since that ordering is defensible too.
      ##
      ## NOT IMPLEMENTABLE IN NHANES -- state as limitations:
      ##   * Subclinical atherosclerosis row: requires CAC. NHANES has none.
      ##   * FH status: cannot be genetically confirmed; severe
      ##     hypercholesterolaemia is identified by LDL-C >=190 only.
      ##   * Diabetes-specific risk modifiers: eGFR and albuminuria are
      ##     available; retinopathy, neuropathy and diabetes duration are not.
      ## ===================================================================

      ## --- Supporting flags -------------------------------------------------
      ## "ASCVD risk factors" is not enumerated in Figure 1; the conventional
      ## set is used and counted so the threshold is visible and auditable.
      ## NOTE: built only from variables that actually exist in derive_vars().
      ## An earlier draft of this block referenced current_smoker,
      ## hypertension, albuminuria and mcq_* -- none of which are defined
      ## anywhere in this script. They were invented. Smoking (SMQ020/SMQ040)
      ## and blood pressure (BPXSY1) are loaded for PREVENT but never derived
      ## into flags, so they are derived inline here.
      smoker_flag = as.integer(coalesce(SMQ040 %in% c(1, 2), FALSE)),
      ## CORRECTED 2026-08-27 (script review, finding 2). Two defects in the
      ## previous definition, which was as.integer(coalesce(BPXSY1 >= 140, FALSE)).
      ##
      ##  1. IT IGNORED TREATMENT. An adult on antihypertensives with a
      ##     controlled pressure of 124 was not counted as hypertensive. The
      ##     guideline's high-risk-condition list says "hypertension" without
      ##     qualification, and on_antihtn is already derived and already passed
      ##     to PREVENT as bp_tx below -- so the risk model was treating these
      ##     people as hypertensive while the guideline's own risk-factor count
      ##     was not. Measured effect: 1,525 flagged (17.1%) versus 3,455
      ##     (38.7%) once treatment is counted.
      ##  2. IT USED THE FIRST READING. avg_sbp averages all available readings
      ##     and is what PREVENT consumes; BPXSY1 alone runs about 0.9 mmHg
      ##     higher and reclassified 246 people in one direction against 72 in
      ##     the other.
      ##
      ## Systolic only: BPXDI1 is not among the columns selected from BPX, so
      ## isolated diastolic hypertension is still missed. Uncommon in this age
      ## range -- retained as a stated limitation.
      ##
      ## Set CONFIG$htn_counts_treatment = FALSE to reproduce the old definition.
      htn_flag = if (isTRUE(CONFIG$htn_counts_treatment)) {
        as.integer(coalesce(avg_sbp >= CONFIG$htn_sbp_threshold, FALSE) |
                   coalesce(on_antihtn == 1, FALSE))
      } else {
        as.integer(coalesce(BPXSY1 >= CONFIG$htn_sbp_threshold, FALSE))
      },
      low_hdl     = as.integer(coalesce(
                       (RIAGENDR == 1 & LBDHDD < 40) |
                       (RIAGENDR == 2 & LBDHDD < 50), FALSE)),
      ## The age term is the only one that moves between the two readings of
      ## "ASCVD risk factor"; smoking, hypertension, low HDL-C and CKD are
      ## worded the same either way. See CONFIG$ascvd_rf_age_sex_specific.
      age_rf = if (isTRUE(CONFIG$ascvd_rf_age_sex_specific)) {
        as.integer(coalesce((RIAGENDR == 1 & RIDAGEYR >= 55) |
                            (RIAGENDR == 2 & RIDAGEYR >= 65), FALSE))
      } else {
        as.integer(coalesce(RIDAGEYR >= 55, FALSE))
      },
      ascvd_rf_count = smoker_flag + htn_flag + low_hdl +
                       as.integer(coalesce(ckd == 1, FALSE)) + age_rf,
      has_ascvd_rf = as.integer(ascvd_rf_count >= 1),

      ## Diabetes-specific risk modifiers, restricted to what NHANES supports.
      ## uacr and egfr are both derived above; retinopathy, neuropathy and
      ## diabetes duration are not collected and cannot be approximated.
      dm_risk_modifier = as.integer(coalesce(egfr < 60, FALSE) |
                                    coalesce(uacr >= 30, FALSE)),

      ## Severe hypercholesterolaemia: LDL-C >=190, or apoB >=140, the
      ## apoB-equivalent the guideline pairs with it (JACC 2026;87:2657).
      severe_hc = as.integer(coalesce(ldl_primary >= 190, FALSE) |
                             coalesce(LBXAPB     >= 140, FALSE)),

      ## Figure 1's apoB goal in primary prevention is conditional on this band.
      tg_in_goal_band = as.integer(coalesce(LBXTR >= 150 & LBXTR <= 499, FALSE)),

      ## --- VERY-HIGH-RISK ASCVD ---------------------------------------------
      ## CORRECTED 2026-08-15 against the guideline's own enumerated definition
      ## (JACC 2026;87, secondary-prevention section, table footnote), verbatim:
      ##
      ##   "Very high risk includes a history of MULTIPLE MAJOR ASCVD EVENTS
      ##    (ACS within past 12 months, history of MI [other than ACS above],
      ##    history of ischemic stroke, symptomatic PAD) or 1 MAJOR ASCVD EVENT
      ##    AND MULTIPLE HIGH-RISK CONDITIONS (age >65 years of age, coronary
      ##    artery revascularization, current smoker, diabetes, history of HF,
      ##    hypertension, LDL-C >100 mg/dL despite maximally tolerated statin +
      ##    ezetimibe)."
      ##
      ## THREE ERRORS IN THE PREVIOUS VERSION:
      ##  1. It counted MCQ160C (CHD) and MCQ160D (angina) as major ASCVD
      ##     events. They are NOT on the guideline's list, which names only
      ##     ACS, MI, ischaemic stroke and symptomatic PAD. NHANES respondents
      ##     frequently report both CHD and angina, so someone who has had
      ##     neither an MI nor a stroke could be classified as having
      ##     "multiple major events". That is a misreading, not a simplification.
      ##  2. It counted CKD as a high-risk condition. CKD is not on the list.
      ##  3. It omitted hypertension, which IS on the list and is already
      ##     computed above as htn_flag.
      ##
      ## MAJOR EVENTS available in NHANES: MI (MCQ160E) and stroke (MCQ160F)
      ## only. ACS timing, symptomatic PAD, revascularisation and heart failure
      ## are not collected, so this remains an approximation -- but now an
      ## approximation of the correct definition rather than of a wrong one.
      ## Because two of the four major events are unobservable, this
      ## UNDER-counts very-high-risk status; state that direction explicitly.
      ascvd_major_event_count = as.integer(coalesce(MCQ160E == 1, FALSE)) +
                                as.integer(coalesce(MCQ160F == 1, FALSE)),

      ## High-risk conditions, restricted to those NHANES supports.
      ## Not available: coronary revascularisation, history of HF, and
      ## "LDL-C >100 despite maximally tolerated statin + ezetimibe" (drug
      ## dose and tolerance are not recorded).
      ascvd_highrisk_cond_count = as.integer(coalesce(RIDAGEYR > 65, FALSE)) +
                                  smoker_flag +
                                  as.integer(coalesce(diabetes == 1, FALSE)) +
                                  htn_flag,

      ## Very-high-risk clinical ASCVD. Two routes, and the second is a
      ## guideline-adjudication decision rather than arithmetic:
      ##   (i)  Figure 10's criteria -- >=2 major events, or 1 major event plus
      ##        >=2 high-risk conditions.
      ##   (ii) Figure 11's entry box -- "Clinical ASCVD at very high risk
      ##        (includes ASCVD with CKD)". Switched by
      ##        CONFIG$ascvd_ckd_very_high; see the note in CONFIG for why the
      ##        two figures disagree and why (ii) is adopted.
      ## Route (ii) is what makes CKD confer the apoB <55 goal, because
      ## g_ascvd_apob keys off ascvd_very_high.
      ascvd_very_high = as.integer(
        coalesce(ascvd_hx == 1, FALSE) &
        (ascvd_major_event_count >= 2 |
         (ascvd_major_event_count >= 1 & ascvd_highrisk_cond_count >= 2) |
         (isTRUE(CONFIG$ascvd_ckd_very_high) & coalesce(ckd == 1, FALSE)))),

      ## Retained so the two routes can be counted separately in the
      ## diagnostics: who qualifies by Figure 10's arithmetic alone.
      ascvd_very_high_fig10 = as.integer(
        coalesce(ascvd_hx == 1, FALSE) &
        (ascvd_major_event_count >= 2 |
         (ascvd_major_event_count >= 1 & ascvd_highrisk_cond_count >= 2))),

      ## Retained only for the diagnostic below; not used in goal assignment.
      ascvd_event_count = as.integer(coalesce(MCQ160C == 1, FALSE)) +
                          as.integer(coalesce(MCQ160D == 1, FALSE)) +
                          as.integer(coalesce(MCQ160E == 1, FALSE)) +
                          as.integer(coalesce(MCQ160F == 1, FALSE)),

      ## =====================================================================
      ## GOAL ASSIGNMENT -- REBUILT 2026-08-16 (reviewer major comment 1)
      ##
      ## WHAT WAS WRONG: goals were assigned by walking the six Figure 1 rows
      ## in a fixed order and stopping at the first match. Two rows state no
      ## apoB goal for some people who match them -- severe hypercholesterol-
      ## aemia (always, absent FH confirmation and therapy-response data) and
      ## hypertriglyceridaemia (outside age 40-75 with a risk factor). Under
      ## first-match those SILENT rows absorbed people and cancelled a goal
      ## that a different matching row stated explicitly.
      ##
      ## The consequence was measured, not assumed: across all 120 orderings
      ## of the five implementable rows the goal-assigned share ran from
      ## 14.7% to 32.3%. The answer depended on our sort order, not on the
      ## patient. That is not publishable.
      ##
      ## THE RULE NOW USED, and it contains no ordering at all:
      ##   1. MEMBERSHIP -- a person has a goal if ANY row that describes them
      ##      states one. A silent row is silent; it does not revoke a goal
      ##      stated elsewhere.
      ##   2. VALUE -- where two rows state different goals, the STRICTER
      ##      (numerically lower) applies. Someone meeting two rows carries
      ##      more risk features, not fewer, and the guideline tightens goals
      ##      as risk accumulates.
      ##
      ## Both steps are order-free: pmin() is commutative, so no permutation
      ## of the rows can change any result. State in the methods that NO
      ## precedence rule is used at any point.
      ##
      ## The first-match assignment is retained in parallel as
      ## *_firstmatch for the labelled sensitivity comparison.
      ## =====================================================================

      ## --- Per-row apoB goals. NA = this row states no apoB goal. ----------
      g_ascvd_apob = case_when(
        ascvd_hx != 1 | is.na(ascvd_hx)                     ~ NA_real_,
        ascvd_very_high == 1                                ~ 55,
        hypertrig == 1                                      ~ 70,
        TRUE                                                ~ NA_real_
      ),
      g_dm_apob = case_when(
        diabetes != 1 | is.na(diabetes)                     ~ NA_real_,
        has_ascvd_rf == 1 | dm_risk_modifier == 1           ~ 70,
        TRUE                                                ~ 90
      ),
      ## Severe hypercholesterolaemia states no apoB goal we can operational-
      ## ise: FH cannot be genetically confirmed in NHANES and the guideline's
      ## therapy-response conditions are unobservable. Silent by construction.
      g_shc_apob = NA_real_,
      g_htg_apob = case_when(
        hypertrig != 1 | is.na(hypertrig)                   ~ NA_real_,
        RIDAGEYR >= 40 & RIDAGEYR <= 75 & has_ascvd_rf == 1 ~ 70,
        TRUE                                                ~ NA_real_
      ),
      g_pp_apob = case_when(
        is.na(prevent_ascvd_10yr)                              ~ NA_real_,
        tg_in_goal_band != 1                                   ~ NA_real_,
        prevent_ascvd_10yr >= CONFIG$prevent_goal_high_risk_cut ~ 70,
        TRUE                                                   ~ 90
      ),

      ## --- Per-row LDL-C goals. Every row states one, so membership is
      ## --- never in question here; only the value can differ.
      ## The CKD branch is Figure 1's "With CKD" bullet in the <55 column of
      ## the Clinical ASCVD row. See CONFIG$ascvd_ckd_ldl55 for why this
      ## changes no apoB goal and no primary-prevention estimate.
      g_ascvd_ldl = case_when(
        ascvd_hx != 1 | is.na(ascvd_hx)           ~ NA_real_,
        ascvd_very_high == 1                      ~ 55,
        isTRUE(CONFIG$ascvd_ckd_ldl55) &
          coalesce(ckd == 1, FALSE)               ~ 55,
        TRUE                                      ~ 70
      ),
      g_dm_ldl = case_when(
        diabetes != 1 | is.na(diabetes)           ~ NA_real_,
        has_ascvd_rf == 1 | dm_risk_modifier == 1 ~ 70,
        TRUE                                      ~ 100
      ),
      g_shc_ldl = case_when(
        severe_hc != 1 | is.na(severe_hc)         ~ NA_real_,
        has_ascvd_rf == 1                         ~ 70,
        TRUE                                      ~ 100
      ),
      g_htg_ldl = case_when(
        hypertrig != 1 | is.na(hypertrig)                   ~ NA_real_,
        RIDAGEYR >= 40 & RIDAGEYR <= 75 & has_ascvd_rf == 1 ~ 70,
        TRUE                                                ~ 100
      ),
      g_pp_ldl = case_when(
        is.na(prevent_ascvd_10yr)                              ~ NA_real_,
        prevent_ascvd_10yr >= CONFIG$prevent_goal_high_risk_cut ~ 70,
        TRUE                                                   ~ 100
      ),

      ## --- Most stringent applicable goal, across all matched rows ---------
      ## pmin(..., na.rm = TRUE) returns Inf when every argument is NA, so the
      ## Inf is mapped back to NA: "no row states a goal for this person".
      apob_goal = {
        v <- pmin(g_ascvd_apob, g_dm_apob, g_shc_apob, g_htg_apob, g_pp_apob,
                  na.rm = TRUE)
        ifelse(is.finite(v), v, NA_real_)
      },
      ldl_goal = {
        v <- pmin(g_ascvd_ldl, g_dm_ldl, g_shc_ldl, g_htg_ldl, g_pp_ldl,
                  na.rm = TRUE)
        ifelse(is.finite(v), v, NA_real_)
      },

      ## --- Descriptive only: which rows a person matches. No longer drives
      ## --- any goal. Retained for tabulation and for Appendix S2.
      n_rows_matched = as.integer(coalesce(ascvd_hx == 1, FALSE)) +
                       as.integer(coalesce(diabetes == 1, FALSE)) +
                       as.integer(coalesce(severe_hc == 1, FALSE)) +
                       as.integer(coalesce(hypertrig == 1, FALSE)) +
                       as.integer(!is.na(prevent_ascvd_10yr)),
      goal_row = case_when(
        ascvd_hx == 1              ~ "clinical_ascvd",
        diabetes == 1              ~ "diabetes",
        severe_hc == 1             ~ "severe_hypercholesterolaemia",
        hypertrig == 1             ~ "hypertriglyceridaemia",
        !is.na(prevent_ascvd_10yr) ~ "primary_prevention",
        TRUE                       ~ NA_character_
      ),

      ## --- SENSITIVITY ARM: the superseded first-match assignment ----------
      apob_goal_firstmatch = case_when(
        goal_row == "clinical_ascvd" & ascvd_very_high == 1 ~ 55,
        goal_row == "clinical_ascvd" & hypertrig == 1       ~ 70,
        goal_row == "diabetes" & (has_ascvd_rf == 1 | dm_risk_modifier == 1) ~ 70,
        goal_row == "diabetes"                              ~ 90,
        goal_row == "hypertriglyceridaemia" &
          RIDAGEYR >= 40 & RIDAGEYR <= 75 & has_ascvd_rf == 1 ~ 70,
        goal_row == "primary_prevention" & tg_in_goal_band == 1 &
          prevent_ascvd_10yr >= CONFIG$prevent_goal_high_risk_cut ~ 70,
        goal_row == "primary_prevention" & tg_in_goal_band == 1 ~ 90,
        TRUE ~ NA_real_
      ),
      ldl_goal_firstmatch = case_when(
        goal_row == "clinical_ascvd" & ascvd_very_high == 1 ~ 55,
        ## Mirrors the union arm so the precedence comparison isolates ORDERING
        ## and nothing else. Without this the two arms would differ in cell
        ## mapping as well, and the sensitivity would measure both at once.
        goal_row == "clinical_ascvd" & isTRUE(CONFIG$ascvd_ckd_ldl55) &
          coalesce(ckd == 1, FALSE)                         ~ 55,
        goal_row == "clinical_ascvd"                        ~ 70,
        goal_row == "diabetes" & (has_ascvd_rf == 1 | dm_risk_modifier == 1) ~ 70,
        goal_row == "diabetes"                              ~ 100,
        goal_row == "severe_hypercholesterolaemia" & has_ascvd_rf == 1 ~ 70,
        goal_row == "severe_hypercholesterolaemia"          ~ 100,
        goal_row == "hypertriglyceridaemia" &
          RIDAGEYR >= 40 & RIDAGEYR <= 75 & has_ascvd_rf == 1 ~ 70,
        goal_row == "hypertriglyceridaemia"                 ~ 100,
        goal_row == "primary_prevention" &
          prevent_ascvd_10yr >= CONFIG$prevent_goal_high_risk_cut ~ 70,
        goal_row == "primary_prevention"                    ~ 100,
        TRUE ~ NA_real_
      ),
      ## Least-stringent variant, for the reviewer-anticipated question of how
      ## much the stringency rule itself matters. Membership is identical to
      ## apob_goal by construction; only the value differs.
      apob_goal_lenient = {
        v <- pmax(g_ascvd_apob, g_dm_apob, g_shc_apob, g_htg_apob, g_pp_apob,
                  na.rm = TRUE)
        ifelse(is.finite(v), v, NA_real_)
      },

      ## Optional/"consider" goals, tracked but never treated as targets.
      apob_goal_optional = case_when(
        goal_row == "clinical_ascvd" & ascvd_very_high == 0 ~ 55,
        TRUE ~ NA_real_
      ),

      has_apob_goal = as.integer(!is.na(apob_goal)),

      ## Retained for backward compatibility with downstream code and for
      ## comparison against the archived single-pathway results.
      risk_tier = case_when(
        ascvd_hx == 1 ~ "secondary_prevention",
        prevent_ascvd_10yr >= CONFIG$prevent_goal_high_risk_cut ~ "high_risk",
        !is.na(prevent_ascvd_10yr) ~ "borderline_intermediate",
        TRUE ~ NA_character_
      ),
      # PRIMARY LDL-C, selected by CONFIG$primary_ldl_equation. Every
      # downstream analysis in this script -- RQ0 through RQ9, the secondary-
      # prevention extension, and the mortality grouping -- reads at_ldl_goal
      # and discordant, so switching the CONFIG value propagates everywhere
      # without any further code changes. The non-selected equation is carried
      # in parallel as the *_alt variables for the Section 11b comparison.
      ## Strict (<) under CONFIG$strict_goal_boundaries; see the note there.
      at_ldl_goal  = as.integer(at_goal(ldl_primary, ldl_goal)),
      # ApoB goal attainment does NOT involve LDL-C, so this variable -- and
      # every result built on it (RQ1b overall ApoB non-attainment, RQ6b, the
      # whole RQ6d treatment-adequacy analysis, RQ6e income/insurance untreated
      # shares) -- is completely unaffected by the equation switch. Those
      # results are identical before and after, by construction.
      above_apob_goal = as.integer(above_goal(LBXAPB, apob_goal)),
      # PRIMARY OUTCOME: discordance = at LDL-C goal but above ApoB goal
      discordant = as.integer(at_ldl_goal == 1 & above_apob_goal == 1),

      # SENSITIVITY PARALLEL: the same two quantities computed from whichever
      # equation is NOT primary, carried alongside so Section 11b can quantify
      # how much of the discordance signal is equation-driven. Not used
      # anywhere else -- any downstream analysis referencing bare at_ldl_goal
      # or discordant is getting the CONFIG-selected primary equation.
      at_ldl_goal_alt = as.integer(at_goal(ldl_alt, ldl_goal)),
      discordant_alt  = as.integer(at_ldl_goal_alt == 1 & above_apob_goal == 1),
      ## The alt equation is not always computable. With Sampson primary, the
      ## cohort now includes adults with TG > 400, for whom NHANES publishes no
      ## Friedewald value -- those people must be labelled explicitly rather
      ## than falling through to "unchanged", which is what a bare TRUE branch
      ## would do with an NA comparison.
      ldl_goal_status_change = case_when(
        is.na(at_ldl_goal_alt)                  ~ "alt_equation_not_computable",
        at_ldl_goal_alt == 1 & at_ldl_goal == 0 ~ "alt_at_goal_primary_above",
        at_ldl_goal_alt == 0 & at_ldl_goal == 1 ~ "alt_above_primary_at_goal",
        TRUE ~ "unchanged"
      ),

      apob_flagged_group = as.integer(diabetes == 1 | ckd == 1 | obesity == 1 |
                                       hypertrig == 1 | at_ldl_goal == 1)
    )
}

## COHORT FILTER -- CORRECTED 2026-08-27 (script review, finding 1).
##
## The requirement used to be !is.na(LBDLDL). NHANES blanks LBDLDL whenever
## triglycerides exceed 400 mg/dL, because Friedewald is invalid there, so that
## filter silently removed every adult with TG > 400 -- 188 people, all of them
## scorable by Sampson, with a median apoB of 120.5 mg/dL against 92.3 in the
## retained sample. It excluded the highest-apoB group in the data, by an
## artefact of an equation the primary analysis does not use.
##
## It also broke the operationalisation of Figure 1: the guideline's primary-
## prevention apoB goal is stated for TG 150-499 mg/dL, and 90 of the excluded
## adults fall inside that band. The maximum triglyceride in the old cohort was
## 398 mg/dL against a stated band of 499.
##
## The filter now requires the CONFIG-selected primary equation to be
## computable. With sampson primary this recovers all 188; with friedewald
## primary it reproduces the old behaviour exactly, which is correct for that
## configuration.
primary_df <- derive_vars(primary_raw) %>%
  filter(
    RIDAGEYR >= CONFIG$age_min, RIDAGEYR <= CONFIG$age_max,
    pregnant == 0,
    ascvd_hx == 0,                 # primary analysis = primary prevention, mirrors source paper
    !is.na(LBXAPB), !is.na(ldl_primary), !is.na(LBXTR),
    WTSAF2YR > 0
  ) %>%
  compute_prevent_risk() %>%
  filter(!is.na(prevent_ascvd_10yr)) %>%   # drop anyone PREVENT couldn't score (missing BP/smoking/etc.)
  assign_goals() %>%
  mutate(n_cycles_pooled = length(primary_suffixes_used),
         wt_pooled = WTSAF2YR / n_cycles_pooled)

## ---------------------------------------------------------------------------
## B1 (round 6): DOES THE PREVENT COMPLETE-CASE EXCLUSION DRIVE THE 33.4%?
##
## 413 otherwise-eligible adults are dropped above because PREVENT could not be
## scored for them. 33.4% is therefore a percentage of the 9,108 who remain,
## and a reviewer will reasonably ask whether the excluded 413 differ.
##
## The union rule makes this answerable without imputing anything. PREVENT is
## needed by exactly ONE Figure 1 row -- primary prevention by estimated risk.
## The diabetes and hypertriglyceridaemia rows state goals without it. So for
## many of the 413 we can still settle MEMBERSHIP (does any row state a goal?)
## even when we cannot settle the VALUE (<70 or <90).
##
## That splits the 413 three ways:
##   KNOWN ASSIGNED   a non-PREVENT row states a goal. Membership certain.
##   KNOWN UNASSIGNED no row states one, and the PREVENT row could not have:
##                    its precondition (triglycerides 150-499) is not met, so
##                    the row is silent whatever the risk estimate.
##   UNDETERMINED     only the PREVENT row could have assigned them.
##
## Reported as bounds over all 9,521: lower treats every UNDETERMINED adult as
## unassigned, upper treats them all as assigned. If 33.4% sits inside a narrow
## band, complete-case exclusion is not driving it.
prevent_missing_bounds <- local({
  elig <- derive_vars(primary_raw) %>%
    filter(
      RIDAGEYR >= CONFIG$age_min, RIDAGEYR <= CONFIG$age_max,
      pregnant == 0, ascvd_hx == 0,
      !is.na(LBXAPB), !is.na(ldl_primary), !is.na(LBXTR),
      WTSAF2YR > 0
    ) %>%
    compute_prevent_risk() %>%
    assign_goals() %>%
    mutate(n_cycles_pooled = length(primary_suffixes_used),
           wt_pooled = WTSAF2YR / n_cycles_pooled,
           prevent_scored = as.integer(!is.na(prevent_ascvd_10yr)),
           ## goal from any row OTHER than primary-prevention-by-risk
           goal_without_prevent = {
             v <- pmin(g_ascvd_apob, g_dm_apob, g_shc_apob, g_htg_apob,
                       na.rm = TRUE)
             as.integer(is.finite(v))
           },
           ## could the PREVENT row ever have spoken for this person?
           prevent_row_possible = as.integer(coalesce(tg_in_goal_band == 1,
                                                      FALSE)),
           status = case_when(
             prevent_scored == 1            ~ "scored",
             goal_without_prevent == 1      ~ "known assigned",
             prevent_row_possible == 0      ~ "known unassigned",
             TRUE                           ~ "undetermined"),
           goal_lower = as.integer(coalesce(has_apob_goal == 1, FALSE) |
                                   (prevent_scored == 0 &
                                    goal_without_prevent == 1)),
           goal_upper = as.integer(goal_lower == 1 |
                                   status == "undetermined"))

  des <- svydesign(ids = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~wt_pooled,
                   data = elig, nest = TRUE)
  lo <- svyciprop(~I(goal_lower == 1), des, method = "beta", na.rm = TRUE)
  hi <- svyciprop(~I(goal_upper == 1), des, method = "beta", na.rm = TRUE)

  tab <- elig %>% count(status, name = "n_unweighted")
  print(as.data.frame(tab))

  out <- tibble(
    bound        = c("Lower (undetermined treated as NOT assigned)",
                     "Upper (undetermined treated as assigned)"),
    n_cohort     = nrow(elig),
    pct_weighted = 100 * c(as.numeric(lo), as.numeric(hi)),
    ci_low       = 100 * c(attr(lo, "ci")[1], attr(hi, "ci")[1]),
    ci_high      = 100 * c(attr(lo, "ci")[2], attr(hi, "ci")[2]),
    ci_method    = "korn-graubard")
  write.csv(out,  "r6_prevent_missing_bounds.csv", row.names = FALSE)
  write.csv(tab,  "r6_prevent_missing_status.csv", row.names = FALSE)
  message("\nB1 PREVENT missingness, over all ", nrow(elig), " eligible adults:")
  message("  complete-case estimate (n = ", sum(elig$prevent_scored), "): 33.4%")
  message("  bounds including the unscored: ",
          sprintf("%.1f%%", out$pct_weighted[1]), " to ",
          sprintf("%.1f%%", out$pct_weighted[2]))
  for (i in seq_len(nrow(tab)))
    message("    ", format(tab$status[i], width = 18), tab$n_unweighted[i])
  out
})

## ---------------------------------------------------------------------------
## B3 (round 6). PREVENT input ranges, and what the systolic clip can do.
##
## Dani asked for a rerun at an SBP ceiling of 200, or an 180-vs-200
## sensitivity. Neither is available: 90-180 is preventr's own validated range,
## not a setting of ours, and outside it estimate_risk() returns NA with
## input_problems set rather than a risk estimate. Raising the ceiling would
## therefore UNSCORE those adults, not score them differently -- a smaller
## denominator, not a corrected one. Scoring an SBP of 195 at all would mean
## hand-coding the published coefficients and overriding the package's
## validation, which is a weaker position than the clip.
##
## So the question is bounded instead, exactly as B1 bounds the unscored 413.
## The argument is one-way: clipping a value DOWN to 180 can only understate
## PREVENT risk, and risk enters Figure 1 only by crossing a threshold, so
## correcting the clip could only ADD a goal, never remove one. The maximum
## possible movement is therefore the weighted share of adults who are both
## above the ceiling AND currently unassigned -- anyone already carrying a
## goal cannot gain one.
##
## The per-variable breakdown is written out at the same time. Until now the
## run reported a single lump ("1344 of 9521 records had >=1 input clipped"),
## which says nothing about WHICH input, and so reads as though blood pressure
## drove it. It does not; body-mass index does, by an order of magnitude.
## Wrapped so a failure here cannot kill an hour-long run. The same
## lesson as kg_cell_summary(): a diagnostic must never be able to cost
## more than it is worth. If it fails, the run continues and the two
## CSVs are absent, which build_supplement.py reports rather than
## silently omitting eTable 18.
prevent_input_ranges <- tryCatch(local({
  elig <- derive_vars(primary_raw) %>%
    filter(
      RIDAGEYR >= CONFIG$age_min, RIDAGEYR <= CONFIG$age_max,
      pregnant == 0, ascvd_hx == 0,
      !is.na(LBXAPB), !is.na(ldl_primary), !is.na(LBXTR),
      WTSAF2YR > 0
    ) %>%
    compute_prevent_risk() %>%
    assign_goals() %>%
    mutate(n_cycles_pooled = length(primary_suffixes_used),
           wt_pooled = WTSAF2YR / n_cycles_pooled)

  ## The bounds below must be the ones compute_prevent_risk() actually applies.
  ## Stated once here and checked against that function's output, so the table
  ## cannot drift away from the code it documents.
  spec <- tibble::tribble(
    ~variable,                  ~raw,        ~clipped,          ~lo,   ~hi,
    "Systolic blood pressure",  "avg_sbp",   "prevent_sbp",     90,    180,
    "Total cholesterol",        "LBXTC",     "prevent_total_c", 130,   320,
    "HDL cholesterol",          "LBDHDD",    "prevent_hdl_c",   20,    100,
    "eGFR",                     "egfr",      "prevent_egfr",    15,    140,
    "Body mass index",          "BMXBMI",    "prevent_bmi",     18.5,  39.9)

  rows <- lapply(seq_len(nrow(spec)), function(i) {
    r <- spec[i, ]
    v <- elig[[r$raw]]
    c_ <- elig[[r$clipped]]
    ## the clipped column must be exactly pmin(pmax(raw, lo), hi)
    want <- pmin(pmax(v, r$lo), r$hi)
    stopifnot(all(is.na(want) == is.na(c_)),
              isTRUE(all.equal(want[!is.na(want)], c_[!is.na(c_)])))
    ok <- !is.na(v)
    tibble(variable = r$variable,
           lower_bound = r$lo, upper_bound = r$hi,
           n_nonmissing = sum(ok),
           n_below = sum(v[ok] < r$lo),
           n_above = sum(v[ok] > r$hi))
  })
  ranges <- dplyr::bind_rows(rows)
  ranges$pct_clipped <- 100 * (ranges$n_below + ranges$n_above) / nrow(elig)
  write.csv(ranges, "r6_prevent_input_ranges.csv", row.names = FALSE)

  ## ---- what the systolic ceiling can do to the headline -------------------
  ## DENOMINATOR. The range counts above are over all eligible adults, which is
  ## the frame the clipping happens in and the one the existing "N of 9,521
  ## records had >=1 input value clipped" message reports. The BOUND must not
  ## use that frame: over 9,521 the goal-assigned share is 33.33%, because the
  ## 413 PREVENT cannot score count as unassigned -- that is the eTable 16
  ## LOWER bound, not the headline. The headline is the complete-case estimate
  ## over the 9,108 scored, and a bound on the ceiling is only readable if it
  ## starts there. Clipping keeps the above-ceiling adults scorable, so none of
  ## them is lost by restricting to the analytic cohort.
  ##
  ## The first version of this block used the 9,521 frame and labelled the
  ## result "as analyzed", which it was not.
  analytic <- elig %>%
    filter(!is.na(prevent_ascvd_10yr)) %>%
    mutate(sbp_above_ceiling = as.integer(coalesce(avg_sbp > 180, FALSE)),
           ## could only GAIN a goal: above the ceiling and not assigned one
           could_gain = as.integer(sbp_above_ceiling == 1 &
                                   coalesce(has_apob_goal, 0L) == 0L),
           goal_upper_sbp = as.integer(coalesce(has_apob_goal == 1, FALSE) |
                                       could_gain == 1))

  des <- svydesign(ids = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~wt_pooled,
                   data = analytic, nest = TRUE)
  lo <- svyciprop(~I(has_apob_goal == 1),  des, method = "beta", na.rm = TRUE)
  hi <- svyciprop(~I(goal_upper_sbp == 1), des, method = "beta", na.rm = TRUE)

  ## The lower row must reproduce the headline exactly; if it ever stops doing
  ## so the denominator has drifted again and the table is not comparable.
  stopifnot(abs(100 * as.numeric(lo) - 33.4436949661563) < 0.01)

  bounds <- tibble(
    bound = c("As analyzed (systolic pressure clipped to 180 mm Hg)",
              "Upper: every adult above the ceiling gains a goal"),
    n_cohort            = nrow(analytic),
    n_eligible          = nrow(elig),
    n_above_ceiling     = sum(analytic$sbp_above_ceiling),
    n_could_gain        = sum(analytic$could_gain),
    pct_weighted        = 100 * c(as.numeric(lo), as.numeric(hi)),
    ci_low              = 100 * c(attr(lo, "ci")[1], attr(hi, "ci")[1]),
    ci_high             = 100 * c(attr(lo, "ci")[2], attr(hi, "ci")[2]),
    preventr_version    = as.character(utils::packageVersion("preventr")),
    ci_method           = "korn-graubard")
  write.csv(bounds, "r6_sbp_clip_bounds.csv", row.names = FALSE)

  message("\nB3 PREVENT input ranges, over ", nrow(elig), " eligible adults:")
  for (i in seq_len(nrow(ranges)))
    message("    ", format(ranges$variable[i], width = 24),
            "below ", format(ranges$n_below[i], width = 5),
            "above ", format(ranges$n_above[i], width = 5),
            sprintf("(%.1f%% clipped)", ranges$pct_clipped[i]))
  message("  above the systolic ceiling (analytic cohort, n = ",
          nrow(analytic), "): ", sum(analytic$sbp_above_ceiling),
          "; of those, unassigned and so able to change: ",
          sum(analytic$could_gain))
  message("  goal-assigned share: ", sprintf("%.2f%%", bounds$pct_weighted[1]),
          " as analyzed, at most ", sprintf("%.2f%%", bounds$pct_weighted[2]),
          " if every one of them gained a goal")
  message("  preventr version: ", bounds$preventr_version[1])
  list(ranges = ranges, bounds = bounds)
}), error = function(e) {
  message("*** B3 PREVENT input-range diagnostic FAILED: ",
          conditionMessage(e),
          "\n    The run continues; r6_prevent_input_ranges.csv and ",
          "r6_sbp_clip_bounds.csv were not written.")
  NULL
})


trend_df <- derive_vars(trend_raw) %>%
  filter(
    RIDAGEYR >= CONFIG$age_min, RIDAGEYR <= CONFIG$age_max,
    pregnant == 0,
    ascvd_hx == 0,
    !is.na(LBXAPB), !is.na(ldl_primary), !is.na(LBXTR),   # see note on primary_df
    WTSAF2YR > 0
  ) %>%
  compute_prevent_risk() %>%
  filter(!is.na(prevent_ascvd_10yr)) %>%
  assign_goals() %>%
  mutate(n_cycles_pooled = length(trend_suffixes),
         wt_pooled = WTSAF2YR / n_cycles_pooled)

## ---- 6. Survey design objects -------------------------------------------------

primary_design <- svydesign(
  id = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~wt_pooled,
  nest = TRUE, data = primary_df
)

trend_design <- svydesign(
  id = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~wt_pooled,
  nest = TRUE, data = trend_df
)

## =============================================================================
## 6b. DOMAIN-ANALYSIS VERIFICATION  (added 2026-09-10, round 4 comment 6)
## =============================================================================
## primary_design above is built on an already-filtered data frame. NCHS
## guidance is the other way round: keep the full design and restrict with
## subset(), so that records outside the analytic cohort still inform the
## variance through the strata and PSUs they belong to.
##
## This block builds the design the recommended way and reports both, rather
## than asserting they agree. It does NOT replace primary_design unless
## CONFIG$use_domain_design is TRUE and the estimates reconcile -- the switch
## is applied at the end of the block, after the comparison has printed, so a
## discrepancy is visible before anything downstream consumes it.
##
## Cost: PREVENT is scored over the whole fasting subsample rather than the
## cohort, which adds a few minutes to the run.

## GUARD: no bare goal comparison anywhere in the project.
## The inclusive convention survived four review rounds because it was written
## out longhand in ten places and only two were ever checked. This scans the
## source of every script in the pipeline for a lipid value compared directly
## against a goal variable, and fails the run if one reappears.
local({
  ## pathway_and_sensitivity_analyses.R was outside this list until 15 September
  ## 2026 -- it was written after the guard and never added, so the newest code
  ## in the pipeline was the only code not checked. It is clean; it is listed
  ## now so it stays that way.
  files <- c("nhanes_apob_analysis.R", "precedence_and_model_analyses.R",
             "precedence_diagnostics.R", "table1_and_figures.R",
             "pathway_and_sensitivity_analyses.R", "compare_blind_coding.R",
             "verify_run.R")
  files <- files[file.exists(files)]
  pat <- paste0("(LBXAPB|ldl_primary|ldl_alt)\\s*(<=|>=|<|>)\\s*",
                "[A-Za-z0-9_$]*(apob_goal|ldl_goal)")
  hits <- unlist(lapply(files, function(f) {
    ln <- readLines(f, warn = FALSE)
    ln <- ln[!grepl("^\\s*#", ln)]          # comments are allowed to show the pattern
    m  <- grep(pat, ln, value = TRUE)
    if (length(m)) paste0(f, ": ", trimws(m)) else character(0)
  }))
  if (length(hits))
    stop("\n\nBare goal comparison found. Use at_goal()/above_goal() instead:\n  ",
         paste(hits, collapse = "\n  "), "\n", call. = FALSE)
  message("  Boundary guard: no bare goal comparison in ", length(files), " scripts.")
})

message("\n=== 6b. DOMAIN-ANALYSIS VERIFICATION ===")

full_fasting_df <- derive_vars(primary_raw) %>%
  filter(!is.na(WTSAF2YR), WTSAF2YR > 0) %>%
  compute_prevent_risk() %>%
  assign_goals() %>%
  mutate(
    n_cycles_pooled = length(primary_suffixes_used),
    wt_pooled = WTSAF2YR / n_cycles_pooled,
    ## The analytic cohort, expressed as a domain indicator rather than as a
    ## filter. Every condition here matches the primary_df filter above; if one
    ## is edited there it must be edited here, and the assertion below is what
    ## catches it if it is not.
    in_primary_cohort = as.integer(
      RIDAGEYR >= CONFIG$age_min & RIDAGEYR <= CONFIG$age_max &
      coalesce(pregnant == 0, FALSE) &
      coalesce(ascvd_hx == 0, FALSE) &
      !is.na(LBXAPB) & !is.na(ldl_primary) & !is.na(LBXTR) &
      !is.na(prevent_ascvd_10yr))
  )

full_fasting_design <- svydesign(
  id = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~wt_pooled,
  nest = TRUE, data = full_fasting_df
)
primary_domain_design <- subset(full_fasting_design, in_primary_cohort == 1)

## The domain must contain exactly the cohort. If these differ, one of the two
## filter specifications has drifted from the other and nothing below is
## comparable -- fail loudly rather than report a mismatch as a finding.
.n_dom <- sum(full_fasting_df$in_primary_cohort == 1, na.rm = TRUE)
message("  Fasting-subsample design frame: ", nrow(full_fasting_df), " records")
message("  Analytic domain within it:      ", .n_dom, " records")
message("  Cohort built by filtering:      ", nrow(primary_df), " records")
if (.n_dom != nrow(primary_df))
  stop("Domain (", .n_dom, ") and filtered cohort (", nrow(primary_df),
       ") differ. The two cohort specifications have diverged; reconcile them ",
       "before reading anything below.", call. = FALSE)

## Headline estimates, both ways. Point estimates should be identical; the
## question Dani raises is whether the standard errors are.
.cmp <- function(lbl, form, d1, d2) {
  e1 <- svyciprop(form, d1, method = "beta", na.rm = TRUE)
  e2 <- svyciprop(form, d2, method = "beta", na.rm = TRUE)
  tibble::tibble(
    quantity   = lbl,
    filtered   = as.numeric(e1),  filt_se = as.numeric(SE(e1)),
    domain     = as.numeric(e2),  dom_se  = as.numeric(SE(e2)),
    est_diff   = as.numeric(e2) - as.numeric(e1),
    se_ratio   = as.numeric(SE(e2)) / as.numeric(SE(e1))
  )
}

.ap_filt <- subset(primary_design,        has_apob_goal == 1)
.ap_dom  <- subset(primary_domain_design, has_apob_goal == 1)

domain_check <- dplyr::bind_rows(
  .cmp("Assigned an apoB goal", ~I(has_apob_goal == 1),
       primary_design, primary_domain_design),
  .cmp("Above apoB goal | assigned", ~I(above_apob_goal == 1), .ap_filt, .ap_dom),
  .cmp("Discordant | assigned",      ~I(discordant == 1),      .ap_filt, .ap_dom)
)
print(as.data.frame(domain_check))
write.csv(domain_check, "domain_design_check.csv", row.names = FALSE)

message("  Design df -- filtered: ", degf(primary_design),
        " | domain: ", degf(primary_domain_design))
message("  Largest absolute difference in a point estimate: ",
        sprintf("%.2e", max(abs(domain_check$est_diff))))
message("  SE ratio (domain / filtered), range: ",
        sprintf("%.4f to %.4f", min(domain_check$se_ratio), max(domain_check$se_ratio)))
message("  Interpretation: point estimates identical to numerical tolerance and ",
        "SE ratios at 1.00 mean the filtering shortcut cost nothing, and that ",
        "is the reportable verification. Any material difference means the ",
        "domain design is the one to use.")

if (isTRUE(CONFIG$use_domain_design)) {
  ## ROW-ALIGNMENT ASSERTION, and it is not ceremonial.
  ## Several downstream steps -- compare_blind_coding.R most importantly --
  ## pass vectors built from primary_df into update() on primary_design, which
  ## matches them BY POSITION. If the domain subset returned the cohort rows in
  ## a different order, or a different set, every one of those values would
  ## attach to the wrong participant and the result would look plausible.
  ## Compare on SEQN, which is the participant identifier.
  .dom_seqn <- primary_domain_design$variables$SEQN
  if (length(.dom_seqn) != nrow(primary_df) || !identical(.dom_seqn, primary_df$SEQN))
    stop("\n\nDomain design rows do not align with primary_df.\n",
         "  domain n = ", length(.dom_seqn), ", primary_df n = ", nrow(primary_df),
         "\n  Downstream steps match these two by position, so the run is ",
         "stopped rather than producing silently misaligned results.\n",
         "  Set CONFIG$use_domain_design = FALSE to proceed on the filtered ",
         "design while this is investigated.\n", call. = FALSE)
  message("  Row alignment verified on SEQN: ", length(.dom_seqn),
          " records, same order as primary_df.")

  primary_design <- primary_domain_design
  message("  CONFIG$use_domain_design is TRUE: primary_design is now the ",
          "domain-subset design. Every estimate below uses it.")
} else {
  message("  CONFIG$use_domain_design is FALSE: primary_design remains the ",
          "filtered design; the comparison above is reported for the record.")
}

## =============================================================================
## ANALYTIC SUBPOPULATION FOR ApoB OUTCOMES  (added 2026-08-15)
## =============================================================================
## Figure 1 of the guideline states an apoB goal only for specific patient
## populations -- in primary prevention, only when TG is 150-499 mg/dL. For
## everyone else the guideline specifies NO apoB goal, so apob_goal is NA and
## above_apob_goal / discordant are NA rather than FALSE.
##
## Leaving those NAs to propagate is what broke the previous run: svytotal
## returned NA for every population count, svyglm silently dropped 1,706
## records, and the mortality cohort collapsed from 8,916 to 2,314 without
## erroring. NA is not a safe default here.
##
## The fix is an explicit analytic subpopulation. Every apoB-dependent estimate
## runs on apob_design; LDL-C-only and descriptive estimates continue to run on
## primary_design. subset() preserves the full design object for variance
## estimation, which is the correct way to analyse a survey subpopulation --
## filtering the data frame first would be wrong.
##
## Both denominators are reported in the manuscript: the guideline-faithful
## restricted estimate as primary, and the applied-to-everyone estimate as a
## sensitivity, so a reader can see what the scope decision costs.
## =============================================================================
apob_design <- subset(primary_design, has_apob_goal == 1)

n_apob_goal <- sum(primary_df$has_apob_goal == 1, na.rm = TRUE)
apob_goal_share <- svyciprop(~I(has_apob_goal == 1), primary_design, method = "beta")
message("ApoB-goal subpopulation: ", n_apob_goal, " of ", nrow(primary_df),
        " primary-prevention records (", round(100 * n_apob_goal / nrow(primary_df), 1),
        "% unweighted). Weighted share of US adults for whom the 2026 guideline ",
        "specifies an apoB goal: ", round(100 * as.numeric(apob_goal_share), 1),
        "% (95% CI ", round(100 * confint(apob_goal_share)[1], 1), "-",
        round(100 * confint(apob_goal_share)[2], 1), "%).")
message("  --> This is a reportable finding in its own right: the guideline ",
        "assigns an apoB goal to a minority of adults, and every apoB estimate ",
        "below is conditioned on membership of that group.")


## ---- 6b. RQ0 -- proxy risk-tier vs real PREVENT risk-tier, head-to-head ---
## Direct before/after comparison: discordant_proxy (the original
## diabetes/CKD-based tiering, Section 5) vs discordant (real PREVENT-ASCVD-
## based tiering, Section 5b/5c). Both are computed on the identical set of
## people (same filters, both now additionally restricted to those PREVENT
## could score), so any difference reflects the tiering method itself, not
## sample composition.

rq0_proxy_prevalence   <- svyciprop(~I(discordant_proxy == 1), apob_design, method = "logit")
rq0_prevent_prevalence <- svyciprop(~I(discordant == 1),       apob_design, method = "logit")
rq0_comparison <- tibble(
  method   = c("diabetes/CKD proxy (original)", "PREVENT-ASCVD-based (this run)"),
  estimate = c(as.numeric(rq0_proxy_prevalence), as.numeric(rq0_prevent_prevalence)),
  ci_low   = c(attr(rq0_proxy_prevalence, "ci")[1], attr(rq0_prevent_prevalence, "ci")[1]),
  ci_high  = c(attr(rq0_proxy_prevalence, "ci")[2], attr(rq0_prevent_prevalence, "ci")[2])
)
print(rq0_comparison)

# How many people's risk tier -- and therefore their LDL-C/ApoB goal --
# actually changed when the proxy was replaced with real PREVENT scoring?
rq0_tier_change <- primary_df %>%
  mutate(tier_changed = risk_tier != risk_tier_proxy) %>%
  summarise(
    n = n(),
    n_tier_changed = sum(tier_changed, na.rm = TRUE),
    pct_tier_changed = mean(tier_changed, na.rm = TRUE)
  )
print(rq0_tier_change)

write.csv(rq0_comparison, "rq0_proxy_vs_prevent_comparison.csv", row.names = FALSE)
write.csv(rq0_tier_change, "rq0_tier_reassignment_summary.csv", row.names = FALSE)

## ---- 7. RQ1 -- overall discordance prevalence & population estimate -----------

rq1_prevalence <- svyciprop(~I(discordant == 1), apob_design, method = "beta")
rq1_summary <- tibble(
  estimate = as.numeric(rq1_prevalence),
  ci_low   = attr(rq1_prevalence, "ci")[1],
  ci_high  = attr(rq1_prevalence, "ci")[2]
)

rq1_population <- svytotal(~I(discordant == 1), apob_design)
# svytotal on a logical/derived indicator returns counts for TRUE/FALSE; extract TRUE row
print(rq1_summary)
print(rq1_population)

## ---- 7b. RQ1b -- overall ApoB goal non-attainment (not just the discordant slice) ----
## RQ1 above only covers the "discordant" cell (at LDL-C goal, ABOVE ApoB
## goal). That understates how many people have not met their ApoB goal at
## all: it excludes anyone who is also above their LDL-C goal (for whom ApoB
## status isn't adding new information, but who are still above ApoB goal).
## This section answers the broader question directly: what fraction of the
## population is above ApoB goal, period, regardless of LDL-C status -- i.e.
## discordant + above_both combined -- plus the full 2x2 breakdown so every
## combination of LDL-C/ApoB goal status is visible in one table.

rq1b_apob_prevalence <- svyciprop(~I(above_apob_goal == 1), apob_design, method = "beta")
rq1b_summary <- tibble(
  estimate = as.numeric(rq1b_apob_prevalence),
  ci_low   = attr(rq1b_apob_prevalence, "ci")[1],
  ci_high  = attr(rq1b_apob_prevalence, "ci")[2]
)
rq1b_population <- svytotal(~I(above_apob_goal == 1), apob_design)
print(rq1b_summary)
print(rq1b_population)

# Full 2x2 breakdown (LDL-C goal x ApoB goal).
#
# CITATION DISCIPLINE -- added 2026-08-27 (script review, finding 4). The
# "discordant" row of this table carries its OWN interval, computed across four
# cells by svyby, which is not the same interval as rq1_prevalence above even
# though the point estimate is identical to seven decimal places. An earlier
# manuscript draft quoted 1.08-2.10 from this table while
# rq1_overall_discordance_prevalence.csv said 1.15-2.20.
#
# rq1_overall_discordance_prevalence.csv IS THE CITABLE SOURCE for the
# discordance estimate and its confidence interval. This table is descriptive:
# cite it for the cell STRUCTURE and the population counts, not for the
# headline interval.
#
# Added to the design via update() rather than rebuilding primary_design from
# scratch, and to primary_df via mutate() for reference/export -- both point at
# the same four mutually exclusive, population-exhaustive cells:
#   concordant_at_goal     = at BOTH LDL-C and ApoB goal
#   discordant              = at LDL-C goal, ABOVE ApoB goal   (RQ1's outcome)
#   ldl_above_apob_at_goal  = ABOVE LDL-C goal, but AT ApoB goal ("reverse" discordance)
#   above_both               = ABOVE both LDL-C and ApoB goal
primary_df <- primary_df %>%
  mutate(goal_cell = case_when(
    at_ldl_goal == 1 & above_apob_goal == 0 ~ "concordant_at_goal",
    at_ldl_goal == 1 & above_apob_goal == 1 ~ "discordant",
    at_ldl_goal == 0 & above_apob_goal == 0 ~ "ldl_above_apob_at_goal",
    at_ldl_goal == 0 & above_apob_goal == 1 ~ "above_both"
  ))

primary_design <- update(
  primary_design,
  goal_cell = factor(
    case_when(
      at_ldl_goal == 1 & above_apob_goal == 0 ~ "concordant_at_goal",
      at_ldl_goal == 1 & above_apob_goal == 1 ~ "discordant",
      at_ldl_goal == 0 & above_apob_goal == 0 ~ "ldl_above_apob_at_goal",
      at_ldl_goal == 0 & above_apob_goal == 1 ~ "above_both"
    ),
    levels = c("concordant_at_goal", "discordant", "ldl_above_apob_at_goal", "above_both")
  )
)

## Refresh: subset() snapshots the design, so columns added by the update()
## above are absent from the subpopulation object unless it is rebuilt.
apob_design <- subset(primary_design, has_apob_goal == 1)

## ---------------------------------------------------------------------------
## KORN-GRAUBARD CELL INTERVALS -- added 2026-09-18 (round 6, "standardize the
## CI methodology").
##
## The paper's Methods state that proportions are reported with Korn-Graubard
## intervals. That was true of every HEADLINE proportion -- all computed with
## svyciprop(method = "beta") -- but NOT of the four-cell goal decompositions,
## which came from svymean() + confint() and are therefore symmetric linearized
## (Wald) intervals. The mismatch was visible in the paper: secondary-prevention
## discordance appeared as 2.94% (1.39-5.42) in the manuscript, from the
## svyciprop route, and as 2.94 (1.13-4.76) in the supplement, from this one.
## Identical point estimate, two intervals, one paper.
##
## This helper recomputes each cell of a decomposition as its own binary
## proportion with method = "beta", so every reported interval in the paper now
## comes from the same estimator and the Methods sentence is true as written.
##
## COST, stated plainly: the cells are no longer drawn from a single joint
## multinomial fit, so they are individually exact rather than jointly
## constrained. This is ordinary practice for reporting cell prevalences and is
## footnoted in the supplement. Point estimates and weighted totals are
## unchanged -- only the interval method moves.
## FAIL-SAFE, and it matters here. This helper is new code sitting a third of
## the way into a run that takes over an hour. If it threw, the whole run would
## die and the hour would be wasted. So it falls back to the previous
## svymean()/confint() computation instead of stopping -- but it does NOT do so
## silently, which would quietly reinstate the very inconsistency round 6 asked
## us to remove. Every row carries a ci_method column recording which estimator
## actually produced it, the message channel says so loudly, and the build
## scripts refuse to typeset a table whose ci_method is not korn-graubard.
kg_cell_summary <- function(design, cell_var = "goal_cell") {
  tot <- svytotal(as.formula(paste0("~", cell_var)), design)
  nm  <- sub(paste0("^", cell_var), "", names(coef(tot)))

  wald <- function() {
    mu <- svymean(as.formula(paste0("~", cell_var)), design)
    tibble(goal_cell  = sub(paste0("^", cell_var), "", names(coef(mu))),
           prevalence = as.numeric(coef(mu)),
           ci_low     = confint(mu)[, 1],
           ci_high    = confint(mu)[, 2],
           population = as.numeric(coef(tot)),
           ci_method  = "linearized-FALLBACK")
  }

  tryCatch({
    lvls <- levels(design$variables[[cell_var]])
    if (!length(lvls)) stop("no factor levels on ", cell_var)
    out <- dplyr::bind_rows(lapply(lvls, function(lv) {
      f  <- as.formula(sprintf("~I(%s == '%s')", cell_var, lv))
      ci <- svyciprop(f, design, method = "beta", na.rm = TRUE)
      tibble(goal_cell  = lv,
             prevalence = as.numeric(ci),
             ci_low     = attr(ci, "ci")[1],
             ci_high    = attr(ci, "ci")[2])
    }))
    out$population <- as.numeric(coef(tot))[match(out$goal_cell, nm)]
    out$ci_method  <- "korn-graubard"
    ## Cheap sanity checks on the result, not just on the absence of an error.
    stopifnot(nrow(out) == length(lvls),
              !any(is.na(out$population)),
              abs(sum(out$prevalence) - 1) < 1e-6,
              all(out$ci_low <= out$prevalence),
              all(out$ci_high >= out$prevalence))
    out
  }, error = function(e) {
    message("\n*** kg_cell_summary FAILED on ", cell_var, ": ",
            conditionMessage(e))
    message("*** Falling back to linearized intervals so the run continues.")
    message("*** The CSV records ci_method = linearized-FALLBACK and the ",
            "document builds will REFUSE to typeset it. Fix before rebuilding.\n")
    wald()
  })
}
## ---------------------------------------------------------------------------

# Built defensively via coef()/confint() rather than relying on svyby's
# column-naming, same lesson as the RQ3 positional-rename fix above --
# exact label formatting can vary by survey-package version.
rq1b_cell_prev  <- svymean(~goal_cell, apob_design)   # retained for the deff/RSE diagnostics below
rq1b_cell_tot   <- svytotal(~goal_cell, apob_design)
rq1b_cell_summary <- kg_cell_summary(apob_design, "goal_cell")
print(rq1b_cell_summary)

message("RQ1b: ", scales::percent(rq1b_summary$estimate, accuracy = 0.01),
        " of the primary-prevention population is above their ApoB goal overall ",
        "(discordant + above_both combined), vs. ",
        scales::percent(rq1_summary$estimate, accuracy = 0.01),
        " who are specifically discordant (at LDL-C goal but above ApoB goal, RQ1's narrower outcome).")

## ---- 8. RQ2 -- discordance among guideline-flagged subgroups ------------------

rq2_by_flagged <- svyby(~I(discordant == 1), ~apob_flagged_group, apob_design,
                         svyciprop, vartype = "ci", method = "logit")

rq2_by_condition <- bind_rows(
  svyby(~I(discordant == 1), ~diabetes,  apob_design, svyciprop, vartype = "ci", method = "logit") %>%
    as_tibble() %>% mutate(condition = "diabetes") %>% rename(level = diabetes),
  svyby(~I(discordant == 1), ~ckd,       apob_design, svyciprop, vartype = "ci", method = "logit") %>%
    as_tibble() %>% mutate(condition = "ckd") %>% rename(level = ckd),
  svyby(~I(discordant == 1), ~obesity,   apob_design, svyciprop, vartype = "ci", method = "logit") %>%
    as_tibble() %>% mutate(condition = "obesity") %>% rename(level = obesity),
  svyby(~I(discordant == 1), ~hypertrig, apob_design, svyciprop, vartype = "ci", method = "logit") %>%
    as_tibble() %>% mutate(condition = "hypertriglyceridemia") %>% rename(level = hypertrig)
)
print(rq2_by_condition)

## ---- 8b. RQ7 -- discordance by income tier and insurance status (health-
## equity extension) -----------------------------------------------------
## Same pattern as RQ2 above, but split by socioeconomic characteristics
## instead of clinical conditions. Tests whether the LDL-C/ApoB mismatch
## itself is more common among lower-income or uninsured adults -- not just
## the secondary-prevention treatment gap examined separately in RQ6e below
## -- which would extend this paper's practical implication from "who
## should get ApoB checked clinically" to "who needs proactive outreach to
## make sure that actually happens." income_tier/insured were computed once
## in derive_vars() and flow through into primary_df/primary_design
## unchanged; NA (refused/don't-know/missing) respondents are excluded from
## each comparison via subset() rather than folded into a tier, since
## silently assuming a tier for them would misrepresent the comparison.

primary_design <- update(
  primary_design,
  income_tier = factor(primary_df$income_tier,
                        levels = c("below_poverty", "low_income", "higher_income")),
  insured = primary_df$insured
)

## Refresh the apoB subpopulation: subset() snapshots the design, so any
## update() to primary_design after this point would otherwise leave
## apob_design without the newly added columns.
apob_design <- subset(primary_design, has_apob_goal == 1)

rq7_income_design <- subset(apob_design, !is.na(income_tier))
rq7_by_income <- svyby(~I(discordant == 1), ~income_tier, rq7_income_design,
                        svyciprop, vartype = "ci", method = "logit")
rq7_by_income_df <- as_tibble(rq7_by_income)
names(rq7_by_income_df)[2] <- "prevalence"
print(rq7_by_income_df)

## Was subset(primary_design, ...) while rq7_income_design used apob_design,
## so the income and insurance rows of RQ7 were computed on different
## denominators and were not comparable to each other.
rq7_insurance_design <- subset(apob_design, !is.na(insured))
rq7_by_insurance <- svyby(~I(discordant == 1), ~insured, rq7_insurance_design,
                           svyciprop, vartype = "ci", method = "logit")
rq7_by_insurance_df <- as_tibble(rq7_by_insurance)
names(rq7_by_insurance_df)[2] <- "prevalence"
print(rq7_by_insurance_df)

write.csv(rq7_by_income_df,    "rq7_discordance_by_income_tier.csv", row.names = FALSE)
write.csv(rq7_by_insurance_df, "rq7_discordance_by_insurance_status.csv", row.names = FALSE)

message("RQ7: discordance prevalence by income tier and insurance status written to ",
        "rq7_discordance_by_income_tier.csv / rq7_discordance_by_insurance_status.csv.")

## ---- 8c. RQ8 -- discordance by statin treatment status (primary prevention) ---
## on_statin already exists on primary_design and has already been used twice:
## as one input to the PREVENT risk equation (Section 5b), and, in the
## secondary-prevention extension, as the basis of the RQ6d treatment-
## adequacy breakdown. Neither of those asks the question this section asks:
## among primary-prevention adults specifically -- the population being
## newly risk-stratified, not the ones who already have ASCVD -- does simply
## being on a statin change how likely someone is to be discordant at all?
##
## Rationale: statins lower LDL-C mainly by upregulating hepatic LDL
## receptors, which increases clearance of LDL particles. In most people
## this lowers cholesterol mass and particle count roughly in parallel, so
## LDL-C and ApoB fall together and discordance risk shouldn't change much.
## But in the small, cholesterol-depleted LDL particle phenotype -- i.e.
## exactly the hypertriglyceridemia/diabetes subgroups already identified as
## the dominant predictors of discordance in RQ4 -- a statin can lower LDL-C
## substantially in mass terms while still leaving a comparatively high
## particle count in circulation, because clearing the same cholesterol mass
## requires clearing more of these smaller particles. That would mean
## treatment can make LDL-C look controlled before it actually resolves the
## particle-count problem -- part of the original clinical rationale for
## ApoB testing in the lipidology literature. A simple prevalence comparison
## (below) tests whether that pattern shows up at the population level; the
## interaction model that follows tests something sharper -- whether statin
## use specifically changes the SIZE of the hypertriglyceridemia effect
## (OR 12.0 in RQ4), rather than just shifting the overall discordance rate.
##
## CAVEAT (same limitation as the rest of this cross-sectional analysis):
## NHANES can only compare different people who happen to be on vs. off a
## statin at one point in time -- it cannot show that starting a statin
## caused discordance to appear in the same person. This is an association
## between treatment status and discordance, not a before/after causal test.

rq8_by_statin <- svyby(~I(discordant == 1), ~on_statin, apob_design,
                        svyciprop, vartype = "ci", method = "logit")
rq8_by_statin_df <- as_tibble(rq8_by_statin)
names(rq8_by_statin_df)[2] <- "prevalence"
print(rq8_by_statin_df)

rq8_by_statin_tot <- svyby(~I(discordant == 1), ~on_statin, apob_design, svytotal)
print(rq8_by_statin_tot)

# Interaction test: does statin use change the SIZE of the
# hypertriglyceridemia effect specifically, rather than just shifting the
# overall discordance rate uniformly? A significant hypertrig:on_statin
# interaction term (OR != 1) would indicate the two variables' joint effect
# is not simply the product of their separate effects -- i.e. statin
# treatment modifies how strongly hypertriglyceridemia predicts discordance.
rq8_interaction_model <- svyglm(
  discordant ~ hypertrig * on_statin + RIDAGEYR + factor(RIAGENDR) + diabetes + ckd,
  design = apob_design, family = quasibinomial()
)
rq8_interaction_or <- broom::tidy(rq8_interaction_model, exponentiate = TRUE, conf.int = TRUE)
print(rq8_interaction_or)

## RQ8b -- THE CONFOUNDING CHECK. This is the analysis that decides whether
## RQ8's headline result above means anything mechanistically.
##
## The problem with the simple comparison: to be classified "discordant" a
## person must FIRST be at LDL-C goal -- that is half the definition. Statins
## move people into the at-LDL-goal state. So a higher discordance rate among
## statin users could reflect nothing more than a larger share of them being
## ELIGIBLE to land in that cell at all, with untreated high-lipid adults
## disproportionately sitting in the "above both goals" cell where they
## cannot be counted as discordant by construction. That is a compositional
## artifact, not evidence that statins leave particle count behind.
##
## The fix: restrict to adults ALREADY at LDL-C goal, then compare statin-
## treated vs statin-naive within that group. Everyone in this comparison is
## equally eligible to be discordant, so the mechanical component is removed
## by design and only the biological question remains -- among people whose
## LDL-C says they're controlled, does being on a statin make it more likely
## that ApoB disagrees?
##
## Interpretation guide for the output:
##   - If the statin gap PERSISTS at similar magnitude here, the biological
##     reading strengthens considerably: statin therapy is associated with
##     residual particle burden that LDL-C is not capturing, which is a
##     direct argument for re-checking ApoB AFTER starting/titrating therapy
##     rather than only at baseline risk stratification.
##   - If the gap COLLAPSES toward null, RQ8's raw ratio was mostly
##     compositional (statins put people at LDL goal, which is the entry
##     condition for the discordant cell) and should be reported as such
##     rather than as a treatment effect.
## Either result is publishable and worth reporting; what is not defensible
## is reporting RQ8 alone without this check.
##
## OUTCOME ON THE FRIEDEWALD-BASED RUN, AND WHY IT MATTERS HERE: under
## NHANES's Friedewald LDL-C this restriction left a 2.7x gap intact (13.7%
## vs 5.1%, non-overlapping CIs), which looked like a strong biological
## finding. It was not. Re-running with the Sampson/NIH equation -- now the
## primary LDL-C definition throughout this script -- shrank it to roughly
## 1.6x with heavily overlapping confidence intervals, because Friedewald
## misclassified statin-treated adults into the at-goal group at nearly
## twice the rate of statin-naive adults (3.7% vs 1.9%). RQ8/RQ8b as
## executed below now run on Sampson and should be reported as a NULL or
## weak association, not as evidence that statins leave particle burden
## behind. Section 11b quantifies the difference explicitly.
##
## NOTE ON CONFOUNDING BY INDICATION: this restriction removes the mechanical
## eligibility artifact but NOT confounding by indication -- statin users
## were prescribed statins because their lipid/risk profile was worse to
## begin with, and residual baseline risk is not fully captured by the
## diabetes/CKD/hypertriglyceridemia adjustment set. Even a persisting gap
## here is an association, not a causal treatment effect.
##
## reclass_design is built independently here rather than reused from
## Section 11 because that section runs LATER in the script (Section 11 vs
## this Section 8c) -- referencing it now would fail with "object
## 'reclass_design' not found". Same subset definition, just constructed
## where it is first needed.

rq8b_at_ldl_goal_design <- subset(apob_design, at_ldl_goal == 1)

rq8b_by_statin <- tryCatch(
  svyby(~I(discordant == 1), ~on_statin, rq8b_at_ldl_goal_design,
        svyciprop, vartype = "ci", method = "logit"),
  error = function(e) {
    message("Could not compute discordance by statin status within the at-LDL-goal ",
            "subgroup (likely a cell too small): ", conditionMessage(e))
    NULL
  }
)

if (!is.null(rq8b_by_statin)) {
  rq8b_by_statin_df <- as_tibble(rq8b_by_statin)
  names(rq8b_by_statin_df)[2] <- "prevalence"
  print(rq8b_by_statin_df)

  # Unweighted cell counts, printed so the sample size behind each estimate is
  # visible alongside it -- the at-LDL-goal subgroup is only ~1,150 records
  # before splitting by statin use, so one of these cells may be thin enough
  # that the CI, not the point estimate, is the honest headline.
  rq8b_cell_counts <- primary_df %>%
    filter(at_ldl_goal == 1, !is.na(apob_goal)) %>%
    group_by(on_statin) %>%
    summarise(
      n = n(),
      n_discordant = sum(discordant, na.rm = TRUE),
      pct_discordant_unweighted = mean(discordant, na.rm = TRUE),
      .groups = "drop"
    )
  print(rq8b_cell_counts)

  write.csv(rq8b_by_statin_df,  "rq8b_discordance_by_statin_among_at_ldl_goal.csv", row.names = FALSE)
  write.csv(rq8b_cell_counts,   "rq8b_at_ldl_goal_statin_cell_counts.csv", row.names = FALSE)
}

write.csv(rq8_by_statin_df,   "rq8_discordance_by_statin_use.csv", row.names = FALSE)
write.csv(rq8_interaction_or, "rq8_hypertrig_statin_interaction_OR.csv", row.names = FALSE)

message("RQ8b: the confounding check. Compare rq8b_* against rq8_* directly. ",
        "If the statin gap persists within the at-LDL-goal subgroup, RQ8 reflects ",
        "residual particle burden on treatment (report as a biological finding, with ",
        "confounding-by-indication caveat). If it collapses, RQ8's 3.9x was mostly a ",
        "compositional artifact of statins moving people into the at-LDL-goal state ",
        "that is the entry condition for being counted discordant at all.")

message("RQ8: discordance prevalence by statin treatment status, and the ",
        "hypertriglyceridemia x statin-use interaction test, written to rq8_*.csv. ",
        "Read the hypertrig:on_statin interaction term as: does statin treatment ",
        "change how strongly hypertriglyceridemia predicts discordance, beyond each ",
        "variable's own separate effect?")

## =============================================================================
## 8d. RQ9 -- PERFORMANCE OF THE GUIDELINE'S ApoB TESTING TRIGGER
##
## THE QUESTION. Every other section of this script asks "who is above goal."
## This one asks something different and arguably more actionable: among adults
## the guideline does NOT recommend for statin therapy -- the group in which no
## one is monitoring lipids at all -- does the guideline's ApoB testing trigger
## successfully identify people with markedly elevated ApoB? And how many are
## missed?
##
## WHY THIS FRAMING. A published analysis (Circulation 2025, VLDbL + NHANES)
## evaluated the 2018 guideline's trigger -- triglycerides >200 mg/dL -- against
## the 2018 risk-enhancer threshold of ApoB >=130 mg/dL, and reported ~12.5%
## yield among trigger-positive vs ~0.6% among trigger-negative, with roughly
## 317,000 US adults missed. The 2026 guideline appears to have BROADENED the
## trigger from triglycerides alone to triglycerides plus diabetes plus low
## achieved LDL-C. This section quantifies whether that broadening reduces the
## missed group -- a guideline-supportive screening-performance question that
## has not been asked.
##
## THRESHOLD PROVENANCE (verified, and important):
##   - ApoB >=130 mg/dL is a 2018 ACC/AHA guideline RISK-ENHANCING FACTOR,
##     stated to correspond to LDL-C >=160 mg/dL. It is NOT a treatment goal.
##     It is retained here as the primary outcome purely for comparability with
##     the published analysis above.
##   - TG >200 mg/dL is the 2018 guideline's stated "relative indication for
##     apoB measurement."
##   - The 2026 guideline's triglyceride threshold is AMBIGUOUS across
##     accessible secondary sources: the AHA "top things to know" summary and
##     the JACC guideline-at-a-glance state >=150 mg/dL, while other
##     professional commentary states >200 mg/dL. The primary guideline text
##     was not obtainable. BOTH are therefore evaluated below so no conclusion
##     depends on resolving it. Settle this against the primary document before
##     reporting.
##
## EXPECTED STRUCTURAL FINDING -- NOT A BUG. Diabetes and CKD make a patient
## automatically statin-eligible under the guideline, so both are ~absent from
## the statin-ineligible denominator by construction. Adding diabetes to the
## trigger set therefore cannot improve detection in THIS population. That is
## itself a reportable result: the diabetes arm of the 2026 trigger is
## redundant among adults not recommended for statins, because those patients
## have already been captured upstream by eligibility. The code below verifies
## this empirically rather than assuming it.
## =============================================================================

trigger_df <- primary_df %>%
  mutate(
    # Statin eligibility per the guideline's own structure: automatic on
    # LDL-C >=190, diabetes, or CKD; otherwise determined by PREVENT risk.
    # This mirrors the framework used by Anderson et al. (JAMA 2026), who
    # reported 8.6% eligible independent of risk estimation and 68.0%
    # requiring risk-based assessment.
    # Uses the ELIGIBILITY cut (therapy "should be considered" at >=5%), not
    # the GOAL cut (stricter goals begin at >=10%). These are different
    # thresholds doing different jobs -- see CONFIG. Using the goal cut here
    # would badly understate how many adults the guideline recommends for
    # therapy.
    statin_eligible_2026 = as.integer(
      coalesce(ldl_primary >= 190, FALSE) |
      diabetes == 1 |
      ckd == 1 |
      coalesce(prevent_ascvd_10yr >= CONFIG$prevent_eligibility_cut, FALSE)
    ),
    # PRIMARY outcome: markedly elevated ApoB (2018 risk-enhancer threshold),
    # chosen for direct comparability with the published analysis.
    apob_markedly_high = as.integer(coalesce(LBXAPB >= CONFIG$apob_high_threshold, FALSE)),
    # Trigger components. Both candidate TG cuts are hard-coded here on
    # purpose: RQ9's whole point is to compare them, so it must not depend on
    # CONFIG$tg_threshold.
    tg_gt200         = as.integer(coalesce(LBXTR > 200, FALSE)),
    tg_ge150         = as.integer(coalesce(LBXTR >= 150, FALSE)),
    low_achieved_ldl = as.integer(coalesce(ldl_primary < 70, FALSE))
  ) %>%
  mutate(
    trigger_2018_tg200      = tg_gt200,
    trigger_tg150           = tg_ge150,
    trigger_tg150_dm        = as.integer(tg_ge150 == 1 | diabetes == 1),
    trigger_2026_full_tg150 = as.integer(tg_ge150 == 1 | diabetes == 1 | low_achieved_ldl == 1),
    trigger_2026_full_tg200 = as.integer(tg_gt200 == 1 | diabetes == 1 | low_achieved_ldl == 1)
  )

# The analytic population: NOT statin-eligible AND not currently taking a
# statin. This is the "no one is watching" group -- the only population in
# which a missed ApoB elevation goes entirely undetected.
untested_df <- trigger_df %>%
  filter(statin_eligible_2026 == 0, on_statin == 0)

message("RQ9 denominator: ", nrow(untested_df), " of ", nrow(trigger_df),
        " primary-prevention adults are neither statin-eligible under the 2026 ",
        "guideline nor currently taking a statin (",
        round(100 * nrow(untested_df) / nrow(trigger_df), 1), "%).")

# Verify the expected structural redundancy rather than assuming it.
message("RQ9 structural check -- within the statin-ineligible group: ",
        sum(untested_df$diabetes, na.rm = TRUE), " have diabetes and ",
        sum(untested_df$ckd, na.rm = TRUE), " have CKD. If these are zero or ",
        "near-zero, the diabetes arm of the 2026 trigger cannot add detection ",
        "in this population, because such patients are already statin-eligible.")

untested_design <- svydesign(
  id = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~wt_pooled,
  nest = TRUE, data = untested_df
)

# Evaluate one trigger: yield among trigger-positive and trigger-negative,
# sensitivity, number needed to test, and the weighted national count of
# missed cases (trigger-negative but markedly elevated ApoB).
evaluate_trigger <- function(trig_name, label, design, df) {
  d <- update(design, .trig = df[[trig_name]])

  yield_pos <- tryCatch(
    svyciprop(~I(apob_markedly_high == 1), subset(d, .trig == 1), method = "logit"),
    error = function(e) NULL)
  yield_neg <- tryCatch(
    svyciprop(~I(apob_markedly_high == 1), subset(d, .trig == 0), method = "logit"),
    error = function(e) NULL)

  # Weighted national count missed = trigger-negative AND markedly high ApoB
  missed_tot <- tryCatch({
    tt <- svytotal(~I(.trig == 0 & apob_markedly_high == 1), d)
    as.numeric(coef(tt)[2])
  }, error = function(e) NA_real_)

  caught_tot <- tryCatch({
    tt <- svytotal(~I(.trig == 1 & apob_markedly_high == 1), d)
    as.numeric(coef(tt)[2])
  }, error = function(e) NA_real_)

  tested_tot <- tryCatch({
    tt <- svytotal(~I(.trig == 1), d)
    as.numeric(coef(tt)[2])
  }, error = function(e) NA_real_)

  tibble(
    trigger              = label,
    n_trigger_pos        = sum(df[[trig_name]] == 1, na.rm = TRUE),
    n_trigger_neg        = sum(df[[trig_name]] == 0, na.rm = TRUE),
    yield_if_tested      = if (is.null(yield_pos)) NA_real_ else as.numeric(yield_pos),
    yield_pos_ci_low     = if (is.null(yield_pos)) NA_real_ else attr(yield_pos, "ci")[1],
    yield_pos_ci_high    = if (is.null(yield_pos)) NA_real_ else attr(yield_pos, "ci")[2],
    yield_if_not_tested  = if (is.null(yield_neg)) NA_real_ else as.numeric(yield_neg),
    yield_neg_ci_low     = if (is.null(yield_neg)) NA_real_ else attr(yield_neg, "ci")[1],
    yield_neg_ci_high    = if (is.null(yield_neg)) NA_real_ else attr(yield_neg, "ci")[2],
    pop_caught           = caught_tot,
    pop_missed           = missed_tot,
    sensitivity          = caught_tot / (caught_tot + missed_tot),
    pop_would_be_tested  = tested_tot,
    number_needed_to_test = ifelse(is.null(yield_pos) || as.numeric(yield_pos) == 0,
                                    NA_real_, 1 / as.numeric(yield_pos))
  )
}

rq9_triggers <- list(
  ## LABELS CORRECTED 2026-08-02. The TG >=150 variant was previously labelled
  ## "2026 full set", which is wrong: the verified 2026 ACC/AHA testing trigger
  ## uses TG >200 mg/dL (AHA "Top Things to Know", updated 13 March 2026,
  ## point 5). The >=150 rows are exploratory alternatives retained to show what
  ## a looser cut would buy; they are NOT the guideline criterion and must not
  ## be reported as such. The guideline-accurate row is flagged [GUIDELINE].
  ## LABELS CORRECTED AGAIN 2026-08-15, and this time from the guideline PDF.
  ## The previous labelling was exactly backwards: it flagged the TG >200 row
  ## as [GUIDELINE] and marked the >=150 rows "NOT guideline", on the authority
  ## of the AHA web summary. The guideline's own Top Ten Take-Home Messages,
  ## point 5, specifies TG >=150 mg/dL. The >=150 full set is the guideline
  ## criterion; >200 is the 2018 risk-enhancer value, retained as sensitivity.
  ##
  ## These five definitions use hard-coded tg_gt200 / tg_ge150 flags and are
  ## deliberately independent of CONFIG$tg_threshold, so that both cut points
  ## are always computed side by side regardless of the CONFIG setting.
  c("trigger_2018_tg200",      "Sensitivity: TG >200 alone (2018 risk-enhancer)"),
  c("trigger_tg150",           "TG >=150 alone"),
  c("trigger_tg150_dm",        "TG >=150 or diabetes"),
  c("trigger_2026_full_tg150", "[GUIDELINE] 2026 full set: TG >=150, diabetes, or LDL-C <70"),
  c("trigger_2026_full_tg200", "Sensitivity: 2026 set with TG cut at >200")
)

rq9_results <- bind_rows(lapply(rq9_triggers, function(x)
  evaluate_trigger(x[1], x[2], untested_design, untested_df)))

print(as.data.frame(rq9_results))

# Secondary outcome: repeat using tier-specific ApoB GOAL non-attainment
# instead of the 2018 risk-enhancer threshold, to show how the picture
# changes when the question shifts from "markedly elevated" to "above
# treatment target." Numbers will be substantially larger.
untested_df2 <- untested_df %>% mutate(apob_markedly_high = above_apob_goal)
untested_design2 <- svydesign(
  id = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~wt_pooled,
  nest = TRUE, data = untested_df2
)
rq9_results_goal <- bind_rows(lapply(rq9_triggers, function(x)
  evaluate_trigger(x[1], x[2], untested_design2, untested_df2)))
print(as.data.frame(rq9_results_goal))

write.csv(rq9_results,      "rq9_apob_trigger_performance_apob130.csv", row.names = FALSE)
write.csv(rq9_results_goal, "rq9_apob_trigger_performance_goalbased.csv", row.names = FALSE)

message("RQ9 complete. rq9_*_apob130.csv uses the 2018 risk-enhancer threshold ",
        "(ApoB >=130) and is directly comparable to the published Circulation ",
        "analysis; rq9_*_goalbased.csv substitutes tier-specific ApoB goal ",
        "non-attainment. Read 'pop_missed' as the weighted national count of ",
        "adults with elevated ApoB who would never be tested under that trigger ",
        "-- the direct analogue of the published ~317,000 figure. Compare the ",
        "2018 row against the 2026 rows to quantify what broadening the trigger buys.")

## =============================================================================
## 8e. RQ10 -- REPLICATION OF THE SELECTIVE-TESTING-YIELD ANALYSIS
##
## Independent reproduction of the cohort and methods used in a companion
## manuscript (Jayakumar A, "Selective Apolipoprotein B Testing Yield and
## Discordantly Elevated ApoB Among U.S. Adults Excluded From 2026 Guideline
## Statin Recommendations"). Target figures to reproduce, from that analysis:
##
##   Cohort         3,138 participants, ~57.09M weighted
##   No trigger     2,943 participants, 53.46M, 18 cases,
##                  0.59% (Korn-Graubard 95% CI 0.30-1.05), ~316,946 adults
##   Trigger        195 participants, 3.63M, 17 cases,
##                  12.51% (95% CI 6.68-20.75), ~454,746 adults
##   Difference     11.92 percentage points (95% CI 5.41-18.43)
##   Ratio          21.11 (95% CI 9.38-47.48); Rao-Scott F(1,94) = 108.17
##
## WHY THIS COHORT IS DIFFERENT FROM RQ9's, AND BETTER FOR THE QUESTION.
## RQ9 excluded only LDL-C >=190 (the automatic statin-eligibility threshold).
## This cohort additionally excludes LDL-C >=160. That matters enormously,
## because ApoB >=130 mg/dL corresponds to roughly LDL-C >=160: retaining the
## 160-189 band means retaining people whose high ApoB is CONCORDANT with a
## high LDL-C that already flags them on a standard panel. Those people are not
## "missed" in any clinically meaningful sense. Excluding them leaves only
## genuine discordance -- high particle count despite acceptable cholesterol
## mass -- which is the construct the question is actually about. RQ9's
## ~2.49 million "missed" figure is inflated for this reason and should not be
## reported; this section supersedes it.
##
## The exclusion of LDL-C <70 removes the other tail, where ApoB >=130 is
## near-impossible, so the retained band is 70-159 mg/dL.
##
## STATISTICAL METHODS mirror the source analysis rather than this script's
## defaults: Korn-Graubard (exact beta) confidence intervals instead of logit,
## appropriate for rare events; relative standard error, design effect, and
## effective sample size for NCHS presentation-standard assessment; Rao-Scott
## adjusted tests; and prevalence difference and ratio from survey-weighted
## identity- and log-link models.
## =============================================================================

rq10_source <- derive_vars(primary_raw) %>%
  compute_prevent_risk()

## =============================================================================
## COHORT FILTER: LDL-C EQUATION  (changed 2026-08-15)
## =============================================================================
## This cohort is defined by an LDL-C window (70 to <160). Which equation draws
## that window decides who is in the analysis, so it is not a technicality.
##
## It previously used LBDLDL (Friedewald), to match the source analysis being
## replicated. That justification no longer holds: the triglyceride trigger moved
## from >200 to >=150 earlier today, which already broke comparability with that
## analysis -- verify_run.R now SKIPS the replication check outright, reporting
## "targets assume TG >200; this run uses >=150 -- comparison not meaningful."
## Keeping Friedewald therefore preserved nothing while contradicting the
## guideline's Class 1 recommendation and the rest of this script.
##
## The filter now reads ldl_primary (Sampson under current CONFIG). Both arms
## are built and reported: the Friedewald-defined cohort runs as a labelled
## sensitivity in Section 10c so the finding can be shown to hold either way.
## With only ~32 cases across both trigger arms, small shifts in cohort
## membership can move the estimates more than their size suggests.
##
## NOTE: the OUTCOME (apoB >=130) is a directly measured value and is unaffected
## by equation choice. Only who enters the cohort changes.
build_rq10_cohort <- function(src, ldl_var, label) {
  ldl <- src[[ldl_var]]
  out <- src %>%
    filter(
      RIDAGEYR >= CONFIG$age_min, RIDAGEYR <= CONFIG$age_max,
      pregnant == 0,
      WTSAF2YR > 0,
      !is.na(LBXAPB), !is.na(.data[[ldl_var]]), !is.na(LBXTR),
      ascvd_hx == 0,                       # no clinical ASCVD
      diabetes == 0,                       # diabetes is an independent indication
      ckd == 0,                            # CKD is an independent indication
      .data[[ldl_var]] >= 70, .data[[ldl_var]] < 160,   # the defining restriction
      on_statin == 0, on_nonstatin_lipid_tx == 0,       # untreated
      !is.na(prevent_ascvd_10yr),
      prevent_ascvd_10yr < 0.03            # below the "can be considered" threshold
    )
  ## BUG FIX 2026-08-15. The 30-year PREVENT criterion was previously applied
  ## only to the primary arm, downstream of this function. The sensitivity arm
  ## called the function and went straight to mutate(), so it never received the
  ## filter -- producing 3,649 vs 2,884 and making the two arms incomparable.
  ## The difference looked like an equation effect and was not. Applying it here
  ## guarantees both arms are built identically in every respect except which
  ## equation draws the LDL-C window, which is the only thing being tested.
  n_pre <- nrow(out)
  if (all(is.na(out$prevent_ascvd_30yr))) {
    message("RQ10 [", label, "]: 30-year PREVENT unavailable -- 10-year ",
            "criterion applied alone. DEVIATION from the replicated cohort ",
            "definition; this cohort will be LARGER than the target of 3,138.")
  } else {
    out <- filter(out, is.na(prevent_ascvd_30yr) | prevent_ascvd_30yr < 0.10)
  }
  message("RQ10 cohort [", label, "]: ", n_pre, " before the 30-year criterion, ",
          nrow(out), " after.")
  out
}

## PRIMARY_EQ_LABEL is defined further down (Section 11b), so the label is
## rebuilt locally here rather than referenced before it exists.
RQ10_PRIMARY_LABEL <- if (CONFIG$primary_ldl_equation == "sampson")
  "Sampson/NIH equation 2" else "Friedewald"
RQ10_ALT_LABEL <- if (CONFIG$primary_ldl_equation == "sampson")
  "Friedewald" else "Sampson/NIH equation 2"

## The 30-year criterion now lives inside build_rq10_cohort(), so both arms
## receive it. Do not reintroduce it here.
rq10_cohort <- build_rq10_cohort(rq10_source, "ldl_primary",
                                 paste0(RQ10_PRIMARY_LABEL, " -- PRIMARY")) %>%
  mutate(
    n_cycles_pooled = length(primary_suffixes_used),
    wt_pooled = WTSAF2YR / n_cycles_pooled,
    apob_high_130 = as.integer(LBXAPB >= 130),
    ## Reads CONFIG rather than hard-coding, so RQ10 cannot drift from RQ9.
    ## Previously fixed at >200; the guideline Take-Home text specifies >=150.
    ## Set CONFIG$tg_threshold = 200 to reproduce the source analysis exactly.
    tg_trigger    = as.integer(LBXTR >= CONFIG$tg_threshold)
  )

message("RQ10 cohort: ", nrow(rq10_cohort), " participants (target 3,138). ",
        "Trigger-positive: ", sum(rq10_cohort$tg_trigger == 1),
        " (target 195); trigger-negative: ", sum(rq10_cohort$tg_trigger == 0),
        " (target 2,943).")
message("RQ10 ApoB >=130 cases: ", sum(rq10_cohort$apob_high_130), " total (target 35); ",
        sum(rq10_cohort$apob_high_130 == 1 & rq10_cohort$tg_trigger == 1), " trigger-positive (target 17); ",
        sum(rq10_cohort$apob_high_130 == 1 & rq10_cohort$tg_trigger == 0), " trigger-negative (target 18).")

rq10_design <- svydesign(id = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~wt_pooled,
                          nest = TRUE, data = rq10_cohort)
message("RQ10 survey degrees of freedom: ", degf(rq10_design), " (target 94).")

## Group-specific yield with Korn-Graubard intervals and NCHS reliability metrics.
rq10_group_stats <- function(design, df, trig_value, label) {
  d  <- subset(design, tg_trigger == trig_value)
  sub <- df %>% filter(tg_trigger == trig_value)

  # method = "beta" is the Korn-Graubard exact interval
  ci <- svyciprop(~I(apob_high_130 == 1), d, method = "beta")
  est <- as.numeric(ci)
  se  <- as.numeric(SE(ci))

  # Design effect and effective sample size from the same estimate
  m    <- svymean(~I(apob_high_130 == 1), d, deff = TRUE)
  deff <- as.numeric(deff(m))[2]
  n_un <- nrow(sub)
  n_eff <- if (is.na(deff) || deff <= 0) NA_real_ else n_un / deff

  pop_group <- as.numeric(coef(svytotal(~I(tg_trigger == trig_value), design))[2])
  pop_cases <- as.numeric(coef(svytotal(~I(apob_high_130 == 1), d))[2])

  tibble(
    group = label,
    n_unweighted = n_un,
    n_cases = sum(sub$apob_high_130, na.rm = TRUE),
    pop_millions = pop_group / 1e6,
    prevalence = est,
    ci_low = attr(ci, "ci")[1],
    ci_high = attr(ci, "ci")[2],
    cases_millions = pop_cases / 1e6,
    rse_pct = 100 * se / est,
    design_effect = deff,
    effective_n = n_eff
  )
}

rq10_results <- bind_rows(
  ## Labels read CONFIG so they cannot contradict the trigger definition, which
  ## they did when tg_trigger moved to >=150 while these strings still said 200.
  rq10_group_stats(rq10_design, rq10_cohort, 0,
                   paste0("No trigger: TG <", CONFIG$tg_threshold, " mg/dL")),
  rq10_group_stats(rq10_design, rq10_cohort, 1,
                   paste0("Trigger: TG >=", CONFIG$tg_threshold, " mg/dL"))
)
print(as.data.frame(rq10_results))

## Absolute prevalence difference -- survey-weighted identity-link model
rq10_pd <- tryCatch({
  m <- svyglm(apob_high_130 ~ tg_trigger, design = rq10_design, family = gaussian())
  ct <- summary(m)$coefficients
  ci <- confint(m)
  tibble(metric = "prevalence difference (percentage points)",
         estimate = 100 * ct["tg_trigger", "Estimate"],
         ci_low = 100 * ci["tg_trigger", 1], ci_high = 100 * ci["tg_trigger", 2])
}, error = function(e) { message("PD model failed: ", conditionMessage(e)); NULL })
if (!is.null(rq10_pd)) print(rq10_pd)

## Prevalence ratio -- survey-weighted quasipoisson with log link
rq10_pr <- tryCatch({
  m <- svyglm(apob_high_130 ~ tg_trigger, design = rq10_design, family = quasipoisson())
  ct <- summary(m)$coefficients
  ci <- confint(m)
  tibble(metric = "prevalence ratio",
         estimate = exp(ct["tg_trigger", "Estimate"]),
         ci_low = exp(ci["tg_trigger", 1]), ci_high = exp(ci["tg_trigger", 2]))
}, error = function(e) { message("PR model failed: ", conditionMessage(e)); NULL })
if (!is.null(rq10_pr)) print(rq10_pr)

## Rao-Scott adjusted test
rq10_rs <- tryCatch(svychisq(~apob_high_130 + tg_trigger, rq10_design, statistic = "F"),
                     error = function(e) { message("Rao-Scott failed: ", conditionMessage(e)); NULL })
if (!is.null(rq10_rs)) print(rq10_rs)

## The derived "missed fraction" -- the share of all ApoB >=130 cases sitting in
## the no-trigger group. Computed for completeness and DELIBERATELY FLAGGED AS
## UNRELIABLE: in the source analysis this estimate had an effective sample size
## of 26 and subpopulation degrees of freedom of 4, and was excluded from the
## headline for that reason. RQ9's "sensitivity" statistic is the same construct
## and was reported there without this scrutiny -- an error this section corrects.
rq10_cases_design <- subset(rq10_design, apob_high_130 == 1)
rq10_missed <- tryCatch(svyciprop(~I(tg_trigger == 0), rq10_cases_design, method = "beta"),
                         error = function(e) { message("Missed-fraction failed: ", conditionMessage(e)); NULL })
if (!is.null(rq10_missed)) {
  message("Missed fraction (share of ApoB>=130 cases lacking the trigger): ",
          round(100 * as.numeric(rq10_missed), 1), "% (95% CI ",
          round(100 * attr(rq10_missed, "ci")[1], 1), "-",
          round(100 * attr(rq10_missed, "ci")[2], 1), "%), on ",
          sum(rq10_cohort$apob_high_130), " unweighted cases and ",
          degf(rq10_cases_design), " subpopulation degrees of freedom. ",
          "DO NOT HEADLINE -- reported only for transparency.")
}

write.csv(rq10_results, "rq10_selective_testing_yield.csv", row.names = FALSE)
if (!is.null(rq10_pd)) write.csv(rq10_pd, "rq10_prevalence_difference.csv", row.names = FALSE)
if (!is.null(rq10_pr)) write.csv(rq10_pr, "rq10_prevalence_ratio.csv", row.names = FALSE)

## =============================================================================
## Section 10c -- RQ10 SENSITIVITY: cohort defined by the OTHER LDL-C equation
## =============================================================================
## The cohort window (LDL-C 70 to <160) decides who is analysed, so the equation
## drawing it is a real analytic choice, not a formatting one. With ~32 cases
## across both trigger arms, a shift of a few dozen participants can move the
## estimates appreciably. This arm rebuilds the whole analysis on the comparison
## equation so the finding can be shown to hold -- or not -- either way.
##
## Report as: "the contrast was X under [primary] and Y under [alternative]".
## If the two disagree materially, the contrast is filter-dependent and must be
## reported as such rather than as a property of the population.
{
  alt_cohort <- build_rq10_cohort(rq10_source, "ldl_alt",
                                  paste0(RQ10_ALT_LABEL, " -- SENSITIVITY")) %>%
    mutate(
      n_cycles_pooled = length(primary_suffixes_used),
      wt_pooled       = WTSAF2YR / n_cycles_pooled,
      apob_high_130   = as.integer(LBXAPB >= 130),
      tg_trigger      = as.integer(LBXTR >= CONFIG$tg_threshold)
    )

  alt_design <- svydesign(id = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~wt_pooled,
                          nest = TRUE, data = alt_cohort)

  alt_tab <- alt_cohort %>%
    group_by(tg_trigger) %>%
    summarise(n = n(), cases = sum(apob_high_130), .groups = "drop")

  alt_pd <- tryCatch({
    m <- svyglm(apob_high_130 ~ tg_trigger, design = alt_design, family = gaussian())
    c(est = 100 * coef(m)[2], lo = 100 * confint(m)[2, 1], hi = 100 * confint(m)[2, 2])
  }, error = function(e) rep(NA_real_, 3))

  alt_pr <- tryCatch({
    m <- svyglm(apob_high_130 ~ tg_trigger, design = alt_design, family = quasipoisson())
    c(est = exp(coef(m)[2]), lo = exp(confint(m)[2, 1]), hi = exp(confint(m)[2, 2]))
  }, error = function(e) rep(NA_real_, 3))

  alt_chi <- tryCatch(
    svychisq(~apob_high_130 + tg_trigger, alt_design, statistic = "F"),
    error = function(e) NULL)

  cat("\n--- RQ10 SENSITIVITY: cohort defined by ", RQ10_ALT_LABEL, " ---\n", sep = "")
  print(alt_tab)
  cat(sprintf("Cohort n: %d (primary arm: %d)\n", nrow(alt_cohort), nrow(rq10_cohort)))
  cat(sprintf("Prevalence difference: %.2f pp (%.2f-%.2f)\n",
              alt_pd[1], alt_pd[2], alt_pd[3]))
  cat(sprintf("Prevalence ratio:      %.2f (%.2f-%.2f)\n",
              alt_pr[1], alt_pr[2], alt_pr[3]))
  if (!is.null(alt_chi))
    cat(sprintf("Rao-Scott F = %.1f, df = %.0f, p = %.3g\n",
                alt_chi$statistic, alt_chi$parameter[2], alt_chi$p.value))

  rq10_equation_sensitivity <- tibble(
    ldl_equation_defining_cohort = c(paste0(RQ10_PRIMARY_LABEL, " (PRIMARY)"),
                                     paste0(RQ10_ALT_LABEL, " (sensitivity)")),
    cohort_n = c(nrow(rq10_cohort), nrow(alt_cohort)),
    cases    = c(sum(rq10_cohort$apob_high_130), sum(alt_cohort$apob_high_130))
  )
  print(rq10_equation_sensitivity)
  write.csv(rq10_equation_sensitivity,
            "rq10_cohort_by_ldl_equation.csv", row.names = FALSE)

  message("Section 10c: if the two arms agree, state it in one line as a ",
          "robustness result. If they diverge, the contrast is dependent on ",
          "which equation drew the cohort window and must be reported that way.")
}

message("RQ10 complete. Compare the printed cohort counts and yields against the ",
        "targets in the section header. Discrepancies most likely arise from the ",
        "30-year PREVENT criterion, the exact definition of 'other independent ",
        "treatment indication', or the age range -- the source analysis reports a ",
        "mean age near 40.8 years, which may indicate a lower age floor than the ",
        "30 used here.")

## ---- 9. RQ3 -- secular trend, 2005-2016 (6 cycles) -----------------------------

rq3_trend <- svyby(~I(discordant == 1), ~cycle, trend_design, svyciprop, vartype = "ci", method = "logit")
rq3_trend_df <- as_tibble(rq3_trend)
# Positional rename: svyby names the 2nd column after the formula (e.g.
# "I(discordant == 1)"), which is awkward to reference -- rename by position
# instead of hardcoding that string, since exact formatting can vary by
# survey-package version.
names(rq3_trend_df)[2] <- "prevalence"
rq3_trend_df <- rq3_trend_df %>% mutate(cycle = factor(cycle, levels = trend_suffixes))
print(rq3_trend_df)

# Mean ApoB by cycle, weighted
rq3_mean_apob <- svyby(~LBXAPB, ~cycle, trend_design, svymean, vartype = "ci")
print(rq3_mean_apob)

trend_plot <- ggplot(rq3_trend_df, aes(x = cycle, y = prevalence, group = 1)) +
  geom_line() + geom_point() +
  geom_errorbar(aes(ymin = ci_l, ymax = ci_u), width = 0.15) +
  scale_y_continuous(labels = scales::percent) +
  labs(title = "Weighted LDL-C/ApoB Discordance Prevalence, NHANES 2005-2016",
       x = "NHANES Cycle", y = "Discordance Prevalence (95% CI)") +
  theme_minimal()
ggsave("apob_discordance_trend.png", trend_plot, width = 7, height = 4.5, dpi = 200)

## ---- 10. RQ4 -- predictors of discordance (weighted logistic regression) ------

## Uses the harmonized race_eth variable (see derive_vars). Previously this
## model referenced factor(RIDRETH3), which silently restricted the regression
## to 2011-2016 while every other section used the full pooled cohort.
rq4_model <- svyglm(
  discordant ~ RIDAGEYR + factor(RIAGENDR) + race_eth +
    BMXBMI + diabetes + ckd + hypertrig,
  design = apob_design, family = quasibinomial()
)
rq4_or_table <- broom::tidy(rq4_model, exponentiate = TRUE, conf.int = TRUE)
print(rq4_or_table)

## Verify the model actually used the whole cohort. If n_used is materially
## below nrow(primary_df), a predictor still carries missingness and the
## regression is describing a subset -- exactly the failure this fix addresses.
message("RQ4 model observations used: ", nobs(rq4_model),
        " of ", nrow(primary_df), " in the primary cohort.")
print(primary_df %>% count(race_eth))

## ---- 11. Reclassification: LDL-C-only tier vs. ApoB-inclusive tier ------------
## Among adults at LDL-C goal, what fraction get "reclassified upward" once
## ApoB is factored in (i.e., are the discordant group)? This is the direct
## analogue of the source paper's "newly eligible" comparison.

reclass_summary <- primary_df %>%
  ## See note below: n and pct must share a denominator.
  filter(at_ldl_goal == 1, !is.na(apob_goal)) %>%
  summarise(
    ## Restricted to !is.na(apob_goal): without it, n() counts everyone at
    ## LDL-C goal while mean(discordant, na.rm=TRUE) silently uses only
    ## those with a guideline-specified apoB goal -- so n and pct in the
    ## same row had different denominators (1,956 vs 250).
    n = n(),
    reclassified_by_apob = sum(discordant, na.rm = TRUE),
    pct_reclassified = mean(discordant, na.rm = TRUE)
  )
print(reclass_summary)

reclass_design <- subset(apob_design, at_ldl_goal == 1)
reclass_weighted <- svyciprop(~I(discordant == 1), reclass_design, method = "beta")
print(reclass_weighted)

## ---------------------------------------------------------------------------
## EXPORT ADDED 2026-08-16. This is the paper's principal clinical estimate --
## among adults who have reached their assigned LDL-C goal, the share still
## above their apoB goal -- and it had NO primary output file. It was
## recoverable only from the console, from the reference row of
## sens_reclassification_by_ldl_equation.csv (a file named for a sensitivity
## analysis), from rev1_precedence_headlines.csv, and from the verification
## log. Every other headline result exports to its own rq* file; this one did
## not. Exported here so it can be cited from a primary results file.
##
## Note the two proportions are NOT interchangeable:
##   pct_unweighted = 48/351 = 13.7%  -- raw sample proportion
##   estimate       = 14.5%           -- survey-weighted population estimate
## Report the weighted one.
## ---------------------------------------------------------------------------
## -----------------------------------------------------------------------------
## NCHS PRESENTATION STANDARD FOR PROPORTIONS
##
## CORRECTED 2026-09-05. The earlier implementation flagged a proportion as
## unreliable when its relative standard error exceeded 30%. That is NOT the
## NCHS rule for proportions. NCHS Data Presentation Standards for Proportions
## (Vital and Health Statistics, Series 2, No. 175) explicitly declines to use
## the RSE for proportions -- "when dividing the SE by very small proportions,
## the RSE can be too conservative, and when dividing the SE by very large
## proportions, the RSE can be too liberal" -- and instead applies, in order:
##
##   1. denominator n >= 30 AND effective n >= 30
##   2. degrees of freedom >= 8
##   3. Korn-Graubard interval width:
##        absolute width <= 0.05                     -> present
##        absolute width >= 0.30                     -> suppress
##        0.05 < absolute width < 0.30 and
##            relative width (width / estimate) > 130% -> suppress
##                                     otherwise      -> present
##   4. numerator > 0 and numerator != denominator
##
## The 30% RSE criterion is retained by NCHS only for means, not proportions.
## RSE is still computed and reported below because it is informative, but it
## no longer drives the flag.
## -----------------------------------------------------------------------------
nchs_present <- function(est, ci_low, ci_high, n, effective_n, df = 94,
                         n_events = NA_integer_) {
  if (is.na(est) || is.na(ci_low) || is.na(ci_high)) return(NA)
  if (n < 30 || is.na(effective_n) || effective_n < 30) return(FALSE)
  if (df < 8) return(FALSE)
  if (!is.na(n_events) && (n_events == 0 || n_events == n)) return(FALSE)
  abs_w <- ci_high - ci_low
  if (abs_w <= 0.05) return(TRUE)
  if (abs_w >= 0.30) return(FALSE)
  rel_w <- abs_w / est
  rel_w <= 1.30
}

reclass_m    <- svymean(~I(discordant == 1), reclass_design, deff = TRUE)
reclass_rse  <- 100 * SE(reclass_m)[2] / coef(reclass_m)[2]
reclass_deff <- as.numeric(deff(reclass_m))[2]

rq1c_reclassification <- tibble(
  population        = "primary_prevention_at_ldl_goal",
  n_unweighted      = reclass_summary$n,
  n_reclassified    = reclass_summary$reclassified_by_apob,
  pct_unweighted    = reclass_summary$pct_reclassified,
  estimate          = as.numeric(reclass_weighted),
  ci_low            = attr(reclass_weighted, "ci")[1],
  ci_high           = attr(reclass_weighted, "ci")[2],
  rse_pct           = round(reclass_rse, 1),
  design_effect     = round(reclass_deff, 2),
  effective_n       = round(reclass_summary$n / reclass_deff),
  nchs_unreliable   = !nchs_present(as.numeric(reclass_weighted),
                                    attr(reclass_weighted, "ci")[1],
                                    attr(reclass_weighted, "ci")[2],
                                    reclass_summary$n,
                                    reclass_summary$n / reclass_deff)
)
write.csv(rq1c_reclassification, "rq1c_reclassification_at_ldl_goal.csv",
          row.names = FALSE)
print(as.data.frame(rq1c_reclassification))
message("RQ1c: reclassification among adults at LDL-C goal written to ",
        "rq1c_reclassification_at_ldl_goal.csv. This is the primary source ",
        "for that estimate -- do not cite it from the sens_* or rev* files.")

## =============================================================================
## 11b. SENSITIVITY ANALYSIS -- HOW MUCH OF "DISCORDANCE" WAS THE LDL-C EQUATION?
##
## ROLE REVERSAL NOTE. In the previous version of this script this section ran
## the OPPOSITE way round: Friedewald was primary and Sampson was the
## sensitivity check. That check found Friedewald materially inflates every
## discordance estimate here, so Sampson is now the primary definition
## (assign_goals(), Section 5c) and this section reports what Friedewald --
## still the NHANES default, and still what most US labs reported during
## 2005-2016 -- would have produced instead. The comparison is retained rather
## than deleted because the DIFFERENCE between the two is itself a finding
## worth publishing.
##
## THE UNDERLYING PROBLEM. Friedewald assumes a FIXED triglyceride-to-VLDL-C
## ratio of 5:1. That assumption degrades in two specific conditions: low
## LDL-C, and elevated triglycerides. Both describe the exact subgroups this
## analysis cares most about -- adults at LDL-C goal (low by definition) and
## adults with hypertriglyceridemia (the dominant discordance predictor,
## OR 12-13 throughout). Where Friedewald underestimates true LDL-C, people
## are misclassified INTO the at-goal group; their ApoB is then correctly
## flagged as above goal, and they get counted as "discordant" when the actual
## error was the LDL-C number.
##
## WHAT THE FIRST RUN OF THIS SECTION FOUND (the reason for the switch):
##   - Friedewald was lower than Sampson in 88% of records, mean bias ~2.2
##     mg/dL -- small in absolute terms, but concentrated at the thresholds.
##   - 104 adults (2.3%) changed LDL-C goal status; 102 of those moved OUT of
##     goal under Sampson, only 2 moved in. Almost entirely one-directional.
##   - The misclassification was DIFFERENTIAL BY TREATMENT: 3.7% of
##     statin-treated adults changed status vs 1.9% of statin-naive.
##   - Overall discordance: 1.86% (Friedewald) -> 1.00% (Sampson).
##   - At-goal reclassification share: 7.16% -> 4.25%.
##   - The RQ8b statin gap collapsed from 2.7x with separated confidence
##     intervals to ~1.6x with heavily overlapping ones -- i.e. the statin
##     finding did NOT survive, and must not be reported as a positive result.
##
## GUIDELINE CONCORDANCE. The 2026 dyslipidemia guideline whose goals this
## script applies preferentially recommends Martin-Hopkins over Friedewald.
## Sampson is much closer to Martin-Hopkins than Friedewald is, so the switch
## also resolves an internal inconsistency: applying the guideline's
## thresholds while using the equation the guideline moved away from.
##
## WHY NOT MARTIN-HOPKINS ITSELF. Martin-Hopkins replaces the fixed /5 divisor
## with a value looked up from a published 180-cell table keyed on non-HDL-C
## and triglyceride strata (ratios spanning roughly 3.1 to 11.9). Those 180
## cell values were NOT independently verified while writing this script, and
## approximating them would silently corrupt every downstream estimate -- far
## worse than omitting the comparison. To add it properly: obtain the table
## from the supplementary material of Martin SS et al., JAMA 2013;310(19):
## 2061-2068, encode it as a lookup, and mirror the sampson_ldl() pattern.
## NHANES itself only began publishing a Martin-Hopkins variable (LBDLDLM)
## from the 2017-2018 cycle onward -- after ApoB measurement had already
## stopped -- so it cannot simply be read from the data for these cycles.
##
## SCOPE LIMITATION. The analytic cohort required non-missing LBDLDL, and
## NHANES blanks that variable when TG > 400 mg/dL. So the most severely
## hypertriglyceridemic adults were excluded upstream and are absent from both
## arms, even though Sampson stays valid to roughly TG 800. This biases the
## comparison toward finding LESS equation-driven difference than truly exists
## in the full population -- the people where Friedewald fails worst are the
## ones missing. Switching the cohort filter to !is.na(ldl_sampson) would
## recover them, at the cost of no longer being sample-comparable with every
## earlier run of this script; deliberately not done here.
## =============================================================================

# at_ldl_goal_alt / discordant_alt / ldl_goal_status_change are assigned in
# assign_goals() alongside their primary counterparts, so both arms exist on
# primary_df and primary_design already -- no separate design object needed.
# "primary" = CONFIG$primary_ldl_equation; "alt" = the other one.

PRIMARY_EQ_LABEL <- if (CONFIG$primary_ldl_equation == "sampson")
  "Sampson/NIH equation 2" else "Friedewald (NHANES default)"
ALT_EQ_LABEL <- if (CONFIG$primary_ldl_equation == "sampson")
  "Friedewald (NHANES default)" else "Sampson/NIH equation 2"
message("Section 11b: primary LDL-C equation is ", PRIMARY_EQ_LABEL,
        "; comparison arm is ", ALT_EQ_LABEL, ".")

# How far apart are the two equations, and in which direction?
# ldl_diff_sampson_minus_fried is always Sampson minus Friedewald regardless of
# which is primary, so its sign is unambiguous: positive means Friedewald read
# lower than Sampson for that person.
sens_ldl_agreement <- primary_df %>%
  summarise(
    n                    = n(),
    mean_friedewald      = mean(LBDLDL, na.rm = TRUE),
    mean_sampson         = mean(ldl_sampson, na.rm = TRUE),
    mean_diff_sampson_minus_fried = mean(ldl_diff_sampson_minus_fried, na.rm = TRUE),
    median_diff          = median(ldl_diff_sampson_minus_fried, na.rm = TRUE),
    pct_friedewald_lower = mean(ldl_diff_sampson_minus_fried > 0, na.rm = TRUE),
    correlation          = cor(LBDLDL, ldl_sampson, use = "complete.obs")
  )
print(sens_ldl_agreement)

# How many adults change LDL-C goal status between equations? This is the entire
# mechanism by which the discordance estimate can move -- if almost nobody
# changes status, the equation cannot be driving the result.
sens_status_change <- primary_df %>%
  count(ldl_goal_status_change) %>%
  mutate(pct = n / sum(n))
print(sens_status_change)

# Same question split by statin use. THIS IS THE KEY DIAGNOSTIC: if
# statin-treated adults disproportionately change status, then any statin-
# stratified discordance comparison built on Friedewald was contaminated by
# differential misclassification rather than measuring biology.
sens_status_change_by_statin <- primary_df %>%
  group_by(on_statin) %>%
  count(ldl_goal_status_change) %>%
  mutate(pct = n / sum(n)) %>%
  ungroup()
print(sens_status_change_by_statin)

# --- Comparison 1: overall discordance prevalence (the RQ1 headline) --------
# The comparison arm is restricted to adults for whom the ALTERNATIVE equation
# is computable. With Sampson primary that excludes TG > 400, where NHANES
# publishes no Friedewald value. This restriction is a property of the
# comparison, not of the primary analysis, and the primary estimate below is
# still the full-cohort one.
sens_alt_design <- subset(apob_design, !is.na(ldl_alt))
sens_alt_n <- sum(!is.na(apob_design$variables$ldl_alt))
message("Section 11b comparison arm restricted to ", sens_alt_n, " of ",
        nrow(apob_design$variables), " goal-assigned records where ",
        ALT_EQ_LABEL, " is computable.")
sens_rq1_alt <- svyciprop(~I(discordant_alt == 1), sens_alt_design, method = "beta")
sens_rq1_comparison <- tibble(
  ldl_equation = c(paste0(PRIMARY_EQ_LABEL, " (PRIMARY, this analysis)"),
                    paste0(ALT_EQ_LABEL, " (for comparison)")),
  estimate = c(as.numeric(rq1_prevalence), as.numeric(sens_rq1_alt)),
  ci_low   = c(attr(rq1_prevalence, "ci")[1], attr(sens_rq1_alt, "ci")[1]),
  ci_high  = c(attr(rq1_prevalence, "ci")[2], attr(sens_rq1_alt, "ci")[2])
)
print(sens_rq1_comparison)

# Population-scale version of the same difference -- how many US adults are
# added to or removed from the "discordant" count purely by choice of equation?
sens_rq1_totals <- tibble(
  ldl_equation = c(PRIMARY_EQ_LABEL, ALT_EQ_LABEL),
  population_discordant = c(
    as.numeric(coef(svytotal(~I(discordant == 1), apob_design))[2]),
    as.numeric(coef(svytotal(~I(discordant_alt == 1), sens_alt_design))[2])
  )
) %>%
  mutate(difference_vs_primary = population_discordant - population_discordant[1])
print(sens_rq1_totals)

# --- Comparison 2: reclassification among adults at LDL-C goal --------------
sens_reclass_alt_design <- subset(apob_design, !is.na(ldl_alt) & at_ldl_goal_alt == 1)
sens_reclass_alt <- tryCatch(
  svyciprop(~I(discordant_alt == 1), sens_reclass_alt_design, method = "logit"),
  error = function(e) {
    message("Could not compute alternative-equation reclassification CI: ", conditionMessage(e))
    NULL
  }
)
if (!is.null(sens_reclass_alt)) {
  sens_reclass_comparison <- tibble(
    ldl_equation = c(PRIMARY_EQ_LABEL, ALT_EQ_LABEL),
    estimate = c(as.numeric(reclass_weighted), as.numeric(sens_reclass_alt)),
    ci_low   = c(attr(reclass_weighted, "ci")[1], attr(sens_reclass_alt, "ci")[1]),
    ci_high  = c(attr(reclass_weighted, "ci")[2], attr(sens_reclass_alt, "ci")[2])
  )
  print(sens_reclass_comparison)
  write.csv(sens_reclass_comparison, "sens_reclassification_by_ldl_equation.csv", row.names = FALSE)
}

# --- Comparison 3: the RQ8b statin comparison under BOTH equations ----------
# RQ8b itself (Section 8c) runs on the CONFIG-selected primary equation. This
# reproduces what the same comparison shows under the alternative, so the
# manuscript can state explicitly how much of any apparent statin effect is
# equation-driven. Under Friedewald this gap looked like a strong biological
# finding (2.7x, separated CIs); under Sampson it collapsed to ~1.6x with
# heavily overlapping CIs.
sens_rq8b_alt <- tryCatch(
  svyby(~I(discordant_alt == 1), ~on_statin, sens_reclass_alt_design,
        svyciprop, vartype = "ci", method = "logit"),
  error = function(e) {
    message("Could not compute alternative-equation statin comparison: ", conditionMessage(e))
    NULL
  }
)
if (!is.null(sens_rq8b_alt)) {
  sens_rq8b_alt_df <- as_tibble(sens_rq8b_alt)
  names(sens_rq8b_alt_df)[2] <- "prevalence_alt_equation"
  print(sens_rq8b_alt_df)
  write.csv(sens_rq8b_alt_df, "sens_rq8b_statin_discordance_alt_equation.csv", row.names = FALSE)
}

write.csv(sens_ldl_agreement,           "sens_ldl_equation_agreement.csv", row.names = FALSE)
write.csv(sens_status_change,           "sens_ldl_goal_status_change.csv", row.names = FALSE)
write.csv(sens_status_change_by_statin, "sens_ldl_goal_status_change_by_statin.csv", row.names = FALSE)
write.csv(sens_rq1_comparison,          "sens_rq1_discordance_by_ldl_equation.csv", row.names = FALSE)
write.csv(sens_rq1_totals,              "sens_rq1_population_totals_by_ldl_equation.csv", row.names = FALSE)

## NOTE: this message previously asserted that Sampson was the primary equation.
## That was a leftover from the period when CONFIG$primary_ldl_equation was set
## to "sampson", and it survived the revert to Friedewald -- so the console text
## contradicted both the CONFIG block and the correctly-labelled table printed
## immediately above it. Rewritten to read CONFIG at runtime so it cannot drift
## out of sync with the analysis again.
message("SENSITIVITY ANALYSIS complete (Section 11b). Primary equation for this ",
        "run: ", PRIMARY_EQ_LABEL, ". Comparison arm: ", ALT_EQ_LABEL, ". ",
        "The tables above show what the comparison arm would have produced. ",
        "The gap between them is a reportable finding in its own ",
        "right: it quantifies how much apparent ApoB discordance is attributable to ",
        "the LDL-C calculation rather than to lipoprotein biology -- recoverable at ",
        "zero marginal cost by changing equations, with no additional test ordered.")

## ---- 12. Export tables ---------------------------------------------------------

write.csv(rq1_summary,       "rq1_overall_discordance_prevalence.csv", row.names = FALSE)
write.csv(rq1b_summary,      "rq1b_overall_apob_nonattainment.csv",    row.names = FALSE)
write.csv(rq1b_cell_summary, "rq1b_goal_cell_breakdown.csv",           row.names = FALSE)
write.csv(rq2_by_condition,  "rq2_discordance_by_subgroup.csv",        row.names = FALSE)
write.csv(rq3_trend_df,      "rq3_discordance_trend_2005_2016.csv",    row.names = FALSE)
write.csv(rq4_or_table,      "rq4_discordance_predictors_OR.csv",      row.names = FALSE)

## =============================================================================
## 12b. SECONDARY PREVENTION EXTENSION (RQ6)
## Everything above (RQ0-RQ4, reclassification) deliberately excludes anyone
## with known ASCVD (ascvd_hx == 1), mirroring the source JAMA paper's
## primary-prevention-only scope. This section asks the identical discordance
## question in the population the source paper didn't touch: adults who
## already HAVE diagnosed ASCVD (self-reported prior MI, CHD, angina, or
## stroke -- MCQ160C/D/E/F) and are therefore already being treated to a
## goal, not newly evaluated for statin eligibility. Arguably this is a
## HIGHER-STAKES version of the same question: secondary prevention carries
## both the tightest LDL-C/ApoB goals (55/55) and the highest absolute event
## rate, so a discordant classification here represents more excess risk per
## person than the same classification in a primary-prevention adult.
##
## PREVENT-EQUATION CAVEAT: the AHA PREVENT-ASCVD equations estimate FIRST-
## event risk and are not validated for people who already have ASCVD --
## applying them here would be scientifically inappropriate, so this section
## does NOT call compute_prevent_risk(). Goal tiering for this population
## does not need PREVENT at all: assign_goals()'s case_when() assigns
## risk_tier = "secondary_prevention" for every ascvd_hx==1 row directly
## (see Section 5c), before it ever looks at prevent_ascvd_10yr. That column
## is set to NA below purely so assign_goals() has something to reference;
## it is never actually consulted for this population.
## =============================================================================

secondary_df <- derive_vars(primary_raw) %>%
  filter(
    RIDAGEYR >= CONFIG$age_min, RIDAGEYR <= CONFIG$age_max,
    pregnant == 0,
    ascvd_hx == 1,                 # SECONDARY prevention: known ASCVD (opposite of primary_df's filter)
    !is.na(LBXAPB), !is.na(ldl_primary), !is.na(LBXTR),   # see note on primary_df
    WTSAF2YR > 0
  ) %>%
  mutate(prevent_ascvd_10yr = NA_real_) %>%
  assign_goals() %>%
  mutate(n_cycles_pooled = length(primary_suffixes_used),
         wt_pooled = WTSAF2YR / n_cycles_pooled)

## LABEL CORRECTED 2026-08-16: this message previously read "(2011-2016)",
## which was a leftover from an earlier cycle selection. secondary_df is built
## from primary_raw with NO cycle filter, so it spans all six pooled cycles.
## The wrong label would have propagated into the manuscript methods.
message("Secondary-prevention analytic sample (2005-2016, all six cycles, ",
        "known ASCVD): ", nrow(secondary_df), " records, cycles present: ",
        paste(sort(unique(secondary_df$cycle)), collapse = ""), ".")

if (nrow(secondary_df) < 100) {
  warning("Secondary-prevention sample is small (n < 100) -- expect wide confidence intervals; ",
          "treat point estimates as exploratory, not precise.")
}

secondary_design <- svydesign(
  id = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~wt_pooled,
  nest = TRUE, data = secondary_df
)

## --- VERY-HIGH-RISK DIAGNOSTIC ------------------------------------------------
## The guideline states, in the same footnote that defines the category: "The
## majority of patients with clinical ASCVD are likely to be at very high risk."
## That is a checkable expectation, so check it rather than assume it. If our
## share falls well below half, the classification is still incomplete -- most
## likely because ACS timing, symptomatic PAD, revascularisation and heart
## failure are all unobservable in NHANES, which biases downward.
##
## Placed AFTER secondary_design is built: an earlier draft referenced the
## design object before it existed. Uses base table() rather than
## tidyr::pivot_wider(), because tidyr is not explicitly attached by this script.
{
  vh_n   <- sum(secondary_df$ascvd_very_high == 1, na.rm = TRUE)
  vh_pct <- 100 * vh_n / nrow(secondary_df)
  vh_w   <- svyciprop(~I(ascvd_very_high == 1), secondary_design, method = "beta")
  message("Very-high-risk ASCVD: ", vh_n, " of ", nrow(secondary_df),
          " (", round(vh_pct, 1), "% unweighted; ",
          round(100 * as.numeric(vh_w), 1), "% weighted, 95% CI ",
          round(100 * confint(vh_w)[1], 1), "-",
          round(100 * confint(vh_w)[2], 1), "%).")
  message("  Guideline expectation: 'the majority of patients with clinical ",
          "ASCVD are likely to be at very high risk'. Observed share is ",
          if (vh_pct >= 50) "consistent with that."
          else paste0("BELOW half -- report as an undercount driven by ",
                      "unobservable major events (ACS timing, symptomatic PAD) ",
                      "and unobservable high-risk conditions (revascularisation, ",
                      "heart failure, LDL-C >100 on maximal therapy)."))
  cat("  Rows = major ASCVD events (MI/stroke), columns = high-risk conditions:\n")
  print(table(major_events   = pmin(secondary_df$ascvd_major_event_count, 2),
              highrisk_conds = pmin(secondary_df$ascvd_highrisk_cond_count, 3)))
}

## --- CKD MAPPING DIAGNOSTIC ---------------------------------------------------
## Reports where CKD enters goal assignment, and proves it is confined to
## secondary prevention. The numbers it prints go into eMethods 1 and 2.
##
## REWRITTEN 2026-09-11 (round 4, comment 9). It previously counted only the
## adults CONFIG$ascvd_ckd_ldl55 moved from LDL-C <70 to <55, computed as
## "has CKD and is NOT very high risk". Once CONFIG$ascvd_ckd_very_high became
## TRUE, every secondary-prevention adult with CKD IS very high risk, so that
## quantity is empty by construction: it printed 0 under a label describing a
## reading the analysis had abandoned, which is the most misleading thing a
## diagnostic can do. It now reports the route that is actually in force.
{
  n_sec     <- nrow(secondary_df)
  n_ckd     <- sum(coalesce(secondary_df$ckd == 1, FALSE))
  ## The adults CKD alone carries into very-high-risk status: they have CKD and
  ## do NOT qualify by the events-and-conditions arithmetic of Figure 10. These
  ## are the adults the Figure 11 reading gains, and the count quoted in the
  ## manuscript Methods.
  via_ckd_only <- with(secondary_df,
                       coalesce(ckd == 1, FALSE) &
                       coalesce(ascvd_very_high == 1, FALSE) &
                       !(coalesce(ascvd_major_event_count >= 2, FALSE) |
                         coalesce(ascvd_major_event_count >= 1 &
                                  ascvd_highrisk_cond_count >= 2, FALSE)))
  n_via_ckd <- sum(via_ckd_only, na.rm = TRUE)

  ## The quantity the manuscript Methods quotes, which is NOT n_via_ckd.
  ## Being carried into very-high-risk status by CKD is not the same as GAINING
  ## a goal: an adult with established ASCVD and triglycerides >=150 already has
  ## <70 from the Clinical ASCVD row, and others take a goal from the diabetes
  ## row, so for them the CKD reading changes the VALUE (<70 to <55) rather than
  ## granting a goal where there was none. Counted here by rebuilding the
  ## counterfactual: very-high-risk WITHOUT the CKD clause, the Clinical ASCVD
  ## row recomputed from it, and the order-free pmin retaken across all five
  ## rows. Adults who hold a goal now and none in the counterfactual are the
  ## ones the Methods sentence is about.
  cf <- secondary_df
  cf$vh_no_ckd <- as.integer(
    coalesce(cf$ascvd_major_event_count >= 2, FALSE) |
    coalesce(cf$ascvd_major_event_count >= 1 &
             cf$ascvd_highrisk_cond_count >= 2, FALSE))
  cf$g_ascvd_cf <- with(cf, ifelse(
    vh_no_ckd == 1, 55,
    ifelse(coalesce(hypertrig == 1, FALSE), 70, NA_real_)))
  cf$apob_goal_cf <- local({
    v <- with(cf, pmin(g_ascvd_cf, g_dm_apob, g_shc_apob, g_htg_apob,
                       g_pp_apob, na.rm = TRUE))
    ifelse(is.finite(v), v, NA_real_)
  })
  gained <- !is.na(secondary_df$apob_goal) & is.na(cf$apob_goal_cf)
  tightened <- !is.na(secondary_df$apob_goal) & !is.na(cf$apob_goal_cf) &
               secondary_df$apob_goal < cf$apob_goal_cf
  n_gained    <- sum(gained, na.rm = TRUE)
  n_tightened <- sum(tightened, na.rm = TRUE)

  message("CKD mapping (secondary prevention only):")
  message("  Switches in force: ascvd_ckd_ldl55 = ",
          isTRUE(CONFIG$ascvd_ckd_ldl55),
          " (LDL-C goal), ascvd_ckd_very_high = ",
          isTRUE(CONFIG$ascvd_ckd_very_high), " (very-high-risk status).")
  message("  Secondary-prevention adults with CKD: ", n_ckd,
          " of ", n_sec,
          " (", round(100 * n_ckd / n_sec, 1), "% unweighted).")
  if (isTRUE(CONFIG$ascvd_ckd_very_high)) {
    message("  Carried into very-high-risk status by CKD alone, i.e. not ",
            "qualifying by the Figure 10 events-and-conditions criteria: ",
            n_via_ckd, ".")
    message("  Of those, GAINED an apolipoprotein B goal they would not ",
            "otherwise hold: ", n_gained,
            "  <-- this is the number quoted in the manuscript Methods.")
    message("  Held a goal already and had it TIGHTENED to <55: ", n_tightened,
            " (mostly triglycerides >=150, which the Clinical ASCVD row ",
            "already answers with <70).")
    if (n_gained + n_tightened != n_via_ckd)
      message("  NOTE: ", n_via_ckd - n_gained - n_tightened,
              " carried in by CKD alone neither gained nor tightened a goal, ",
              "their goal coming from another row at <55 or stricter.")
  } else {
    message("  ascvd_ckd_very_high is FALSE: CKD carries the LDL-C and ",
            "non-HDL-C goals only and no apoB goal follows from it.")
  }

  ## Assertions, not messages: these must hold or the run is wrong.
  ## Both are about PRIMARY prevention. Neither claims the CKD reading leaves
  ## apoB alone everywhere -- since 2026-09-10 it deliberately does not, in
  ## secondary prevention. What must never happen is a primary-prevention
  ## estimate moving, because ascvd_hx == 1 excludes those adults from that
  ## cohort before any goal is assigned.
  stopifnot(
    "CKD branch leaked into primary prevention" =
      !any(coalesce(primary_df$ascvd_hx == 1, FALSE)),
    "CKD branch altered a PRIMARY-PREVENTION apoB goal" =
      identical(
        primary_df$apob_goal,
        {
          v <- with(primary_df, pmin(g_ascvd_apob, g_dm_apob, g_shc_apob,
                                     g_htg_apob, g_pp_apob, na.rm = TRUE))
          ifelse(is.finite(v), v, NA_real_)
        })
  )
  message("  Assertions passed: no adult with established ASCVD sits in the ",
          "primary-prevention cohort, and no primary-prevention apoB goal was ",
          "touched by the CKD reading.")
}


## Refresh: subset() snapshots the design, so columns added by the update()
## above are absent from the subpopulation object unless it is rebuilt.
secondary_apob_design <- subset(secondary_design, has_apob_goal == 1)

## RQ6a -- overall discordance prevalence among secondary-prevention adults,
## directly comparable to RQ1's primary-prevention estimate.
rq6_prevalence <- svyciprop(~I(discordant == 1), secondary_apob_design, method = "beta")
rq6_summary <- tibble(
  population = "secondary_prevention",
  estimate = as.numeric(rq6_prevalence),
  ci_low   = attr(rq6_prevalence, "ci")[1],
  ci_high  = attr(rq6_prevalence, "ci")[2]
)
rq6_population <- svytotal(~I(discordant == 1), secondary_apob_design)
print(rq6_summary)
print(rq6_population)

## --- NCHS PRESENTATION-STANDARD CHECK -----------------------------------------
## Added 2026-08-15 as a 30% relative-standard-error rule; CORRECTED 2026-09-05
## to the criteria NCHS actually applies to proportions (Series 2, No. 175 --
## see nchs_present() above, defined near RQ1c). The RSE is still computed and
## written to the CSV because it is informative about precision, but it is not
## the standard and does not drive the flag. Computed and attached here so the
## verdict travels with the number into the CSV, rather than living only in a
## verification report someone may not read.
{
  m6  <- svymean(~I(discordant == 1), secondary_apob_design, deff = TRUE)
  e6  <- as.numeric(coef(m6))[2]; s6 <- as.numeric(SE(m6))[2]
  d6  <- as.numeric(deff(m6))[2]
  rse6 <- 100 * s6 / e6
  n6  <- nrow(secondary_apob_design$variables)
  eff6 <- n6 / d6
  ci6  <- c(attr(rq6_prevalence, "ci")[1], attr(rq6_prevalence, "ci")[2])
  flag <- !nchs_present(e6, ci6[1], ci6[2], n6, eff6)
  rq6_summary$n_unweighted   <- n6
  rq6_summary$rse_pct        <- round(rse6, 1)
  rq6_summary$effective_n    <- round(eff6)
  rq6_summary$nchs_unreliable <- flag
  print(rq6_summary)
  if (flag) {
    message("*** NCHS PRESENTATION STANDARD: secondary-prevention discordance ***")
    message("    Fails Series 2 No. 175. CI ", round(100*ci6[1], 2), "%-",
            round(100*ci6[2], 2), "% (absolute width ",
            round(ci6[2] - ci6[1], 4), "), effective n = ", round(eff6),
            ", unweighted n = ", n6, ", RSE = ", round(rse6, 1), "%.")
    message("    Suppress it, or carry an explicit unreliability marker ",
            "wherever it appears. Do not report it as a headline figure.")
  } else {
    message("Secondary-prevention discordance MEETS NCHS presentation ",
            "standards for proportions (Series 2 No. 175): CI ",
            round(100*ci6[1], 2), "%-", round(100*ci6[2], 2),
            "% (absolute width ", round(ci6[2] - ci6[1], 4),
            " <= 0.05), effective n = ", round(eff6), ", n = ", n6,
            ". RSE is ", round(rse6, 1), "%, which is NOT an NCHS criterion ",
            "for proportions; the estimate is imprecise but presentable.")
  }
}

## RQ6b -- overall ApoB non-attainment + full 2x2 breakdown, same structure as RQ1b.
rq6b_apob_prevalence <- svyciprop(~I(above_apob_goal == 1), secondary_apob_design, method = "beta")
rq6b_summary <- tibble(
  population = "secondary_prevention",
  estimate = as.numeric(rq6b_apob_prevalence),
  ci_low   = attr(rq6b_apob_prevalence, "ci")[1],
  ci_high  = attr(rq6b_apob_prevalence, "ci")[2]
)
print(rq6b_summary)

secondary_df <- secondary_df %>%
  mutate(goal_cell = case_when(
    at_ldl_goal == 1 & above_apob_goal == 0 ~ "concordant_at_goal",
    at_ldl_goal == 1 & above_apob_goal == 1 ~ "discordant",
    at_ldl_goal == 0 & above_apob_goal == 0 ~ "ldl_above_apob_at_goal",
    at_ldl_goal == 0 & above_apob_goal == 1 ~ "above_both"
  ))

secondary_design <- update(
  secondary_design,
  goal_cell = factor(
    case_when(
      at_ldl_goal == 1 & above_apob_goal == 0 ~ "concordant_at_goal",
      at_ldl_goal == 1 & above_apob_goal == 1 ~ "discordant",
      at_ldl_goal == 0 & above_apob_goal == 0 ~ "ldl_above_apob_at_goal",
      at_ldl_goal == 0 & above_apob_goal == 1 ~ "above_both"
    ),
    levels = c("concordant_at_goal", "discordant", "ldl_above_apob_at_goal", "above_both")
  )
)


## Rebuild the subpopulation: the update() immediately above added goal_cell to
## secondary_design, and subset() snapshots, so the existing
## secondary_apob_design predates that column. Without this the next two lines
## fail with "object 'goal_cell' not found".
secondary_apob_design <- subset(secondary_design, has_apob_goal == 1)

rq6_cell_prev <- svymean(~goal_cell, secondary_apob_design)   # retained for diagnostics
rq6_cell_tot  <- svytotal(~goal_cell, secondary_apob_design)
rq6_cell_summary <- kg_cell_summary(secondary_apob_design, "goal_cell")
print(rq6_cell_summary)

## RQ6c -- direct head-to-head: is discordance MORE common in secondary
## prevention than primary prevention? This is the single number that most
## directly answers "does checking ApoB have value in secondary prevention"
## -- if secondary-prevention discordance is similar to or higher than
## primary-prevention discordance, that's a strong argument the guideline's
## ApoB requirement matters just as much (or more) in patients already being
## treated, not just those newly being risk-stratified.
rq6_vs_rq1_comparison <- bind_rows(
  rq1_summary %>% mutate(population = "primary_prevention", n = nrow(primary_df)),
  rq6_summary %>% mutate(n = nrow(secondary_df))
) %>%
  select(population, n, estimate, ci_low, ci_high)
print(rq6_vs_rq1_comparison)

## Reclassification among secondary-prevention adults already at LDL-C goal
## -- direct analogue of Section 11's primary-prevention reclassification.
secondary_reclass_summary <- secondary_df %>%
  filter(at_ldl_goal == 1, !is.na(apob_goal)) %>%
  summarise(
    n = n(),
    reclassified_by_apob = sum(discordant, na.rm = TRUE),
    pct_reclassified = mean(discordant, na.rm = TRUE)
  )
print(secondary_reclass_summary)

## Secondary-prevention apoB subpopulation. Same rationale as apob_design:
## Figure 1 assigns an apoB goal to clinical-ASCVD patients at very high risk,
## and to those with hypertriglyceridaemia; for others it offers only an
## optional "consider apoB <55". Restricting here keeps the secondary-prevention
## estimates on the same footing as the primary ones.
secondary_apob_design <- subset(secondary_design, has_apob_goal == 1)
message("Secondary-prevention apoB-goal subpopulation: ",
        sum(secondary_df$has_apob_goal == 1, na.rm = TRUE), " of ",
        nrow(secondary_df), " records.")

secondary_reclass_design <- subset(secondary_apob_design, at_ldl_goal == 1)
secondary_reclass_weighted <- tryCatch(
  svyciprop(~I(discordant == 1), secondary_reclass_design, method = "logit"),
  error = function(e) {
    message("Could not compute weighted reclassification CI for secondary prevention ",
            "(likely too few at-goal secondary-prevention records in this cell): ", conditionMessage(e))
    NULL
  }
)
if (!is.null(secondary_reclass_weighted)) print(secondary_reclass_weighted)

## RQ6d -- statin/lipid-therapy use and treatment adequacy among secondary-
## prevention adults. "Adequate treatment" is operationalized pragmatically
## here as ANY guideline-directed lipid-lowering therapy (statin and/or
## non-statin add-on) combined with being AT ApoB goal. NHANES's public
## RXQ_RX file does not reliably expose statin DOSE/INTENSITY in a
## structured, cross-cycle-consistent way, so the standard clinical
## high/moderate/low-intensity statin classification is NOT attempted here
## -- that is a genuine granularity limit of this data source, disclosed
## rather than approximated with a guess. What this section CAN answer
## cleanly is the more basic and arguably more actionable split: is someone
## above ApoB goal because they're on NO lipid therapy at all (a
## diagnosis/access/adherence gap), or because they're already on therapy
## and still not at goal (a treatment-intensification gap -- the group most
## directly relevant to a "consider add-on therapy" clinical message).

secondary_df <- secondary_df %>%
  mutate(
    on_any_lipid_tx = as.integer(on_statin == 1 | on_nonstatin_lipid_tx == 1),
    treatment_status = case_when(
      on_any_lipid_tx == 0 & above_apob_goal == 1 ~ "untreated_above_goal",
      on_any_lipid_tx == 1 & above_apob_goal == 1 ~ "treated_above_goal",
      on_any_lipid_tx == 0 & above_apob_goal == 0 ~ "untreated_at_goal",
      on_any_lipid_tx == 1 & above_apob_goal == 0 ~ "treated_at_goal"
    )
  )

secondary_design <- update(
  secondary_design,
  on_statin = secondary_df$on_statin,
  on_nonstatin_lipid_tx = secondary_df$on_nonstatin_lipid_tx,
  on_any_lipid_tx = secondary_df$on_any_lipid_tx,
  treatment_status = factor(
    secondary_df$treatment_status,
    levels = c("untreated_above_goal", "treated_above_goal", "untreated_at_goal", "treated_at_goal")
  ),
  # income_tier/insured are needed here (not just later in RQ6e) because
  # svydesign objects only expose variables that have been attached via
  # update() -- referencing them in a subset()/svyby() call below without
  # this step would fail with "object 'income_tier' not found" even though
  # the columns exist in secondary_df itself.
  income_tier = factor(secondary_df$income_tier,
                        levels = c("below_poverty", "low_income", "higher_income")),
  insured = secondary_df$insured
)

## Refresh: subset() snapshots the design, so columns added by the update()
## above are absent from the subpopulation object unless it is rebuilt.
secondary_apob_design <- subset(secondary_design, has_apob_goal == 1)

rq6d_statin_prevalence <- svyciprop(~I(on_statin == 1), secondary_design, method = "logit")
rq6d_statin_summary <- tibble(
  metric   = "on_any_statin",
  estimate = as.numeric(rq6d_statin_prevalence),
  ci_low   = attr(rq6d_statin_prevalence, "ci")[1],
  ci_high  = attr(rq6d_statin_prevalence, "ci")[2]
)
print(rq6d_statin_summary)

rq6d_anylipidtx_prevalence <- svyciprop(~I(on_any_lipid_tx == 1), secondary_design, method = "logit")
rq6d_anylipidtx_summary <- tibble(
  metric   = "on_any_lipid_lowering_therapy",
  estimate = as.numeric(rq6d_anylipidtx_prevalence),
  ci_low   = attr(rq6d_anylipidtx_prevalence, "ci")[1],
  ci_high  = attr(rq6d_anylipidtx_prevalence, "ci")[2]
)
print(rq6d_anylipidtx_summary)

rq6d_treatment_prev <- svymean(~treatment_status, secondary_apob_design)
rq6d_treatment_tot  <- svytotal(~treatment_status, secondary_apob_design)
rq6d_treatment_summary <- tibble(
  treatment_status = sub("^treatment_status", "", names(coef(rq6d_treatment_prev))),
  prevalence = as.numeric(coef(rq6d_treatment_prev)),
  ci_low     = confint(rq6d_treatment_prev)[, 1],
  ci_high    = confint(rq6d_treatment_prev)[, 2],
  population = as.numeric(coef(rq6d_treatment_tot))
)
print(rq6d_treatment_summary)

# Among only those above ApoB goal (the clinically actionable slice): what
# share are undertreated (no lipid therapy at all, a diagnosis/access gap)
# vs already on therapy but still not at goal (a treatment-intensity gap)?
rq6d_among_above_goal <- subset(secondary_apob_design, above_apob_goal == 1)
rq6d_untreated_share <- tryCatch(
  svyciprop(~I(on_any_lipid_tx == 0), rq6d_among_above_goal, method = "logit"),
  error = function(e) {
    message("Could not compute weighted untreated-share CI (likely too few records in this cell): ",
            conditionMessage(e))
    NULL
  }
)
if (!is.null(rq6d_untreated_share)) {
  rq6d_untreated_summary <- tibble(
    metric   = "untreated_share_among_above_apob_goal",
    estimate = as.numeric(rq6d_untreated_share),
    ci_low   = attr(rq6d_untreated_share, "ci")[1],
    ci_high  = attr(rq6d_untreated_share, "ci")[2]
  )
  print(rq6d_untreated_summary)
  write.csv(rq6d_untreated_summary, "rq6d_secondary_prevention_untreated_share.csv", row.names = FALSE)
}

## =============================================================================
## RQ6d-ALL -- the DELIVERY question, on the full secondary-prevention cohort
## =============================================================================
## Added 2026-08-15. The estimate above is conditioned on having a
## guideline-specified apoB goal, which after the very-high-risk correction
## leaves n = 562 of 1,021. That is the right denominator for the
## guideline-specific question ("among adults the guideline assigns an apoB goal
## and who exceed it, how many are untreated?"), but it is the WRONG denominator
## for the clinical question most readers will have in mind: "how many adults
## with established ASCVD are on no lipid-lowering therapy at all?"
##
## That second question needs no apoB goal, so restricting to the subpopulation
## discards nearly half the cohort for nothing. Across three denominators the
## same headline moved 36.0% -> 32.3% -> 38.8%; those are three different
## questions, not an unstable estimate. Report both, labelled, so the movement
## is legible rather than alarming.
rq6d_untreated_all <- svyciprop(~I(on_any_lipid_tx == 0), secondary_design,
                                method = "logit")
rq6d_untreated_all_summary <- tibble(
  metric   = "untreated_share_ALL_secondary_prevention",
  n_unweighted = nrow(secondary_df),
  estimate = as.numeric(rq6d_untreated_all),
  ci_low   = attr(rq6d_untreated_all, "ci")[1],
  ci_high  = attr(rq6d_untreated_all, "ci")[2]
)
print(rq6d_untreated_all_summary)
write.csv(rq6d_untreated_all_summary,
          "rq6d_secondary_prevention_untreated_share_ALL.csv", row.names = FALSE)

## Same contrast by income and insurance, on the full cohort.
rq6e_all_income <- svyby(~I(on_any_lipid_tx == 0), ~income_tier,
                         subset(secondary_design, !is.na(income_tier)),
                         svyciprop, vartype = "ci", method = "logit")
rq6e_all_income_df <- as_tibble(rq6e_all_income)
names(rq6e_all_income_df)[2] <- "untreated_share"
cat("\nUntreated share by income -- ALL secondary prevention (n=",
    nrow(secondary_df), "):\n", sep = "")
print(rq6e_all_income_df)
write.csv(rq6e_all_income_df,
          "rq6e_untreated_by_income_ALL.csv", row.names = FALSE)

rq6e_all_insurance <- svyby(~I(on_any_lipid_tx == 0), ~insured,
                            subset(secondary_design, !is.na(insured)),
                            svyciprop, vartype = "ci", method = "logit")
rq6e_all_insurance_df <- as_tibble(rq6e_all_insurance)
names(rq6e_all_insurance_df)[2] <- "untreated_share"
cat("Untreated share by insurance -- ALL secondary prevention:\n")
print(rq6e_all_insurance_df)
## Unweighted denominators, so the triple-stratification concern the reviewer
## raised can be judged directly rather than inferred.
cat("Unweighted cell counts (insurance x untreated), full cohort:\n")
print(table(insured   = secondary_df$insured,
            untreated = secondary_df$on_any_lipid_tx == 0, useNA = "ifany"))
write.csv(rq6e_all_insurance_df,
          "rq6e_untreated_by_insurance_ALL.csv", row.names = FALSE)

message("RQ6d-ALL: the two untreated-share estimates answer different questions. ",
        "The restricted one (apoB-goal subpopulation, above goal) is the ",
        "guideline-specific measure; the full-cohort one is the delivery ",
        "measure. Report both with their denominators stated.")

write.csv(rq6d_statin_summary,      "rq6d_secondary_prevention_statin_use.csv", row.names = FALSE)
write.csv(rq6d_anylipidtx_summary,  "rq6d_secondary_prevention_any_lipid_tx_use.csv", row.names = FALSE)
write.csv(rq6d_treatment_summary,   "rq6d_secondary_prevention_treatment_status.csv", row.names = FALSE)

message("RQ6d: treatment-adequacy breakdown written to rq6d_*.csv. Read treatment_status as: ",
        "untreated_above_goal = diagnosis/access gap (no lipid therapy at all despite established ASCVD); ",
        "treated_above_goal = treatment-intensity gap (already on therapy, still not at ApoB goal -- ",
        "candidates for add-on therapy); treated_at_goal / untreated_at_goal = at goal either way.")

## RQ6e -- does the treatment gap identified in RQ6d concentrate among
## uninsured or lower-income secondary-prevention adults? This is the
## disparities extension of the treatment-adequacy question: it tests
## whether the "untreated_above_goal" share specifically (no lipid therapy
## at all, despite established ASCVD) is higher in the uninsured/lower-
## income slice, which would reframe RQ6d from a purely clinical finding
## into a health-access finding -- i.e., whether the gap looks more like
## "doctors aren't escalating therapy" or "patients aren't getting to a
## doctor/pharmacy in the first place." Sample sizes here are small (the
## whole secondary-prevention group is ~500 records before further
## splitting by income/insurance), so these are flagged explicitly as
## exploratory regardless of what they show.

rq6e_income_design <- subset(secondary_apob_design, !is.na(income_tier))
rq6e_untreated_by_income <- tryCatch(
  svyby(~I(on_any_lipid_tx == 0), ~income_tier, rq6e_income_design,
        svyciprop, vartype = "ci", method = "logit"),
  error = function(e) {
    message("Could not compute untreated share by income tier -- likely a cell too small ",
            "in this subgroup (", conditionMessage(e), "). Skipping this comparison.")
    NULL
  }
)
if (!is.null(rq6e_untreated_by_income)) {
  rq6e_untreated_by_income_df <- as_tibble(rq6e_untreated_by_income)
  names(rq6e_untreated_by_income_df)[2] <- "untreated_share"
  print(rq6e_untreated_by_income_df)
  write.csv(rq6e_untreated_by_income_df, "rq6e_secondary_prevention_untreated_by_income.csv", row.names = FALSE)
}

rq6e_insurance_design <- subset(secondary_apob_design, !is.na(insured))
rq6e_untreated_by_insurance <- tryCatch(
  svyby(~I(on_any_lipid_tx == 0), ~insured, rq6e_insurance_design,
        svyciprop, vartype = "ci", method = "logit"),
  error = function(e) {
    message("Could not compute untreated share by insurance status -- likely a cell too small ",
            "in this subgroup (", conditionMessage(e), "). Skipping this comparison.")
    NULL
  }
)
if (!is.null(rq6e_untreated_by_insurance)) {
  rq6e_untreated_by_insurance_df <- as_tibble(rq6e_untreated_by_insurance)
  names(rq6e_untreated_by_insurance_df)[2] <- "untreated_share"
  print(rq6e_untreated_by_insurance_df)
  write.csv(rq6e_untreated_by_insurance_df, "rq6e_secondary_prevention_untreated_by_insurance.csv", row.names = FALSE)
}

message("RQ6e: untreated-share breakdown by income tier and insurance status written to ",
        "rq6e_*.csv where the sample supported it. insured=0 means confirmed uninsured; ",
        "refused/don't-know insurance and income responses are coded NA and excluded from ",
        "these comparisons entirely, not assumed either way.")

write.csv(rq6_summary,               "rq6_secondary_prevention_discordance.csv", row.names = FALSE)
write.csv(rq6b_summary,              "rq6b_secondary_prevention_apob_nonattainment.csv", row.names = FALSE)
write.csv(rq6_cell_summary,          "rq6_secondary_prevention_goal_cell_breakdown.csv", row.names = FALSE)
write.csv(rq6_vs_rq1_comparison,     "rq6_primary_vs_secondary_comparison.csv", row.names = FALSE)
write.csv(secondary_reclass_summary, "rq6_secondary_prevention_reclassification.csv", row.names = FALSE)

message("Secondary-prevention extension complete. See rq6_*.csv and rq6d_*.csv files. n=", nrow(secondary_df),
        " secondary-prevention records with complete lipid panels, 2011-2016 pooled cycles. ",
        "Compare rq6_primary_vs_secondary_comparison.csv directly against RQ1's estimate to see ",
        "whether ApoB checking matters more, less, or about the same in already-treated patients.")

## =============================================================================
## 13-16. MORTALITY LINKAGE EXTENSION
## Adds NCHS's 2019 public-use Linked Mortality Files (NDI follow-up through
## Dec 31, 2019) so discordance/reclassification isn't just a cross-sectional
## snapshot -- this tests whether it actually predicts worse survival, which
## is exactly the "prognostic validity" gap flagged as a limitation earlier,
## and a step neither this analysis nor the source JAMA paper had taken yet.
##
## Source: NCHS Data Linkage, public-use LMF, 1999-2018 NHANES.
##   https://www.cdc.gov/nchs/data-linkage/mortality-public.htm
## File naming/columns confirmed against three independent sources: (1) the
## CDC file listing at ftp.cdc.gov/.../linked_mortality/, (2) the CDC 2019
## public-use LMF data dictionary/codebook (dated 2022-04-28), and (3) an
## NCHS-authored R read-in script (R_ReadInProgramAllSurveys.R) mirrored on
## GitHub, whose variable names/codings match the codebook exactly. NCHS
## keeps this ASCII layout stable release-over-release specifically so old
## read-in code keeps working -- only the follow-up cutoff and re-linkage
## results change between vintages, not the column positions -- but this
## was not independently re-verified against a raw 2019-vintage .dat file
## byte-for-byte (ftp.cdc.gov was not reachable from this environment to do
## that final check). The sanity-check print statements below exist for
## exactly this reason: if mortstat or ucod_leading show impossible values
## after a download, the fixed-width offsets have drifted and must be
## re-checked against the current codebook before trusting anything
## downstream.
##
## KNOWN CAVEAT: for 2015-2016 NHANES (cycle I) specifically, NCHS only
## releases UCOD_LEADING categories 001 (heart disease), 002 (cancer), and
## 010 (all other) -- cerebrovascular deaths are folded into "010" rather
## than broken out separately, unlike cycles D-H. This makes the
## "cardiovascular death" definition below (heart disease + cerebrovascular)
## slightly less complete for cycle I than for earlier cycles. This is a
## genuine cross-cycle comparability limitation to disclose if reporting the
## CVD-specific results, not a bug in this script.
## =============================================================================

mort_base_url <- "https://ftp.cdc.gov/pub/health_statistics/nchs/datalinkage/linked_mortality"

read_mortality <- function(cycle_label) {
  fname <- paste0("NHANES_", gsub("-", "_", cycle_label), "_MORT_2019_PUBLIC.dat")
  local_path <- file.path(data_dir, fname)
  if (!file.exists(local_path)) {
    url <- paste0(mort_base_url, "/", fname)
    message("Downloading ", url)
    ok <- tryCatch({ download.file(url, local_path, mode = "wb", quiet = TRUE); TRUE },
                    error = function(e) { message("  -> FAILED: ", e$message); FALSE })
    if (!ok || !file.exists(local_path) || file.size(local_path) == 0) {
      warning("Could not download mortality file for ", cycle_label)
      return(NULL)
    }
  }

  ## Fixed-width layout: NCHS public-use Linked Mortality File, 2019 release.
  ## The offsets are asserted by the value tables printed at the end of this
  ## function -- mortstat must be 0/1/NA and ucod_leading 1-10/NA, which they
  ## are for all 60,936 records across the six cycles.
  ##
  ## Two properties of these files govern how they must be read.
  ##
  ## 1. dodqtr, dodyear, wgt_new and sa_wgt_new are RESERVED BUT EMPTY in the
  ##    public-use release: every record carries spaces in those positions.
  ##    They are read and then dropped, rather than omitted from the spec, so
  ##    that the columns after them stay anchored to the published layout.
  ##
  ## 2. Trailing whitespace is stripped in the source files, so a record ends
  ##    at column 46, 47 or 48 depending on how wide its last value is. 43,483
  ##    of the 60,936 records -- 71%, and every record in the 2013-2014 and
  ##    2015-2016 files -- stop before column 48. Declaring permth_exm as
  ##    ending at a fixed column 48 therefore made readr report "one or more
  ##    parsing issues" once per cycle, six warnings per run, on a file that
  ##    was being read correctly. Ending the field at NA means "to the end of
  ##    the record", which is what these files actually are. Verified
  ##    2026-09-07: the parsed values are identical under both specifications
  ##    for all 60,936 records.
  mort <- readr::read_fwf(
    local_path,
    readr::fwf_cols(
      publicid     = c(1, 14),
      eligstat     = c(15, 15),
      mortstat     = c(16, 16),
      ucod_leading = c(17, 19),
      diabetes     = c(20, 20),
      hyperten     = c(21, 21),
      dodqtr       = c(22, 22),
      dodyear      = c(23, 26),
      wgt_new      = c(27, 34),
      sa_wgt_new   = c(35, 42),
      permth_int   = c(43, 45),
      permth_exm   = c(46, NA)
    ),
    col_types = readr::cols(
      publicid     = readr::col_character(),
      eligstat     = readr::col_integer(),
      mortstat     = readr::col_integer(),
      ucod_leading = readr::col_integer(),
      diabetes     = readr::col_integer(),
      hyperten     = readr::col_integer(),
      dodqtr       = readr::col_integer(),
      dodyear      = readr::col_integer(),
      wgt_new      = readr::col_double(),
      sa_wgt_new   = readr::col_double(),
      permth_int   = readr::col_integer(),
      permth_exm   = readr::col_integer()
    ),
    na = c(".", "")
  ) %>%
    mutate(SEQN = as.numeric(substr(publicid, 1, 5))) %>%
    # NOTE: renaming diabetes/hyperten -> mcod_diabetes/mcod_hyperten here is
    # load-bearing. trend_df already has its own clinical `diabetes` column
    # (from derive_vars(), lab/self-report based); the mortality file's
    # `diabetes` is a completely different thing -- a "was diabetes listed as
    # a contributing cause of death" flag. Keeping the collision would make
    # inner_join() silently rename both to diabetes.x/diabetes.y, which is
    # exactly what broke the Cox models earlier (bare `diabetes` never
    # existed in surv_df post-join, hence "object 'diabetes' not found").
    select(SEQN, eligstat, mortstat, ucod_leading,
           mcod_diabetes = diabetes, mcod_hyperten = hyperten,
           permth_int, permth_exm)

  # ---- STOP AND CHECK THIS OUTPUT ----
  # mortstat must be only 0 / 1 / NA. ucod_leading must be only 1-10 / NA.
  # If either shows anything else, the fixed-width column offsets above have
  # drifted for this release and must be re-verified against the current
  # codebook PDF before trusting the survival results downstream.
  message("  [", cycle_label, "] mortstat table: ",
          paste(capture.output(table(mort$mortstat, useNA = "ifany")), collapse = " | "))
  message("  [", cycle_label, "] ucod_leading table: ",
          paste(capture.output(table(mort$ucod_leading, useNA = "ifany")), collapse = " | "))

  mort
}

mortality_list <- purrr::map(trend_suffixes, function(sfx) {
  cyc_label <- cycles$cycle[cycles$suffix == sfx]
  m <- read_mortality(cyc_label)
  if (!is.null(m)) m$cycle <- sfx
  m
})
mortality_raw <- bind_rows(mortality_list)

## ---- 14. Merge mortality onto the pooled cohort & derive survival vars ----
## Uses the FULL trend cohort (2005-2016, 6 cycles) rather than just the
## 3-cycle primary window: for survival analysis, more person-time and more
## deaths matter more than ApoB-assay-vintage consistency, and 2005-2006 in
## particular contributes up to ~14 years of follow-up (vs. ~4 years for
## 2015-2016), which is valuable for event accrual given how rare the
## discordant group is (~2% of the cohort).

surv_df <- trend_df %>%
  inner_join(mortality_raw, by = c("SEQN", "cycle")) %>%
  filter(eligstat == 1, !is.na(mortstat), !is.na(permth_int)) %>%
  mutate(
    time_yrs = permth_int / 12,
    # Age-as-timescale fields for the Cox models below (see Section 16 note):
    # age_entry = age at NHANES exam, age_exit = age at death/censoring.
    # Left-truncating at age_entry and following to age_exit puts age itself
    # in the risk-set machinery rather than treating it as an ordinary linear
    # covariate -- the standard fix when age is by far the dominant mortality
    # predictor, as it is here.
    age_entry = RIDAGEYR,
    age_exit  = RIDAGEYR + time_yrs,
    death_allcause = mortstat,
    death_cvd = as.integer(mortstat == 1 & ucod_leading %in% c(1, 5)),  # heart disease + cerebrovascular
    # Three-group comparison for a dose-response-style check on whether ApoB
    # is adding real prognostic information beyond LDL-C:
    #   concordant_at_goal:     at BOTH goals                  (reference)
    #   discordant:             at LDL-C goal, ABOVE ApoB goal (group of interest)
    #   ldl_above_apob_at_goal: above LDL-C goal, AT ApoB goal (reverse discordance)
    #   above_both:             above BOTH goals
    #
    # CORRECTED 2026-08-02. The previous definition collapsed everyone above
    # LDL-C goal into "above_both" regardless of ApoB status, so that stratum
    # silently pooled two distinct lipid phenotypes and its label asserted
    # something untrue of roughly a quarter of its members. The four cells now
    # match the RQ1 goal-cell cross-tabulation exactly, which also makes the
    # two sections mutually checkable.
    mort_group = case_when(
      at_ldl_goal == 1 & above_apob_goal == 0 ~ "concordant_at_goal",
      at_ldl_goal == 1 & above_apob_goal == 1 ~ "discordant",
      at_ldl_goal == 0 & above_apob_goal == 0 ~ "ldl_above_apob_at_goal",
      at_ldl_goal == 0 & above_apob_goal == 1 ~ "above_both",
      TRUE ~ NA_character_
    ),
    mort_group = factor(mort_group, levels = c("concordant_at_goal", "discordant",
                                               "ldl_above_apob_at_goal", "above_both")),
    # Pre-computed factor for sex, referenced as sex_f in the Cox models below
    # (kept as a bare column regardless of the fix below, since it's good
    # practice for svycoxph formulas either way).
    sex_f = factor(RIAGENDR)
  ) %>%
  filter(!is.na(mort_group), time_yrs > 0, age_exit > age_entry)

message("Mortality-linked analytic sample: ", nrow(surv_df), " person-records, ",
        sum(surv_df$death_allcause), " all-cause deaths, ",
        sum(surv_df$death_cvd), " CVD deaths, follow-up through Dec 31, 2019.")

## NOTE: svycoxph (unlike svyglm/svyby/svyciprop, which all worked fine above
## on tibble-backed designs) has a known incompatibility with tibbles -- it
## can throw "all variables must be in design= argument" even when every
## formula variable genuinely is a column, because its internal model-frame
## construction doesn't handle tbl_df the same way it handles a base
## data.frame. Coercing to a plain data.frame before building this specific
## design object is the documented fix.
surv_design <- svydesign(
  id = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~wt_pooled,
  nest = TRUE, data = as.data.frame(surv_df)
)

## ---- 15. Weighted crude mortality rate by group (deaths per 1,000 PY) -----

rate_table <- surv_df %>%
  group_by(mort_group) %>%
  summarise(
    n = n(),
    deaths_allcause = sum(death_allcause),
    deaths_cvd = sum(death_cvd),
    person_years = sum(time_yrs),
    crude_rate_allcause_per1000py = 1000 * deaths_allcause / person_years,
    crude_rate_cvd_per1000py = 1000 * deaths_cvd / person_years,
    .groups = "drop"
  )
print(rate_table)
write.csv(rate_table, "rq5_mortality_rate_by_group.csv", row.names = FALSE)

## ---- 16. Survey-weighted Cox models ----------------------------------------
## Cause-specific hazard approach for the CVD model (competing non-CVD
## deaths are censored at their own death time) -- the simpler standard
## default. A full Fine-Gray subdistribution-hazard model (cmprsk package)
## would be a defensible upgrade if competing-risk bias is a concern; not
## implemented here.
##
## TWO FIXES APPLIED HERE following the first successful run, which produced
## non-significant, clinically implausible protective-direction (HR<1)
## estimates for the discordant/above_both groups vs. concordant_at_goal:
##   (1) Age-as-timescale: Surv(age_entry, age_exit, event) replaces
##       Surv(time_yrs, event) ~ ... + RIDAGEYR. Age is by far the dominant
##       mortality predictor here; modeling it as the timescale (left
##       truncation at exam age, follow-up to age at death/censoring) is a
##       more standard and typically better-powered way to handle it than a
##       linear covariate, especially with a rare exposure and few events.
##   (2) on_statin added as a covariate. The leading hypothesis for the
##       original counterintuitive result was confounding by indication /
##       the "healthy adherer" effect: concordant_at_goal likely
##       over-represents statin-treated adults, who tend to look protected
##       in observational data for reasons beyond their lipid numbers
##       (closer monitoring, better overall adherence to preventive care).
##       Without adjusting for treatment, that can drag the comparison
##       groups' hazard ratios toward or below 1 regardless of the true
##       biological relationship between discordance and risk. on_statin was
##       already computed in Section 3 (from RXQ_RX drug-name matching) and
##       carried through trend_df/surv_df untouched -- it just wasn't in the
##       model formula before.
## Both changes are additive to the original model structure; mort_group
## remains the exposure of interest and concordant_at_goal remains the
## reference level.

## Each model below is written as a direct, un-wrapped tryCatch() at the top
## level of the script (not inside a helper function). This is deliberate:
## coxph()/svycoxph() use non-standard evaluation to find their formula's
## variables, and passing formula/data through a custom wrapper function is a
## well-known way to break that lookup (R ends up searching the wrong
## environment for the columns). Writing each call out directly, even
## repetitively, avoids that trap entirely.
##
## Primary attempt: svycoxph (correct Taylor-series survey design variance).
## Fallback: weighted coxph() with a robust sandwich variance clustered on
## SDMVPSU. The fallback is an approximation -- it accounts for clustering
## and unequal weights but NOT for stratification (SDMVSTRA), so its
## standard errors will differ somewhat from a true survey-design Cox model.
## Treat fallback results as directionally reliable but not final-manuscript
## precise; if this script consistently lands on the fallback path, that's
## worth investigating separately (e.g. check packageVersion("survey")).

surv_flat <- as.data.frame(surv_df)

cox_allcause_method <- "svycoxph (survey-weighted, correct design variance), age-as-timescale + on_statin"
cox_allcause <- tryCatch(
  svycoxph(
    Surv(age_entry, age_exit, death_allcause) ~ mort_group + sex_f + diabetes + ckd + hypertrig + on_statin,
    design = surv_design
  ),
  error = function(e) {
    message("svycoxph failed for the all-cause model (", conditionMessage(e), ") -- ",
            "falling back to weighted coxph() with robust clustering on SDMVPSU.")
    cox_allcause_method <<- "coxph() fallback (weighted, clustered SE, no stratification), age-as-timescale + on_statin"
    survival::coxph(
      Surv(age_entry, age_exit, death_allcause) ~ mort_group + sex_f + diabetes + ckd + hypertrig + on_statin,
      data = surv_flat, weights = surv_flat$wt_pooled, cluster = surv_flat$SDMVPSU
    )
  }
)
message("All-cause Cox model fit via: ", cox_allcause_method)
cox_allcause_hr <- broom::tidy(cox_allcause, exponentiate = TRUE, conf.int = TRUE)
print(cox_allcause_hr)

cox_cvd_method <- "svycoxph (survey-weighted, correct design variance), age-as-timescale + on_statin"
cox_cvd <- tryCatch(
  svycoxph(
    Surv(age_entry, age_exit, death_cvd) ~ mort_group + sex_f + diabetes + ckd + hypertrig + on_statin,
    design = surv_design
  ),
  error = function(e) {
    message("svycoxph failed for the CVD model (", conditionMessage(e), ") -- ",
            "falling back to weighted coxph() with robust clustering on SDMVPSU.")
    cox_cvd_method <<- "coxph() fallback (weighted, clustered SE, no stratification), age-as-timescale + on_statin"
    survival::coxph(
      Surv(age_entry, age_exit, death_cvd) ~ mort_group + sex_f + diabetes + ckd + hypertrig + on_statin,
      data = surv_flat, weights = surv_flat$wt_pooled, cluster = surv_flat$SDMVPSU
    )
  }
)
message("CVD Cox model fit via: ", cox_cvd_method)
cox_cvd_hr <- broom::tidy(cox_cvd, exponentiate = TRUE, conf.int = TRUE)
print(cox_cvd_hr)

write.csv(cox_allcause_hr, "rq5_cox_allcause_HR.csv", row.names = FALSE)
write.csv(cox_cvd_hr,      "rq5_cox_cvd_HR.csv",      row.names = FALSE)

## ---- Forest-style plot of adjusted HRs for the mort_group term ------------

hr_plot_df <- bind_rows(
  cox_allcause_hr %>% filter(str_detect(term, "mort_group")) %>% mutate(outcome = "All-cause mortality"),
  cox_cvd_hr      %>% filter(str_detect(term, "mort_group")) %>% mutate(outcome = "CVD mortality")
) %>%
  mutate(term = str_remove(term, "mort_group"))

hr_plot <- ggplot(hr_plot_df, aes(x = term, y = estimate, ymin = conf.low, ymax = conf.high, color = outcome)) +
  geom_pointrange(position = position_dodge(width = 0.4)) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  coord_flip() +
  labs(title = "Adjusted Hazard Ratios vs. Concordant-at-Goal Reference",
       subtitle = "NHANES 2005-2016, linked mortality through Dec 31, 2019",
       x = NULL, y = "Hazard Ratio (95% CI)", color = NULL) +
  theme_minimal()
ggsave("rq5_mortality_hazard_ratios.png", hr_plot, width = 7, height = 4, dpi = 200)

message("Mortality linkage complete. See rq5_mortality_rate_by_group.csv, ",
        "rq5_cox_allcause_HR.csv, rq5_cox_cvd_HR.csv, and rq5_mortality_hazard_ratios.png")

message("Done. Tables written to CSV; trend figure written to apob_discordance_trend.png")
## Reads CONFIG at runtime rather than hard-coding the window, which had gone
## stale when the primary cohort was extended from three cycles to six.
message("Reminder: primary cohort for this run = ", CONFIG$primary_cycle_set,
        " (", length(primary_suffixes_used), " cycles, ",
        if (CONFIG$primary_cycle_set == "all6") "2005-2016" else "2011-2016", "); ",
        "trend analysis covers 2005-2016. ApoB was measured only in cycles D-I; ",
        "no ApoB data exist past 2015-2016 as of this writing. ",
        "Mortality follow-up runs through Dec 31, 2019 regardless of NHANES exam year.")
