# NHANES apolipoprotein B implementation analysis — code

R code implementing Figure 1 of the 2026 ACC/AHA dyslipidemia guideline across
NHANES 2005–2016, and computing goal assignment, goal attainment and the
associated sensitivity analyses.

Code state as of 15 September 2026.

No data are distributed here. Everything the analysis needs is downloaded at
run time from CDC and NCHS.

---

## Run order

The whole pipeline is one command:

```r
setwd("<this folder>")
source("run_all.R")
```

`run_all.R` runs the seven scripts below in order, in a single R session, and
writes a complete transcript to `run_log_<date>_<time>.txt`. Use it rather than
sourcing by hand: RStudio keeps only the last 1,000 console lines and this
pipeline prints far more, so without the log the beginning of every run is
lost. The log is also what makes two runs comparable — when a number moves,
diffing this run's log against the previous one shows which step it moved at.

| Order | Script | What it does | Time |
|---|---|---|---|
| 1 | `nhanes_apob_analysis.R` | Main analysis: download, merge, derive, implement Figure 1, every research question | 45-95 min |
| 2 | `precedence_and_model_analyses.R` | Precedence sensitivity across all 120 orderings, denominator variants, EPV, Firth | <1 min |
| 3 | `precedence_diagnostics.R` | Two diagnostics: precedence driver, stringency sensitivity | fast |
| 4 | `verify_run.R` | 29 automated checks; writes `verification_report.csv` | fast |
| 5 | `table1_and_figures.R` | Baseline characteristics table and the figures | fast |
| 6 | `pathway_and_sensitivity_analyses.R` | Goal-assignment pathways, secondary-prevention routes, phenotype at LDL-C goal, hypertension threshold, stroke-only ASCVD, no-age-limit secondary prevention, ASCVD risk-factor age rule | 33-45 min |
| 7 | `compare_blind_coding.R` | Agreement, Cohen's kappa, and the effect of coding differences on the headline estimates | fast |

Each later script reads objects the earlier ones leave in the global
environment, and each carries a preflight check that stops with instructions
rather than cascading errors if run out of sequence. Do not restart R or clear
the workspace part-way through.

`pathway_and_sensitivity_analyses.R` runs last of the article-producing
scripts deliberately: it is the most recently written, so anything it gets
wrong costs only itself — every other output is already on disk by the time
it starts.

`compare_blind_coding.R` **reads the independent coder's returned file, which is
not distributed here**, so it is marked optional in `run_all.R`: it is skipped
with a message rather than aborting the run. It is included so the comparison
method is inspectable.

Two further scripts sit outside the sequence and are run on their own, after the
main analysis, in the same session:

| Script | Runs how |
|---|---|
| `cohort_definition_diagnostics.R` | Diagnostics on cohort eligibility and the hypertension definition; changes nothing. |
| `export_blind_dataset.R` | Rebuilds the de-identified dataset that was given to an independent coder, guarding against 24 forbidden column patterns. |

## Data

Downloaded automatically on first run:

- NHANES `.xpt` files and the six NCHS Linked Mortality Files (2019 public-use
  release) into `nhanes_raw/`
- parsed, merged cohorts cached as `.rds` into `nhanes_cache/`

Subsequent runs read the cache. Delete `nhanes_cache/` to force a clean
re-parse, `nhanes_raw/` to force a fresh download. Neither directory is
committed; both are in `.gitignore`.

Outputs — roughly 57 CSVs named by research question (`rq*`, `rev*`, `sens_*`),
`verification_report.csv`, `session_info.txt` and the figures — are written into
the working directory when the scripts run. They are not committed either.

## Packages, and the versions the published results were produced under

**R 4.6.1 (2026-06-24 ucrt), x86_64-w64-mingw32/x64, Windows 11 x64.**
Numeric locale `LC_NUMERIC=C`.

| Package | Version | | Package | Version |
|---|---|---|---|---|
| haven | 2.5.5 | | broom | 1.0.13 |
| dplyr | 1.2.1 | | readr | 2.2.0 |
| tidyr | 1.3.2 | | stringr | 1.6.0 |
| tibble | 3.3.1 | | survival | 3.8-6 |
| purrr | 1.2.2 | | Matrix | 1.7-5 |
| survey | 4.5 | | preventr | 0.12.0 |
| srvyr | 1.3.1 | | logistf | 1.26.1 |
| ggplot2 | 4.0.3 | | mitools | 2.4 |
| scales | 1.4.0 | | | |

These are the versions recorded by `session_info.txt` on the run that produced
the results in the paper, not a wish-list. Every run writes that file again, so
if your numbers differ from the published ones, diffing it against this table is
the first place to look.

`survey` is the one to match most carefully: variance estimation and the
Korn-Graubard intervals come from it, and its defaults have changed across major
versions.

To install these exact versions:

```r
install.packages("remotes")
remotes::install_version("survey",  "4.5")
remotes::install_version("srvyr",   "1.3.1")
remotes::install_version("dplyr",   "1.2.1")
# ... and so on for the table above
```

Or, to capture your own environment as a lockfile once the packages are in
place:

```r
install.packages("renv")
renv::init(bare = TRUE)
renv::snapshot()      # writes renv.lock
```

**Nothing is installed without consent.** A missing package stops the run with
the exact `install.packages()` call to fix it. To opt in to automatic
installation for a session:

```r
options(apob_auto_install = TRUE)
```

`preventr` is not on the default CRAN mirror and comes from its r-universe
build; the error message gives the full command.

---

## What each script does

| Script | Lines | Purpose |
|---|---|---|
| `run_all.R` | 154 | Runs the pipeline in order in one session and writes a full transcript to `run_log_<date>_<time>.txt` |
| `nhanes_apob_analysis.R` | 3864 | Loads and merges six NHANES cycles, derives all variables, implements Figure 1 of the guideline, runs every research question, writes ~40 CSVs |
| `precedence_and_model_analyses.R` | 711 | Precedence sensitivity across all 120 orderings, denominator variants, events-per-variable, unadjusted correlates, Firth penalised regression, and the Sampson validity-ceiling sensitivity |
| `precedence_diagnostics.R` | 142 | Two diagnostics: which Figure 1 row drives the precedence span, and the stringency sensitivity |
| `verify_run.R` | 429 | 29 independent checks: convergence, separation, goal-cell consistency, survey design sanity, NCHS presentation standards, cycle influence, mortality event counts, config provenance |
| `table1_and_figures.R` | 398 | Baseline characteristics table and the figures |
| `pathway_and_sensitivity_analyses.R` | 547 | Goal assignment by Figure 1 pathway and their overlap, secondary-prevention routes, phenotype of adults at their LDL-C goal, hypertension threshold sensitivity, stroke-only ASCVD, secondary prevention without the PREVENT age restriction, and the ASCVD risk-factor age rule |
| `compare_blind_coding.R` | 224 | Agreement, Cohen's kappa, weighted kappa, and the effect of coding differences on the headline estimates |
| `export_blind_dataset.R` | 115 | Builds the de-identified dataset released to the independent coder, guarding against 24 forbidden column patterns |
| `cohort_definition_diagnostics.R` | 184 | What the LDL-C eligibility filter excludes, and the hypertension definition under first-reading versus averaged blood pressure |

---

## Analytic choices, all in one place

Every tunable decision lives in the `CONFIG` block at the top of
`nhanes_apob_analysis.R`. Nothing downstream hard-codes these; change a value,
re-source, and every result updates.

| Setting | Value |
|---|---|
| `primary_cycle_set` | `"all6"` (2005–2016) |
| `primary_ldl_equation` | `"sampson"` |
| `apply_apob_calibration` | `TRUE` |
| `htn_counts_treatment` | `TRUE` |
| `ascvd_ckd_ldl55` | `TRUE` |
| `tg_threshold` | `150` mg/dL |
| `apob_high_threshold` | `130` mg/dL |
| `prevent_goal_high_risk_cut` | `0.10` |
| `prevent_eligibility_cut` | `0.05` |
| `age_min` / `age_max` | `30` / `79` |

Each carries a comment recording its provenance — which guideline page, table or
CDC document the value came from, and where an earlier version got it wrong.
Line numbers are deliberately not given here: they drift with every edit, and
stale ones are worse than none.

---

## The three things worth reading first

**1. Goal assignment — `assign_goals()`.**

Figure 1's rows are not mutually exclusive and the guideline states no
precedence rule, so goals are assigned by a rule containing no ordering: an
apoB goal applies if **any** matched row states one (union), and where matched
rows differ the **most stringent** value applies (`pmin`). Both operations are
order-independent, so no permutation of the rows can change any result.

Per-row goals are held in explicit columns (`g_ascvd_apob`, `g_dm_apob`,
`g_shc_apob`, `g_htg_apob`, `g_pp_apob`), each `NA` where that row states no
apoB goal, so a silent row cannot cancel a goal another row states. The
superseded first-match assignment is retained in parallel as `*_firstmatch` for
the sensitivity comparison, and a most-lenient variant as `apob_goal_lenient`.

One branch sits outside the per-row columns. Figure 1 places a "With CKD"
bullet in the <55 mg/dL column of the Clinical ASCVD row, so chronic kidney
disease is a second route into that column for adults with established ASCVD.
`g_ascvd_ldl` implements it under `CONFIG$ascvd_ckd_ldl55`; `g_ascvd_apob` is
deliberately **not** modified, because the apoB <55 value is stated under "At
very high risk" alone, and reading it otherwise would require inferring
very-high-risk status from CKD, which the figure does not say. Setting the flag
to `FALSE` reproduces the superseded mapping exactly. A diagnostic block prints
how many adults moved, and `stopifnot()` asserts that no primary-prevention row
and no apoB goal was touched.

**2. ApoB cross-instrument calibration — `derive_vars()`.**

NHANES measured apoB on three instrument and laboratory combinations across the
six cycles. CDC publishes bridging equations for both transitions; these are
applied before anything reads `LBXAPB`, with the uncalibrated value preserved as
`LBXAPB_raw`. The 2005–2006 equation is **inverted** (× 0.923); the 2015–2016
equation already maps Cobas onto ProSpec and is applied **as published**
(× 0.9988 − 2.880).

**3. Cohort construction — where `primary_df` and `secondary_df` are built from
`derive_vars()`.**

Requires a calculable Sampson LDL-C rather than a published Friedewald value.
NHANES omits the latter above 400 mg/dL of triglycerides, so conditioning on it
would exclude the most severely hypertriglyceridaemic adults — the group the
guideline's hypertriglyceridaemia provisions address, and the highest-apoB group
in the data.

---

## Verification

`verify_run.R` runs 29 checks and writes `verification_report.csv` plus
`session_info.txt` — R version, platform, the version of every attached and
loaded package, and the `CONFIG` values actually used. If your numbers differ
from someone else's, that file is where the difference is diagnosable.

Expected state on a clean run: **27 pass, 1 fail, 1 skipped.**

- The **failure** is a prespecified events-per-variable gate closing as
  designed: 10 events against 70 required in the discordant mortality stratum.
  It is meant to fail, and it means no survival model should be fitted on that
  stratum.
- The **skipped** check compares against external targets that assume a
  triglyceride threshold of >200 mg/dL, which is not comparable to the
  guideline's ≥150. See the note at `tg_threshold` in `nhanes_apob_analysis.R`.

Leave-one-cycle-out influence is worth reading even though it now passes. The
discordant cell holds 9 to 14 unweighted cases per cycle, so dropping any single
cycle moves the pooled estimate appreciably: −14.2% dropping 2005–2006, +13.7%
dropping 2015–2016, against a 15% gate. The direction is unsystematic, and this
is a precision limit of a small cell rather than a cycle effect. Full numbers in
`rev4_cycle_influence_on_share.csv`.

This check carried a **warning** at −16.5% before the round-4 corrections
(strict goal boundaries and domain estimation); it now sits inside the gate.

Separately, the goal-assigned population produced by `assign_goals()` and by an
independently written implementation of the union rule agree exactly (n = 3,324
by both routes). That verifies the coding of the rule, **not** the reading of
Figure 1, which was checked by independent double-coding — see
`compare_blind_coding.R` for the method.

---

## Known limitations of the code

Cohort filters are applied before `svydesign()` rather than by subsetting a
full-file design. Subpopulation analyses within the cohort **do** use
`subset()`, which is what matters for the apoB estimates, but the
cohort-construction step departs from strict NCHS practice and has not been
tested against the alternative.

Firth penalisation is not implemented for complex survey designs. Its intervals
are not design-based and are point-estimate sensitivities only.

`options(survey.lonely.psu = "adjust")` is set, needed for the
secondary-prevention subgroup where some strata contribute a single PSU.

---

## Licence and data

Code released under the MIT Licence (`LICENSE`). The licence covers the code
only. NHANES and the NCHS Linked Mortality Files are not redistributed here;
they are US Government public-use datasets downloaded at run time from CDC, and
their use is governed by the NCHS Data User Agreement, which prohibits any
attempt to identify individual respondents. No participant-level file is
committed to this repository.
