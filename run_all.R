## =============================================================================
## run_all.R -- run the whole pipeline and keep a complete transcript
##
##   setwd("<the folder holding these scripts>")
##   source("run_all.R")
##
## Runs the seven scripts in order, in one session, exactly as running them by
## hand does -- every object still lands in the global environment, so nothing
## downstream behaves differently.
##
## What it adds is the log. Everything R prints, on BOTH channels (printed
## results, and the message/warning channel that carries the "[PASS]" lines and
## the "wrote ..." lines), goes to run_log_<date>_<time>.txt as well as to the
## screen. RStudio keeps only the last 1,000 console lines by default and this
## pipeline prints far more than that, so without a log the beginning of every
## run is lost by the time it finishes.
##
## The log is also what makes two runs comparable. When a number moves, diffing
## this run's log against the previous one shows exactly which step it moved at,
## which is faster and more reliable than reasoning about what should have
## changed.
##
## The transcript is closed by on.exit(), so it closes even if a script stops
## with an error -- otherwise a failed run would leave the console silently
## redirected into a file and look like R had hung.
## =============================================================================

run_all <- function(
  ## pathway_and_sensitivity_analyses.R runs LAST of the article-producing
  ## scripts deliberately: it is the most recently written, so anything it
  ## gets wrong costs only itself -- every other output is already on disk by
  ## the time it starts.
  ## compare_blind_coding.R also runs here. It was previously run by hand, which
  ## meant its outputs (blind_agreement_summary.csv and friends) were left on
  ## disk from a run predating the strict-goal-boundary correction while every
  ## other output had moved on. eTable 10 of the supplement compares the two
  ## readings' headline estimates, so a stale blind file puts a 2.25% next to a
  ## 2.54% that were never computed under the same rule. It runs after
  ## pathway_and_sensitivity_analyses.R because it needs nothing from it and,
  ## for the same reason, any failure then costs only itself. It is skipped, not failed, when the blind
  ## coder's returned file is absent -- that file is not in the repository.
  scripts = c("nhanes_apob_analysis.R",
              "precedence_and_model_analyses.R",
              "precedence_diagnostics.R",
              "verify_run.R",
              "table1_and_figures.R",
              "pathway_and_sensitivity_analyses.R",
              "compare_blind_coding.R"),
  ## Scripts that need an input the repository does not carry. A stop() from one
  ## of these is reported and the run continues; a stop() from any other script
  ## still aborts, which is what you want from the ones that produce the article.
  optional = "compare_blind_coding.R",
  log_file = NULL
) {
  missing <- scripts[!file.exists(scripts)]
  if (length(missing))
    stop("\n\nNot found in ", getwd(), ":\n  ",
         paste(missing, collapse = "\n  "),
         "\n\nIs the working directory right?\n", call. = FALSE)

  ## Any sink left open by an earlier failed run would nest inside this one and
  ## the console would stay dead afterwards. Clear them first.
  ##
  ## NOTE the asymmetry, which is a trap: sink.number() returns 0 when nothing
  ## is diverted, but sink.number(type = "message") returns 2 -- the connection
  ## number of stderr -- because messages are always going somewhere. Testing
  ## it with "> 0" loops forever.
  if (sink.number(type = "message") != 2) sink(type = "message")
  while (sink.number() > 0) sink()

  if (is.null(log_file))
    log_file <- paste0("run_log_", format(Sys.time(), "%Y%m%d_%H%M"), ".txt")
  con <- file(log_file, open = "wt")

  ## split = TRUE keeps output visible while it is being written.
  sink(con, split = TRUE)
  sink(con, type = "message")

  ## Runs on normal completion AND on error, which is the point of it.
  on.exit({
    sink(type = "message")
    sink()
    close(con)
    cat("\nTranscript written to ", file.path(getwd(), log_file), "\n", sep = "")
  }, add = TRUE)

  cat("=================================================================\n")
  cat(" PIPELINE RUN  ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n", sep = "")
  cat(" Working directory: ", getwd(), "\n", sep = "")
  cat(" R: ", R.version.string, " on ", R.version$platform, "\n", sep = "")
  cat("=================================================================\n\n")

  t0 <- Sys.time()
  timings <- data.frame(script = character(0), minutes = numeric(0),
                        stringsAsFactors = FALSE)

  for (f in scripts) {
    cat("\n-----------------------------------------------------------------\n")
    cat(">>> ", f, "   started ", format(Sys.time(), "%H:%M:%S"), "\n", sep = "")
    cat("-----------------------------------------------------------------\n")
    ts <- Sys.time()
    ## local = FALSE is source()'s default and means the global environment,
    ## so calling it from inside this function changes nothing about where the
    ## objects end up.
    if (f %in% optional) {
      ok <- tryCatch({ source(f, local = FALSE, echo = FALSE); TRUE },
                     error = function(e) {
                       cat("\n*** SKIPPED ", f, ": ", conditionMessage(e), "\n",
                           sep = "")
                       FALSE
                     })
      if (!ok) {
        timings <- rbind(timings,
                         data.frame(script = paste0(f, "  [SKIPPED]"),
                                    minutes = NA_real_,
                                    stringsAsFactors = FALSE))
        next
      }
    } else {
      source(f, local = FALSE, echo = FALSE)
    }
    mins <- as.numeric(difftime(Sys.time(), ts, units = "mins"))
    timings <- rbind(timings,
                     data.frame(script = f, minutes = round(mins, 2),
                                stringsAsFactors = FALSE))
    cat("\n<<< ", f, " finished in ", round(mins, 2), " min\n", sep = "")
  }

  cat("\n=================================================================\n")
  cat(" RUN COMPLETE in ",
      round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1),
      " minutes\n", sep = "")
  cat("=================================================================\n")
  print(timings, row.names = FALSE)

  ## A one-line summary at the very end, so the state of the run is readable
  ## without scrolling back through it.
  if (exists("verification_report", envir = .GlobalEnv)) {
    vr <- get("verification_report", envir = .GlobalEnv)
    cat("\nVerification: ",
        paste(paste0(names(table(vr$status)), "=", as.integer(table(vr$status))),
              collapse = "  "),
        "   (", nrow(vr), " checks)\n", sep = "")
    bad <- vr[vr$status %in% c("FAIL", "WARN"), c("description", "detail")]
    if (nrow(bad)) {
      cat("\nNot passing:\n")
      for (i in seq_len(nrow(bad)))
        cat("  - ", bad$description[i], ": ", bad$detail[i], "\n", sep = "")
    }
  }

  ## ---- session and package provenance -----------------------------------
  ## Round-6 (Dani): "preserve package/session information". The README lists
  ## the versions this analysis was developed against; this records the ones
  ## that actually produced THIS run, next to its outputs, so a reader can see
  ## whether their environment differs from the submitted one rather than
  ## having to trust a table written by hand.
  si <- utils::sessionInfo()
  writeLines(utils::capture.output(print(si)), "session_info.txt")

  pk <- c("survey", "srvyr", "dplyr", "tidyr", "purrr", "haven", "ggplot2",
          "readr", "preventr", "tibble", "scales")
  pk <- pk[vapply(pk, requireNamespace, logical(1), quietly = TRUE)]
  vers <- data.frame(
    package = pk,
    version = vapply(pk, function(x) as.character(utils::packageVersion(x)), ""),
    stringsAsFactors = FALSE)
  vers <- rbind(data.frame(package = "R",
                           version = paste(R.version$major, R.version$minor,
                                           sep = "."),
                           stringsAsFactors = FALSE),
                vers)
  write.csv(vers, "session_packages.csv", row.names = FALSE)
  cat("\nEnvironment written to session_info.txt and session_packages.csv\n")
  print(vers, row.names = FALSE)

  invisible(timings)
}

run_all()
