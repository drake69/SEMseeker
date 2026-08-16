# ---------------------------------------------------------------------------
# Run the test suite against the INSTALLED package and print explicit totals.
#
#   Rscript dev/run-suite.R                       # whole suite
#   Rscript dev/run-suite.R 7-association_analysis  # one file, to iterate
#
# Why a file and not a -e one-liner: in PowerShell a double-quoted string
# interpolates $passed, $failed, $error, so `sum(df$passed)` reaches R as
# `sum(df)` and the run ends in a parse error after the tests have already
# taken their time. A script has no quoting layer to get wrong.
#
# Why load_package = "installed" and not devtools::test(): under load_all a
# helper file's closure is package:SEMseeker, whose environment chain reaches
# neither the global environment nor setup.R's, so the two harnesses do not
# measure the same thing - and the one that skips tests silently is the one that
# hides defects. See AI-254.
#
# Why explicit totals and not the default reporter: the totals come from the
# returned data frame, so they cannot be lost to a truncated console.
# ---------------------------------------------------------------------------

filter <- commandArgs(trailingOnly = TRUE)
filter <- if (length(filter) && nzchar(filter[1])) filter[1] else NULL

res <- testthat::test_dir("tests/testthat",
                          package        = "SEMseeker",
                          load_package   = "installed",
                          filter         = filter,
                          reporter       = "silent",
                          stop_on_failure = FALSE)

df <- as.data.frame(res)

## The report goes to a FILE as well as to the console, under dev/test-reports/,
## named after the platform. On a VM the repository is a shared folder, so the
## file lands on the host and can be read there directly - no copying a console
## buffer between two machines, and nothing lost to a window that scrolled.
report_dir <- file.path("dev", "test-reports")
dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)
report <- file.path(report_dir,
                    paste0("suite-", gsub("[^A-Za-z0-9]+", "-", R.version$platform),
                           if (!is.null(filter)) paste0("-", gsub("[^A-Za-z0-9]+", "-", filter)) else "",
                           ".txt"))

emit <- function(...) {
  line <- paste0(...)
  cat(line, "\n", sep = "")
  cat(line, "\n", sep = "", file = report, append = TRUE)
}

if (file.exists(report)) unlink(report)

emit("================================================")
emit(R.version.string, " / ", R.version$platform)
emit("run at        : ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
if (!is.null(filter)) emit("filter        : ", filter)
emit("PASS          : ", sum(df$passed))
emit("FAIL          : ", sum(df$failed))
emit("ERROR         : ", sum(df$error))
emit("SKIP          : ", sum(df$skipped))
emit("================================================")

bad_idx <- which(df$failed > 0 | df$error)
if (length(bad_idx)) {
  emit("")
  emit("Not green:")
  for (i in bad_idx)
    emit("  ", df$file[i], " :: ", df$test[i],
         "  (failed=", df$failed[i], " error=", df$error[i], ")")

  ## The counts say WHERE, the messages say WHAT. Without them a report from
  ## another machine only tells you that something is wrong, which is the
  ## position we were in with a CI job that died leaving no log.
  emit("")
  emit("================ messages ================")
  for (i in bad_idx) {
    emit("")
    emit("--- ", df$file[i], " :: ", df$test[i])
    for (r in res[[i]]$results) {
      if (inherits(r, c("expectation_failure", "expectation_error"))) {
        emit("  [", class(r)[1], "] ",
             gsub("\n", "\n      ", conditionMessage(r)))
        srcref <- r$srcref
        if (!is.null(srcref))
          emit("      at line ", srcref[1])
      }
    }
  }
  emit("")
  emit("==========================================")
} else {
  emit("")
  emit("Nothing failed.")
}

cat("\nReport written to: ", report, "\n", sep = "")
