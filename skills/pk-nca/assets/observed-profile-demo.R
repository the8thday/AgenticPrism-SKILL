# Educational synthetic profile, not an experimental dataset or validated run.
# Run from any directory: Rscript --vanilla observed-profile-demo.R NEW_OUTPUT_DIR
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Supply one new output directory")
if (file.exists(args[[1L]])) stop("Output path already exists")
if (!requireNamespace("PKNCA", quietly = TRUE)) {
  stop("PKNCA is missing; use a project-local R library. No packages were installed.")
}

profile <- data.frame(subject = "synthetic-1",
                      time = c(0, 1, 2, 4, 8),
                      concentration = c(1, 7, 5, 3, 1))
intervals <- data.frame(start = 0, end = 8,
                        cmax = TRUE, tmax = TRUE, auclast = TRUE)
choices <- list(auc.method = "lin up/log down")
input <- PKNCA::PKNCAdata(
  PKNCA::PKNCAconc(profile, concentration ~ time | subject),
  intervals = intervals, options = choices
)
if (!dir.create(args[[1L]], recursive = TRUE)) stop("Cannot create output directory")
out <- function(name) file.path(args[[1L]], name)
write.csv(profile, out("synthetic-profile.csv"), row.names = FALSE)
write.csv(intervals, out("intervals.csv"), row.names = FALSE)
saveRDS(input, out("input.rds"))
dput(PKNCA::PKNCA.options(), file = out("package-default-options.R"))
dput(choices, file = out("explicit-options.R"))
writeLines(c("Educational demo; not a study report or validation benchmark.",
             "Synthetic units: hours and mg/L. No dose or terminal endpoints.",
             "Package warnings and reportability still require review."),
           out("scope.txt"))
capture.output(sessionInfo(), file = out("sessionInfo.txt"))
capture.output(citation("PKNCA"), file = out("package-citation.txt"))
warnings_seen <- character()
result <- tryCatch(withCallingHandlers(PKNCA::pk.nca(input), warning = function(w) {
  warnings_seen <<- c(warnings_seen, conditionMessage(w))
  invokeRestart("muffleWarning")
}), error = function(e) {
  writeLines(conditionMessage(e), out("error.txt"))
  stop(e)
}, finally = writeLines(warnings_seen, out("warnings.txt")))
saveRDS(result, out("result.rds"))
write.csv(as.data.frame(result), out("package-results.csv"), row.names = FALSE)
capture.output(summary(result), file = out("package-summary.txt"))
