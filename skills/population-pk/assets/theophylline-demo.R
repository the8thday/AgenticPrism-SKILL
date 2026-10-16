# Educational adaptation of https://nlmixr2.org/articles/running_nlmixr.html
# Uses the public teaching dataset distributed by nlmixr2data; not user data.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Supply one new output directory")
if (file.exists(args[[1L]])) stop("Output path already exists")
required <- c("nlmixr2", "nlmixr2est", "rxode2", "nlmixr2data")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing packages: ", paste(missing, collapse = ", "),
                          ". No packages were installed.")
library(nlmixr2)

# Apparent clearance and volume for oral-only data, with F fixed to one.
# These initial values and variability structure are educational choices.
oral_model <- function() {
  ini({
    log_ka <- 0.4
    log_cl <- 1.1
    log_v <- 3.4
    eta_ka ~ 0.5
    eta_cl ~ 0.2
    eta_v ~ 0.1
    error_sd <- 0.8
  })
  model({
    ka <- exp(log_ka + eta_ka)
    cl <- exp(log_cl + eta_cl)
    v <- exp(log_v + eta_v)
    linCmt() ~ add(error_sd)
  })
}
events <- nlmixr2data::theo_sd
seed <- 20261013L
set.seed(seed)
control <- nlmixr2est::saemControl(seed = seed, nBurn = 200, nEm = 300, print = 0)
if (!dir.create(args[[1L]], recursive = TRUE)) stop("Cannot create output directory")
out <- function(name) file.path(args[[1L]], name)
write.csv(events, out("packaged-events.csv"), row.names = FALSE)
dput(oral_model, file = out("model.R"))
saveRDS(control, out("control.rds"))
writeLines(as.character(seed), out("seed.txt"))
capture.output(sessionInfo(), file = out("sessionInfo.txt"))
capture.output(lapply(required, citation), file = out("package-citations.txt"))
writeLines(c("Educational template; no claim of a validated fit.",
             "Data: nlmixr2data::theo_sd; original event coding retained.",
             "Review dataset units, fit diagnostics and CL/F, V/F interpretation.",
             "VPC, bootstrap and interpretation facts are not generated here."),
           out("scope.txt"))
warnings_seen <- character()
fit <- tryCatch(withCallingHandlers(
  nlmixr2::nlmixr2(oral_model, events, est = "saem", control = control),
  warning = function(w) {
    warnings_seen <<- c(warnings_seen, conditionMessage(w))
    invokeRestart("muffleWarning")
  }
), error = function(e) {
  writeLines(conditionMessage(e), out("error.txt"))
  stop(e)
}, finally = writeLines(warnings_seen, out("warnings.txt")))
saveRDS(fit, out("fit.rds"))
capture.output(print(fit), file = out("fit-summary.txt"))
