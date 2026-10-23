# Fixed PK bridge entry point (0.15.0). No eval/parse or formula text from user data.
# PKNCA performs every NCA calculation; this script only maps JSON to PKNCA calls.
args <- commandArgs(TRUE)
if (length(args) >= 3 && nzchar(args[3])) .libPaths(c(args[3], .libPaths()))
suppressPackageStartupMessages(library(jsonlite))
request <- fromJSON(args[1], simplifyVector=TRUE)
for (p in names(request$pins)) {
  if (length(find.package(p, quiet=TRUE)) == 0 || as.character(packageVersion(p)) != request$pins[[p]])
    stop(paste("Required optional package", p, "version", request$pins[[p]]))
}
suppressPackageStartupMessages(library(PKNCA))
versions <- list(R=R.version.string, packages=lapply(names(request$pins), function(p) as.character(packageVersion(p))))
names(versions$packages) <- names(request$pins)
versions$loaded_packages <- as.list(vapply(loadedNamespaces(), function(p) as.character(packageVersion(p)), character(1)))

nca_pass <- function(conc, dose, p, params) {
  pc <- if (isTRUE(p$declared_points)) PKNCAconc(conc, conc~time|subject, include_half.life="include_hl")
        else PKNCAconc(conc, conc~time|subject)
  pd <- PKNCAdose(dose, dose~time|subject, route=p$route, duration="duration")
  intervals <- data.frame(start=0, end=Inf)
  for (k in params) intervals[[k]] <- TRUE
  opts <- list(auc.method=p$auc_method,
               conc.blq=list(first=p$blq$first, middle=p$blq$middle, last=p$blq$last),
               conc.na="drop", min.hl.points=as.integer(p$min_points),
               allow.tmax.in.half.life=FALSE, adj.r.squared.factor=as.numeric(p$adj_r2_factor),
               first.tmax=TRUE)
  warnings <- character()
  res <- withCallingHandlers(pk.nca(PKNCAdata(pc, pd, intervals=intervals, options=opts)),
    warning=function(w) {warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")})
  df <- as.data.frame(res)
  list(rows=list(subject=as.character(df$subject), parameter=as.character(df$PPTESTCD),
                 value=as.numeric(df$PPORRES), exclude=as.character(df$exclude)),
       warnings=unique(warnings), options=opts)
}

run_pknca <- function(p) {
  conc <- data.frame(subject=as.character(p$subject), time=as.numeric(p$time),
                     conc=as.numeric(p$concentration), stringsAsFactors=FALSE)
  if (isTRUE(p$declared_points)) conc$include_hl <- as.logical(p$include_half_life)
  dose <- data.frame(subject=as.character(p$dose_subject), time=0, dose=as.numeric(p$dose),
                     duration=as.numeric(p$duration), stringsAsFactors=FALSE)
  first <- nca_pass(conc, dose, p, p$parameters_pass1)
  out <- list(pass1=first, r_environment=versions,
              pknca_c0_method=eval(formals(PKNCA::pk.calc.c0)$method))
  if (isTRUE(p$bolus)) {
    r <- first$rows
    c0 <- r$value[r$parameter == "c0"]; names(c0) <- r$subject[r$parameter == "c0"]
    second <- conc
    for (s in names(c0)) {
      if (!is.finite(c0[[s]])) next
      # The predose placeholder at time 0 is replaced by PKNCA's back-extrapolated C0.
      second$conc[second$subject == s & second$time == 0] <- c0[[s]]
      if (isTRUE(p$declared_points)) second$include_hl[second$subject == s & second$time == 0] <- FALSE
    }
    out$pass2 <- nca_pass(second, dose, p, p$parameters_pass2)
    out$c0_inserted <- as.list(c0)
  }
  out
}

result <- tryCatch(if (request$method == "pk_doctor") list(r_environment=versions)
                   else if (request$method == "pknca") run_pknca(request$data)
                   else stop("Unsupported R method"),
                   error=function(e) list(error=conditionMessage(e), r_environment=versions))
write_json(result, args[2], digits=I(17), auto_unbox=TRUE, na="null", null="null")
