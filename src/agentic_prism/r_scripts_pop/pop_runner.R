# Fixed population-PK bridge entry point (0.15.1). nlmixr2 (FOCEi) estimates every value.
# Models come from the fixed templates below; parameter names arriving in JSON are matched
# against a whitelist and only ever become symbols through bquote, never parsed text.
args <- commandArgs(TRUE)
if (length(args) >= 3 && nzchar(args[3])) .libPaths(c(args[3], .libPaths()))
suppressPackageStartupMessages(library(jsonlite))
request <- fromJSON(args[1], simplifyVector=TRUE)
for (p in names(request$pins)) {
  # package_version comparison: R normalises "15.6.0-1" to 15.6.0.1.
  if (length(find.package(p, quiet=TRUE)) == 0 || packageVersion(p) != package_version(request$pins[[p]]))
    stop(paste("Required optional package", p, "version", request$pins[[p]]))
}
suppressPackageStartupMessages({library(rxode2); library(nlmixr2est)})
versions <- list(R=R.version.string, packages=lapply(names(request$pins), function(p) as.character(packageVersion(p))))
names(versions$packages) <- names(request$pins)
versions$loaded_packages <- as.list(vapply(loadedNamespaces(), function(p) as.character(packageVersion(p)), character(1)))

PARAMS <- list(linear_2cmt=c("cl","v1","q","v2"),
               linear_mm_2cmt=c("cl","v1","q","v2","vmax","km"),
               qss_tmdd_2cmt=c("cl","v1","q","v2","kss","kint","kdeg","r0"))

template <- function(model, route, species) {
  # Every structural parameter carries an eta here; undeclared etas are removed afterwards.
  if (model == "linear_2cmt") f <- function() {
    ini({ lcl <- 0; lv1 <- 0; lq <- 0; lv2 <- 0; lka <- 0; lfbio <- fix(0)
          eta_cl ~ 0.1; eta_v1 ~ 0.1; eta_q ~ 0.1; eta_v2 ~ 0.1; eta_ka ~ 0.1; prop_sd <- 0.2; add_sd <- 0.1 })
    model({ cl <- exp(lcl + eta_cl); v1 <- exp(lv1 + eta_v1); q <- exp(lq + eta_q); v2 <- exp(lv2 + eta_v2)
            ka <- exp(lka + eta_ka); fbio <- exp(lfbio)
            cp <- central/v1
            d/dt(depot) <- -ka*depot
            d/dt(central) <- ka*depot - cl*cp - q*cp + q*peripheral/v2
            d/dt(peripheral) <- q*cp - q*peripheral/v2
            f(depot) <- fbio
            cp ~ prop(prop_sd) + add(add_sd) })
  }
  if (model == "linear_mm_2cmt") f <- function() {
    ini({ lcl <- 0; lv1 <- 0; lq <- 0; lv2 <- 0; lvmax <- 0; lkm <- 0; lka <- 0; lfbio <- fix(0)
          eta_cl ~ 0.1; eta_v1 ~ 0.1; eta_q ~ 0.1; eta_v2 ~ 0.1; eta_vmax ~ 0.1; eta_km ~ 0.1; eta_ka ~ 0.1
          prop_sd <- 0.2; add_sd <- 0.1 })
    model({ cl <- exp(lcl + eta_cl); v1 <- exp(lv1 + eta_v1); q <- exp(lq + eta_q); v2 <- exp(lv2 + eta_v2)
            vmax <- exp(lvmax + eta_vmax); km <- exp(lkm + eta_km); ka <- exp(lka + eta_ka); fbio <- exp(lfbio)
            cp <- central/v1
            d/dt(depot) <- -ka*depot
            d/dt(central) <- ka*depot - cl*cp - vmax*cp/(km + cp) - q*cp + q*peripheral/v2
            d/dt(peripheral) <- q*cp - q*peripheral/v2
            f(depot) <- fbio
            cp ~ prop(prop_sd) + add(add_sd) })
  }
  if (model == "qss_tmdd_2cmt") f <- function() {
    ini({ lcl <- 0; lv1 <- 0; lq <- 0; lv2 <- 0; lkss <- 0; lkint <- 0; lkdeg <- 0; lr0 <- 0; lka <- 0; lfbio <- fix(0)
          eta_cl ~ 0.1; eta_v1 ~ 0.1; eta_q ~ 0.1; eta_v2 ~ 0.1; eta_kss ~ 0.1; eta_kint ~ 0.1; eta_kdeg ~ 0.1; eta_r0 ~ 0.1; eta_ka ~ 0.1
          prop_sd <- 0.2; add_sd <- 0.1 })
    model({ cl <- exp(lcl + eta_cl); v1 <- exp(lv1 + eta_v1); q <- exp(lq + eta_q); v2 <- exp(lv2 + eta_v2)
            kss <- exp(lkss + eta_kss); kint <- exp(lkint + eta_kint); kdeg <- exp(lkdeg + eta_kdeg); r0 <- exp(lr0 + eta_r0)
            ka <- exp(lka + eta_ka); fbio <- exp(lfbio)
            ctot <- central/v1
            cfree <- 0.5*((ctot - rtot - kss) + sqrt((ctot - rtot - kss)^2 + 4*kss*ctot))
            d/dt(depot) <- -ka*depot
            d/dt(central) <- ka*depot - cl*cfree - q*cfree + q*peripheral/v2 - kint*rtot*cfree/(kss + cfree)*v1
            d/dt(peripheral) <- q*cfree - q*peripheral/v2
            d/dt(rtot) <- kdeg*r0 - kdeg*rtot - (kint - kdeg)*rtot*cfree/(kss + cfree)
            rtot(0) <- r0
            f(depot) <- fbio
            cp <- ctot
            cp ~ prop(prop_sd) + add(add_sd) })
  }
  f
}

build <- function(p) {
  allowed <- PARAMS[[p$model]]; if (is.null(allowed)) stop("Unsupported structural model")
  m <- template(p$model, p$route, p$species)()
  if (p$model == "qss_tmdd_2cmt" && identical(p$species, "free"))
    m <- model(m, cp <- cfree)
  use <- c(allowed, if (p$route == "sc") "ka")
  # Initial estimates and fixed values (log scale), by whitelisted name only.
  for (k in names(p$init)) {
    if (!(k %in% use)) stop("Initial estimate for an unused parameter: ", k)
    th <- as.name(paste0("l", k)); v <- log(p$init[[k]])
    m <- if (k %in% p$fixed) eval(bquote(ini(m, .(th) <- fix(.(v))))) else eval(bquote(ini(m, .(th) <- .(v))))
  }
  if (p$route == "sc") m <- eval(bquote(ini(m, lfbio <- fix(.(log(p$bioavailability))))))
  # Remove random effects that were not declared, and every eta of an unused or fixed parameter.
  all_eta <- c(PARAMS[[p$model]], "ka")
  for (k in all_eta) {
    if (k %in% p$eta && k %in% use && !(k %in% p$fixed)) {
      m <- eval(bquote(ini(m, .(as.name(paste0("eta_", k))) ~ .(p$eta_init))))
    } else {
      m <- eval(bquote(model(m, .(as.name(k)) <- exp(.(as.name(paste0("l", k)))))))
    }
  }
  # IV routes keep the (never dosed) depot with ka fixed; it contributes nothing.
  if (p$route != "sc") m <- ini(m, lka <- fix(0))
  if (p$residual == "proportional") {
    m <- model(m, cp ~ prop(prop_sd)); m <- eval(bquote(ini(m, prop_sd <- .(p$residual_init$proportional))))
  } else {
    m <- eval(bquote(ini(m, prop_sd <- .(p$residual_init$proportional))))
    m <- eval(bquote(ini(m, add_sd <- .(p$residual_init$additive))))
  }
  m
}

run_fit <- function(p) {
  d <- data.frame(ID=as.integer(p$data$id), TIME=as.numeric(p$data$time), DV=as.numeric(p$data$dv),
                  AMT=as.numeric(p$data$amt), EVID=as.integer(p$data$evid), RATE=as.numeric(p$data$rate),
                  CMT=as.character(p$data$cmt), stringsAsFactors=FALSE)
  m <- build(p)
  warnings <- character()
  set.seed(p$seed)
  ctl <- foceiControl(print=0, maxOuterIterations=p$max_outer_iterations, covMethod="r,s", calcTables=TRUE)
  f <- withCallingHandlers(nlmixr2(m, d, est="focei", control=ctl),
    warning=function(w) {warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")},
    message=function(x) invokeRestart("muffleMessage"))
  pf <- f$parFixedDf
  theta <- data.frame(name=rownames(pf), estimate=pf$Estimate, se=pf$SE, rse=pf$`%RSE`, back=pf$`Back-transformed`,
                      ci_low=pf$`CI Lower`, ci_high=pf$`CI Upper`,
                      fixed=as.logical(m$iniDf$fix[match(rownames(pf), m$iniDf$name)]))
  om <- f$omega; cor <- tryCatch(f$cor, error=function(e) NULL)
  sh <- f$shrink
  etas <- grep("^eta_", colnames(om), value=TRUE)
  tab <- as.data.frame(f)
  ind <- intersect(c("ID","TIME","DV","PRED","IPRED","IWRES","RES","IRES", etas, PARAMS[[p$model]], "ka", "fbio"), names(tab))
  list(theta=theta, omega_diagonal=as.list(setNames(diag(om), colnames(om))),
       # Shrinkage is parFixedDf's Shrink(SD)% on the parameter's row (f$shrink holds EBE moments).
       shrinkage_sd_percent=as.list(setNames(vapply(etas, function(e) as.numeric(pf[paste0("l", sub("^eta_", "", e)), "Shrink(SD)%"]), numeric(1)), etas)),
       objective=f$objf, aic=AIC(f), bic=BIC(f), n_obs=nrow(tab), n_subjects=length(unique(tab$ID)),
       cov_method=as.character(f$covMethod), theta_correlation=if (is.null(cor)) NULL else list(names=rownames(cor), matrix=unname(cor)),
       table=as.list(tab[, ind]), message=paste(f$message, collapse="; "), warnings=unique(warnings),
       model_code=paste(deparse(as.function(m)), collapse="\n"), r_environment=versions)
}

result <- tryCatch(if (request$method == "pop_doctor") list(r_environment=versions)
                   else if (request$method == "nlmixr2_fit") run_fit(request$data)
                   else stop("Unsupported R method"),
                   error=function(e) list(error=conditionMessage(e), r_environment=versions))
write_json(result, args[2], digits=I(17), auto_unbox=TRUE, na="null", null="null", dataframe="columns")
