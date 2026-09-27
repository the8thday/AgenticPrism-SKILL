# Fixed bridge entry point. No eval/parse or formula text from user data.
args <- commandArgs(TRUE)
if (length(args) >= 3 && nzchar(args[3])) .libPaths(c(args[3], .libPaths()))
suppressPackageStartupMessages(library(jsonlite))
request <- fromJSON(args[1], simplifyVector=TRUE)
for (p in names(request$pins)) {
  if (length(find.package(p, quiet=TRUE)) == 0 || as.character(packageVersion(p)) != request$pins[[p]])
    stop(paste("Required optional package",p,"version",request$pins[[p]]))
}
suppressPackageStartupMessages(library(mmrm))
versions <- list(R=R.version.string, packages=lapply(names(request$pins), function(p) as.character(packageVersion(p))))
names(versions$packages) <- names(request$pins)
# Include all loaded package versions, including TMB/Matrix that affect numerics.
versions$loaded_packages <- as.list(vapply(loadedNamespaces(), function(p) as.character(packageVersion(p)), character(1)))
run_mmrm <- function(p) {
  d <- data.frame(y=as.numeric(p$y), id=factor(p$ids), visit=factor(p$visits, levels=0:(p$k-1)))
  d$X <- I(as.matrix(p$x))
  form <- if (p$covariance == "unstructured") y ~ 0 + X + us(visit | id) else if (p$covariance == "ar1") y ~ 0 + X + ar1(visit | id) else stop("Unsupported covariance")
  warnings <- character()
  f <- withCallingHandlers(mmrm(form, data=d, reml=TRUE, control=mmrm_control(method="Kenward-Roger",optimizer="BFGS",optimizer_control=list(reltol=1e-12,maxit=3000),accept_singular=FALSE,drop_visit_levels=FALSE)), warning=function(w) {warnings <<- c(warnings,conditionMessage(w));invokeRestart("muffleWarning")})
  if (component(f,"convergence") != 0) stop("MMRM did not converge")
  joint <- df_md(f, as.matrix(p$joint))
  tests <- list(f_statistic=joint$f_stat,df_numerator=joint$num_df,df_denominator=joint$denom_df,p_value=joint$p_val)
  rows <- list()
  if (length(p$contrasts)>0) for(i in seq_len(nrow(p$contrasts))) {
    r <- df_1d(f,p$contrasts[i,]);margin<-qt(1-(1-p$level)/(2*nrow(p$contrasts)),r$df)*r$se
    rows[[i]]<-list(estimate=r$est,standard_error=r$se,df=r$df,statistic=r$t_stat,p_unadjusted=r$p_val,ci_low=r$est-margin,ci_high=r$est+margin)
  }
  list(beta=unname(coef(f)),beta_covariance=unname(component(f,"beta_vcov")),covariance_matrix=unname(component(f,"varcor")),tests=tests,contrasts=rows,
       warnings=warnings,r_environment=versions)
}
result <- tryCatch(if(request$method=="doctor") list(r_environment=versions) else if(request$method=="mmrm") run_mmrm(request$data) else stop("Unsupported R method"),error=function(e)list(error=conditionMessage(e),r_environment=versions))
write_json(result,args[2],digits=17,auto_unbox=TRUE,na="null",null="null")
