# Win Ratio rho sensitivity: common DGM and locked independent analysis.
wr_root <- normalizePath(Sys.getenv("RHO_SENSITIVITY_VERSION_DIR"), mustWork=TRUE)
source(file.path(wr_root, "rho_sensitivity_parameter_helper.R"))
wr_cases <- c("IMpower131", "KEYNOTE598", "OAK", "survival_colon")
wr_scenarios <- c("H0", "H1-standard", "H1-observed")
wr_grid <- list(IMpower131=c(low=.3,high=.7), KEYNOTE598=c(low=.3,high=.7),
                OAK=c(low=.5,high=.9), survival_colon=c(low=.7,high=.98))
wr_hash <- function(x) digest::digest(x, algo="sha256", serialize=TRUE)
wr_file_hash <- function(p) digest::digest(file=p, algo="sha256")
wr_seed_add <- function(x, y) as.integer((as.double(x)+as.double(y)) %% 2147483646)
wr_model_base <- function(case,scenario,level,mode) {
  # Disjoint one-million-seed blocks isolate each case, dependence level and
  # scenario while accommodating the replicate-by-look seed increments.
  unit <- (match(case,wr_cases)-1L)*6L+(match(level,c("low","high"))-1L)*3L+match(scenario,wr_scenarios)-1L
  as.integer(switch(mode,formal=1800000000,smoke=1900000000,pilot=2000000000)+unit*1000000)
}
wr_key <- function(s) c(H0="H0", `H1-standard`="H1_calibration", `H1-observed`="H1_observed")[[s]]
wr_cfg <- function(case) {
  stopifnot(case %in% wr_cases)
  e <- new.env(parent=baseenv())
  sys.source(file.path(wr_root,"vendor_snapshot",case,"config.R"), e)
  cfg <- as.list(e)
  manifest <- read.delim(file.path(wr_root,"vendor_snapshot",case,"LOCKED_PARAMETER_MANIFEST.tsv"))
  m <- manifest[manifest$case_key == case,,drop=FALSE]
  expected <- list(IMpower131=c(.875,3.05),KEYNOTE598=c(.865,1.85),OAK=c(.875,3.4),survival_colon=c(.850,1.862))[[case]]
  stopifnot(nrow(m)==1, identical(as.numeric(c(m$c_F,m$gamma_E)),expected),
            m$N_max==cfg$N.max, m$interim_counts==paste(cfg$interim_cohort,collapse="/"),
            cfg$design==1, !cfg$enable_superiority, cfg$eta_interim==1, cfg$eta_final==1,
            cfg$Niter==10000, cfg$thin_MCMC==5)
  keep <- c("N.max","interim_cohort","recruit.int","cens_upper","final_analysis_time",
            "m0","L0","v0","S0","time_max","thin_MCMC","Niter","eta_interim","eta_final")
  out <- cfg[keep]
  out$c_F <- m$c_F; out$gamma_E <- m$gamma_E
  out$balanced_randomization <- case %in% c("OAK","KEYNOTE598")
  out$look_cohort <- c(cfg$interim_cohort,cfg$N.max)
  out$analysis_row_rule <- "interim_n_minus_1_final_all_n_v2"
  out
}
wr_parameters <- function(case, scenario, level) {
  stopifnot(case %in% wr_cases, scenario %in% wr_scenarios,level %in% c("low","high"))
  arm <- lapply(c("control","treatment"), function(a) rho_get_arm(case,scenario,a,level))
  names(arm) <- c("control","treatment")
  for (a in arm) {
    stopifnot(all(is.finite(unlist(a))), all(a$weight>0),
              a$rho == wr_grid[[case]][[level]],
              length(a$component)==if(case %in% c("OAK","KEYNOTE598")) 2L else 1L)
    for (k in seq_along(a$component)) {
      S <- matrix(c(a$sigma_log_R[k]^2,a$covariance[k],a$covariance[k],a$sigma_log_D[k]^2),2)
      stopifnot(min(eigen(S,symmetric=TRUE,only.values=TRUE)$values)>0)
    }
    if (case=="OAK") stopifnot(length(unique(a$mu_log_D))==1L,length(unique(a$sigma_log_D))==1L)
  }
  if(scenario=="H0") stopifnot(identical(arm$control,arm$treatment))
  # Total log-time covariance includes between-component covariance.
  # Used as the locked Win Ratio sampler's fallback initialization.
  a <- arm$control
  means <- cbind(a$mu_log_R,a$mu_log_D)
  mu <- colSums(means*a$weight)
  S <- matrix(0,2,2)
  for(k in seq_along(a$weight)) {
    delta <- means[k,]-mu
    S <- S+a$weight[k]*(matrix(c(a$sigma_log_R[k]^2,a$covariance[k],a$covariance[k],a$sigma_log_D[k]^2),2)+tcrossprod(delta))
  }
  stopifnot(min(eigen(S,symmetric=TRUE,only.values=TRUE)$values)>0)
  list(arms=arm, bundle=rho_make_generator_parameters_by_scenario(case,level)[[wr_key(scenario)]],
       Sigma_init=S, marginal_log_correlation=cov2cor(S)[1,2],rho=a$rho)
}
wr_generator <- function(case) {
  f <- c(IMpower131="cox_IMpower.R",KEYNOTE598="cox_KEYNOTE598.R",OAK="cox_OAK.R",survival_colon="cox_survival_colon.R")[[case]]
  e <- new.env(parent=globalenv())
  sys.source(file.path(wr_root,"reference_dgm","data_simulation",f),e)
  e$data.simulation
}
wr_generate <- function(gen, case, scenario, p, cfg, n) {
  args <- list(N.sim=as.integer(n),N.max=cfg$N.max,cens_upper=cfg$cens_upper)
  if(case=="IMpower131") {
    args <- c(args,list(mu.trt=p$bundle$treatment$mu,Sigma.trt=p$bundle$treatment$Sigma,
                       mu.ctrl=p$bundle$control$mu,Sigma.ctrl=p$bundle$control$Sigma))
  } else args <- c(args,list(scenario=wr_key(scenario), parameters=p$bundle,
                            balanced_randomization=cfg$balanced_randomization))
  do.call(gen,args)
}
wr_check_data <- function(d,cfg) {
  stopifnot(is.array(d),dim(d)[1]==cfg$N.max,dim(d)[2]==7L,
            all(is.finite(d[,1:5,])),all(d[,1:3,]>0),all(d[,3,]<=cfg$cens_upper),
            all(d[,4,] %in% 0:1), all(d[,5,]==seq_len(cfg$N.max)))
  if(cfg$balanced_randomization) for(i in seq_len(dim(d)[3])) stopifnot(sum(d[,4,i])==cfg$N.max/2)
  invisible(TRUE)
}
wr_model <- function(case) {
  requireNamespace("Rcpp"); requireNamespace("MASS"); requireNamespace("survival")
  e <- new.env(parent=globalenv())
  directory <- file.path(wr_root,"vendor_snapshot",case)
  dll <- dyn.load(file.path(directory,paste0("BayesianWinRatio",.Platform$dynlib.ext)))
  routines <- getDLLRegisteredRoutines(dll)[[".Call"]]
  for(nm in names(routines)) assign(nm,routines[[nm]],e)
  sys.source(file.path(directory,"RcppExports.R"),e)
  sys.source(file.path(directory,"win_ratio_final.R"),e)
  e$winratio
}
wr_compute <- function(d,cfg,Sigma_init,model,model_base,id,Niter) {
  pp <- rep(NA_real_,4); wr <- rep(NA_real_,4); errors <- character(); warnings <- character()
  seeds <- vapply(c(1L,2L,3L,9L),function(j) wr_seed_add(model_base,id*10+j),integer(1))
  for(j in 1:4) {
    n <- cfg$look_cohort[j]; final <- j==4
    rows <- seq_len(if(final) n else n-1L)
    entry <- cfg$recruit.int*(rows-1)
    now <- if(final) cfg$final_analysis_time else cfg$recruit.int*(n-1)
    stopifnot(now>=max(entry),length(rows)==if(final) cfg$N.max else n-1)
    set.seed(seeds[j],kind="L'Ecuyer-CMRG",normal.kind="Inversion",sample.kind="Rejection")
    ans <- tryCatch(withCallingHandlers(model(
      currentdd=d[rows,,drop=FALSE], n_current=n, Time_current=now, Time_entry=entry,
      m0=cfg$m0,L0=cfg$L0,v0=cfg$v0,S0=cfg$S0,time_max=cfg$time_max,
      eta=if(final) cfg$eta_final else cfg$eta_interim,lambda=0,N.max=cfg$N.max,
      thin_MCMC=cfg$thin_MCMC,Niter=Niter,Sigma.ctrl=Sigma_init),
      warning=function(w) {warnings <<- c(warnings,paste(j,conditionMessage(w),sep=":"));invokeRestart("muffleWarning")}),
      error=function(e) {errors <<- c(errors,paste(j,conditionMessage(e),sep=":"));NULL})
    if(!is.null(ans) && length(ans$probs)==1 && length(ans$WR)==1) {pp[j]<-ans$probs;wr[j]<-ans$WR}
  }
  valid <- all(is.finite(pp) & pp>=0 & pp<=1) && all(is.finite(wr)) && !length(errors)
  list(PP=pp,WR=wr,valid_fit=valid,error=paste(errors,collapse=" | "),warning=paste(unique(warnings),collapse=" | "),model_seeds=seeds)
}
wr_replay <- function(pp,cfg) {
  stopifnot(length(pp)==4,all(is.finite(pp)),all(pp>=0 & pp<=1))
  for(j in 1:3) if(pp[j] <= cfg$c_F*(cfg$look_cohort[j]/cfg$N.max)^cfg$gamma_E)
    return(list(decision="futility",go=FALSE,stop_stage=j,n_enrolled=cfg$look_cohort[j],n_analyzed=cfg$look_cohort[j]-1L))
  list(decision=if(pp[4]>cfg$c_F) "go" else "no_go",go=pp[4]>cfg$c_F,
       stop_stage=4L,n_enrolled=cfg$N.max,n_analyzed=cfg$N.max)
}
wr_atomic_rds <- function(x,p) {
  tmp <- paste0(p,".tmp_",Sys.getpid());saveRDS(x,tmp)
  if(!file.rename(tmp,p)) stop("Atomic checkpoint rename failed: ",p)
}
wr_wilson <- function(k,n) {
  if(n==0) return(c(NA_real_,NA_real_))
  z <- qnorm(.975);p <- k/n;den <- 1+z*z/n
  center <- (p+z*z/(2*n))/den; half <- z*sqrt(p*(1-p)/n+z*z/(4*n*n))/den
  c(center-half,center+half)
}
wr_verify_sources <- function() {
  helper <- new.env(parent=baseenv())
  sys.source(file.path(wr_root,"verify_runtime_sources.R"),helper)
  helper$verify_runtime_sources(wr_root)
  invisible(TRUE)
}
wr_code_hash <- function() {
  files <- c("winratio_rho_common.R","validate_winratio_rho.R","preflight_winratio_rho.R",
             "audit_winratio_rho.R","RUNTIME_SHA256SUMS.tsv","verify_runtime_sources.R")
  wr_hash(setNames(vapply(file.path(wr_root,files),wr_file_hash,character(1)),files))
}
