args <- commandArgs(TRUE)
arg <- function(key,default=NULL) {
  hit <- args[startsWith(args,paste0("--",key,"="))]
  if(length(hit)>1L) stop("Duplicate argument: ",key)
  if(!length(hit)) return(default)
  sub(paste0("^--",key,"="),"",hit)
}
script <- sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[1])
Sys.setenv(RHO_SENSITIVITY_VERSION_DIR=normalizePath(dirname(script),mustWork=TRUE))
source(file.path(Sys.getenv("RHO_SENSITIVITY_VERSION_DIR"),"winratio_rho_common.R"))
if(any(!grepl("^--(case|level|mode|n|cores|chunk)=",args))) stop("Unknown argument")
case <- arg("case");level <- arg("level");mode <- arg("mode")
stopifnot(case %in% wr_cases,level %in% c("low","high"),mode %in% c("smoke","pilot","formal"))
integer_arg <- function(x,lo,hi) {
  y <- suppressWarnings(as.numeric(x))
  stopifnot(length(y)==1L,is.finite(y),y==floor(y),y>=lo,y<=hi)
  as.integer(y)
}
n <- integer_arg(arg("n",switch(mode,smoke="1",pilot="20",formal="10000")),1,10000)
if(mode=="formal" && n!=10000L || mode=="smoke" && n!=1L) stop("Mode replicate count is locked")
if(mode=="pilot" && n>100L) stop("Pilot n must be <= 100")
cores <- integer_arg(arg("cores",Sys.getenv("NUM_CORES","1")),1,70)
chunk <- integer_arg(arg("chunk","20"),1,1000)
master <- 8252026L
if(nzchar(Sys.getenv("RHO_SENSITIVITY_MASTER_SEED")) && Sys.getenv("RHO_SENSITIVITY_MASTER_SEED")!=as.character(master))
  stop("Master seed is locked to the formal value 8252026")
wr_verify_sources()
cfg <- wr_cfg(case)
Niter <- if(mode=="smoke") 1000L else cfg$Niter
code_hash <- wr_code_hash();parameter_hash <- wr_file_hash(rho_parameter_file)
unit <- file.path(wr_root,"outputs",mode,case,level)
dir.create(unit,recursive=TRUE,showWarnings=FALSE)
if(!dir.create(file.path(unit,".active_lock"),showWarnings=FALSE)) stop("Unit already active or stale lock requires inspection: ",unit)
writeLines(as.character(Sys.getpid()),file.path(unit,".active_lock","pid"))
unit_run <- function() {
  on.exit(unlink(file.path(unit,".active_lock"),recursive=TRUE),add=TRUE)
  marker <- file.path(unit,"RUN_COMPLETE.txt")
  if(file.exists(marker)) unlink(marker)
  gen <- wr_generator(case);model <- wr_model(case)
  data_offset <- switch(mode,formal=0,smoke=500000000,pilot=600000000)
  model_namespace <- switch(mode,formal=1800000000,smoke=1900000000,pilot=2000000000)
  summaries <- list()
  for(scenario in wr_scenarios) {
    p <- wr_parameters(case,scenario,level)
    data_seed <- wr_seed_add(rho_data_seed(master,rho_case_name(case),scenario,level),data_offset)
    model_base <- wr_model_base(case,scenario,level,mode)
    signature <- list(version="WR_RHO_KM_TARGET",case=case,scenario=scenario,
      level=level,rho=p$rho,mode=mode,n=n,Niter=Niter,cfg=cfg,parameters=p,
      code_sha256=code_hash,parameter_table_sha256=parameter_hash,
      master_seed=master,data_seed=data_seed,data_offset=data_offset,model_namespace=model_namespace,model_base=model_base,
      data_rng=c("Mersenne-Twister","Inversion","Rejection"),model_rng=c("L'Ecuyer-CMRG","Inversion","Rejection"),
      R_version=R.version.string,package_versions=vapply(c("Rcpp","MASS","survival"),function(x) as.character(packageVersion(x)),character(1)))
    signature_hash <- wr_hash(signature)
    scdir <- file.path(unit,scenario);dir.create(scdir,showWarnings=FALSE)
    metadata_path <- file.path(scdir,"metadata.rds")
    if(file.exists(metadata_path) && !identical(readRDS(metadata_path),signature)) stop("Incompatible run configuration; choose a new output directory: ",scdir)
    wr_atomic_rds(signature,metadata_path)
    jsonlite::write_json(signature,file.path(scdir,"RUN_METADATA.json"),pretty=TRUE,auto_unbox=TRUE,digits=17)
    set.seed(data_seed,kind="Mersenne-Twister",normal.kind="Inversion",sample.kind="Rejection")
    data <- wr_generate(gen,case,scenario,p,cfg,n)
    wr_check_data(data,cfg)
    # A common bulk RNG stream follows the locked formal sequence; model seeds never alter it.
    # First-five-column hashes compare latent times, censoring, allocation and ids.
    data_hashes <- vapply(seq_len(n),function(i) wr_hash(unname(data[,1:5,i])),character(1))
    chunks <- split(seq_len(n),ceiling(seq_len(n)/chunk))
    fits <- vector("list",n)
    for(ids in chunks) {
      cp <- file.path(scdir,sprintf("pp_%06d_%06d.rds",min(ids),max(ids)))
      sig <- list(run=signature_hash,ids=as.integer(ids),data_hashes=data_hashes[ids])
      if(file.exists(cp)) {
        saved <- readRDS(cp)
        if(!identical(saved$signature,sig) || length(saved$fits)!=length(ids)) stop("Incompatible checkpoint: ",cp)
        fit <- saved$fits
        cat("REUSE",scenario,min(ids),max(ids),"\n")
      } else {
        fit <- parallel::mclapply(ids,function(i) wr_compute(data[,,i],cfg,p$Sigma_init,model,model_base,i,Niter),
                                 mc.cores=min(cores,length(ids)),mc.set.seed=FALSE)
        if(any(vapply(fit,function(x) !is.list(x)||inherits(x,"try-error"),logical(1)))) stop("Worker failed; no complete checkpoint saved")
        wr_atomic_rds(list(signature=sig,fits=fit),cp)
        cat("SAVED",scenario,min(ids),max(ids),"\n")
      }
      fits[ids] <- fit
    }
    records <- lapply(seq_len(n),function(i) {
      f <- fits[[i]]
      # Invalid fits are never silently treated as no-go or replaced by new seeds.
      expected_valid <- all(is.finite(f$PP)&f$PP>=0&f$PP<=1) && all(is.finite(f$WR)) && !nzchar(f$error)
      stopifnot(identical(f$valid_fit,expected_valid),length(f$PP)==4L)
      decision <- if(f$valid_fit) wr_replay(f$PP,cfg) else list(decision="invalid",go=NA,stop_stage=NA_integer_,n_enrolled=NA_integer_,n_analyzed=NA_integer_)
      data.frame(case=rho_case_name(case),scenario=scenario,rho_level=level,rho=p$rho,mode=mode,
        replicate_id=i,valid_fit=f$valid_fit,PP_40=f$PP[1],PP_60=f$PP[2],PP_80=f$PP[3],PP_final=f$PP[4],
        WR_40=f$WR[1],WR_60=f$WR[2],WR_80=f$WR[3],WR_final=f$WR[4],
        decision=decision$decision,go=decision$go,stop_stage=decision$stop_stage,
        n_enrolled=decision$n_enrolled,n_analyzed=decision$n_analyzed,
        data_seed=data_seed,data_stream_index=i,data_sha256=data_hashes[i],
        model_seed_40=f$model_seeds[1],model_seed_60=f$model_seeds[2],model_seed_80=f$model_seeds[3],model_seed_final=f$model_seeds[4],
        parameter_table_sha256=parameter_hash,code_sha256=code_hash,dgm_sha256=wr_hash(p),run_signature_sha256=signature_hash,
        Niter=Niter,thin_MCMC=cfg$thin_MCMC,warning=f$warning,error=f$error,stringsAsFactors=FALSE)
    })
    rows <- do.call(rbind,records)
    write.csv(rows,file.path(scdir,"replicates.csv"),row.names=FALSE,na="NA")
    nv <- sum(rows$valid_fit);ng <- sum(rows$go,na.rm=TRUE)
    pg <- if(nv) ng/nv else NA_real_;ci <- wr_wilson(ng,nv)
    fut <- if(nv) sum(rows$decision=="futility")/nv else NA_real_
    no <- if(nv) sum(rows$decision=="no_go")/nv else NA_real_
    if(nv) stopifnot(abs(pg+fut+no-1)<1e-12)
    summaries[[scenario]] <- data.frame(case=rho_case_name(case),scenario=scenario,rho_level=level,rho=p$rho,
      metric=if(scenario=="H0") "Type1_error" else if(scenario=="H1-standard") "Power" else "Go_probability",
      n_requested=n,n_valid=nv,n_invalid=n-nv,valid_fit_rate=nv/n,
      P_go=pg,Type1_error=if(scenario=="H0") pg else NA_real_,Power=if(scenario=="H1-standard") pg else NA_real_,
      go_lower=ci[1],go_upper=ci[2],go_MCSE=if(nv) sqrt(pg*(1-pg)/nv) else NA_real_,
      PET_futility=fut,P_final_nogo=no,decision_path_sum=pg+fut+no,
      ESS=if(nv) mean(rows$n_enrolled[rows$valid_fit]) else NA_real_,
      ESS_over_N=if(nv) mean(rows$n_enrolled[rows$valid_fit])/cfg$N.max else NA_real_,
      mean_n_analyzed=if(nv) mean(rows$n_analyzed[rows$valid_fit]) else NA_real_,
      N_max=cfg$N.max,c_F=cfg$c_F,gamma_E=cfg$gamma_E,eta_interim=cfg$eta_interim,eta_final=cfg$eta_final,
      Niter=Niter,thin_MCMC=cfg$thin_MCMC,mode=mode,parameter_table_sha256=parameter_hash,code_sha256=code_hash,
      valid_fit_definition="all_4_PP_finite_in_0_1_and_WR_finite_no_error;not_MCMC_convergence",stringsAsFactors=FALSE)
    write.csv(summaries[[scenario]],file.path(scdir,"summary.csv"),row.names=FALSE)
    rm(data);gc()
    if(nv!=n) stop("Invalid replicate(s) recorded; no completion marker: ",scenario)
  }
  write.csv(do.call(rbind,summaries),file.path(unit,"summary.csv"),row.names=FALSE)
  writeLines(c("status=PASS",paste0("case=",case),paste0("rho_level=",level),paste0("mode=",mode),
    paste0("n_per_scenario=",n),"completed_scenarios=3",paste0("code_sha256=",code_hash),
    paste0("parameter_table_sha256=",parameter_hash),paste0("completed_at=",Sys.time())),marker)
  capture.output(sessionInfo(),file=file.path(unit,"sessionInfo.txt"))
  cat("UNIT_COMPLETE",case,level,mode,"\n")
}
unit_run()
