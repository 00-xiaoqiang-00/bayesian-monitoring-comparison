script <- sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[1])
Sys.setenv(RHO_SENSITIVITY_VERSION_DIR=normalizePath(dirname(script),mustWork=TRUE))
source(file.path(Sys.getenv("RHO_SENSITIVITY_VERSION_DIR"),"winratio_rho_common.R"))
wr_verify_sources()
for(f in list.files(wr_root,pattern="\\.R$",full.names=TRUE)) parse(f)
accept <- read.csv(file.path(wr_root,"rho_sensitivity_km_target_acceptance.csv"))
stopifnot(all(tolower(accept$eligible_for_sensitivity_validation)=="true"))
rows <- list();baseline_rows <- list();seed_rows <- list()
for(case in wr_cases) {
  cfg <- wr_cfg(case);gen <- wr_generator(case)
  for(level in c("low","high")) for(scenario in wr_scenarios) {
    p <- wr_parameters(case,scenario,level)
    seed <- rho_data_seed(8252026L,rho_case_name(case),scenario,level)
    set.seed(seed,kind="Mersenne-Twister",normal.kind="Inversion",sample.kind="Rejection")
    d <- wr_generate(gen,case,scenario,p,cfg,2L)
    wr_check_data(d,cfg)
    # Replay/boundary tests verify exact <= for futility and strict > at final.
    b <- cfg$c_F*(cfg$look_cohort[1:3]/cfg$N.max)^cfg$gamma_E
    for(j in 1:3) {
      pp <- rep(1,4);pp[j] <- b[j]
      stopifnot(wr_replay(pp,cfg)$stop_stage==j,wr_replay(pp,cfg)$decision=="futility")
    }
    stopifnot(wr_replay(c(1,1,1,cfg$c_F),cfg)$decision=="no_go",
              wr_replay(c(1,1,1,cfg$c_F+1e-12),cfg)$decision=="go",
              wr_replay(rep(1,4),cfg)$n_analyzed==cfg$N.max)
    # Invalid values must be rejected rather than counted as no-go.
    stopifnot(inherits(try(wr_replay(c(1,1,NA,1),cfg),silent=TRUE),"try-error"))
    idx <- length(rows)+1L
    rows[[idx]] <- data.frame(case=case,scenario=scenario,level=level,rho=p$rho,
      formal_data_seed=seed,first_data_sha256=wr_hash(unname(d[,1:5,1])),
      second_data_sha256=wr_hash(unname(d[,1:5,2])),
      dgm_sha256=wr_hash(p),control_marginal_log_correlation=p$marginal_log_correlation,
      parameter_components=length(p$arms$control$component),status="PASS")
    for(mode in c("formal","smoke","pilot")) {
      offd <- switch(mode,formal=0,smoke=500000000,pilot=600000000)
      model_base <- wr_model_base(case,scenario,level,mode)
      model_seeds <- as.vector(outer(seq_len(if(mode=="formal") 10000 else if(mode=="pilot") 100 else 1),c(1L,2L,3L,9L),function(i,j) wr_seed_add(model_base,i*10+j)))
      seed_rows[[length(seed_rows)+1L]] <- list(data=wr_seed_add(seed,offd),model=model_seeds)
    }
    cat("PREFLIGHT_UNIT_PASS",case,level,scenario,"\n")
  }
}
ds <- vapply(seed_rows,`[[`,integer(1),"data")
ms <- unlist(lapply(seed_rows,`[[`,"model"),use.names=FALSE)
stopifnot(!anyDuplicated(ds),!anyDuplicated(ms),!any(ds %in% ms),all(ds>0),all(ms>0))
write.csv(do.call(rbind,rows),file.path(wr_root,"PREFLIGHT_UNITS.csv"),row.names=FALSE)
writeLines(c("status=PASS","scenario_units=24","generated_trials=48","decision_boundary_tests=PASS",
             "invalid_PP_rejection=PASS","formal_smoke_pilot_seed_collision_check=PASS",
             paste0("code_sha256=",wr_code_hash()),paste0("parameter_table_sha256=",wr_file_hash(rho_parameter_file))),
           file.path(wr_root,"PREFLIGHT_STATUS.txt"))
cat("PREFLIGHT_PASS 24 scenario units; 48 generated trials; boundaries and seed namespaces\n")
