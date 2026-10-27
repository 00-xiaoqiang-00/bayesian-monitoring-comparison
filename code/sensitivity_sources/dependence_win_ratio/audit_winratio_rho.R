script <- sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[1])
Sys.setenv(RHO_SENSITIVITY_VERSION_DIR=normalizePath(dirname(script),mustWork=TRUE))
source(file.path(Sys.getenv("RHO_SENSITIVITY_VERSION_DIR"),"winratio_rho_common.R"))
mode <- commandArgs(TRUE)
stopifnot(length(mode)==1L,mode %in% c("smoke","pilot","formal"))
wr_verify_sources();all_summaries <- list();units <- 0L
for(case in wr_cases) for(level in c("low","high")) {
  unit <- file.path(wr_root,"outputs",mode,case,level);cfg <- wr_cfg(case)
  stopifnot(file.exists(file.path(unit,"RUN_COMPLETE.txt")))
  marker <- readLines(file.path(unit,"RUN_COMPLETE.txt"))
  stopifnot(paste0("code_sha256=",wr_code_hash()) %in% marker)
  summary <- read.csv(file.path(unit,"summary.csv"))
  stopifnot(nrow(summary)==3L,setequal(summary$scenario,wr_scenarios))
  for(sc in wr_scenarios) {
    d <- read.csv(file.path(unit,sc,"replicates.csv")); s <- summary[summary$scenario==sc,]
    meta <- readRDS(file.path(unit,sc,"metadata.rds"))
    stopifnot(nrow(d)==s$n_requested,nrow(d)==meta$n,all(d$valid_fit),
      identical(d$replicate_id,seq_len(nrow(d))),all(d$mode==mode),
      all(d$code_sha256==wr_code_hash()),all(d$parameter_table_sha256==wr_file_hash(rho_parameter_file)),
      all(d$run_signature_sha256==wr_hash(meta)),all(d$dgm_sha256==wr_hash(wr_parameters(case,sc,level))))
    for(i in seq_len(nrow(d))) {
      pp <- unlist(d[i,c("PP_40","PP_60","PP_80","PP_final")],use.names=FALSE)
      replay <- wr_replay(pp,cfg)
      stopifnot(d$decision[i]==replay$decision,d$go[i]==replay$go,
                d$stop_stage[i]==replay$stop_stage,d$n_enrolled[i]==replay$n_enrolled,d$n_analyzed[i]==replay$n_analyzed)
    }
    stopifnot(s$n_valid==nrow(d),s$n_invalid==0,s$valid_fit_rate==1,
              abs(s$P_go-mean(d$go))<1e-12,abs(s$PET_futility-mean(d$decision=="futility"))<1e-12,
              abs(s$P_final_nogo-mean(d$decision=="no_go"))<1e-12,
              abs(s$ESS-mean(d$n_enrolled))<1e-9,abs(s$decision_path_sum-1)<1e-12)
    if(sc=="H1-observed") stopifnot(is.na(s$Power),s$metric=="Go_probability")
    units <- units+1L
  }
  all_summaries[[length(all_summaries)+1L]] <- summary
}
stopifnot(units==24L)
write.csv(do.call(rbind,all_summaries),file.path(wr_root,paste0(toupper(mode),"_SUMMARY.csv")),row.names=FALSE)
writeLines(c("status=PASS",paste0("mode=",mode),"completed_case_rho_units=8","completed_scenario_units=24",
             "replicate_replay_and_summary_checks=PASS",paste0("code_sha256=",wr_code_hash()),
             "finite_fit_check_does_not_establish_MCMC_convergence=true"),
           file.path(wr_root,paste0(toupper(mode),"_AUDIT_STATUS.txt")))
cat("AUDIT_PASS",mode,"8 case-rho units; 24 scenario units\n")
