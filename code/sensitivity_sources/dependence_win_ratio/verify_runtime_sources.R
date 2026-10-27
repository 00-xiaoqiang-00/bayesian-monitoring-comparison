# Check the independently built runtime; source-component provenance is recorded separately.
runtime_required_sources <- function() {
  cases <- c("IMpower131","KEYNOTE598","OAK","survival_colon")
  top <- c("winratio_rho_common.R","validate_winratio_rho.R","preflight_winratio_rho.R",
           "audit_winratio_rho.R","verify_runtime_sources.R","prepare_runtime.R",
           "rho_sensitivity_parameter_helper.R","rho_sensitivity_km_target_parameters.csv",
           "rho_sensitivity_km_target_acceptance.csv","run_winratio_rho.sh",
           "SOURCE_PROVENANCE.json")
  model <- c("win_ratio_final.R","RcppExports.R","RcppExports.cpp","update_theta_V9.cpp",
             "update_theta_univariate_v3.cpp","winratio.cpp","config.R","DESCRIPTION",
             "LOCKED_PARAMETER_MANIFEST.tsv","LICENSE","LICENSE.md",
             paste0("BayesianWinRatio",.Platform$dynlib.ext))
  c(top, unlist(lapply(cases,function(x) file.path("vendor_snapshot",x,model))),
    file.path("reference_dgm","data_simulation",
              c("cox_IMpower.R","cox_KEYNOTE598.R","cox_OAK.R","cox_survival_colon.R")))
}

verify_runtime_sources <- function(root) {
  root <- normalizePath(root,mustWork=TRUE)
  path <- file.path(root,"RUNTIME_SHA256SUMS.tsv")
  if(!file.exists(path)) stop("Runtime not prepared: run prepare_runtime.R with an empty output directory")
  m <- utils::read.delim(path,stringsAsFactors=FALSE,check.names=FALSE)
  if(!identical(names(m),c("path","sha256")) || !nrow(m) || anyNA(m) ||
     anyDuplicated(m$path) || any(!grepl("^[0-9a-f]{64}$",m$sha256)) ||
     any(grepl("(^[/\\\\]|^[A-Za-z]:|(^|[/\\\\])\\.\\.([/\\\\]|$))",m$path))) {
    stop("Malformed runtime source manifest")
  }
  canonical <- function(x) gsub("\\\\","/",x)
  if(anyDuplicated(canonical(m$path)) || !setequal(canonical(m$path),canonical(runtime_required_sources()))) stop("Runtime manifest has missing, duplicate or unexpected entries")
  files <- file.path(root,m$path)
  if(any(!file.exists(files)) || any(file.info(files)$isdir)) stop("Runtime source file missing")
  hashes <- vapply(files,function(p) digest::digest(file=p,algo="sha256"),character(1))
  if(!identical(unname(hashes),m$sha256)) stop("Runtime source hash mismatch")
  invisible(TRUE)
}

if(sys.nframe()==0L) {
  args <- commandArgs(TRUE)
  if(length(args)!=1L) stop("Usage: Rscript verify_runtime_sources.R RUNTIME_DIRECTORY")
  verify_runtime_sources(args[1])
  cat("RUNTIME_SOURCE_HASHES_PASS\n")
}
