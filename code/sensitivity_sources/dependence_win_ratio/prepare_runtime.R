# Build the registered Win Ratio library only. No data generation or fitting.
args <- commandArgs(TRUE)
if(length(args)!=1L || !startsWith(args[1],"--output=")) {
  stop("Usage: Rscript prepare_runtime.R --output=NEW_EMPTY_RUNTIME_DIRECTORY")
}
script <- sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[1])
source_root <- normalizePath(dirname(script),mustWork=TRUE)
output <- sub("^--output=","",args[1])
if(!nzchar(output)) stop("Output directory is empty")
parent <- normalizePath(dirname(output),mustWork=TRUE)
output <- file.path(parent,basename(output))
if(file.exists(output) || dir.exists(output)) stop("Output directory must not already exist")
canonical <- function(x) gsub("\\\\","/",x)
if(startsWith(paste0(canonical(output),"/"),paste0(canonical(source_root),"/"))) stop("Output must be outside the source tree")
dependencies <- c("Rcpp","RcppArmadillo","RcppDist","digest","jsonlite","MASS","survival")
for(p in dependencies) if(!requireNamespace(p,quietly=TRUE)) stop("Dependency missing: ",p)
if(.Platform$OS.type=="windows" && nzchar(Sys.getenv("RTOOLS44_HOME"))) {
  rtools <- normalizePath(Sys.getenv("RTOOLS44_HOME"),mustWork=TRUE)
  tool_paths <- file.path(rtools,c("usr/bin","x86_64-w64-mingw32.static.posix/bin"))
  if(any(!dir.exists(tool_paths))) stop("RTOOLS44_HOME lacks the Rtools44 compiler directories")
  Sys.setenv(PATH=paste(c(tool_paths,Sys.getenv("PATH")),collapse=.Platform$path.sep))
}
if(!nzchar(Sys.which("make"))) stop("make not found; configure the C++ toolchain (Windows: RTOOLS44_HOME)")
includes <- vapply(c("Rcpp","RcppArmadillo","RcppDist"),function(p) {
  normalizePath(system.file("include",package=p),mustWork=TRUE)
},character(1))
compiler_config <- function(x) system2(file.path(R.home("bin"),"R"),c("CMD","config",x),stdout=TRUE,stderr=TRUE)
compiler <- compiler_config("CXX")
receipt <- list(runtime_version="WR_RHO_KM_TARGET_portable_runtime",
                R=R.version.string,platform=R.version$platform,compiler=compiler,
                package_versions=setNames(lapply(dependencies,function(x) as.character(utils::packageVersion(x))),dependencies),
                analysis_executed=FALSE,builds=list())
files <- list.files(source_root,recursive=TRUE,full.names=TRUE,all.files=TRUE,no..=TRUE)
files <- files[!file.info(files)$isdir]
relative <- substring(files,nchar(source_root)+2L)
dir.create(output,recursive=FALSE)
for(i in seq_along(files)) {
  target <- file.path(output,relative[i]);dir.create(dirname(target),recursive=TRUE,showWarnings=FALSE)
  if(!file.copy(files[i],target,overwrite=FALSE)) stop("Cannot copy source: ",files[i])
}
# Materialize the required runtime layout from canonical, byte-identical sources.
rho_shared <- normalizePath(file.path(source_root,"..","..","common","rho"),mustWork=TRUE)
for(name in c("rho_sensitivity_parameter_helper.R","rho_sensitivity_km_target_parameters.csv")) {
  source_file <- file.path(rho_shared,name)
  target <- file.path(output,name)
  if(!file.copy(source_file,target,overwrite=FALSE) ||
     !identical(digest::digest(file=source_file,algo="sha256"),digest::digest(file=target,algo="sha256"))) {
    stop("Cannot materialize shared dependence source: ",name)
  }
}
winratio_root <- normalizePath(file.path(source_root,"..","..","vendor","bayesian_win_ratio"),mustWork=TRUE)
model_package <- normalizePath(file.path(winratio_root,"engine"),mustWork=TRUE)
model_cpp <- c("RcppExports.cpp","update_theta_V9.cpp","update_theta_univariate_v3.cpp","winratio.cpp")
case_ids <- c(IMpower131="IMpower131",KEYNOTE598="keynote598",OAK="OAK_A_ITT",survival_colon="survival_colon")
for(case in c("IMpower131","KEYNOTE598","OAK","survival_colon")) {
  target_root <- file.path(output,"vendor_snapshot",case)
  dir.create(target_root,recursive=TRUE,showWarnings=FALSE)
  source_files <- c(file.path(model_package,"src",model_cpp),
                    file.path(model_package,"R","RcppExports.R"),
                    file.path(source_root,"shared_model",c("win_ratio_final.R","LOCKED_PARAMETER_MANIFEST.tsv")),
                    file.path(model_package,c("DESCRIPTION","LICENSE","LICENSE.md")),
                    file.path(winratio_root,"cases",case_ids[[case]],"config.R"))
  targets <- file.path(target_root,basename(source_files))
  if(any(!file.exists(source_files)) || !all(file.copy(source_files,targets,overwrite=FALSE))) {
    stop("Cannot materialize canonical Win Ratio runtime sources for ",case)
  }
  source_hashes <- vapply(source_files,function(p) digest::digest(file=p,algo="sha256"),character(1))
  target_hashes <- vapply(targets,function(p) digest::digest(file=p,algo="sha256"),character(1))
  if(!identical(unname(source_hashes),unname(target_hashes))) stop("Runtime source copy differs for ",case)
}
dgm_sources <- c(cox_IMpower.R="IMpower_shared.R",cox_KEYNOTE598.R="KEYNOTE598_shared.R",
                 cox_OAK.R="OAK_shared.R",cox_survival_colon.R="survival_colon_shared.R")
dgm_root <- normalizePath(file.path(source_root,"..","..","common","dgm"),mustWork=TRUE)
dgm_target <- file.path(output,"reference_dgm","data_simulation")
dir.create(dgm_target,recursive=TRUE,showWarnings=FALSE)
for(name in names(dgm_sources)) {
  source_file <- file.path(dgm_root,dgm_sources[[name]])
  target <- file.path(dgm_target,name)
  if(!file.copy(source_file,target,overwrite=FALSE) ||
     !identical(digest::digest(file=source_file,algo="sha256"),digest::digest(file=target,algo="sha256"))) {
    stop("Cannot materialize canonical data generator: ",name)
  }
}
source(file.path(output,"verify_runtime_sources.R"))
cpp <- c("RcppExports.cpp","update_theta_V9.cpp","update_theta_univariate_v3.cpp","winratio.cpp")
build_case <- function(case) {
  stage <- tempfile(paste0("wr_build_",case,"_"));dir.create(stage)
  # Compile in a separate temporary directory.
  if(!all(file.copy(file.path(output,"vendor_snapshot",case,cpp),stage))) stop("Cannot copy C++ sources")
  makevars <- file.path(stage,"Makevars.runtime")
  include_flags <- paste(paste0('-I"',canonical(includes),'"'),collapse=" ")
  writeLines(c(paste0("PKG_CPPFLAGS = ",include_flags),
               "PKG_CXXFLAGS = $(SHLIB_OPENMP_CXXFLAGS)",
               "PKG_LIBS = $(SHLIB_OPENMP_CXXFLAGS) $(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS)"),makevars)
  old_dir <- getwd();old_makevars <- Sys.getenv("R_MAKEVARS_USER",unset=NA_character_)
  on.exit({setwd(old_dir);if(is.na(old_makevars)) Sys.unsetenv("R_MAKEVARS_USER") else Sys.setenv(R_MAKEVARS_USER=old_makevars)},add=TRUE)
  setwd(stage);Sys.setenv(R_MAKEVARS_USER=makevars)
  library_name <- paste0("BayesianWinRatio",.Platform$dynlib.ext)
  log_path <- file.path(output,paste0("BUILD_",case,".log"))
  status <- system2(file.path(R.home("bin"),"R"),
                    c("CMD","SHLIB","-o",shQuote(library_name),shQuote(cpp)),
                    stdout=log_path,stderr=log_path)
  if(status!=0L || !file.exists(library_name)) stop("Build failed; inspect ",log_path,"; stage ",stage)
  dll <- dyn.load(normalizePath(library_name,mustWork=TRUE))
  on.exit(try(dyn.unload(dll[["path"]]),silent=TRUE),add=TRUE)
  calls <- getDLLRegisteredRoutines(dll)[[".Call"]]
  expected <- c(`_BayesianWinRatio_update_theta`=9L,
                `_BayesianWinRatio_update_theta_univariate`=9L,
                `_BayesianWinRatio_compare`=8L)
  if(!setequal(names(calls),names(expected)) ||
     any(vapply(calls,function(x) as.integer(x$numParameters),integer(1))[names(expected)]!=expected)) {
    stop("Compiled routine registration mismatch")
  }
  dyn.unload(dll[["path"]])
  target <- file.path(output,"vendor_snapshot",case,library_name)
  if(!file.copy(library_name,target,overwrite=FALSE)) stop("Cannot copy built library")
  list(case=case,binary=canonical(substring(target,nchar(output)+2L)),
       binary_sha256=digest::digest(file=target,algo="sha256"),registered_routines=as.list(expected))
}
for(case in c("IMpower131","KEYNOTE598","OAK","survival_colon")) receipt$builds[[case]] <- build_case(case)
paths <- canonical(runtime_required_sources())
if(any(!file.exists(file.path(output,paths)))) stop("Required runtime file missing")
manifest <- data.frame(path=paths,sha256=vapply(file.path(output,paths),function(p) digest::digest(file=p,algo="sha256"),character(1)),stringsAsFactors=FALSE)
utils::write.table(manifest,file.path(output,"RUNTIME_SHA256SUMS.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
verify_runtime_sources(output)
receipt$runtime_manifest_sha256 <- digest::digest(file=file.path(output,"RUNTIME_SHA256SUMS.tsv"),algo="sha256")
jsonlite::write_json(receipt,file.path(output,"RUNTIME_BUILD_RECEIPT.json"),pretty=TRUE,auto_unbox=TRUE)
cat("BUILD_ONLY_COMPLETE; no analysis executed\n",normalizePath(output),"\n")
