# Load a case config and prepare its output root.
#
# Usage:
#   source("code/vendor/bayesian_win_ratio/shared/load_case_config.R")
#   cfg <- load_case_config("survival_colon", repo_root = WR_ROOT)
#   # or parse --case=: cfg <- load_case_config(repo_root = WR_ROOT)

load_case_config <- function(case_id = NULL,
                             repo_root = getwd(),
                             args = commandArgs(TRUE)) {
  repo_root <- normalizePath(repo_root, winslash = "/", mustWork = TRUE)
  if (is.null(case_id) || !nzchar(case_id)) {
    # Support --case=CASE_ID, case=CASE_ID, or the first bare argument.
    m <- regexec("(?:--)?case=([^[:space:]]+)", paste(args, collapse = " "))
    hit <- regmatches(paste(args, collapse = " "), m)[[1]]
    if (length(hit) >= 2) {
      case_id <- hit[2]
    } else if (length(args) >= 1 && !grepl("=", args[1])) {
      case_id <- args[1]
    } else {
      stop("Pass --case=IMpower131|keynote598|OAK_A_ITT|survival_colon.")
    }
  }

  case_id <- switch(case_id, IMpower = "IMpower131", KEYNOTE598 = "keynote598",
                    OAK = "OAK_A_ITT", case_id)
  if (!case_id %in% c("IMpower131", "keynote598", "OAK_A_ITT", "survival_colon")) {
    stop("Unknown Win Ratio case configuration: ", case_id)
  }
  config_path <- file.path(repo_root, "cases", case_id, "config.R")
  if (!file.exists(config_path)) {
    stop("Missing case config: ", config_path)
  }

  env <- new.env(parent = baseenv())
  sys.source(config_path, envir = env)

  if (is.null(env$case_id) || !identical(as.character(env$case_id), as.character(case_id))) {
    warning("config case_id=\"", env$case_id, "\" differs from requested \"", case_id, "\". Using folder name.")
    env$case_id <- case_id
  }

  required_fields <- c(
    "N.max", "interim_cohort", "design", "recruit.int", "cens_upper",
    "final_analysis_time", "m0", "L0", "v0", "S0", "time_max",
    "thin_MCMC", "Niter", "Sigma.ctrl", "Sigma_h0_trt", "Sigma_h1_trt",
    "alpha", "power_target", "L_sim", "mc.cores", "enable_superiority",
    "rebuild_pp_cache", "eta_interim", "eta_final", "c_F_grid", "gamma_E_grid",
    "mu_h0_trt", "mu_h0_ctrl", "mu_h1_trt", "mu_h1_ctrl", "seed",
    "scenarios_extra", "design_version"
  )
  missing_fields <- required_fields[!vapply(required_fields, exists, logical(1),
                                           envir = env, inherits = FALSE)]
  if (length(missing_fields)) {
    stop("Missing required case configuration fields: ", paste(missing_fields, collapse = ", "))
  }
  for (field in setdiff(required_fields, "final_analysis_time")) {
    value <- env[[field]]
    if (is.null(value) || !length(value) || (is.numeric(value) && any(!is.finite(value)))) {
      stop("Invalid required case configuration field: ", field)
    }
  }
  if (!is.list(env$scenarios_extra) || !length(env$scenarios_extra)) {
    stop("Case configuration must explicitly provide its observed scenario.")
  }
  case_root <- file.path(repo_root, "output", "cases", case_id)
  env$repo_root <- repo_root
  env$config_path <- config_path
  env$out_dir <- case_root
  if (!dir.exists(env$out_dir)) {
    dir.create(env$out_dir, recursive = TRUE, showWarnings = FALSE)
  }

  cat("Loaded case:", case_id, "\n")
  cat("  config:", config_path, "\n")
  cat("  output:", case_root, "\n")
  env
}

# Convenience: list available cases
list_cases <- function(repo_root = getwd()) {
  repo_root <- normalizePath(repo_root, winslash = "/", mustWork = TRUE)
  cases_dir <- file.path(repo_root, "cases")
  if (!dir.exists(cases_dir)) return(character(0))
  names <- list.dirs(cases_dir, full.names = FALSE, recursive = FALSE)
  names[file.exists(file.path(cases_dir, names, "config.R"))]
}
