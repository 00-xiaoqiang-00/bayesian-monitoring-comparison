# Resolve shared RMUT sources for a case-specific analysis entry.
rmut_dependency_paths <- function(directory) {
  directory <- normalizePath(directory, winslash = "/", mustWork = TRUE)
  scope <- basename(dirname(directory))
  case <- basename(directory)
  if (!scope %in% c("calibration", "validation", "S8", "rho") ||
      !case %in% c("IMpower", "KEYNOTE598", "OAK", "survival_colon")) {
    stop("Unknown RMUT case directory: ", directory)
  }
  common <- normalizePath(file.path(directory, "..", "..", "common"),
                          winslash = "/", mustWork = TRUE)
  fit <- if (scope == "calibration") "rmut_calibration_fit_helper.R" else "rmut_stable_fit_helper.R"
  runtime <- if (scope %in% c("calibration", "validation")) "production_runtime.R" else "rmut_sensitivity_runtime.R"
  paths <- c(
    "rmut_pairing.R" = file.path(common, "rmut_pairing.R"),
    "rmut_stable_fit_helper.R" = file.path(common, fit),
    "data_simulation_final.R" = file.path(directory, "..", "..", "..", "common", "dgm", paste0(case, "_shared.R")),
    "runtime" = file.path(common, runtime),
    "rmut_source_paths.R" = file.path(common, "rmut_source_paths.R"),
    "source_path.R" = file.path(directory, "..", "..", "..", "common", "source_path.R"),
    "analysis_helpers.R" = file.path(directory, "..", "..", "..", "common", "analysis_helpers.R")
  )
  names(paths)[names(paths) == "runtime"] <- runtime
  if (any(!file.exists(paths))) stop("Missing RMUT dependency: ", paste(paths[!file.exists(paths)], collapse = ", "))
  setNames(normalizePath(unname(paths), winslash = "/", mustWork = TRUE), names(paths))
}
