#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
mode_hit <- args[grepl("^--mode=", args)]
run_mode <- if (length(mode_hit) == 1L) {
  sub("^--mode=", "", mode_hit[1L])
} else {
  "formal_10000"
}

script_hit <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_hit) != 1L) stop("Run this entry with Rscript.")
sens_root <- dirname(normalizePath(sub("^--file=", "", script_hit[1L]), winslash = "/", mustWork = TRUE))
source(file.path(sens_root, "R", "interim_sensitivity_core.R"))

cases <- c("IMpower131", "KEYNOTE598", "OAK", "survival_colon")
oc_list <- list()
delta_list <- list()
for (case_key in cases) {
  oc_path <- file.path(sens_root, "outputs", run_mode, case_key, "results", "oc_latest.csv")
  delta_path <- file.path(sens_root, "outputs", run_mode, case_key, "results",
                          "paired_delta_latest.csv")
  if (file.exists(oc_path)) {
    oc_list[[case_key]] <- utils::read.csv(oc_path, stringsAsFactors = FALSE)
  }
  if (file.exists(delta_path)) {
    delta_list[[case_key]] <- utils::read.csv(delta_path, stringsAsFactors = FALSE)
  }
}

out_dir <- file.path(sens_root, "outputs", run_mode, "_summary")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

if (length(oc_list)) {
  oc <- do.call(rbind, oc_list)
  utils::write.csv(oc, file.path(out_dir, "oc_all_cases.csv"), row.names = FALSE)
  print(oc)
}
if (length(delta_list)) {
  delta <- do.call(rbind, delta_list)
  delta$primary_analysis <- delta$scenario != "H1_observed"
  delta$note <- ifelse(
    delta$scenario == "H1_observed",
    "descriptive; do not retune c_F/gamma_E",
    "primary locked-design sensitivity"
  )
  utils::write.csv(delta, file.path(out_dir, "paired_delta_all_cases.csv"),
                   row.names = FALSE)
  print(delta)
}

cat("Summary written to ", out_dir, "\n", sep = "")
