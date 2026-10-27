# ANALYSIS ROW RULE: interim n-1; final all N.max.
# ============================================================================
# Grid calibration using posterior-probability paths and the LG decision rule.
#
# Decision rules:
#   Interim futility: P(WR > eta_interim | data) <= c_F * (n/N.max)^gamma_E
#   Final go/no-go:   P(WR > eta_final   | data) >  c_F
#
# MCMC posterior probabilities are computed at every look and saved to
# calib_stats_H0.rds / calib_stats_H1.rds. The (c_F, gamma_E) grid replays
# the decision rule on these cached probabilities.
# This is equivalent when eta is fixed and boundaries depend only on (c_F, gamma_E).
#
# Candidate tables report H0 go, H1-standard power, PET and ESS.
# Formal validation reads the supplied locked parameter manifest.
# ============================================================================

repo_override <- Sys.getenv("BWR_REPO_ROOT", unset = "")
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (nzchar(repo_override)) {
  repo_root <- normalizePath(repo_override, winslash = "/", mustWork = TRUE)
} else if (length(script_arg)) {
  script_file <- normalizePath(sub("^--file=", "", script_arg[1]),
                               winslash = "/", mustWork = TRUE)
  repo_root <- normalizePath(file.path(dirname(script_file), ".."),
                             winslash = "/", mustWork = TRUE)
} else {
  repo_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}

pkgload::load_all(file.path(repo_root, "engine"), quiet = TRUE)
cat("load_all success\n")

library(parallel)
RNGkind("L'Ecuyer-CMRG")
source(file.path(repo_root, "shared", "checkpoint_pp_cache.R"))

# ---------------------------------------------------------------------------
# 1) Case configuration
# ---------------------------------------------------------------------------
case_args <- commandArgs(trailingOnly = TRUE)
case_hit <- case_args[grepl("^--case=", case_args)]
if (length(case_hit) != 1L) stop("Pass exactly one --case=CASE_ID.")
case_id <- sub("^--case=", "", case_hit[1L])
source(file.path(repo_root, "shared", "load_case_config.R"))
cfg <- load_case_config(case_id, repo_root = repo_root)
case_id <- cfg$case_id
is_impower_case <- identical(case_id, "IMpower131")

N.max <- cfg$N.max
interim_cohort <- cfg$interim_cohort
design <- cfg$design
recruit.int <- cfg$recruit.int
cens_upper <- cfg$cens_upper
final_analysis_time <- cfg$final_analysis_time
m0 <- cfg$m0
L0 <- cfg$L0
v0 <- cfg$v0
S0 <- cfg$S0
time_max <- cfg$time_max
thin_MCMC <- cfg$thin_MCMC
Niter <- cfg$Niter
Sigma.ctrl <- cfg$Sigma.ctrl
Sigma_h0_trt <- cfg$Sigma_h0_trt
Sigma_h1_trt <- cfg$Sigma_h1_trt
alpha <- cfg$alpha
power_target <- cfg$power_target
L_sim <- cfg$L_sim
mc.cores <- cfg$mc.cores
enable_superiority <- cfg$enable_superiority
rebuild_pp_cache <- cfg$rebuild_pp_cache
eta_interim <- cfg$eta_interim
eta_final_grid <- c(cfg$eta_final)
c_F_grid <- cfg$c_F_grid
gamma_E_grid <- cfg$gamma_E_grid
mu_h0_trt <- cfg$mu_h0_trt
mu_h0_ctrl <- cfg$mu_h0_ctrl
mu_h1_trt <- cfg$mu_h1_trt
mu_h1_ctrl <- cfg$mu_h1_ctrl
seed <- cfg$seed
if (design != 1) {
  stop("PP-cache calibration currently supports design = 1 only.")
}
if (length(eta_final_grid) != 1) {
  stop("PP-cache path currently supports a single eta_final; got: ",
       paste(eta_final_grid, collapse = ", "))
}
eta_final <- eta_final_grid[1]

set.seed(seed)
run_id <- format(Sys.time(), "%Y%m%d_%H%M%S")

output_version <- as.character(cfg$design_version)
out_dir <- file.path(cfg$out_dir, "design_versions", output_version)
cache_dir <- file.path(out_dir, "cache")
latest_dir <- file.path(out_dir, "latest")
run_dir <- file.path(out_dir, "runs", run_id)
for (d in c(out_dir, cache_dir, latest_dir, run_dir)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

look_cohort <- sort(unique(as.integer(c(interim_cohort, N.max))))
look_cohort <- look_cohort[look_cohort > 0 & look_cohort <= N.max]
n_looks <- length(look_cohort)
analysis_row_rule <- "interim_n_minus_1_final_all_n_v2"
cat("Looks:", paste(look_cohort, collapse = ", "), "\n")

generation_fp_h0 <- NULL
generation_fp_h1 <- NULL
if (is_impower_case) {
  data_h0 <- data.simulation(
    N.sim = L_sim, N.max = N.max, mu.trt = mu_h0_trt, Sigma.trt = Sigma_h0_trt,
    mu.ctrl = mu_h0_ctrl, Sigma.ctrl = Sigma.ctrl, cens_upper = cens_upper
  )
  data_h1 <- data.simulation(
    N.sim = L_sim, N.max = N.max, mu.trt = mu_h1_trt, Sigma.trt = Sigma_h1_trt,
    mu.ctrl = mu_h1_ctrl, Sigma.ctrl = Sigma.ctrl, cens_upper = cens_upper
  )
} else {
  generation_fp_h0 <- case.generator.fingerprint(
    cfg = cfg, scenario = "H0",
    mu.trt = mu_h0_trt, Sigma.trt = Sigma_h0_trt,
    mu.ctrl = mu_h0_ctrl, Sigma.ctrl = Sigma.ctrl,
    cens_upper = cens_upper
  )
  generation_fp_h1 <- case.generator.fingerprint(
    cfg = cfg, scenario = "H1_calibration",
    mu.trt = mu_h1_trt, Sigma.trt = Sigma_h1_trt,
    mu.ctrl = mu_h1_ctrl, Sigma.ctrl = Sigma.ctrl,
    cens_upper = cens_upper
  )

  data_h0 <- case.data.simulation(
    N.sim = L_sim, N.max = N.max, mu.trt = mu_h0_trt, Sigma.trt = Sigma_h0_trt,
    mu.ctrl = mu_h0_ctrl, Sigma.ctrl = Sigma.ctrl, cens_upper = cens_upper,
    cfg = cfg, scenario = "H0"
  )
  data_h1 <- case.data.simulation(
    N.sim = L_sim, N.max = N.max, mu.trt = mu_h1_trt, Sigma.trt = Sigma_h1_trt,
    mu.ctrl = mu_h1_ctrl, Sigma.ctrl = Sigma.ctrl, cens_upper = cens_upper,
    cfg = cfg, scenario = "H1_calibration"
  )
}

# ---------------------------------------------------------------------------
# 2) Build / load PP caches for each scenario
# ---------------------------------------------------------------------------
# At interim looks the newest enrollee has zero follow-up and is omitted.
# At the prespecified final calendar cutoff all N.max participants are used.
compute_pp_one_trial <- function(myData) {
  pp <- numeric(n_looks)
  names(pp) <- as.character(look_cohort)
  for (j in seq_len(n_looks)) {
    n.current <- look_cohort[j]
    is_final_look <- (n.current == N.max)
    Time.entry <- recruit.int * (seq_len(n.current) - 1)
    enrollment_time <- recruit.int * (n.current - 1)
    Time.current <- if (is_final_look && !is.null(final_analysis_time)) {
      final_analysis_time
    } else {
      enrollment_time
    }
    if (Time.current < enrollment_time) {
      stop("Analysis time precedes current-cohort enrollment.")
    }
    analysis_rows <- if (is_final_look) {
      seq_len(n.current)
    } else {
      seq_len(n.current - 1L)
    }
    temp <- myData[analysis_rows, , drop = FALSE]
    Time.entry.analysis <- Time.entry[analysis_rows]
    eta_use <- if (is_final_look) eta_final else eta_interim

    result <- winratio(
      currentdd = temp,
      n_current = n.current,
      Time_current = Time.current,
      Time_entry = Time.entry.analysis,
      m0 = m0,
      L0 = L0,
      v0 = v0,
      S0 = S0,
      time_max = time_max,
      eta = eta_use,
      lambda = 0,
      N.max = N.max,
      thin_MCMC = thin_MCMC,
      Niter = Niter,
      Sigma.ctrl = Sigma.ctrl
    )
    pp[j] <- as.numeric(result$probs)
  }
  pp
}

build_pp_cache <- function(data_arr, nsim, label, generation_fingerprint = NULL) {
  cat(sprintf("\nBuilding PP cache for %s (%d trials x %d looks)...\n",
              label, nsim, n_looks))
  t0 <- Sys.time()
  pp_mat <- if (is_impower_case) {
    build_pp_matrix_checkpointed(data_arr, nsim, label)
  } else {
    build_pp_matrix_checkpointed(data_arr, nsim, label, generation_fingerprint)
  }
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  cat(sprintf("  %s cache done in %.1f min\n", label, elapsed))
  cache <- list(
    PP = pp_mat,
    look_cohort = look_cohort,
    N.max = N.max,
    eta_interim = eta_interim,
    eta_final = eta_final,
    Niter = Niter,
    thin_MCMC = thin_MCMC,
    recruit.int = recruit.int,
    cens_upper = cens_upper,
    final_analysis_time = final_analysis_time,
    analysis_row_rule = analysis_row_rule,
    L_sim = nsim,
    label = label,
    generation_fingerprint = generation_fingerprint,
    created_at = as.character(Sys.time())
  )
  if (is_impower_case) cache$generation_fingerprint <- NULL
  cache
}

cache_h0_path <- file.path(cache_dir, "calib_stats_H0.rds")
cache_h1_path <- file.path(cache_dir, "calib_stats_H1.rds")

cache_basic_compatible <- function(cache, expected_generation_fingerprint = NULL) {
  if (is.null(cache$PP)) return(FALSE)
  if (!identical(as.integer(cache$look_cohort), as.integer(look_cohort))) return(FALSE)
  if (!isTRUE(cache$N.max == N.max)) return(FALSE)
  if (!isTRUE(cache$eta_interim == eta_interim) || !isTRUE(cache$eta_final == eta_final)) return(FALSE)
  if (!isTRUE(cache$Niter == Niter) || !isTRUE(cache$thin_MCMC == thin_MCMC)) return(FALSE)
  if (!isTRUE(all.equal(cache$recruit.int, recruit.int))) return(FALSE)
  if (!isTRUE(nrow(cache$PP) == L_sim)) return(FALSE)
  if (!isTRUE(all.equal(cache$final_analysis_time, final_analysis_time))) return(FALSE)
  if (!identical(cache$analysis_row_rule, analysis_row_rule)) return(FALSE)
  # When cens_upper metadata is absent, use the fixed value 41.
  cache_cens <- if (is.null(cache$cens_upper)) 41 else cache$cens_upper
  if (!isTRUE(all.equal(cache_cens, cens_upper))) return(FALSE)
  if (is_impower_case) return(TRUE)
  cache_gen <- cache$generation_fingerprint
  if (identical(case_id, "survival_colon") && is.null(cache_gen) &&
      !is.null(cache$fingerprint$generation)) {
    cache_gen <- cache$fingerprint$generation
  }
  if (is.null(cache_gen)) {
    if (identical(case_id, "OAK_A_ITT")) {
      if (identical(expected_generation_fingerprint$generator_type,
                    "oak_lognormal_mixture_v1")) return(FALSE)
    } else if (!identical(expected_generation_fingerprint$generator_type,
                          "bivariate_lognormal_v1")) return(FALSE)
  } else if (!isTRUE(all.equal(
    cache_gen, expected_generation_fingerprint,
    tolerance = 1e-10, check.attributes = FALSE
  ))) {
    return(FALSE)
  }
  TRUE
}

# Per-scenario: reuse compatible cache; build only if missing/incompatible
# (or if rebuild_pp_cache=TRUE forces overwrite).
load_or_build_cache <- function(path, data_arr, label, generation_fingerprint) {
  run_snap <- file.path(run_dir, paste0("calib_stats_", label, "_", run_id, ".rds"))
  if (!rebuild_pp_cache && file.exists(path)) {
    cache <- readRDS(path)
    if (cache_basic_compatible(cache, generation_fingerprint)) {
      cat(sprintf("Reusing %s cache: %s\n", label, path))
      saveRDS(cache, run_snap)
      return(cache)
    }
    cat(sprintf("%s cache incompatible, rebuilding: %s\n", label, path))
  } else if (file.exists(path) && rebuild_pp_cache) {
    cat(sprintf("rebuild_pp_cache=TRUE, rebuilding %s: %s\n", label, path))
  } else {
    cat(sprintf("No %s cache, building: %s\n", label, path))
  }
  cache <- build_pp_cache(data_arr, L_sim, label, generation_fingerprint)
  atomic_save_rds(cache, path)
  atomic_save_rds(cache, run_snap)
  cat("Wrote:", path, "\n")
  cache
}

cache_h0 <- load_or_build_cache(
  cache_h0_path, data_h0, "H0", generation_fp_h0
)
cache_h1 <- load_or_build_cache(
  cache_h1_path, data_h1, "H1", generation_fp_h1
)

pp_h0 <- cache_h0$PP
pp_h1 <- cache_h1$PP

# ---------------------------------------------------------------------------
# 3) Replay the decision rule on cached PP without MCMC
# ---------------------------------------------------------------------------
replay_one <- function(pp_row, c_F_value, gamma_E) {
  decision <- NA_character_
  go <- NA
  n_final <- NA_real_
  for (j in seq_len(n_looks)) {
    n.current <- look_cohort[j]
    is_final_look <- (n.current == N.max)
    probs <- as.numeric(pp_row[j])
    ratio <- n.current / N.max
    fut_bound <- c_F_value * (ratio^gamma_E)
    sup_bound <- if (isTRUE(enable_superiority)) 1 - ratio^gamma_E else Inf

    if (!is_final_look) {
      if (!is.na(probs) && probs <= fut_bound) {
        decision <- "futility"
        go <- FALSE
        n_final <- n.current
        break
      } else if (!is.na(probs) && probs >= sup_bound) {
        decision <- "superiority"
        go <- TRUE
        n_final <- n.current
        break
      }
    } else {
      n_final <- n.current
      if (!is.na(probs) && probs > c_F_value) {
        decision <- "go"
        go <- TRUE
      } else {
        decision <- "no_go"
        go <- FALSE
      }
      break
    }
  }
  list(go = isTRUE(go), decision = decision, pts.stop = n_final)
}

oc_from_cache <- function(pp_mat, c_F_value, gamma_E) {
  nsim <- nrow(pp_mat)
  go_vec <- logical(nsim)
  fut_vec <- logical(nsim)
  sup_vec <- logical(nsim)
  ess_vec <- numeric(nsim)
  for (i in seq_len(nsim)) {
    r <- replay_one(pp_mat[i, ], c_F_value, gamma_E)
    go_vec[i] <- r$go
    fut_vec[i] <- identical(r$decision, "futility")
    sup_vec[i] <- identical(r$decision, "superiority")
    ess_vec[i] <- r$pts.stop
  }
  list(
    P_go            = mean(go_vec, na.rm = TRUE),
    PET_futility    = mean(fut_vec, na.rm = TRUE),
    P_early_success = mean(sup_vec, na.rm = TRUE),
    ESS             = mean(ess_vec, na.rm = TRUE)
  )
}

eval_cell <- function(c_F_value, gamma_E, phase = "grid") {
  h0_oc <- oc_from_cache(pp_h0, c_F_value, gamma_E)
  h1_oc <- oc_from_cache(pp_h1, c_F_value, gamma_E)
  data.frame(
    c_F              = c_F_value,
    eta_interim      = eta_interim,
    eta_final        = eta_final,
    gamma_E          = gamma_E,
    P_H0_go          = h0_oc$P_go,
    PET_H0           = h0_oc$PET_futility,
    ESS_H0           = h0_oc$ESS,
    Power_H1         = h1_oc$P_go,
    ESS_H1           = h1_oc$ESS,
    calibration_phase = phase,
    stringsAsFactors = FALSE
  )
}

cat(sprintf(
  "\nReplaying grid: %d c_F x %d gamma_E = %d cells (no MCMC)...\n",
  length(c_F_grid), length(gamma_E_grid), length(c_F_grid) * length(gamma_E_grid)
))
t_grid0 <- Sys.time()
all_rows <- list()
k <- 0L
for (cf in c_F_grid) {
  cat(sprintf("\n=== c_F=%.3f ===\n", cf))
  for (g in gamma_E_grid) {
    k <- k + 1L
    if (k %% 20L == 1L || g == max(gamma_E_grid)) {
      cat(sprintf("  [grid][c_F=%.3f][gamma_E=%.3f]\n", cf, g))
    }
    all_rows[[k]] <- eval_cell(cf, g, "grid")
  }
}
all_results_df <- do.call(rbind, all_rows)
cat(sprintf(
  "Grid replay done in %.1f sec\n",
  as.numeric(difftime(Sys.time(), t_grid0, units = "secs"))
))

grid_results_df <- all_results_df
cat("\nFull calibration table (cached-PP grid):\n")
print(all_results_df)

# ---------------------------------------------------------------------------
# 4) Candidate summaries
# ---------------------------------------------------------------------------
all_results_df$run_id <- rep(run_id, nrow(all_results_df))
all_results_df$seed <- rep(seed, nrow(all_results_df))
type1_feasible <- all_results_df[all_results_df$P_H0_go <= alpha, , drop = FALSE]
power_feasible <- type1_feasible[type1_feasible$Power_H1 >= power_target, , drop = FALSE]

write.csv(all_results_df, file.path(run_dir, paste0("calibration_LG_table_", run_id, ".csv")), row.names = FALSE)
write.csv(type1_feasible, file.path(run_dir, paste0("calibration_LG_type1_feasible_", run_id, ".csv")), row.names = FALSE)
write.csv(power_feasible, file.path(run_dir, paste0("calibration_LG_power_feasible_", run_id, ".csv")), row.names = FALSE)
write.csv(all_results_df, file.path(latest_dir, "calibration_LG_table.csv"), row.names = FALSE)
write.csv(type1_feasible, file.path(latest_dir, "calibration_LG_type1_feasible.csv"), row.names = FALSE)
write.csv(power_feasible, file.path(latest_dir, "calibration_LG_power_feasible.csv"), row.names = FALSE)

result_obj <- list(
  results = all_results_df,
  type1_feasible = type1_feasible,
  power_feasible = power_feasible,
  cache_paths = list(H0 = cache_h0_path, H1 = cache_h1_path),
  config = list(
    case_id = cfg$case_id,
    N.max = N.max, interim_cohort = interim_cohort, look_cohort = look_cohort,
    design = design, eta_interim = eta_interim, eta_final = eta_final,
    recruit.int = recruit.int, cens_upper = cens_upper,
    final_analysis_time = final_analysis_time,
    analysis_row_rule = analysis_row_rule,
    c_F_grid = c_F_grid, gamma_E_grid = gamma_E_grid,
    L_sim = L_sim, Niter = Niter, alpha = alpha, power_target = power_target,
    enable_superiority = enable_superiority, rebuild_pp_cache = rebuild_pp_cache,
    Sigma.ctrl = Sigma.ctrl, Sigma_h0_trt = Sigma_h0_trt, Sigma_h1_trt = Sigma_h1_trt,
    mu_h0_trt = mu_h0_trt, mu_h0_ctrl = mu_h0_ctrl,
    mu_h1_trt = mu_h1_trt, mu_h1_ctrl = mu_h1_ctrl,
    clinical_targets = list(recurrence_hr = cfg$calibration_target_recurrence_hr,
                            death_hr = cfg$calibration_target_death_hr),
    generator_version = cfg$generator_version
  )
)
saveRDS(result_obj, file.path(run_dir, paste0("calibration_LG_result_", run_id, ".rds")))
saveRDS(result_obj, file.path(latest_dir, "calibration_LG_result_latest.rds"))

cat("\nCandidate evaluation completed. run_id =", run_id, "\n")
cat("Candidate rows:", nrow(all_results_df), "| H0-filtered rows:", nrow(type1_feasible),
    "| H0/H1-filtered rows:", nrow(power_feasible), "\n")
cat("PP caches:\n", cache_h0_path, "\n", cache_h1_path, "\n")
cat("Candidate table:", file.path(latest_dir, "calibration_LG_table.csv"), "\n")
