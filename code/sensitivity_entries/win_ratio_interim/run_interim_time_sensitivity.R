#!/usr/bin/env Rscript

# Locked-design interim-time sensitivity.
# Formal 40/60/80 caches are read-only. New 50% PP is written only under
# this directory. look_code is mapped by sample size, never by new-look index.

args <- commandArgs(trailingOnly = TRUE)
case_hit <- args[grepl("^--case=", args)]
if (length(case_hit) != 1L) {
  stop("Pass exactly one --case=IMpower131|KEYNOTE598|OAK|survival_colon")
}
case_key <- sub("^--case=", "", case_hit[1L])
run_mode <- if ("--smoke" %in% args) "smoke" else "formal_10000"
dry_run <- "--dry-run" %in% args
replay_only <- "--replay-only" %in% args
force_rebuild_look50 <- "--rebuild-look50" %in% args
force_run <- "--force" %in% args

script_hit <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_hit) != 1L) stop("Run this entry with Rscript.")
sens_root <- dirname(normalizePath(sub("^--file=", "", script_hit[1L]), winslash = "/", mustWork = TRUE))
project_root <- Sys.getenv("PROJECT_ROOT", unset = "")
winratio_root <- if (nzchar(project_root)) {
  normalizePath(file.path(project_root, "code", "vendor", "bayesian_win_ratio"), winslash = "/", mustWork = TRUE)
} else {
  normalizePath(file.path(sens_root, "..", "..", "vendor", "bayesian_win_ratio"), winslash = "/", mustWork = TRUE)
}
source(file.path(sens_root, "R", "interim_sensitivity_core.R"))

schedule_all <- read_tsv(file.path(sens_root, "SCHEDULE_MANIFEST.tsv"))
seed_all <- read_tsv(file.path(sens_root, "SEED_MANIFEST.tsv"))
sched <- schedule_all[schedule_all$case_key == case_key, , drop = FALSE]
if (nrow(sched) != 1L) stop("SCHEDULE_MANIFEST must have one row for ", case_key)
seed_case <- seed_all[seed_all$case_key == case_key, , drop = FALSE]
if (nrow(seed_case) != 3L) stop("SEED_MANIFEST must have 3 scenarios for ", case_key)

locked_registry <- list(
  IMpower131 = list(
    repo_root = winratio_root,
    case_id = "IMpower131"
  ),
  KEYNOTE598 = list(
    repo_root = winratio_root,
    case_id = "keynote598"
  ),
  OAK = list(
    repo_root = winratio_root,
    case_id = "OAK_A_ITT"
  ),
  survival_colon = list(
    repo_root = winratio_root,
    case_id = "survival_colon"
  )
)
if (!case_key %in% names(locked_registry)) stop("Unknown case key: ", case_key)
locked <- locked_registry[[case_key]]
repo_root <- normalizePath(locked$repo_root, winslash = "/", mustWork = TRUE)

pkgload::load_all(file.path(repo_root, "engine"), quiet = TRUE)
cat("load_all success:", repo_root, "\n")
library(parallel)
RNGkind("L'Ecuyer-CMRG")

source(file.path(repo_root, "shared", "load_case_config.R"))
cfg <- load_case_config(locked$case_id, repo_root = repo_root)

N.max <- as.integer(sched$N_max)
if (!identical(as.integer(cfg$N.max), N.max)) {
  stop("config N.max=", cfg$N.max, " differs from schedule N_max=", N.max)
}
n_40 <- as.integer(sched$n_40)
n_50 <- as.integer(sched$n_50)
n_60 <- as.integer(sched$n_60)
n_80 <- as.integer(sched$n_80)
n_final <- as.integer(sched$n_final)
c_F <- as.numeric(sched$c_F)
gamma_E <- as.numeric(sched$gamma_E)

looks_primary <- as.integer(c(n_40, n_60, n_80, n_final))
looks_sens <- as.integer(c(n_50, n_80, n_final))
n_keys <- as.character(c(n_40, n_50, n_60, n_80, n_final))
if (any(duplicated(n_keys))) {
  stop("Look sample sizes must be unique within a case: ", paste(n_keys, collapse = "/"))
}
look_codes <- c(
  setNames(as.integer(sched$look_code_40), as.character(n_40)),
  setNames(as.integer(sched$look_code_50), as.character(n_50)),
  setNames(as.integer(sched$look_code_60), as.character(n_60)),
  setNames(as.integer(sched$look_code_80), as.character(n_80)),
  final = as.integer(sched$look_code_final)
)
if (!identical(as.integer(unique(seed_case$look_code_50)), as.integer(sched$look_code_50)) ||
    !identical(as.integer(unique(seed_case$look_code_80)), as.integer(sched$look_code_80)) ||
    !identical(as.integer(unique(seed_case$look_code_final)), as.integer(sched$look_code_final))) {
  stop("SEED_MANIFEST look_code columns do not match SCHEDULE_MANIFEST")
}

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
L_sim <- 10000L
mc.cores <- 66L
enable_superiority <- isTRUE(cfg$enable_superiority)
eta_interim <- cfg$eta_interim
eta_final <- cfg$eta_final
analysis_row_rule <- "interim_n_minus_1_final_all_n_v2"
Sigma.ctrl_default <- cfg$Sigma.ctrl
chunk_size <- 140L

if (identical(run_mode, "smoke")) {
  L_sim <- 1L
  Niter <- 1000L
  mc.cores <- 1L
  chunk_size <- 1L
}
if (design != 1) stop("PP-cache validation currently supports design = 1 only.")

if (identical(run_mode, "formal_10000") && !replay_only && !dry_run && !force_run) {
  pid_file <- Sys.getenv("PID_FILE", unset = file.path(tempdir(), "formal_run.pid"))
  if (file.exists(pid_file)) {
    formal_pid <- as.integer(readLines(pid_file, warn = FALSE)[1L])
    if (!is.na(formal_pid) &&
        file.exists(file.path("/proc", as.character(formal_pid)))) {
      stop("Formal independent validation still running (pid ", formal_pid,
           "). Wait, or pass --force.")
    }
  }
}

out_root <- file.path(sens_root, "outputs", run_mode, case_key)
cache_dir <- file.path(out_root, "cache")
val_dir <- file.path(out_root, "results")
for (d in c(out_root, cache_dir, val_dir)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

best_path <- file.path(repo_root, "validation",
                       "LOCKED_PARAMETER_MANIFEST.tsv")
best <- read_tsv(best_path)
best <- best[best$case_key == case_key, , drop = FALSE]
if (nrow(best) != 1L) stop("Locked manifest must contain exactly one row for ", case_key)
if (!isTRUE(all.equal(as.numeric(best$c_F), c_F)) ||
    !isTRUE(all.equal(as.numeric(best$gamma_E), gamma_E))) {
  stop("Schedule c_F/gamma_E do not match locked manifest for ", case_key)
}
if ("eta_interim" %in% names(best) && !is.na(best$eta_interim)) {
  eta_interim <- as.numeric(best$eta_interim)
}
if ("eta_final" %in% names(best) && !is.na(best$eta_final)) {
  eta_final <- as.numeric(best$eta_final)
}

scenarios <- c(
  list(
    list(
      name = "H0_null",
      role = "null",
      generator_scenario = "H0",
      note = "Calibration null: identical arms (Type I)",
      cache_file = "calib_stats_H0.rds",
      mu.trt = cfg$mu_h0_trt,
      mu.ctrl = cfg$mu_h0_ctrl,
      Sigma.trt = cfg$Sigma_h0_trt,
      Sigma.ctrl = cfg$Sigma.ctrl
    ),
    list(
      name = "H1_calibration",
      role = "alternative",
      generator_scenario = "H1_calibration",
      note = "Calibration alternative (Power)",
      cache_file = "calib_stats_H1.rds",
      mu.trt = cfg$mu_h1_trt,
      mu.ctrl = cfg$mu_h1_ctrl,
      Sigma.trt = cfg$Sigma_h1_trt,
      Sigma.ctrl = cfg$Sigma.ctrl
    )
  ),
  cfg$scenarios_extra
)
if (length(scenarios) != 3L ||
    !identical(vapply(scenarios, `[[`, "", "name"),
               c("H0_null", "H1_calibration", "H1_observed"))) {
  stop("Expected scenarios H0_null, H1_calibration, H1_observed")
}

seed_offset <- if (identical(run_mode, "smoke")) {
  as.integer(seed_case$smoke_offset[1L])
} else {
  0L
}
for (i in seq_along(scenarios)) {
  nm <- scenarios[[i]]$name
  seed_row <- seed_case[seed_case$scenario == nm, , drop = FALSE]
  if (nrow(seed_row) != 1L) stop("No seed row for ", nm)
  scenarios[[i]]$data_seed_base <- as.integer(seed_row$data_seed_base + seed_offset)
  scenarios[[i]]$model_seed_base <- as.integer(seed_row$model_seed_base + seed_offset)
}

generation_fingerprint <- function(sc, Sigma.ctrl_use, generator_scenario) {
  if (exists("case.generator.fingerprint", mode = "function", inherits = TRUE)) {
    return(case.generator.fingerprint(
      cfg = cfg, scenario = generator_scenario,
      mu.trt = sc$mu.trt, Sigma.trt = sc$Sigma.trt,
      mu.ctrl = sc$mu.ctrl, Sigma.ctrl = Sigma.ctrl_use,
      cens_upper = cens_upper
    ))
  }
  list(
    generator_type = "bivariate_lognormal_v1",
    generator_version = 1L,
    scenario = generator_scenario,
    mu.trt = as.numeric(sc$mu.trt),
    Sigma.trt = as.numeric(sc$Sigma.trt),
    mu.ctrl = as.numeric(sc$mu.ctrl),
    Sigma.ctrl = as.numeric(Sigma.ctrl_use),
    cens_upper = cens_upper
  )
}

formal_cache_ok <- function(cache, sc, Sigma.ctrl_use, generation) {
  if (is.null(cache$PP) || !isTRUE(nrow(cache$PP) == L_sim)) return(FALSE)
  if (any(!is.finite(cache$PP))) return(FALSE)
  if (!identical(as.integer(cache$look_cohort), as.integer(looks_primary))) {
    return(FALSE)
  }
  if (!isTRUE(as.integer(cache$N.max) == N.max)) return(FALSE)
  if (!isTRUE(as.integer(cache$Niter) == as.integer(Niter)) ||
      !isTRUE(as.integer(cache$thin_MCMC) == as.integer(thin_MCMC))) {
    return(FALSE)
  }
  if (!isTRUE(all.equal(cache$eta_interim, eta_interim)) ||
      !isTRUE(all.equal(cache$eta_final, eta_final))) {
    return(FALSE)
  }
  if (!isTRUE(all.equal(cache$recruit.int, recruit.int))) return(FALSE)
  cache_cens <- if (is.null(cache$cens_upper)) 41 else cache$cens_upper
  if (!isTRUE(all.equal(cache_cens, cens_upper))) return(FALSE)
  if (!isTRUE(all.equal(cache$final_analysis_time, final_analysis_time))) return(FALSE)
  if (!identical(cache$analysis_row_rule, analysis_row_rule)) return(FALSE)
  if (!is.null(cache$seed_namespace) &&
      !identical(cache$seed_namespace, "SEED-IV-V1-20260731")) return(FALSE)
  if (is.null(cache$data_seed_base) ||
      !isTRUE(cache$data_seed_base == sc$data_seed_base)) return(FALSE)
  if (is.null(cache$model_seed_base) ||
      !isTRUE(cache$model_seed_base == sc$model_seed_base)) return(FALSE)
  cache_generation <- if (!is.null(cache$fingerprint$generation)) {
    cache$fingerprint$generation
  } else {
    cache$generation_fingerprint
  }
  if (!is.null(cache_generation) && !is.null(generation)) {
    if (!isTRUE(all.equal(
      cache_generation, generation,
      tolerance = 1e-10, check.attributes = FALSE
    ))) return(FALSE)
  }
  TRUE
}

generate_independent_data <- function(sc, Sigma.ctrl_use, generator_scenario) {
  cat(sprintf("\nGenerating %d independent trials for %s...\n", L_sim, sc$name))
  data_list <- parallel::mclapply(seq_len(L_sim), function(replicate_id) {
    set.seed(data_seed(sc$data_seed_base, replicate_id))
    one <- if (exists("case.data.simulation", mode = "function", inherits = TRUE)) {
      case.data.simulation(
        N.sim = 1L, N.max = N.max,
        mu.trt = sc$mu.trt, Sigma.trt = sc$Sigma.trt,
        mu.ctrl = sc$mu.ctrl, Sigma.ctrl = Sigma.ctrl_use,
        cens_upper = cens_upper, cfg = cfg, scenario = generator_scenario
      )
    } else {
      data.simulation(
        N.sim = 1L, N.max = N.max,
        mu.trt = sc$mu.trt, Sigma.trt = sc$Sigma.trt,
        mu.ctrl = sc$mu.ctrl, Sigma.ctrl = Sigma.ctrl_use,
        cens_upper = cens_upper
      )
    }
    if (length(dim(one)) == 3L) one[, , 1L, drop = FALSE][, , 1L] else as.matrix(one)
  }, mc.cores = mc.cores, mc.set.seed = FALSE)
  first <- data_list[[1L]]
  data_arr <- array(
    NA_real_, dim = c(nrow(first), ncol(first), L_sim),
    dimnames = list(rownames(first), colnames(first), paste0("sim", seq_len(L_sim)))
  )
  for (i in seq_len(L_sim)) data_arr[, , i] <- data_list[[i]]
  data_arr
}

build_look50_pp <- function(sc, data_arr, Sigma.ctrl_use, generation) {
  label <- sc$name
  checkpoint_dir <- file.path(cache_dir, "checkpoints", label)
  if (!dir.exists(checkpoint_dir)) dir.create(checkpoint_dir, recursive = TRUE)
  indices <- seq_len(L_sim)
  chunk_id <- ceiling(indices / chunk_size)
  chunks <- split(indices, chunk_id)
  pp_chunks <- vector("list", length(chunks))
  look_code_50 <- as.integer(look_codes[[as.character(n_50)]])
  cat(sprintf(
    "  Building 50%% PP n=%d look_code=%d; %d chunks of %d\n",
    n_50, look_code_50, length(chunks), chunk_size
  ))
  t0 <- Sys.time()
  for (i in seq_along(chunks)) {
    idx <- as.integer(chunks[[i]])
    checkpoint_path <- file.path(
      checkpoint_dir,
      sprintf("chunk_%06d_%06d.rds", min(idx), max(idx))
    )
    expected <- look50_checkpoint_signature(
      idx, case_key, label, n_50, look_code_50, N.max, Niter, thin_MCMC,
      eta_interim, eta_final, recruit.int, final_analysis_time, cens_upper,
      analysis_row_rule, sc$model_seed_base, sc$data_seed_base, generation
    )
    checkpoint <- NULL
    if (!force_rebuild_look50 && file.exists(checkpoint_path)) {
      checkpoint <- tryCatch(readRDS(checkpoint_path), error = function(e) NULL)
    }
    if (pp_checkpoint_compatible(checkpoint, expected)) {
      pp_chunks[[i]] <- checkpoint$PP
      cat(sprintf("  [%s] reused look50 chunk %d/%d (%d-%d)\n",
                  label, i, length(chunks), min(idx), max(idx)))
      next
    }
    chunk_started <- Sys.time()
    pp_list <- parallel::mclapply(idx, function(x) {
      compute_pp_one_look(
        data_arr[, , x], n_50, N.max, look_code_50, Sigma.ctrl_use,
        sc$model_seed_base, x, recruit.int, final_analysis_time,
        m0, L0, v0, S0, time_max, eta_interim, eta_final, thin_MCMC, Niter
      )
    }, mc.cores = min(as.integer(mc.cores), length(idx)), mc.set.seed = FALSE)
    bad <- vapply(pp_list, function(x) {
      inherits(x, "try-error") || !is.numeric(x) || length(x) != 1L || !is.finite(x)
    }, logical(1))
    if (any(bad)) {
      stop("Invalid 50% PP in ", label, " trials: ", paste(idx[bad], collapse = ", "))
    }
    pp_mat <- matrix(unlist(pp_list), ncol = 1L)
    colnames(pp_mat) <- as.character(n_50)
    rownames(pp_mat) <- paste0("sim", idx)
    checkpoint <- list(
      signature = expected,
      PP = pp_mat,
      completed_at = as.character(Sys.time())
    )
    atomic_save_rds(checkpoint, checkpoint_path)
    pp_chunks[[i]] <- pp_mat
    cat(sprintf(
      "  [%s] saved look50 chunk %d/%d (%d-%d) in %.1f min\n",
      label, i, length(chunks), min(idx), max(idx),
      as.numeric(difftime(Sys.time(), chunk_started, units = "mins"))
    ))
  }
  pp50 <- do.call(rbind, pp_chunks)
  rownames(pp50) <- paste0("sim", seq_len(L_sim))
  colnames(pp50) <- as.character(n_50)
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  assembled <- list(
    PP = pp50,
    look_cohort = n_50,
    look_code = look_code_50,
    N.max = N.max,
    L_sim = L_sim,
    scenario = label,
    build_minutes = elapsed,
    created_at = as.character(Sys.time())
  )
  atomic_save_rds(assembled, file.path(cache_dir, paste0("pp_look50_", label, ".rds")))
  assembled
}

build_smoke_pp_matrix <- function(sc, data_arr, Sigma.ctrl_use, look_ns) {
  pp <- matrix(NA_real_, nrow = L_sim, ncol = length(look_ns))
  colnames(pp) <- as.character(look_ns)
  for (j in seq_along(look_ns)) {
    n.current <- look_ns[j]
    code <- look_code_for_n(n.current, N.max, look_codes)
    for (x in seq_len(L_sim)) {
      pp[x, j] <- compute_pp_one_look(
        data_arr[, , x], n.current, N.max, code, Sigma.ctrl_use,
        sc$model_seed_base, x, recruit.int, final_analysis_time,
        m0, L0, v0, S0, time_max, eta_interim, eta_final, thin_MCMC, Niter
      )
    }
  }
  rownames(pp) <- paste0("sim", seq_len(L_sim))
  pp
}

cat(sprintf(
  "Using locked design: c_F=%.4f, gamma_E=%.4f\n", c_F, gamma_E
))
cat("Mode:", run_mode, "| L_sim:", L_sim, "| Niter:", Niter,
    "| cores:", mc.cores, "\n")
cat("Primary looks:", paste(looks_primary, collapse = ", "), "\n")
cat("Sensitivity looks:", paste(looks_sens, collapse = ", "), "\n")
cat("look_code 50/80/final:",
    look_codes[[as.character(n_50)]], "/",
    look_codes[[as.character(n_80)]], "/",
    look_codes[["final"]], "\n")
cat("Output:", out_root, "\n")

if (dry_run) {
  cat("DRY RUN\n")
  for (sc in scenarios) {
    fpath <- formal_cache_path(repo_root, case_key, sc$cache_file)
    cat(sprintf(
      "  %s  data_seed=%d  model_seed=%d  formal_cache=%s exists=%s\n",
      sc$name, sc$data_seed_base, sc$model_seed_base, fpath, file.exists(fpath)
    ))
  }
  cat("DRY RUN PASSED\n")
  quit(save = "no", status = 0L)
}

official_oc <- NULL
official_path <- official_oc_path(repo_root, case_key)
if (file.exists(official_path) && identical(run_mode, "formal_10000")) {
  official_oc <- utils::read.csv(official_path, stringsAsFactors = FALSE)
}

run_id <- format(Sys.time(), "%Y%m%d_%H%M%S")
oc_rows <- list()
replicate_store <- list()
replay_checks <- list()

for (sc in scenarios) {
  role <- sc$role
  Sigma.ctrl_use <- if (!is.null(sc$Sigma.ctrl)) sc$Sigma.ctrl else Sigma.ctrl_default
  generator_scenario <- resolve_generator_scenario(sc)
  generation <- generation_fingerprint(sc, Sigma.ctrl_use, generator_scenario)
  cat(sprintf("\n=== Scenario: %s (role=%s) ===\n", sc$name, role))

  pp_primary <- NULL
  pp_sens <- NULL
  reused_80_final <- FALSE

  if (identical(run_mode, "formal_10000")) {
    fpath <- formal_cache_path(repo_root, case_key, sc$cache_file)
    if (!file.exists(fpath)) {
      msg <- paste0("Formal cache missing: ", fpath)
      if (replay_only || dry_run) {
        cat("  SKIP:", msg, "\n")
        next
      }
      stop(msg)
    }
    cache <- readRDS(fpath)
    if (!formal_cache_ok(cache, sc, Sigma.ctrl_use, generation)) {
      stop("Formal cache incompatible: ", fpath)
    }
    pp_primary <- as.matrix(cache$PP)
    colnames(pp_primary) <- as.character(as.integer(cache$look_cohort))
    rownames(pp_primary) <- paste0("sim", seq_len(L_sim))
    col80 <- extract_pp_column(cache, n_80)
    col_final <- extract_pp_column(cache, n_final)
    reused_80_final <- TRUE

    if (replay_only) {
      look50_path <- file.path(cache_dir, paste0("pp_look50_", sc$name, ".rds"))
      if (!file.exists(look50_path)) {
        cat("  replay-only: 50% PP not built yet; skipping 50/80 replay\n")
      } else {
        look50 <- readRDS(look50_path)
        pp_sens <- cbind(
          look50$PP[, 1L],
          col80,
          col_final
        )
        colnames(pp_sens) <- as.character(looks_sens)
        rownames(pp_sens) <- rownames(pp_primary)
      }
    } else {
      data_arr <- generate_independent_data(sc, Sigma.ctrl_use, generator_scenario)
      look50 <- build_look50_pp(sc, data_arr, Sigma.ctrl_use, generation)
      if (!identical(nrow(look50$PP), L_sim)) {
        stop("50% PP row count ", nrow(look50$PP), " != L_sim ", L_sim)
      }
      pp_sens <- cbind(look50$PP[, 1L], col80, col_final)
      colnames(pp_sens) <- as.character(looks_sens)
      rownames(pp_sens) <- rownames(pp_primary)
    }
  } else {
    data_arr <- generate_independent_data(sc, Sigma.ctrl_use, generator_scenario)
    pp_sens <- build_smoke_pp_matrix(sc, data_arr, Sigma.ctrl_use, looks_sens)
    pp_primary <- build_smoke_pp_matrix(sc, data_arr, Sigma.ctrl_use, looks_primary)
  }

  replay_pri <- replay_schedule(pp_primary, looks_primary, N.max, c_F, gamma_E,
                                enable_superiority)
  oc_pri <- oc_row(role, replay_pri$oc, L_sim)
  oc_pri$case_key <- case_key
  oc_pri$scenario <- sc$name
  oc_pri$schedule <- "40_60_80"
  oc_pri$looks <- paste(looks_primary, collapse = "/")
  oc_pri$c_F <- c_F
  oc_pri$gamma_E <- gamma_E
  oc_pri$reused_80_final <- reused_80_final
  oc_pri$run_id <- run_id
  oc_rows[[length(oc_rows) + 1L]] <- oc_pri

  if (identical(run_mode, "formal_10000") && !is.null(official_oc)) {
    off <- official_oc[official_oc$scenario == sc$name, , drop = FALSE]
    if (nrow(off) == 1L) {
      chk <- compare_oc_to_official(replay_pri$oc, off)
      replay_checks[[sc$name]] <- chk
      cat(sprintf(
        "  40/60/80 replay vs official: %s  max|diff|=%.3e\n",
        if (chk$ok) "PASS" else "FAIL",
        max(abs(chk$diffs))
      ))
      if (!chk$ok) {
        stop("40/60/80 replay does not match official OC for ", sc$name)
      }
    }
  }

  cat(sprintf(
    "  40/60/80  P(go)=%.4f  PET=%.4f  ESS=%.2f\n",
    replay_pri$oc$P_go, replay_pri$oc$PET_futility, replay_pri$oc$ESS
  ))

  if (!is.null(pp_sens)) {
    replay_sens <- replay_schedule(pp_sens, looks_sens, N.max, c_F, gamma_E,
                                   enable_superiority)
    oc_sens <- oc_row(role, replay_sens$oc, L_sim)
    oc_sens$case_key <- case_key
    oc_sens$scenario <- sc$name
    oc_sens$schedule <- "50_80"
    oc_sens$looks <- paste(looks_sens, collapse = "/")
    oc_sens$c_F <- c_F
    oc_sens$gamma_E <- gamma_E
    oc_sens$reused_80_final <- reused_80_final
    oc_sens$run_id <- run_id
    oc_rows[[length(oc_rows) + 1L]] <- oc_sens
    cat(sprintf(
      "  50/80     P(go)=%.4f  PET=%.4f  ESS=%.2f\n",
      replay_sens$oc$P_go, replay_sens$oc$PET_futility, replay_sens$oc$ESS
    ))
    replicate_store[[sc$name]] <- data.frame(
      replicate_id = seq_len(L_sim),
      case_key = case_key,
      scenario = sc$name,
      go_40_60_80 = replay_pri$go,
      fut_40_60_80 = replay_pri$futility,
      nogo_40_60_80 = replay_pri$final_nogo,
      n_40_60_80 = replay_pri$n_used,
      stage_40_60_80 = replay_pri$stop_stage,
      go_50_80 = replay_sens$go,
      fut_50_80 = replay_sens$futility,
      nogo_50_80 = replay_sens$final_nogo,
      n_50_80 = replay_sens$n_used,
      stage_50_80 = replay_sens$stop_stage,
      stringsAsFactors = FALSE
    )
  }
}

if (!length(oc_rows)) {
  cat("No scenarios completed (formal caches missing?).\n")
  quit(save = "no", status = 0L)
}

oc_df <- do.call(rbind, oc_rows)
oc_csv <- file.path(val_dir, paste0("oc_", run_id, ".csv"))
oc_latest <- file.path(val_dir, "oc_latest.csv")
utils::write.csv(oc_df, oc_csv, row.names = FALSE)
utils::write.csv(oc_df, oc_latest, row.names = FALSE)

if (length(replicate_store)) {
  rep_df <- do.call(rbind, replicate_store)
  atomic_save_rds(rep_df, file.path(val_dir, paste0("replicates_", run_id, ".rds")))
  atomic_save_rds(rep_df, file.path(val_dir, "replicates_latest.rds"))
  delta_rows <- lapply(split(rep_df, rep_df$scenario), function(d) {
    go <- paired_mean_se(d$go_40_60_80, d$go_50_80)
    fut <- paired_mean_se(d$fut_40_60_80, d$fut_50_80)
    ess <- paired_mean_se(d$n_40_60_80, d$n_50_80)
    data.frame(
      case_key = case_key,
      scenario = d$scenario[1L],
      role = scenarios[[match(d$scenario[1L], vapply(scenarios, `[[`, "", "name"))]]$role,
      L_sim = nrow(d),
      P_go_40_60_80 = go$mean_x,
      P_go_50_80 = go$mean_y,
      delta_P_go = go$delta,
      se_delta_P_go = go$se,
      PET_40_60_80 = fut$mean_x,
      PET_50_80 = fut$mean_y,
      delta_PET = fut$delta,
      se_delta_PET = fut$se,
      ESS_40_60_80 = ess$mean_x,
      ESS_50_80 = ess$mean_y,
      delta_ESS = ess$delta,
      se_delta_ESS = ess$se,
      stringsAsFactors = FALSE
    )
  })
  delta_df <- do.call(rbind, delta_rows)
  utils::write.csv(delta_df, file.path(val_dir, paste0("paired_delta_", run_id, ".csv")),
                   row.names = FALSE)
  utils::write.csv(delta_df, file.path(val_dir, "paired_delta_latest.csv"),
                   row.names = FALSE)
}

saveRDS(
  list(
    oc = oc_df,
    replay_checks = replay_checks,
    design = list(
      c_F = c_F, gamma_E = gamma_E,
      looks_primary = looks_primary, looks_sens = looks_sens,
      look_codes = look_codes, L_sim = L_sim, Niter = Niter
    ),
    run_id = run_id,
    mode = run_mode
  ),
  file.path(val_dir, paste0("sensitivity_", run_id, ".rds"))
)

cat("\nWrote:\n  ", oc_latest, "\n")
cat("run_id =", run_id, "\n")
