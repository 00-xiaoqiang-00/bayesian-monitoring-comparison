# Independent simulation evaluation with locked decision parameters.
# Analysis row rule: interim n-1; final all N.max.
# Data and model seeds are independent of parameter calibration.

args <- commandArgs(trailingOnly = TRUE)
case_hit <- args[grepl("^--case=", args)]
if (length(case_hit) != 1L) {
  stop("Pass exactly one --case=IMpower131|KEYNOTE598|OAK|survival_colon")
}
case_key <- sub("^--case=", "", case_hit[1L])

project_root <- Sys.getenv("PROJECT_ROOT", unset = "")
script_hit <- commandArgs()[grepl("^--file=", commandArgs())]
if (!nzchar(project_root) && length(script_hit) == 1L) {
  script_path <- normalizePath(sub("^--file=", "", script_hit[1L]), winslash = "/", mustWork = TRUE)
  script_dir <- dirname(script_path)
  up_levels <- 4L
  project_root <- script_dir
  for (level in seq_len(up_levels)) project_root <- dirname(project_root)
  project_root <- normalizePath(project_root, winslash = "/", mustWork = TRUE)
}
if (!nzchar(project_root)) project_root <- getwd()
winratio_root <- normalizePath(file.path(project_root, "code", "vendor", "bayesian_win_ratio"), winslash = "/", mustWork = TRUE)
run_mode <- if ("--smoke" %in% args) "smoke" else "formal_10000"
dry_run <- "--dry-run" %in% args

locked_registry <- list(
  IMpower131 = list(
    repo_root = winratio_root,
    case_id = "IMpower131", c_F = 0.875, gamma_E = 3.05,
    data_seed = c(H0_null = 111000000L, H1_calibration = 112000000L, H1_observed = 113000000L),
    model_seed = c(H0_null = 561000000L, H1_calibration = 562000000L, H1_observed = 563000000L)
  ),
  KEYNOTE598 = list(
    repo_root = winratio_root,
    case_id = "keynote598", c_F = 0.865, gamma_E = 1.85,
    data_seed = c(H0_null = 121000000L, H1_calibration = 122000000L, H1_observed = 123000000L),
    model_seed = c(H0_null = 611000000L, H1_calibration = 612000000L, H1_observed = 613000000L)
  ),
  OAK = list(
    repo_root = winratio_root,
    case_id = "OAK_A_ITT", c_F = 0.875, gamma_E = 3.40,
    data_seed = c(H0_null = 131000000L, H1_calibration = 132000000L, H1_observed = 133000000L),
    model_seed = c(H0_null = 661000000L, H1_calibration = 662000000L, H1_observed = 663000000L)
  ),
  survival_colon = list(
    repo_root = winratio_root,
    case_id = "survival_colon", c_F = 0.850, gamma_E = 1.862,
    data_seed = c(H0_null = 141000000L, H1_calibration = 142000000L, H1_observed = 143000000L),
    model_seed = c(H0_null = 711000000L, H1_calibration = 712000000L, H1_observed = 713000000L)
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
L_sim <- 10000L
mc.cores <- 70L
enable_superiority <- cfg$enable_superiority
eta_interim <- cfg$eta_interim
eta_final <- cfg$eta_final
force_rebuild_cache <- FALSE
Sigma.ctrl_default <- cfg$Sigma.ctrl

if (identical(run_mode, "smoke")) {
  L_sim <- 1L
  Niter <- 1000L
  mc.cores <- 1L
}
if (design != 1) stop("PP-cache validation currently supports design = 1 only.")

# 2) Use immutable locked parameters and isolated outputs
lock_root <- file.path(winratio_root, "validation")
out_dir <- file.path(repo_root, "output", "validation", case_key, run_mode)
cache_dir <- file.path(out_dir, "cache")
val_dir <- file.path(out_dir, "results")
for (d in c(lock_root, out_dir, cache_dir, val_dir)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

best_path <- file.path(lock_root, "LOCKED_PARAMETER_MANIFEST.tsv")
if (!file.exists(best_path)) {
  stop("Missing locked parameter manifest: ", best_path)
}
best <- read.delim(best_path, stringsAsFactors = FALSE)
best <- best[best$case_key == case_key, , drop = FALSE]
if (nrow(best) != 1L) stop("Locked manifest must contain exactly one row for ", case_key)

# Use the locked registry values for the decision boundaries.
c_F_override     <- locked$c_F
gamma_E_override <- locked$gamma_E

c_F     <- if (is.null(c_F_override)) as.numeric(best$c_F) else c_F_override
gamma_E <- if (is.null(gamma_E_override)) as.numeric(best$gamma_E) else gamma_E_override
if ("eta_interim" %in% names(best) && !is.na(best$eta_interim)) {
  eta_interim <- as.numeric(best$eta_interim)
}
if ("eta_final" %in% names(best) && !is.na(best$eta_final)) {
  eta_final <- as.numeric(best$eta_final)
}

look_cohort <- sort(unique(as.integer(c(interim_cohort, N.max))))
look_cohort <- look_cohort[look_cohort > 0 & look_cohort <= N.max]
n_looks <- length(look_cohort)
analysis_row_rule <- "interim_n_minus_1_final_all_n_v2"

cat(sprintf(
  "Using design: c_F=%.4f, gamma_E=%.4f, eta_interim=%.2f, eta_final=%.2f\n",
  c_F, gamma_E, eta_interim, eta_final
))
cat(sprintf("Locked source: %s\n", best_path))
cat("Mode:", run_mode, "| L_sim:", L_sim, "| Niter:", Niter,
    "| cores:", mc.cores, "\n")
cat("Looks:", paste(look_cohort, collapse = ", "), "\n")
cat("Cache dir:", cache_dir, "\n")

# 3) Validation scenarios
#    Required: name, role ("null" | "alternative"), mu.trt, mu.ctrl, Sigma.trt
#    Optional: Sigma.ctrl, note, cache_file (override cache path basename)
#    role="null"       -> Type1_error := P(go)
#    role="alternative"-> Power      := P(go)
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
        note = if (!is.null(cfg$calibration_target_recurrence_hr) &&
                   !is.null(cfg$calibration_target_death_hr)) {
          sprintf(
            "Clinical alternative: recurrence HR=%.2f, death HR=%.2f (Power)",
            cfg$calibration_target_recurrence_hr,
            cfg$calibration_target_death_hr
          )
        } else {
          "Calibration alternative (Power)"
        },
        cache_file = "calib_stats_H1.rds",
        mu.trt = cfg$mu_h1_trt,
        mu.ctrl = cfg$mu_h1_ctrl,
        Sigma.trt = cfg$Sigma_h1_trt,
        Sigma.ctrl = cfg$Sigma.ctrl
      )
    ),
    cfg$scenarios_extra
)

seed_offset <- if (identical(run_mode, "smoke")) 800000000L else 0L
for (i in seq_along(scenarios)) {
  nm <- scenarios[[i]]$name
  if (!nm %in% names(locked$data_seed) || !nm %in% names(locked$model_seed)) {
    stop("No locked seed allocation for scenario: ", nm)
  }
  scenarios[[i]]$data_seed_base <- as.integer(locked$data_seed[[nm]] + seed_offset)
  scenarios[[i]]$model_seed_base <- as.integer(locked$model_seed[[nm]] + seed_offset)
}

if (dry_run) {
  cat("DRY RUN PASSED\n")
  cat("Case:", case_key, "| project:", repo_root, "\n")
  cat("Parameters: c_F=", c_F, ", gamma_E=", gamma_E, "\n", sep = "")
  cat("Scenarios and seeds:\n")
  print(do.call(rbind, lapply(scenarios, function(sc) data.frame(
    scenario = sc$name,
    data_seed_base = sc$data_seed_base,
    model_seed_base = sc$model_seed_base
  ))))
  quit(save = "no", status = 0L)
}

# 4) PP cache + replay helpers
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

scenario_fingerprint <- function(sc) {
  Sigma.ctrl_use <- if (!is.null(sc$Sigma.ctrl)) sc$Sigma.ctrl else Sigma.ctrl_default
  generator_scenario <- if (!is.null(sc$generator_scenario)) {
    sc$generator_scenario
  } else if (identical(sc$role, "null")) {
    "H0"
  } else {
    "H1"
  }
  list(
    name = sc$name,
    generation = generation_fingerprint(sc, Sigma.ctrl_use, generator_scenario),
    mu.trt = as.numeric(sc$mu.trt),
    mu.ctrl = as.numeric(sc$mu.ctrl),
    Sigma.trt = as.numeric(sc$Sigma.trt),
    Sigma.ctrl = as.numeric(Sigma.ctrl_use),
    N.max = N.max,
    look_cohort = as.integer(look_cohort),
    recruit.int = recruit.int,
    final_analysis_time = final_analysis_time,
    analysis_row_rule = analysis_row_rule,
    eta_interim = eta_interim,
    eta_final = eta_final,
    Niter = Niter,
    thin_MCMC = thin_MCMC,
    L_sim = L_sim,
    cens_upper = cens_upper,
    seed_namespace = "SEED-IV-V1-20260731",
    data_seed_base = sc$data_seed_base,
    data_seed_formula = "base + replicate_id",
    model_seed_base = sc$model_seed_base,
    model_seed_formula = "base + replicate_id*10 + look_code"
  )
}

cache_path_for <- function(sc) {
  fname <- if (!is.null(sc$cache_file)) sc$cache_file else paste0("calib_stats_", sc$name, ".rds")
  file.path(cache_dir, fname)
}

compute_pp_one_trial <- function(myData, Sigma.ctrl_use,
                                 model_seed_base, replicate_id) {
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
    look_code <- if (is_final_look) 9L else as.integer(j)
    set.seed(as.integer(model_seed_base + replicate_id * 10L + look_code))

    result <- winratio(
      currentdd = temp,
      n_current = n.current,
      Time_current = Time.current,
      Time_entry = Time.entry.analysis,
      m0 = m0, L0 = L0, v0 = v0, S0 = S0,
      time_max = time_max,
      eta = eta_use,
      lambda = 0,
      N.max = N.max,
      thin_MCMC = thin_MCMC,
      Niter = Niter,
      Sigma.ctrl = Sigma.ctrl_use
    )
    pp[j] <- as.numeric(result$probs)
  }
  pp
}

build_pp_cache <- function(data_arr, nsim, label, Sigma.ctrl_use, fingerprint,
                           model_seed_base) {
  cat(sprintf("\nBuilding PP cache for %s (%d trials x %d looks)...\n",
              label, nsim, n_looks))
  t0 <- Sys.time()
  pp_list <- parallel::mclapply(seq_len(nsim), function(x) {
    compute_pp_one_trial(
      data_arr[, , x], Sigma.ctrl_use,
      model_seed_base = model_seed_base, replicate_id = x
    )
  }, mc.cores = mc.cores, mc.set.seed = FALSE)
  pp_mat <- do.call(rbind, pp_list)
  colnames(pp_mat) <- as.character(look_cohort)
  rownames(pp_mat) <- paste0("sim", seq_len(nsim))
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  cat(sprintf("  %s cache done in %.1f min\n", label, elapsed))
  list(
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
    seed_namespace = "SEED-IV-V1-20260731",
    data_seed_base = fingerprint$data_seed_base,
    data_seed_formula = fingerprint$data_seed_formula,
    model_seed_base = fingerprint$model_seed_base,
    model_seed_formula = fingerprint$model_seed_formula,
    fingerprint = fingerprint,
    created_at = as.character(Sys.time()),
    build_minutes = elapsed
  )
}

cache_compatible <- function(cache, fingerprint) {
  if (is.null(cache$PP)) return(FALSE)
  if (!identical(as.integer(cache$look_cohort), as.integer(look_cohort))) return(FALSE)
  if (!isTRUE(cache$N.max == N.max)) return(FALSE)
  if (!isTRUE(cache$eta_interim == eta_interim) || !isTRUE(cache$eta_final == eta_final)) return(FALSE)
  if (!isTRUE(cache$Niter == Niter) || !isTRUE(cache$thin_MCMC == thin_MCMC)) return(FALSE)
  if (!isTRUE(nrow(cache$PP) == L_sim)) return(FALSE)
  if (!isTRUE(all.equal(cache$final_analysis_time, final_analysis_time))) return(FALSE)
  if (!identical(cache$analysis_row_rule, analysis_row_rule)) return(FALSE)
  # When cens_upper metadata is absent, use the fixed value 41.
  cache_cens <- if (is.null(cache$cens_upper)) 41 else cache$cens_upper
  if (!isTRUE(all.equal(cache_cens, cens_upper))) return(FALSE)
  # Calibration caches store the generator fingerprint directly, whereas
  # validation-built caches nest it inside the full scenario fingerprint.
  cache_generation <- if (!is.null(cache$fingerprint$generation)) {
    cache$fingerprint$generation
  } else {
    cache$generation_fingerprint
  }
  if (!is.null(cache_generation)) {
    if (!isTRUE(all.equal(
      cache_generation, fingerprint$generation,
      tolerance = 1e-10, check.attributes = FALSE
    ))) return(FALSE)
  } else if (!identical(fingerprint$generation$generator_type,
                        "bivariate_lognormal_v1")) {
    return(FALSE)
  }
  # If fingerprint exists, require matching generative settings.
  # Case-specific caches must include the complete generator fingerprint.
  if (!is.null(cache$fingerprint)) {
    keys <- c("generation", "mu.trt", "mu.ctrl", "Sigma.trt", "Sigma.ctrl",
              "N.max", "look_cohort", "recruit.int", "final_analysis_time",
              "analysis_row_rule",
              "eta_interim", "eta_final",
              "Niter", "thin_MCMC", "L_sim", "cens_upper",
              "seed_namespace", "data_seed_base", "data_seed_formula",
              "model_seed_base", "model_seed_formula")
    for (k in keys) {
      if (!is.null(cache$fingerprint[[k]]) && !is.null(fingerprint[[k]])) {
        if (!isTRUE(all.equal(cache$fingerprint[[k]], fingerprint[[k]], tolerance = 1e-10))) {
          return(FALSE)
        }
      }
    }
  }
  TRUE
}

generate_independent_data <- function(sc, Sigma.ctrl_use, generator_scenario) {
  cat(sprintf("\nGenerating %d independent trials for %s...\n", L_sim, sc$name))
  data_list <- parallel::mclapply(seq_len(L_sim), function(replicate_id) {
    set.seed(as.integer(sc$data_seed_base + replicate_id))
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

get_or_build_cache <- function(sc) {
  Sigma.ctrl_use <- if (!is.null(sc$Sigma.ctrl)) sc$Sigma.ctrl else Sigma.ctrl_default
  fp <- scenario_fingerprint(sc)
  path <- cache_path_for(sc)

  if (!force_rebuild_cache && file.exists(path)) {
    cache <- readRDS(path)
    if (cache_compatible(cache, fp)) {
      cat(sprintf("  Reusing cache: %s\n", path))
      return(list(cache = cache, path = path, built = FALSE, Sigma.ctrl_use = Sigma.ctrl_use))
    }
    cat(sprintf("  Cache exists but incompatible, rebuilding: %s\n", path))
  } else if (file.exists(path) && force_rebuild_cache) {
    cat(sprintf("  force_rebuild_cache=TRUE, rebuilding: %s\n", path))
  } else {
    cat(sprintf("  No cache found, building: %s\n", path))
  }

  generator_scenario <- if (!is.null(sc$generator_scenario)) {
    sc$generator_scenario
  } else if (identical(sc$role, "null")) {
    "H0"
  } else {
    "H1"
  }
  data_arr <- generate_independent_data(sc, Sigma.ctrl_use, generator_scenario)
  cache <- build_pp_cache(
    data_arr, L_sim, sc$name, Sigma.ctrl_use, fp,
    model_seed_base = sc$model_seed_base
  )
  saveRDS(cache, path)
  cat("  Saved cache:", path, "\n")
  list(cache = cache, path = path, built = TRUE, Sigma.ctrl_use = Sigma.ctrl_use)
}

replay_one <- function(pp_row, c_F_value, gamma_E_value) {
  decision <- NA_character_
  go <- NA
  n_final <- NA_real_
  for (j in seq_len(n_looks)) {
    n.current <- look_cohort[j]
    is_final_look <- (n.current == N.max)
    probs <- as.numeric(pp_row[j])
    ratio <- n.current / N.max
    fut_bound <- c_F_value * (ratio^gamma_E_value)
    sup_bound <- if (isTRUE(enable_superiority)) 1 - ratio^gamma_E_value else Inf

    if (!is_final_look) {
      if (!is.na(probs) && probs <= fut_bound) {
        decision <- "futility"; go <- FALSE; n_final <- n.current; break
      } else if (!is.na(probs) && probs >= sup_bound) {
        decision <- "superiority"; go <- TRUE; n_final <- n.current; break
      }
    } else {
      n_final <- n.current
      if (!is.na(probs) && probs > c_F_value) {
        decision <- "go"; go <- TRUE
      } else {
        decision <- "no_go"; go <- FALSE
      }
      break
    }
  }
  list(go = isTRUE(go), decision = decision, pts.stop = n_final)
}

oc_from_cache <- function(pp_mat, c_F_value, gamma_E_value) {
  nsim <- nrow(pp_mat)
  go_vec <- logical(nsim)
  fut_vec <- logical(nsim)
  sup_vec <- logical(nsim)
  nogo_vec <- logical(nsim)
  ess_vec <- numeric(nsim)
  for (i in seq_len(nsim)) {
    r <- replay_one(pp_mat[i, ], c_F_value, gamma_E_value)
    go_vec[i] <- r$go
    fut_vec[i] <- identical(r$decision, "futility")
    sup_vec[i] <- identical(r$decision, "superiority")
    nogo_vec[i] <- identical(r$decision, "no_go")
    ess_vec[i] <- r$pts.stop
  }
  list(
    P_go = mean(go_vec, na.rm = TRUE),
    PET_futility = mean(fut_vec, na.rm = TRUE),
    P_early_success = mean(sup_vec, na.rm = TRUE),
    P_final_nogo = mean(nogo_vec, na.rm = TRUE),
    ESS = mean(ess_vec, na.rm = TRUE)
  )
}

validate_scenario <- function(sc) {
  role <- if (is.null(sc$role)) "alternative" else sc$role
  if (!role %in% c("null", "alternative")) {
    stop("scenario '", sc$name, "': role must be 'null' or 'alternative', got: ", role)
  }
  cat(sprintf("\n=== Scenario: %s (role=%s) ===\n", sc$name, role))
  if (!is.null(sc$note)) cat("  note:", sc$note, "\n")
  cat("  mu.trt =", paste(round(sc$mu.trt, 4), collapse = ", "),
      " | mu.ctrl =", paste(round(sc$mu.ctrl, 4), collapse = ", "), "\n")

  t0 <- Sys.time()
  got <- get_or_build_cache(sc)
  oc <- oc_from_cache(got$cache$PP, c_F, gamma_E)
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "mins"))

  Type1_error <- if (identical(role, "null")) oc$P_go else NA_real_
  Power <- if (identical(role, "alternative")) oc$P_go else NA_real_

  cat(sprintf(
    "  P(go)=%.4f  Type1=%s  Power=%s  PET_fut=%.4f  ESS=%.2f  cache_%s  (%.1f min)\n",
    oc$P_go,
    if (is.na(Type1_error)) "NA" else sprintf("%.4f", Type1_error),
    if (is.na(Power)) "NA" else sprintf("%.4f", Power),
    oc$PET_futility, oc$ESS,
    if (got$built) "built" else "reused",
    elapsed
  ))

  data.frame(
    scenario = sc$name,
    role = role,
    note = if (is.null(sc$note)) NA_character_ else sc$note,
    cache_path = got$path,
    cache_built = got$built,
    c_F = c_F,
    gamma_E = gamma_E,
    eta_interim = eta_interim,
    eta_final = eta_final,
    Type1_error = Type1_error,
    Power = Power,
    P_go = oc$P_go,
    PET_futility = oc$PET_futility,
    P_early_success = oc$P_early_success,
    P_final_nogo = oc$P_final_nogo,
    ESS = oc$ESS,
    L_sim = L_sim,
    Niter = Niter,
    N.max = N.max,
    seed_namespace = "SEED-IV-V1-20260731",
    data_seed_base = sc$data_seed_base,
    model_seed_base = sc$model_seed_base,
    elapsed_min = elapsed,
    mu_trt_1 = sc$mu.trt[1],
    mu_trt_2 = sc$mu.trt[2],
    mu_ctrl_1 = sc$mu.ctrl[1],
    mu_ctrl_2 = sc$mu.ctrl[2],
    stringsAsFactors = FALSE
  )
}

# 5) Run all scenarios
run_id <- format(Sys.time(), "%Y%m%d_%H%M%S")
cat("run_id =", run_id, "\n")
cat("Scenarios:", paste(vapply(scenarios, `[[`, "", "name"), collapse = ", "), "\n")
cat("PP-cache: reuse if compatible; new/changed scenarios are built and saved.\n")

rows <- lapply(scenarios, validate_scenario)
results_df <- do.call(rbind, rows)
results_df$run_id <- run_id
results_df$seed_namespace <- "SEED-IV-V1-20260731"
results_df$calib_source <- best_path

oc_summary <- results_df[, c(
  "scenario", "role", "Type1_error", "Power",
  "P_go", "PET_futility", "P_final_nogo", "ESS",
  "cache_built", "elapsed_min"
)]

type1_rows <- results_df[results_df$role == "null", , drop = FALSE]
power_rows <- results_df[results_df$role == "alternative", , drop = FALSE]
key_metrics <- rbind(
  if (nrow(type1_rows) > 0) {
    data.frame(
      metric = "Type1_error",
      scenario = type1_rows$scenario,
      value = type1_rows$Type1_error,
      PET_futility = type1_rows$PET_futility,
      ESS = type1_rows$ESS,
      stringsAsFactors = FALSE
    )
  } else NULL,
  if (nrow(power_rows) > 0) {
    data.frame(
      metric = "Power",
      scenario = power_rows$scenario,
      value = power_rows$Power,
      PET_futility = power_rows$PET_futility,
      ESS = power_rows$ESS,
      stringsAsFactors = FALSE
    )
  } else NULL
)

out_csv <- file.path(val_dir, paste0("validation_LG_", run_id, ".csv"))
latest_csv <- file.path(val_dir, "validation_LG_latest.csv")
oc_csv <- file.path(val_dir, paste0("validation_LG_oc_", run_id, ".csv"))
oc_latest <- file.path(val_dir, "validation_LG_oc_latest.csv")
key_csv <- file.path(val_dir, paste0("validation_LG_key_metrics_", run_id, ".csv"))
key_latest <- file.path(val_dir, "validation_LG_key_metrics_latest.csv")

write.csv(results_df, out_csv, row.names = FALSE)
write.csv(results_df, latest_csv, row.names = FALSE)
write.csv(oc_summary, oc_csv, row.names = FALSE)
write.csv(oc_summary, oc_latest, row.names = FALSE)
write.csv(key_metrics, key_csv, row.names = FALSE)
write.csv(key_metrics, key_latest, row.names = FALSE)

saveRDS(
  list(
    results = results_df,
    oc_summary = oc_summary,
    key_metrics = key_metrics,
    design = list(
      c_F = c_F, gamma_E = gamma_E,
      eta_interim = eta_interim, eta_final = eta_final,
      N.max = N.max, interim_cohort = interim_cohort, look_cohort = look_cohort,
      recruit.int = recruit.int, cens_upper = cens_upper,
      final_analysis_time = final_analysis_time,
      analysis_row_rule = analysis_row_rule,
      L_sim = L_sim, Niter = Niter,
      seed_namespace = "SEED-IV-V1-20260731",
      data_seed_formula = "base + replicate_id",
      model_seed_formula = "base + replicate_id*10 + look_code",
      force_rebuild_cache = force_rebuild_cache
    ),
    clinical_targets = if (is.null(cfg)) NULL else list(
      recurrence_hr = cfg$calibration_target_recurrence_hr,
      death_hr = cfg$calibration_target_death_hr
    ),
    scenarios = scenarios,
    locked_parameter_manifest_row = best,
    cache_dir = cache_dir
  ),
  file.path(val_dir, paste0("validation_LG_", run_id, ".rds"))
)

cat("\n===== Key metrics (Type I / Power) =====\n")
print(key_metrics)
cat("\n===== Full OC summary =====\n")
print(oc_summary)
cat("\nCaches:\n")
print(results_df[, c("scenario", "cache_path", "cache_built")])
cat("\nWrote:\n")
cat("  ", out_csv, "\n")
cat("  ", oc_latest, "\n")
cat("  ", key_latest, "\n")
cat("  ", latest_csv, "\n")
