# Shared helpers for locked-design interim-time sensitivity.
# Formal PP caches are read-only inputs.

checkpoint_format_version <- "winratio_interim_sens_look50_v1"

atomic_save_rds <- function(object, path) {
  parent <- dirname(path)
  if (!dir.exists(parent)) dir.create(parent, recursive = TRUE)
  tmp <- paste0(
    path, ".tmp_", Sys.getpid(), "_",
    format(Sys.time(), "%Y%m%d%H%M%OS6")
  )
  saveRDS(object, tmp)
  if (!file.rename(tmp, path)) {
    unlink(tmp)
    stop("Atomic RDS rename failed: ", path)
  }
  invisible(path)
}

read_tsv <- function(path) {
  if (!file.exists(path)) stop("Missing file: ", path)
  utils::read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
}

mcse_prop <- function(p, n) {
  p <- as.numeric(p)
  n <- as.numeric(n)
  sqrt(pmax(p * (1 - p), 0) / n)
}

paired_mean_se <- function(x, y) {
  d <- as.numeric(y) - as.numeric(x)
  n <- length(d)
  list(
    mean_x = mean(x),
    mean_y = mean(y),
    delta = mean(d),
    se = if (n < 2L) NA_real_ else stats::sd(d) / sqrt(n),
    n = n
  )
}

# Match formal validation: base + replicate_id, and
# base + replicate_id*10 + look_code, with look_code in 0..9.
integer_seed <- function(value, label) {
  value <- as.double(value)
  if (!is.finite(value) || value < 1 || value > .Machine$integer.max) {
    stop("RNG seed overflow or non-finite for ", label, ": ", value)
  }
  as.integer(value)
}

data_seed <- function(base, replicate_id) {
  integer_seed(as.double(base) + as.double(replicate_id), "data_seed")
}

model_seed <- function(base, replicate_id, look_code) {
  look_code <- as.integer(look_code)
  if (is.na(look_code) || look_code < 0L || look_code > 9L) {
    stop("look_code must be a single digit 0-9 to avoid seed collisions; got ",
         look_code)
  }
  integer_seed(
    as.double(base) + as.double(replicate_id) * 10 + as.double(look_code),
    "model_seed"
  )
}

resolve_generator_scenario <- function(sc) {
  if (!is.null(sc$generator_scenario)) {
    sc$generator_scenario
  } else if (identical(sc$role, "null")) {
    "H0"
  } else {
    "H1"
  }
}

look_code_for_n <- function(n.current, N.max, codes) {
  n.current <- as.integer(n.current)
  N.max <- as.integer(N.max)
  if (n.current == N.max) return(as.integer(codes[["final"]]))
  key <- as.character(n.current)
  if (!key %in% names(codes)) {
    stop("No look_code mapped for n=", n.current)
  }
  as.integer(codes[[key]])
}

extract_pp_column <- function(cache, n) {
  pp <- as.matrix(cache$PP)
  col_name <- as.character(as.integer(n))
  if (col_name %in% colnames(pp)) {
    return(as.numeric(pp[, col_name]))
  }
  looks <- as.integer(cache$look_cohort)
  hit <- which(looks == as.integer(n))
  if (length(hit) != 1L) {
    stop("Cannot find PP column for n=", n, " in look_cohort=",
         paste(looks, collapse = "/"))
  }
  as.numeric(pp[, hit])
}

replay_one <- function(pp_row, look_cohort, N.max, c_F, gamma_E,
                       enable_superiority = FALSE) {
  look_cohort <- as.integer(look_cohort)
  N.max <- as.integer(N.max)
  decision <- NA_character_
  go <- NA
  n_final <- NA_real_
  stop_stage <- NA_character_
  for (j in seq_along(look_cohort)) {
    n.current <- look_cohort[j]
    is_final_look <- (n.current == N.max)
    probs <- as.numeric(pp_row[[j]])
    ratio <- n.current / N.max
    fut_bound <- c_F * (ratio^gamma_E)
    sup_bound <- if (isTRUE(enable_superiority)) 1 - ratio^gamma_E else Inf
    stage_name <- if (is_final_look) {
      "FINAL"
    } else {
      paste0("INTERIM_", j)
    }
    if (!is_final_look) {
      if (!is.na(probs) && probs <= fut_bound) {
        decision <- "futility"
        go <- FALSE
        n_final <- n.current
        stop_stage <- stage_name
        break
      } else if (!is.na(probs) && probs >= sup_bound) {
        decision <- "superiority"
        go <- TRUE
        n_final <- n.current
        stop_stage <- stage_name
        break
      }
    } else {
      n_final <- n.current
      stop_stage <- "FINAL"
      if (!is.na(probs) && probs > c_F) {
        decision <- "go"
        go <- TRUE
      } else {
        decision <- "no_go"
        go <- FALSE
      }
      break
    }
  }
  list(
    go = isTRUE(go),
    decision = decision,
    pts.stop = n_final,
    stop_stage = stop_stage
  )
}

replay_schedule <- function(pp_mat, look_cohort, N.max, c_F, gamma_E,
                            enable_superiority = FALSE) {
  nsim <- nrow(pp_mat)
  go_vec <- logical(nsim)
  fut_vec <- logical(nsim)
  sup_vec <- logical(nsim)
  nogo_vec <- logical(nsim)
  ess_vec <- numeric(nsim)
  stage_vec <- character(nsim)
  for (i in seq_len(nsim)) {
    r <- replay_one(
      pp_mat[i, ], look_cohort, N.max, c_F, gamma_E, enable_superiority
    )
    go_vec[i] <- r$go
    fut_vec[i] <- identical(r$decision, "futility")
    sup_vec[i] <- identical(r$decision, "superiority")
    nogo_vec[i] <- identical(r$decision, "no_go")
    ess_vec[i] <- r$pts.stop
    stage_vec[i] <- r$stop_stage
  }
  list(
    go = go_vec,
    futility = fut_vec,
    superiority = sup_vec,
    final_nogo = nogo_vec,
    n_used = ess_vec,
    stop_stage = stage_vec,
    oc = list(
      P_go = mean(go_vec),
      PET_futility = mean(fut_vec),
      P_early_success = mean(sup_vec),
      P_final_nogo = mean(nogo_vec),
      ESS = mean(ess_vec)
    )
  )
}

oc_row <- function(role, oc, L_sim) {
  Type1_error <- if (identical(role, "null")) oc$P_go else NA_real_
  Power <- if (identical(role, "alternative")) oc$P_go else NA_real_
  data.frame(
    role = role,
    Type1_error = Type1_error,
    Power = Power,
    P_go = oc$P_go,
    PET_futility = oc$PET_futility,
    P_early_success = oc$P_early_success,
    P_final_nogo = oc$P_final_nogo,
    ESS = oc$ESS,
    L_sim = L_sim,
    stringsAsFactors = FALSE
  )
}

compare_oc_to_official <- function(got, official, tol = 1e-10) {
  keys <- c("P_go", "PET_futility", "P_final_nogo", "ESS")
  diffs <- vapply(keys, function(k) {
    as.numeric(got[[k]]) - as.numeric(official[[k]])
  }, numeric(1))
  list(
    ok = all(abs(diffs) <= tol),
    diffs = diffs,
    tol = tol
  )
}

analysis_time <- function(n.current, N.max, recruit.int, final_analysis_time) {
  enrollment_time <- recruit.int * (n.current - 1)
  Time.current <- if (n.current == N.max && !is.null(final_analysis_time)) {
    final_analysis_time
  } else {
    enrollment_time
  }
  if (Time.current < enrollment_time) {
    stop("Analysis time precedes current-cohort enrollment.")
  }
  list(
    Time.entry = recruit.int * (seq_len(n.current) - 1),
    Time.current = Time.current,
    enrollment_time = enrollment_time,
    analysis_rows = if (n.current == N.max) {
      seq_len(n.current)
    } else {
      seq_len(n.current - 1L)
    }
  )
}

compute_pp_one_look <- function(myData, n.current, N.max, look_code,
                                Sigma.ctrl_use, model_seed_base, replicate_id,
                                recruit.int, final_analysis_time,
                                m0, L0, v0, S0, time_max,
                                eta_interim, eta_final,
                                thin_MCMC, Niter) {
  timing <- analysis_time(n.current, N.max, recruit.int, final_analysis_time)
  is_final_look <- (as.integer(n.current) == as.integer(N.max))
  temp <- myData[timing$analysis_rows, , drop = FALSE]
  Time.entry.analysis <- timing$Time.entry[timing$analysis_rows]
  eta_use <- if (is_final_look) eta_final else eta_interim
  set.seed(model_seed(model_seed_base, replicate_id, look_code))
  result <- winratio(
    currentdd = temp,
    n_current = n.current,
    Time_current = timing$Time.current,
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
  as.numeric(result$probs)
}

look50_checkpoint_signature <- function(indices, case_key, scenario,
                                        n_50, look_code_50, N.max, Niter,
                                        thin_MCMC, eta_interim, eta_final,
                                        recruit.int, final_analysis_time,
                                        cens_upper, analysis_row_rule,
                                        model_seed_base, data_seed_base,
                                        generation) {
  list(
    checkpoint_format_version = checkpoint_format_version,
    indices = as.integer(indices),
    case_key = as.character(case_key),
    scenario = as.character(scenario),
    n_50 = as.integer(n_50),
    look_code_50 = as.integer(look_code_50),
    N.max = as.integer(N.max),
    Niter = as.integer(Niter),
    thin_MCMC = as.integer(thin_MCMC),
    eta_interim = as.numeric(eta_interim),
    eta_final = as.numeric(eta_final),
    recruit.int = as.numeric(recruit.int),
    final_analysis_time = final_analysis_time,
    cens_upper = cens_upper,
    analysis_row_rule = analysis_row_rule,
    seed_namespace = "SEED-IV-V1-20260731",
    data_seed_base = as.integer(data_seed_base),
    model_seed_base = as.integer(model_seed_base),
    generation = generation
  )
}

pp_checkpoint_compatible <- function(checkpoint, expected_signature) {
  if (!is.list(checkpoint) || is.null(checkpoint$signature) ||
      is.null(checkpoint$PP)) {
    return(FALSE)
  }
  if (!isTRUE(all.equal(
    checkpoint$signature, expected_signature,
    tolerance = 1e-12, check.attributes = FALSE
  ))) {
    return(FALSE)
  }
  is.matrix(checkpoint$PP) &&
    nrow(checkpoint$PP) == length(expected_signature$indices) &&
    ncol(checkpoint$PP) == 1L &&
    all(is.finite(checkpoint$PP))
}

formal_cache_path <- function(repo_root, case_key, scenario_cache_file) {
  file.path(
    repo_root,
    "output", "validation", case_key, "formal_10000", "cache",
    scenario_cache_file
  )
}

official_oc_path <- function(repo_root, case_key) {
  file.path(
    repo_root,
    "output", "validation", case_key, "formal_10000", "results",
    "validation_LG_oc_latest.csv"
  )
}
