# Chunked, reproducible PP-cache construction for long Win Ratio calibration.
#
# The calling calibration script provides the design objects, the simulated
# data array, and compute_pp_one_trial(). Each completed chunk is written
# atomically. Re-running the same script reuses compatible chunks and computes
# only missing ones.

checkpoint_format_version <- "winratio_pp_chunk_v2_final_all_n"

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

checkpoint_global <- function(name, default = NULL) {
  if (exists(name, envir = .GlobalEnv, inherits = FALSE)) {
    get(name, envir = .GlobalEnv, inherits = FALSE)
  } else {
    default
  }
}

scenario_data_fingerprint <- function(label) {
  suffix <- if (identical(label, "H0")) "h0" else "h1"
  list(
    mu_trt = checkpoint_global(paste0("mu_", suffix, "_trt")),
    mu_ctrl = checkpoint_global(paste0("mu_", suffix, "_ctrl")),
    Sigma_trt = checkpoint_global(paste0("Sigma_", suffix, "_trt")),
    Sigma_ctrl = checkpoint_global("Sigma.ctrl"),
    cens_upper = checkpoint_global("cens_upper")
  )
}

pp_checkpoint_signature <- function(
    indices, label, generation_fingerprint, scenario_seed_offset) {
  list(
    checkpoint_format_version = checkpoint_format_version,
    indices = as.integer(indices),
    label = as.character(label),
    calibration_seed = as.integer(seed),
    scenario_seed_offset = as.integer(scenario_seed_offset),
    rng_scheme = "fixed_seed_per_trial_lecuyer_v1",
    look_cohort = as.integer(look_cohort),
    N.max = as.integer(N.max),
    Niter = as.integer(Niter),
    thin_MCMC = as.integer(thin_MCMC),
    eta_interim = as.numeric(eta_interim),
    eta_final = as.numeric(eta_final),
    recruit.int = as.numeric(recruit.int),
    final_analysis_time = final_analysis_time,
    analysis_row_rule = checkpoint_global(
      "analysis_row_rule", "interim_n_minus_1_final_all_n_v2"
    ),
    prior = list(
      m0 = checkpoint_global("m0"),
      L0 = checkpoint_global("L0"),
      v0 = checkpoint_global("v0"),
      S0 = checkpoint_global("S0"),
      time_max = checkpoint_global("time_max")
    ),
    scenario_data = scenario_data_fingerprint(label),
    generation_fingerprint = generation_fingerprint
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
    ncol(checkpoint$PP) == length(expected_signature$look_cohort) &&
    all(is.finite(checkpoint$PP))
}

trial_mcmc_seed <- function(index, scenario_seed_offset) {
  # R integer seeds must stay below .Machine$integer.max.
  value <- (
    as.double(seed) + as.double(scenario_seed_offset) + as.double(index)
  ) %% 2147483646
  as.integer(value + 1)
}

build_pp_matrix_checkpointed <- function(
    data_arr, nsim, label, generation_fingerprint = NULL) {
  chunk_size <- suppressWarnings(as.integer(
    Sys.getenv("CALIBRATION_CHUNK_SIZE", unset = "140")
  ))
  if (is.na(chunk_size) || chunk_size < 1L) {
    stop("CALIBRATION_CHUNK_SIZE must be a positive integer.")
  }

  scenario_seed_offset <- if (identical(label, "H0")) {
    100000000L
  } else if (identical(label, "H1")) {
    200000000L
  } else {
    stop("Unsupported calibration scenario label: ", label)
  }

  checkpoint_dir <- file.path(cache_dir, "checkpoints", label)
  if (!dir.exists(checkpoint_dir)) {
    dir.create(checkpoint_dir, recursive = TRUE)
  }

  indices <- seq_len(nsim)
  chunk_id <- ceiling(indices / chunk_size)
  chunks <- split(indices, chunk_id)
  pp_chunks <- vector("list", length(chunks))

  cat(sprintf(
    "  Checkpointing every %d trials (%d chunks); directory: %s\n",
    chunk_size, length(chunks), checkpoint_dir
  ))

  for (i in seq_along(chunks)) {
    idx <- as.integer(chunks[[i]])
    checkpoint_path <- file.path(
      checkpoint_dir,
      sprintf("chunk_%06d_%06d.rds", min(idx), max(idx))
    )
    expected_signature <- pp_checkpoint_signature(
      idx, label, generation_fingerprint, scenario_seed_offset
    )

    checkpoint <- NULL
    if (!isTRUE(rebuild_pp_cache) && file.exists(checkpoint_path)) {
      checkpoint <- tryCatch(readRDS(checkpoint_path), error = function(e) NULL)
    }

    if (pp_checkpoint_compatible(checkpoint, expected_signature)) {
      pp_chunks[[i]] <- checkpoint$PP
      cat(sprintf(
        "  [%s] reused chunk %d/%d (%d-%d)\n",
        label, i, length(chunks), min(idx), max(idx)
      ))
      next
    }

    chunk_started <- Sys.time()
    pp_list <- parallel::mclapply(
      idx,
      function(x) {
        set.seed(
          trial_mcmc_seed(x, scenario_seed_offset),
          kind = "L'Ecuyer-CMRG"
        )
        compute_pp_one_trial(data_arr[, , x])
      },
      mc.cores = min(as.integer(mc.cores), length(idx)),
      mc.set.seed = FALSE
    )

    bad <- vapply(pp_list, function(x) {
      inherits(x, "try-error") || !is.numeric(x) ||
        length(x) != n_looks || any(!is.finite(x))
    }, logical(1))
    if (any(bad)) {
      stop(
        "Invalid PP result(s) in ", label, " chunk ",
        min(idx), "-", max(idx), "; trial indices: ",
        paste(idx[bad], collapse = ", ")
      )
    }

    pp_mat <- do.call(rbind, pp_list)
    colnames(pp_mat) <- as.character(look_cohort)
    rownames(pp_mat) <- paste0("sim", idx)
    checkpoint <- list(
      signature = expected_signature,
      PP = pp_mat,
      completed_at = as.character(Sys.time())
    )
    atomic_save_rds(checkpoint, checkpoint_path)
    pp_chunks[[i]] <- pp_mat
    elapsed <- as.numeric(difftime(
      Sys.time(), chunk_started, units = "mins"
    ))
    cat(sprintf(
      "  [%s] saved chunk %d/%d (%d-%d) in %.1f min\n",
      label, i, length(chunks), min(idx), max(idx), elapsed
    ))
  }

  pp_mat <- do.call(rbind, pp_chunks)
  rownames(pp_mat) <- paste0("sim", seq_len(nsim))
  colnames(pp_mat) <- as.character(look_cohort)
  pp_mat
}
