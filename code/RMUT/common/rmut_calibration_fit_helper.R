# Two-stage Stan fitting for RMUT calibration-cache generation.
#
# The same simulated dataset is retained across attempts. A retry changes only
# sampler settings and the model seed. Every component fit returns a structured
# result containing its attempt count, diagnostics and final failure reason.

rmut_chain_diagnostics <- function(fit) {
  raw_draws <- try(
    rstan::extract(
      fit$fit,
      permuted = FALSE,
      inc_warmup = FALSE
    ),
    silent = TRUE
  )
  if (inherits(raw_draws, "try-error")) {
    return(list(
      ok = FALSE,
      max_rhat = NA_real_,
      min_ess_bulk = NA_real_,
      min_ess_tail = NA_real_
    ))
  }

  draws_array <- posterior::as_draws_array(raw_draws)
  relevant_variables <- grep(
    "^(b_|b\\[|shape$|sigma$)",
    posterior::variables(draws_array),
    value = TRUE
  )
  if (length(relevant_variables) == 0L) {
    return(list(
      ok = FALSE,
      max_rhat = NA_real_,
      min_ess_bulk = NA_real_,
      min_ess_tail = NA_real_
    ))
  }

  diagnostics <- try(
    posterior::summarise_draws(
      posterior::subset_draws(
        draws_array,
        variable = relevant_variables
      ),
      rhat = posterior::rhat,
      ess_bulk = posterior::ess_bulk,
      ess_tail = posterior::ess_tail
    ),
    silent = TRUE
  )
  if (inherits(diagnostics, "try-error") ||
      !all(c("rhat", "ess_bulk", "ess_tail") %in% names(diagnostics))) {
    return(list(
      ok = FALSE,
      max_rhat = NA_real_,
      min_ess_bulk = NA_real_,
      min_ess_tail = NA_real_
    ))
  }

  rhat_values <- diagnostics$rhat[is.finite(diagnostics$rhat)]
  ess_bulk_values <- diagnostics$ess_bulk[
    is.finite(diagnostics$ess_bulk)
  ]
  ess_tail_values <- diagnostics$ess_tail[
    is.finite(diagnostics$ess_tail)
  ]
  max_rhat <- if (length(rhat_values)) max(rhat_values) else NA_real_
  min_ess_bulk <- if (length(ess_bulk_values)) {
    min(ess_bulk_values)
  } else {
    NA_real_
  }
  min_ess_tail <- if (length(ess_tail_values)) {
    min(ess_tail_values)
  } else {
    NA_real_
  }

  list(
    ok = is.finite(max_rhat) &&
      is.finite(min_ess_bulk) &&
      is.finite(min_ess_tail) &&
      max_rhat <= 1.01 &&
      min_ess_bulk >= 400 &&
      min_ess_tail >= 400,
    max_rhat = max_rhat,
    min_ess_bulk = min_ess_bulk,
    min_ess_tail = min_ess_tail
  )
}

rmut_sampler_diagnostics <- function(fit, max_treedepth) {
  sampler_params <- try(
    rstan::get_sampler_params(fit$fit, inc_warmup = FALSE),
    silent = TRUE
  )
  if (inherits(sampler_params, "try-error") ||
      length(sampler_params) == 0L) {
    return(list(
      ok = FALSE,
      divergent_count = NA_integer_,
      treedepth_hit_count = NA_integer_,
      max_rhat = NA_real_,
      min_ess_bulk = NA_real_,
      min_ess_tail = NA_real_
    ))
  }

  required_columns <- c("divergent__", "treedepth__")
  columns_ok <- all(vapply(
    sampler_params,
    function(x) all(required_columns %in% colnames(x)),
    logical(1)
  ))
  if (!columns_ok) {
    return(list(
      ok = FALSE,
      divergent_count = NA_integer_,
      treedepth_hit_count = NA_integer_,
      max_rhat = NA_real_,
      min_ess_bulk = NA_real_,
      min_ess_tail = NA_real_
    ))
  }

  divergent_count <- sum(vapply(
    sampler_params,
    function(x) sum(x[, "divergent__"] > 0),
    numeric(1)
  ))
  treedepth_hit_count <- sum(vapply(
    sampler_params,
    function(x) sum(x[, "treedepth__"] >= max_treedepth),
    numeric(1)
  ))
  chain_diagnostics <- rmut_chain_diagnostics(fit)

  list(
    ok = divergent_count == 0L &&
      treedepth_hit_count == 0L &&
      isTRUE(chain_diagnostics$ok),
    divergent_count = as.integer(divergent_count),
    treedepth_hit_count = as.integer(treedepth_hit_count),
    max_rhat = chain_diagnostics$max_rhat,
    min_ess_bulk = chain_diagnostics$min_ess_bulk,
    min_ess_tail = chain_diagnostics$min_ess_tail
  )
}

rmut_format_diagnostic_failure <- function(label, diagnostics) {
  sprintf(
    paste0(
      "%s_sampler_diagnostics_failed",
      "_divergent_%s_treedepth_%s",
      "_rhat_%s_essbulk_%s_esstail_%s"
    ),
    label,
    diagnostics$divergent_count,
    diagnostics$treedepth_hit_count,
    format(diagnostics$max_rhat, digits = 4),
    format(diagnostics$min_ess_bulk, digits = 6),
    format(diagnostics$min_ess_tail, digits = 6)
  )
}

rmut_fit_with_stable_retry <- function(base_model, newdata, seed) {
  attempts <- list(
    list(
      label = "base",
      iter = bayes_iter,
      warmup = bayes_warmup,
      adapt_delta = bayes_adapt,
      max_treedepth = bayes_treedepth,
      seed_offset = 0L
    ),
    list(
      label = "retry",
      iter = bayes_retry_iter,
      warmup = bayes_retry_warmup,
      adapt_delta = bayes_retry_adapt,
      max_treedepth = bayes_retry_treedepth,
      seed_offset = 500000L
    )
  )

  last_reason <- "not_attempted"
  retry_trigger_reason <- ""
  last_diagnostics <- list(
    divergent_count = NA_integer_,
    treedepth_hit_count = NA_integer_,
    max_rhat = NA_real_,
    min_ess_bulk = NA_real_,
    min_ess_tail = NA_real_
  )

  for (attempt_index in seq_along(attempts)) {
    attempt <- attempts[[attempt_index]]
    attempt_seed <- as.integer(
      (as.double(seed) + attempt$seed_offset) %%
        .Machine$integer.max
    )
    fit <- try(
      update(
        base_model,
        newdata = newdata,
        chains = bayes_chains,
        cores = 1,
        iter = attempt$iter,
        warmup = attempt$warmup,
        seed = attempt_seed,
        refresh = 0,
        silent = 2,
        control = list(
          adapt_delta = attempt$adapt_delta,
          max_treedepth = attempt$max_treedepth
        )
      ),
      silent = TRUE
    )

    if (inherits(fit, "try-error")) {
      last_reason <- paste0(attempt$label, "_fit_error")
      if (attempt_index == 1L) {
        retry_trigger_reason <- last_reason
      }
      next
    }

    diagnostics <- rmut_sampler_diagnostics(
      fit,
      attempt$max_treedepth
    )
    last_diagnostics <- diagnostics
    if (!isTRUE(diagnostics$ok)) {
      last_reason <- rmut_format_diagnostic_failure(
        attempt$label,
        diagnostics
      )
      if (attempt_index == 1L) {
        retry_trigger_reason <- last_reason
      }
      next
    }

    draws <- extract_valid_aft_draws(fit)
    if (is.null(draws)) {
      last_reason <- paste0(attempt$label, "_posterior_extraction_failed")
      if (attempt_index == 1L) {
        retry_trigger_reason <- last_reason
      }
      next
    }

    return(list(
      ok = TRUE,
      fit = fit,
      draws = draws,
      attempt = attempt$label,
      attempts_used = as.integer(attempt_index),
      failure_reason = "",
      retry_trigger_reason = retry_trigger_reason,
      diagnostics = diagnostics
    ))
  }

  list(
    ok = FALSE,
    fit = NULL,
    draws = NULL,
    attempt = "failed_after_retry",
    attempts_used = as.integer(length(attempts)),
    failure_reason = last_reason,
    retry_trigger_reason = retry_trigger_reason,
    diagnostics = last_diagnostics
  )
}

rmut_fit_record <- function(result, stage, component) {
  diagnostics <- result$diagnostics
  data.frame(
    stage = as.character(stage),
    component = as.character(component),
    ok = isTRUE(result$ok),
    n_draws = if (is.null(result$draws)) 0L else nrow(result$draws),
    attempt = as.character(result$attempt),
    attempts_used = as.integer(result$attempts_used),
    failure_reason = as.character(result$failure_reason),
    retry_trigger_reason = as.character(result$retry_trigger_reason),
    divergent_count = as.integer(diagnostics$divergent_count),
    treedepth_hit_count = as.integer(
      diagnostics$treedepth_hit_count
    ),
    max_rhat = as.numeric(diagnostics$max_rhat),
    min_ess_bulk = as.numeric(diagnostics$min_ess_bulk),
    min_ess_tail = as.numeric(diagnostics$min_ess_tail),
    stringsAsFactors = FALSE
  )
}

rmut_empty_fit_records <- function() {
  data.frame(
    stage = character(),
    component = character(),
    ok = logical(),
    attempt = character(),
    attempts_used = integer(),
    failure_reason = character(),
    retry_trigger_reason = character(),
    divergent_count = integer(),
    treedepth_hit_count = integer(),
    max_rhat = numeric(),
    min_ess_bulk = numeric(),
    min_ess_tail = numeric(),
    stringsAsFactors = FALSE
  )
}

rmut_make_replicate_result <- function(
  attempt_id,
  interim_pp,
  pp_final,
  fit_records
) {
  records <- if (length(fit_records)) {
    do.call(rbind, fit_records)
  } else {
    rmut_empty_fit_records()
  }
  pp_ok <- length(interim_pp) > 0L &&
    all(is.finite(interim_pp)) &&
    length(pp_final) == 1L &&
    is.finite(pp_final)
  records_ok <- nrow(records) > 0L && all(records$ok)
  fit_ok <- pp_ok && records_ok

  failure_parts <- if (nrow(records)) {
    failed_records <- records[!records$ok, , drop = FALSE]
    if (nrow(failed_records)) {
      sprintf(
        "%s:%s:%s",
        failed_records$stage,
        failed_records$component,
        failed_records$failure_reason
      )
    } else {
      character()
    }
  } else {
    character()
  }
  if (!pp_ok && !length(failure_parts)) {
    failure_parts <- "posterior_probability_nonfinite"
  }

  list(
    attempt_id = as.integer(attempt_id),
    fit_ok = fit_ok,
    failure_reason = paste(failure_parts, collapse = " | "),
    fit_attempts = if (nrow(records)) {
      sum(records$attempts_used)
    } else {
      0L
    },
    n_fit_calls = nrow(records),
    n_base_success = if (nrow(records)) {
      sum(records$ok & records$attempt == "base")
    } else {
      0L
    },
    n_retry_success = if (nrow(records)) {
      sum(records$ok & records$attempt == "retry")
    } else {
      0L
    },
    n_failed_fit_calls = if (nrow(records)) sum(!records$ok) else 0L,
    PP = as.numeric(interim_pp),
    PP_final = as.numeric(pp_final),
    fit_records = records
  )
}

rmut_pp_from_fit_results <- function(
  death_result,
  composite_result,
  tau,
  weight,
  pair_seed = NULL
) {
  if (!isTRUE(death_result$ok) || !isTRUE(composite_result$ok)) {
    return(NA_real_)
  }
  draws_d <- death_result$draws
  draws_c <- composite_result$draws
  required_columns <- c("b_Intercept", "b_group", "sigma")
  if (is.null(draws_d) ||
      is.null(draws_c) ||
      !all(required_columns %in% colnames(draws_d)) ||
      !all(required_columns %in% colnames(draws_c))) {
    return(NA_real_)
  }

  rmut_ctrl <- compute_rmut_samples(
    mu_D = draws_d$b_Intercept,
    sigma_D = draws_d$sigma,
    mu_C = draws_c$b_Intercept,
    sigma_C = draws_c$sigma,
    tau = tau,
    w_R = weight,
    pair_seed = pair_seed
  )
  rmut_trt <- compute_rmut_samples(
    mu_D = draws_d$b_Intercept + draws_d$b_group,
    sigma_D = draws_d$sigma,
    mu_C = draws_c$b_Intercept + draws_c$b_group,
    sigma_C = draws_c$sigma,
    tau = tau,
    w_R = weight,
    pair_seed = pair_seed
  )
  pp <- mean((rmut_trt - rmut_ctrl) > 0)
  if (length(pp) == 1L && is.finite(pp)) pp else NA_real_
}

rmut_cache_is_current <- function(stats_list, n_expected, n_interims) {
  if (!is.list(stats_list) || length(stats_list) != n_expected) {
    return(FALSE)
  }
  required_fields <- c(
    "attempt_id",
    "fit_ok",
    "failure_reason",
    "fit_attempts",
    "n_fit_calls",
    "n_base_success",
    "n_retry_success",
    "n_failed_fit_calls",
    "PP",
    "PP_final",
    "fit_records"
  )
  all(vapply(
    stats_list,
    function(item) {
      is.list(item) &&
        all(required_fields %in% names(item)) &&
        length(item$PP) == n_interims &&
        is.data.frame(item$fit_records)
    },
    logical(1)
  ))
}

rmut_write_scenario_audit <- function(
  stats_list,
  output_dir,
  scenario_id,
  scenario_label,
  metrics
) {
  replicate_rows <- lapply(stats_list, function(item) {
    pp_values <- as.list(item$PP)
    names(pp_values) <- paste0("PP_interim_", seq_along(item$PP))
    data.frame(
      attempt_id = item$attempt_id,
      fit_ok = item$fit_ok,
      failure_reason = item$failure_reason,
      fit_attempts = item$fit_attempts,
      n_fit_calls = item$n_fit_calls,
      n_base_success = item$n_base_success,
      n_retry_success = item$n_retry_success,
      n_failed_fit_calls = item$n_failed_fit_calls,
      as.data.frame(pp_values, check.names = FALSE),
      PP_final = item$PP_final,
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
  })
  replicate_table <- do.call(rbind, replicate_rows)
  write.csv(
    replicate_table,
    file.path(
      output_dir,
      sprintf("calibration_replicates_%s.csv", scenario_id)
    ),
    row.names = FALSE,
    na = ""
  )

  failure_rows <- lapply(stats_list, function(item) {
    records <- item$fit_records
    if (!nrow(records)) {
      return(NULL)
    }
    records <- records[!records$ok, , drop = FALSE]
    if (!nrow(records)) {
      return(NULL)
    }
    data.frame(
      attempt_id = item$attempt_id,
      records,
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
  })
  failure_rows <- Filter(Negate(is.null), failure_rows)
  failure_table <- if (length(failure_rows)) {
    do.call(rbind, failure_rows)
  } else {
    data.frame(
      attempt_id = integer(),
      rmut_empty_fit_records(),
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
  }
  write.csv(
    failure_table,
    file.path(
      output_dir,
      sprintf("calibration_failures_%s.csv", scenario_id)
    ),
    row.names = FALSE,
    na = ""
  )

  n_attempted <- length(stats_list)
  n_valid <- sum(vapply(
    stats_list,
    function(item) isTRUE(item$fit_ok),
    logical(1)
  ))
  summary_row <- data.frame(
    scenario_id = scenario_id,
    scenario_label = scenario_label,
    n_attempted = n_attempted,
    n_valid = n_valid,
    n_invalid = n_attempted - n_valid,
    valid_rate = n_valid / n_attempted,
    n_component_fit_calls = sum(vapply(
      stats_list,
      function(item) item$n_fit_calls,
      numeric(1)
    )),
    n_mcmc_fit_attempts = sum(vapply(
      stats_list,
      function(item) item$fit_attempts,
      numeric(1)
    )),
    n_base_success = sum(vapply(
      stats_list,
      function(item) item$n_base_success,
      numeric(1)
    )),
    n_retry_success = sum(vapply(
      stats_list,
      function(item) item$n_retry_success,
      numeric(1)
    )),
    n_failed_after_retry = sum(vapply(
      stats_list,
      function(item) item$n_failed_fit_calls,
      numeric(1)
    )),
    n_base_failures_triggering_retry = sum(vapply(
      stats_list,
      function(item) {
        if (nrow(item$fit_records)) {
          sum(item$fit_records$attempts_used > 1L)
        } else {
          0L
        }
      },
      numeric(1)
    )),
    n_replicates_with_retry = sum(vapply(
      stats_list,
      function(item) {
        nrow(item$fit_records) > 0L &&
          any(item$fit_records$attempts_used > 1L)
      },
      logical(1)
    )),
    PET = metrics$PET,
    ESS = metrics$ESS,
    reference_success_rate = metrics$success_rate,
    stringsAsFactors = FALSE
  )
  write.csv(
    summary_row,
    file.path(
      output_dir,
      sprintf("calibration_summary_%s.csv", scenario_id)
    ),
    row.names = FALSE,
    na = ""
  )
  invisible(summary_row)
}
