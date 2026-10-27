# Stable two-stage Stan fitting for RMUT independent validation.
# The same simulated dataset is retained across attempts. Only the numerical
# sampler settings and retry seed change.

rmut_sampler_diagnostics_ok <- function(fit, max_treedepth) {
  sampler_params <- try(
    rstan::get_sampler_params(fit$fit, inc_warmup = FALSE),
    silent = TRUE
  )
  if (inherits(sampler_params, "try-error") ||
      length(sampler_params) == 0L) {
    return(FALSE)
  }

  required_columns <- c("divergent__", "treedepth__")
  if (any(vapply(
    sampler_params,
    function(x) !all(required_columns %in% colnames(x)),
    logical(1)
  ))) {
    return(FALSE)
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

  divergent_count == 0L && treedepth_hit_count == 0L
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

  audit_attempts <- list()
  for (attempt_index in seq_along(attempts)) {
    attempt <- attempts[[attempt_index]]
    attempt_seed <- as.integer(
      (as.double(seed) + attempt$seed_offset) %% .Machine$integer.max
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

    sampler_ok <- !inherits(fit, "try-error") && rmut_sampler_diagnostics_ok(fit, attempt$max_treedepth)
    accepted_draws <- if (sampler_ok) extract_valid_aft_draws(fit) else NULL
    passed <- sampler_ok && !is.null(accepted_draws)
    audit_attempts[[attempt_index]] <- list(
      attempt = attempt$label, seed = attempt_seed, chains = bayes_chains,
      iter = attempt$iter, warmup = attempt$warmup,
      adapt_delta = attempt$adapt_delta, max_treedepth = attempt$max_treedepth,
      sampler_ok = sampler_ok, posterior_diagnostics_passed = !is.null(accepted_draws),
      posterior_rows = if (!is.null(accepted_draws)) nrow(accepted_draws) else NA_integer_,
      diagnostics = if (!is.null(accepted_draws)) attr(accepted_draws, "rmut_posterior_diagnostics") else NULL
    )
    if (passed) {
      attr(fit, "iv_attempt") <- attempt$label
      attr(fit, "iv_attempt_index") <- attempt_index
      attr(fit, "rmut_fit_audit") <- list(
        passed = TRUE, attempt = attempt$label, posterior_rows = nrow(accepted_draws),
        attempts = audit_attempts)
      return(fit)
    }
  }

  failed <- structure(
    "RMUT AFT fit failed convergence diagnostics after stable retry",
    class = "try-error"
  )
  attr(failed, "rmut_fit_audit") <- list(
    passed = FALSE, attempt = "failed", posterior_rows = NA_integer_, attempts = audit_attempts)
  failed
}
