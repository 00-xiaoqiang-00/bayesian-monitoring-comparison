# Stable two-stage Weibull fitting for RMST calibration and validation
#
# The simulated dataset is unchanged across attempts. A retry changes only
# sampler settings and the model seed. A fit is accepted only when:
#   - Stan reports no divergent transitions;
#   - no transition reaches max_treedepth;
#   - relevant model parameters have R-hat <= 1.01 and bulk ESS >= 400;
#   - posterior RMST input draws can be extracted and are finite.

rmst_chain_diagnostics <- function(fit) {
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
      min_ess_bulk = NA_real_
    ))
  }
  draws_array <- posterior::as_draws_array(raw_draws)

  variable_names <- posterior::variables(draws_array)
  relevant_variables <- grep(
    "^(b_|b\\[|shape$|sigma$)",
    variable_names,
    value = TRUE
  )
  if (length(relevant_variables) == 0L) {
    return(list(
      ok = FALSE,
      max_rhat = NA_real_,
      min_ess_bulk = NA_real_
    ))
  }

  diagnostics <- try(
    posterior::summarise_draws(
      posterior::subset_draws(
        draws_array,
        variable = relevant_variables
      ),
      rhat = posterior::rhat,
      ess_bulk = posterior::ess_bulk
    ),
    silent = TRUE
  )
  if (inherits(diagnostics, "try-error") ||
      !all(c("rhat", "ess_bulk") %in% names(diagnostics))) {
    return(list(
      ok = FALSE,
      max_rhat = NA_real_,
      min_ess_bulk = NA_real_
    ))
  }

  rhat_values <- diagnostics$rhat[is.finite(diagnostics$rhat)]
  ess_values <- diagnostics$ess_bulk[is.finite(diagnostics$ess_bulk)]
  max_rhat <- if (length(rhat_values)) max(rhat_values) else NA_real_
  min_ess_bulk <- if (length(ess_values)) min(ess_values) else NA_real_

  list(
    ok = is.finite(max_rhat) &&
      is.finite(min_ess_bulk) &&
      max_rhat <= 1.01 &&
      min_ess_bulk >= 400,
    max_rhat = max_rhat,
    min_ess_bulk = min_ess_bulk
  )
}

rmst_sampler_diagnostics <- function(fit, max_treedepth) {
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
      min_ess_bulk = NA_real_
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
      min_ess_bulk = NA_real_
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

  chain_diagnostics <- rmst_chain_diagnostics(fit)
  list(
    ok = divergent_count == 0L &&
      treedepth_hit_count == 0L &&
      isTRUE(chain_diagnostics$ok),
    divergent_count = as.integer(divergent_count),
    treedepth_hit_count = as.integer(treedepth_hit_count),
    max_rhat = chain_diagnostics$max_rhat,
    min_ess_bulk = chain_diagnostics$min_ess_bulk
  )
}

rmst_extract_valid_draws <- function(fit) {
  newdata <- data.frame(group = c(0, 1))
  mu_draws <- try(
    posterior_linpred(
      fit,
      newdata = newdata,
      transform = TRUE,
      re_formula = NA
    ),
    silent = TRUE
  )
  alpha_draws <- try(
    posterior_linpred(
      fit,
      dpar = "shape",
      transform = TRUE
    ),
    silent = TRUE
  )

  if (inherits(mu_draws, "try-error") ||
      inherits(alpha_draws, "try-error") ||
      !is.matrix(mu_draws) ||
      ncol(mu_draws) < 2L ||
      !is.matrix(alpha_draws) ||
      ncol(alpha_draws) < 1L) {
    return(NULL)
  }

  draws <- list(
    mu0 = as.numeric(mu_draws[, 1]),
    mu1 = as.numeric(mu_draws[, 2]),
    alpha = as.numeric(alpha_draws[, 1])
  )
  draws_ok <- length(draws$mu0) > 0L &&
    length(draws$mu0) == length(draws$mu1) &&
    length(draws$mu0) == length(draws$alpha) &&
    all(is.finite(draws$mu0)) &&
    all(is.finite(draws$mu1)) &&
    all(is.finite(draws$alpha)) &&
    all(draws$mu0 > 0) &&
    all(draws$mu1 > 0) &&
    all(draws$alpha > 0)

  if (!draws_ok) NULL else draws
}

rmst_fit_with_stable_retry <- function(base_model, newdata, seed) {
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
  for (attempt in attempts) {
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
      next
    }

    diagnostics <- rmst_sampler_diagnostics(
      fit,
      attempt$max_treedepth
    )
    if (!isTRUE(diagnostics$ok)) {
      last_reason <- sprintf(
        paste0(
          "%s_sampler_diagnostics_failed_divergent_%s",
          "_treedepth_%s_rhat_%s_essbulk_%s"
        ),
        attempt$label,
        diagnostics$divergent_count,
        diagnostics$treedepth_hit_count,
        format(diagnostics$max_rhat, digits = 4),
        format(diagnostics$min_ess_bulk, digits = 6)
      )
      next
    }

    draws <- rmst_extract_valid_draws(fit)
    if (is.null(draws)) {
      last_reason <- paste0(attempt$label, "_posterior_extraction_failed")
      next
    }

    return(list(
      ok = TRUE,
      attempt = attempt$label,
      failure_reason = "",
      draws = draws
    ))
  }

  list(
    ok = FALSE,
    attempt = "failed_after_retry",
    failure_reason = last_reason,
    draws = NULL
  )
}
