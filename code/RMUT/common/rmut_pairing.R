# Scalar checks run before integer conversion, which otherwise truncates values
# or turns overflow into NA.
rmut_check_integer_scalar <- function(value, name, minimum = 0,
                                      maximum = .Machine$integer.max) {
  if (!is.numeric(value) || length(value) != 1L || !is.null(dim(value)) ||
      is.na(value) || !is.finite(value) || value < minimum ||
      value > maximum || value != floor(value)) {
    stop(name, " must be one finite integer in [", minimum, ", ", maximum, "]",
         call. = FALSE)
  }
  as.integer(value)
}

# Construct a deterministic pairing seed from replicate, look, scenario and case identifiers.

rmut_pair_seed <- function(replicate_id, look_id, scenario_id=0L,
                           case_id=0L, base_seed=720260922L) {
  replicate_id <- rmut_check_integer_scalar(replicate_id, "replicate_id")
  look_id <- rmut_check_integer_scalar(look_id, "look_id")
  scenario_id <- rmut_check_integer_scalar(scenario_id, "scenario_id")
  case_id <- rmut_check_integer_scalar(case_id, "case_id")
  base_seed <- rmut_check_integer_scalar(base_seed, "base_seed")
  x <- as.double(base_seed) + 10000*as.double(replicate_id) +
    101*as.double(look_id) + 1000000*as.double(scenario_id) +
    100000000*as.double(case_id)
  if (!is.finite(x) || x > 2^53 - 1) {
    stop("Seed identifiers exceed the range of exact integer arithmetic", call. = FALSE)
  }
  x <- x %% (.Machine$integer.max - 1)
  # Map zero at the modulo boundary to positive seed 1.
  if (x == 0) x <- 1
  as.integer(x)
}

# RMUT posterior-row alignment.
# Used by calibration, validation and sensitivity analyses.

rmut_pair_indices <- function(n_death, n_composite, pair_seed = NULL) {
  n_death <- rmut_check_integer_scalar(n_death, "n_death", minimum = 1)
  n_composite <- rmut_check_integer_scalar(n_composite, "n_composite", minimum = 1)
  if (!is.null(pair_seed)) {
    pair_seed <- rmut_check_integer_scalar(pair_seed, "pair_seed")
  }
  if (n_death == n_composite) {
    return(list(death = seq_len(n_death), composite = seq_len(n_composite)))
  }
  if (is.null(pair_seed)) {
    stop("pair_seed is required when posterior lengths differ")
  }

  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  old_kind <- RNGkind()
  on.exit({
    do.call(RNGkind, as.list(old_kind))
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)

  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  set.seed(pair_seed)
  n_pair <- min(n_death, n_composite)
  i_death <- if (n_death == n_pair) seq_len(n_death) else sample.int(n_death, n_pair)
  i_composite <- if (n_composite == n_pair) seq_len(n_composite) else sample.int(n_composite, n_pair)
  list(death = i_death, composite = i_composite)
}

compute_rmut_samples_corrected <- function(
  mu_D, sigma_D, mu_C, sigma_C, tau, w_R, M = 100L, pair_seed = NULL
) {
  parameters <- list(mu_D = mu_D, sigma_D = sigma_D,
                     mu_C = mu_C, sigma_C = sigma_C)
  valid_parameter <- vapply(parameters, function(value) {
    is.numeric(value) && is.null(dim(value)) && length(value) > 0L &&
      all(is.finite(value))
  }, logical(1L))
  if (!all(valid_parameter)) {
    stop("Posterior parameters must be non-empty finite numeric vectors", call. = FALSE)
  }
  n_death <- c(length(mu_D), length(sigma_D))
  n_composite <- c(length(mu_C), length(sigma_C))
  if (any(n_death < 1L) || any(n_composite < 1L) ||
      length(unique(n_death)) != 1L || length(unique(n_composite)) != 1L) {
    stop("Each model's intercept and sigma vectors must have matching non-empty lengths")
  }
  if (n_death[[1L]] != n_composite[[1L]] && is.null(pair_seed)) {
    stop("pair_seed is required when posterior model lengths differ")
  }
  if (!is.numeric(tau) || length(tau) != 1L || !is.null(dim(tau)) ||
      is.na(tau) || !is.finite(tau) || tau <= 0 ||
      !is.numeric(w_R) || length(w_R) != 1L || !is.null(dim(w_R)) ||
      is.na(w_R) || !is.finite(w_R) || w_R < 0 || w_R > 1) {
    stop("Invalid RMUT integration settings")
  }
  M <- rmut_check_integer_scalar(M, "M", minimum = 1)
  if (any(c(sigma_D, sigma_C) <= 0)) {
    stop("Posterior parameters must be finite with positive sigma")
  }

  # Align the two complete posterior rows before integration.
  pair <- rmut_pair_indices(length(mu_D), length(mu_C), pair_seed)
  mu_D <- mu_D[pair$death]
  sigma_D <- sigma_D[pair$death]
  mu_C <- mu_C[pair$composite]
  sigma_C <- sigma_C[pair$composite]
  stopifnot(length(mu_D) == length(mu_C), length(sigma_D) == length(sigma_C))

  n <- length(mu_D)
  t_seq <- tau * (seq_len(M) - 0.5) / M
  log_t_mat <- matrix(log(t_seq), nrow = n, ncol = M, byrow = TRUE)
  mu_D_mat <- matrix(mu_D, nrow = n, ncol = M)
  sigma_D_mat <- matrix(sigma_D, nrow = n, ncol = M)
  mu_C_mat <- matrix(mu_C, nrow = n, ncol = M)
  sigma_C_mat <- matrix(sigma_C, nrow = n, ncol = M)

  S_D <- pnorm(log_t_mat, mean = mu_D_mat, sd = sigma_D_mat, lower.tail = FALSE)
  S_C <- pnorm(log_t_mat, mean = mu_C_mat, sd = sigma_C_mat, lower.tail = FALSE)
  rowSums((1 - w_R) * S_D + w_R * S_C) * tau / M
}


compute_rmut_samples <- compute_rmut_samples_corrected

