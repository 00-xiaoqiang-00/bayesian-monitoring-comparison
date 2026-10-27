# OAK A-ITT arm-specific two-component lognormal data generator.
# H1-standard uses the specified treatment locations and mixture weight
# with the control component log-scale dispersions.
# H1-observed uses the arm-specific parameters specified below.
#
# Returns an N.max x 7 x N.sim array. OAK event times use the specified
# component parameters; optional mu/Sigma arguments support the common API.

oak_trial_settings <- function() {
  list(
    N.max = 850L,
    interim_cohort = c(340L, 510L, 680L),
    interim_fracs = c(0.40, 0.60, 0.80),
    recruitment_duration = 8.6,
    recruit.int = 8.6 / 850,
    final_analysis_time = 27.9,
    cens_upper = 109.3,
    rho = 0.7,
    generator_version = 5L
  )
}

oak_generator_parameters <- function() {
  list(
    rho = 0.7,
    control = list(
      pi = 0.13123216948698738,
      mu_R1 = 0.2666273699841142,
      sigma_R1 = 0.3000000000000000,
      mu_R2 = 1.6274656470707884,
      sigma_R2 = 0.8489991311979127,
      mu_D = 2.3088355183235585,
      sigma_D = 1.0177
    ),
    H1_calibration = list(
      pi = 0.10183833564436579,
      mu_R1 = 0.510611455725049,
      sigma_R1 = 0.3000000000000000,
      mu_R2 = 1.9615174106704887,
      sigma_R2 = 0.8489991311979127,
      mu_D = 2.636754680552583,
      sigma_D = 1.0177
    ),
    H1_observed = list(
      pi = 0.2623242204350347,
      mu_R1 = 0.3645597370078664,
      sigma_R1 = 0.0806101279283579,
      mu_R2 = 1.677883933719392,
      sigma_R2 = 1.2693769377027608,
      mu_D = 2.636754680552583,
      sigma_D = 1.1709
    ),
    targets = list(
      H1_calibration_pfs_hr = 0.70,
      H1_calibration_pfs_median = 5.7,
      H1_observed_pfs_hr = 0.95,
      os_hr = 0.73,
      H1_standard_variance_rule = "component log-scale sigmas and rho equal control; location and mixture weight retained"
    )
  )
}

oak_validate_mixture <- function(x, label) {
  required <- c("pi", "mu_R1", "sigma_R1", "mu_R2", "sigma_R2",
                "mu_D", "sigma_D")
  missing <- required[!vapply(required, function(nm) {
    !is.null(x[[nm]]) && length(x[[nm]]) == 1L && is.finite(x[[nm]])
  }, logical(1))]
  if (length(missing)) {
    stop(label, " is missing valid parameters: ", paste(missing, collapse = ", "))
  }
  if (x$pi < 0.05 || x$pi > 0.95) {
    stop(label, " pi must lie in [0.05, 0.95].")
  }
  if (x$mu_R1 >= x$mu_R2) {
    stop(label, " must satisfy mu_R1 < mu_R2.")
  }
  if (any(c(x$sigma_R1, x$sigma_R2, x$sigma_D) <= 0)) {
    stop(label, " log-scale standard deviations must be positive.")
  }
  invisible(TRUE)
}

oak_mixture_moments <- function(x, rho = oak_trial_settings()$rho) {
  oak_validate_mixture(x, "mixture")
  mean_R <- x$pi * x$mu_R1 + (1 - x$pi) * x$mu_R2
  var_R <- x$pi * (x$sigma_R1^2 + (x$mu_R1 - mean_R)^2) +
    (1 - x$pi) * (x$sigma_R2^2 + (x$mu_R2 - mean_R)^2)
  cov_RD <- rho * x$sigma_D *
    (x$pi * x$sigma_R1 + (1 - x$pi) * x$sigma_R2)
  list(
    mu = c(mean_R, x$mu_D),
    Sigma = matrix(
      c(var_R, cov_RD, cov_RD, x$sigma_D^2),
      nrow = 2L, byrow = TRUE,
      dimnames = list(c("log_R", "log_D"), c("log_R", "log_D"))
    )
  )
}

oak_scenario_parameters <- function(
    scenario = c("H1_observed", "H1_calibration", "H0", "H1"),
    parameters = oak_generator_parameters()) {
  scenario <- match.arg(scenario)
  if (scenario == "H1") scenario <- "H1_calibration"
  treatment <- if (scenario == "H0") {
    parameters$control
  } else {
    parameters[[scenario]]
  }
  list(
    scenario = scenario,
    treatment = treatment,
    control = parameters$control,
    rho = parameters$rho
  )
}

oak_draw_mixture <- function(n, x, rho) {
  if (n == 0L) return(matrix(numeric(0), nrow = 0L, ncol = 2L))
  oak_validate_mixture(x, "mixture")
  component_1 <- stats::runif(n) < x$pi
  mu_R <- ifelse(component_1, x$mu_R1, x$mu_R2)
  sigma_R <- ifelse(component_1, x$sigma_R1, x$sigma_R2)
  z_R <- stats::rnorm(n)
  z_D_independent <- stats::rnorm(n)
  z_D <- rho * z_R + sqrt(1 - rho^2) * z_D_independent
  cbind(
    recurrence_t = exp(mu_R + sigma_R * z_R),
    death_t = exp(x$mu_D + x$sigma_D * z_D)
  )
}

data.simulation <- function(
    N.sim,
    N.max = oak_trial_settings()$N.max,
    mu.trt = NULL,
    Sigma.trt = NULL,
    mu.ctrl = NULL,
    Sigma.ctrl = NULL,
    cens_upper = oak_trial_settings()$cens_upper,
    scenario = c("H1_observed", "H1_calibration", "H0", "H1"),
    balanced_randomization = TRUE,
    parameters = oak_generator_parameters()) {

  scenario <- match.arg(scenario)
  selected <- oak_scenario_parameters(scenario, parameters)
  N.sim <- as.integer(N.sim)
  N.max <- as.integer(N.max)
  if (N.sim < 1L || N.max < 2L) stop("N.sim >= 1 and N.max >= 2 are required.")
  if (!is.finite(cens_upper) || cens_upper <= 0) {
    stop("cens_upper must be a positive finite number.")
  }
  if (selected$rho <= -1 || selected$rho >= 1) {
    stop("rho must lie strictly between -1 and 1.")
  }
  oak_validate_mixture(selected$treatment, "treatment")
  oak_validate_mixture(selected$control, "control")

  variables <- c("recurrence_t", "death_t", "censor_t", "group", "id",
                 "delta1", "delta2")
  out <- array(
    NA_real_, dim = c(N.max, 7L, N.sim),
    dimnames = list(NULL, Variables = variables,
                    Simulations = paste("Simulation", seq_len(N.sim)))
  )

  make_groups <- function(n) {
    if (!isTRUE(balanced_randomization)) {
      return(stats::rbinom(n, size = 1L, prob = 0.5))
    }
    n_trt <- n %/% 2L
    if ((n %% 2L) == 1L) n_trt <- n_trt + stats::rbinom(1L, 1L, 0.5)
    sample(c(rep.int(1, n_trt), rep.int(0, n - n_trt)), replace = FALSE)
  }

  for (j in seq_len(N.sim)) {
    group <- make_groups(N.max)
    idx_trt <- which(group == 1)
    idx_ctrl <- which(group == 0)
    event_times <- matrix(NA_real_, nrow = N.max, ncol = 2L)
    event_times[idx_trt, ] <- oak_draw_mixture(
      length(idx_trt), selected$treatment, selected$rho
    )
    event_times[idx_ctrl, ] <- oak_draw_mixture(
      length(idx_ctrl), selected$control, selected$rho
    )
    out[, "recurrence_t", j] <- event_times[, 1L]
    out[, "death_t", j] <- event_times[, 2L]
    out[, "censor_t", j] <- stats::runif(N.max, 0, cens_upper)
    out[, "group", j] <- group
    out[, "id", j] <- seq_len(N.max)
  }

  attr(out, "scenario") <- selected$scenario
  attr(out, "generator_version") <- oak_trial_settings()$generator_version
  out
}

oak_entry_times <- function(N.max = oak_trial_settings()$N.max) {
  oak_trial_settings()$recruit.int * (seq_len(as.integer(N.max)) - 1L)
}

oak_admin_followup <- function(
    entry_times = oak_entry_times(),
    analysis_time = oak_trial_settings()$final_analysis_time) {
  pmax(analysis_time - entry_times, 0)
}
