# survival::colon Obs versus Lev+5FU data generator.
#
# The latent recurrence and death times are generated jointly on the log scale.
# Returns an N.max x 7 x N.sim array for Cox, RMST and RMUT simulations.

survival_colon_trial_settings <- function() {
  list(
    N.max = 620L,
    interim_cohort = c(248L, 372L, 496L),
    interim_fracs = c(0.40, 0.60, 0.80),
    recruitment_duration = 43,
    recruit.int = 43 / 620,
    final_analysis_time = 109,
    cens_upper = 1900,
    tau_fixed = c(17.1, 25.7, 34.3, 60),
    generator_version = 1L
  )
}

survival_colon_generator_parameters <- function() {
  Sigma_ctrl <- matrix(
    c(3.29767175, 2.19701400,
      2.19701400, 1.77436755),
    nrow = 2L, byrow = TRUE,
    dimnames = list(c("log_R", "log_D"), c("log_R", "log_D"))
  )

  list(
    control = list(
      mu = c(log(38.1268), log(64.5907)),
      Sigma = Sigma_ctrl
    ),
    H1_calibration = list(
      mu = c(log(92.6930), log(116.2198)),
      Sigma = Sigma_ctrl
    ),
    H1_observed = list(
      mu = c(log(87.1562), log(111.4653)),
      Sigma = matrix(
        c(4.46376802, 3.41770993,
          3.41770993, 2.92813972),
        nrow = 2L, byrow = TRUE,
        dimnames = list(c("log_R", "log_D"), c("log_R", "log_D"))
      )
    ),
    targets = list(
      H1_calibration_recurrence_hr = 0.59,
      H1_calibration_death_hr = 0.67,
      H1_observed_recurrence_hr = 0.5990,
      H1_observed_death_hr = 0.6240
    )
  )
}

survival_colon_validate_arm <- function(x, label) {
  if (!is.list(x) || length(x$mu) != 2L || any(!is.finite(x$mu))) {
    stop(label, " mu must contain two finite log-time means.")
  }
  if (!is.matrix(x$Sigma) || !identical(dim(x$Sigma), c(2L, 2L)) ||
      any(!is.finite(x$Sigma)) ||
      max(abs(x$Sigma - t(x$Sigma))) > 1e-10 ||
      min(eigen(x$Sigma, symmetric = TRUE, only.values = TRUE)$values) <= 0) {
    stop(label, " Sigma must be a finite symmetric positive-definite 2x2 matrix.")
  }
  invisible(TRUE)
}

survival_colon_scenario_parameters <- function(
    scenario = c("H1_observed", "H1_calibration", "H0", "H1"),
    parameters = survival_colon_generator_parameters()) {
  scenario <- match.arg(scenario)
  if (scenario == "H1") scenario <- "H1_calibration"
  treatment <- if (scenario == "H0") parameters$control else parameters[[scenario]]
  survival_colon_validate_arm(parameters$control, "control")
  survival_colon_validate_arm(treatment, "treatment")
  list(
    scenario = scenario,
    treatment = treatment,
    control = parameters$control
  )
}

survival_colon_draw_times <- function(n, arm) {
  if (n == 0L) return(matrix(numeric(0), nrow = 0L, ncol = 2L))
  exp(MASS::mvrnorm(n = n, mu = arm$mu, Sigma = arm$Sigma))
}

data.simulation <- function(
    N.sim,
    N.max = survival_colon_trial_settings()$N.max,
    mu.trt = NULL,
    Sigma.trt = NULL,
    mu.ctrl = NULL,
    Sigma.ctrl = NULL,
    cens_upper = survival_colon_trial_settings()$cens_upper,
    scenario = c("H1_observed", "H1_calibration", "H0", "H1"),
    balanced_randomization = FALSE,
    parameters = survival_colon_generator_parameters()) {

  scenario <- match.arg(scenario)
  selected <- survival_colon_scenario_parameters(scenario, parameters)
  N.sim <- as.integer(N.sim)
  N.max <- as.integer(N.max)
  if (N.sim < 1L || N.max < 2L) stop("N.sim >= 1 and N.max >= 2 are required.")
  if (!is.finite(cens_upper) || cens_upper <= 0) {
    stop("cens_upper must be a positive finite number.")
  }

  variables <- c("recurrence_t", "death_t", "censor_t", "group", "id",
                 "delta1", "delta2")
  out <- array(
    NA_real_, dim = c(N.max, 7L, N.sim),
    dimnames = list(NULL, Variables = variables,
                    Simulations = paste("Simulation", seq_len(N.sim)))
  )

  make_groups <- function(n) {
    if (!isTRUE(balanced_randomization)) {
      return(sample(c(0, 1), n, replace = TRUE))
    }
    n_trt <- n %/% 2L
    if ((n %% 2L) == 1L) n_trt <- n_trt + stats::rbinom(1L, 1L, 0.5)
    sample(c(rep.int(1, n_trt), rep.int(0, n - n_trt)), replace = FALSE)
  }

  for (j in seq_len(N.sim)) {
    group <- make_groups(N.max)
    idx_trt <- which(group == 1L)
    idx_ctrl <- which(group == 0L)
    event_times <- matrix(NA_real_, nrow = N.max, ncol = 2L)
    event_times[idx_trt, ] <- survival_colon_draw_times(
      length(idx_trt), selected$treatment
    )
    event_times[idx_ctrl, ] <- survival_colon_draw_times(
      length(idx_ctrl), selected$control
    )

    out[, "recurrence_t", j] <- event_times[, 1L]
    out[, "death_t", j] <- event_times[, 2L]
    out[, "censor_t", j] <- stats::runif(N.max, 0, cens_upper)
    out[, "group", j] <- group
    out[, "id", j] <- seq_len(N.max)
  }

  attr(out, "scenario") <- selected$scenario
  attr(out, "generator_version") <- survival_colon_trial_settings()$generator_version
  out
}

survival_colon_entry_times <- function(
    N.max = survival_colon_trial_settings()$N.max) {
  survival_colon_trial_settings()$recruit.int *
    (seq_len(as.integer(N.max)) - 1L)
}

survival_colon_admin_followup <- function(
    entry_times = survival_colon_entry_times(),
    analysis_time = survival_colon_trial_settings()$final_analysis_time) {
  pmax(analysis_time - entry_times, 0)
}
