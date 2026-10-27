# KEYNOTE-598 shared-component progression/death data generator.
#
# R is latent progression time, D is latent death time, and PFS=min(R,D).
# Each patient draws one shared mixture label K that selects both progression
# and death parameters. The two numerical components are not clinical
# subgroups.

keynote598_trial_settings <- function(tau_fixed = c(4, 6, 9, 18)) {
  list(
    N.max = 568L,
    interim_cohort = c(228L, 340L, 454L),
    interim_fracs = c(228, 340, 454) / 568,
    recruitment_duration = 19.3,
    recruit.int = 19.3 / 568,
    final_analysis_time = 31.7,
    cens_upper = 141.8559446075,
    tau_fixed = tau_fixed,
    time_max = 120,
    balanced_randomization = TRUE,
    generator_version = 1L
  )
}

keynote598_component_sigma <- function(sigma_R, sigma_D, rho = 0.5) {
  matrix(
    c(
      sigma_R^2, rho * sigma_R * sigma_D,
      rho * sigma_R * sigma_D, sigma_D^2
    ),
    nrow = 2L, byrow = TRUE,
    dimnames = list(c("log_R", "log_D"), c("log_R", "log_D"))
  )
}

keynote598_generator_parameters <- function() {
  control <- list(
    pi = c(0.2070545742, 0.7929454258),
    mu_R = c(log(2.4837219024), log(16.4600817241)),
    sigma_R = c(0.4340057961, 1.5317036442),
    mu_D = c(log(2.7290953670), log(28.6890202294)),
    sigma_D = c(1.5000000000, 0.9670082477),
    rho = 0.5
  )
  control$Sigma <- list(
    keynote598_component_sigma(
      control$sigma_R[1L], control$sigma_D[1L], control$rho
    ),
    keynote598_component_sigma(
      control$sigma_R[2L], control$sigma_D[2L], control$rho
    )
  )

  observed <- list(
    pi = c(0.1936579096, 0.8063420904),
    mu_R = c(log(2.0212342835), log(15.2003022720)),
    sigma_R = c(0.1500000000, 1.3810016053),
    mu_D = c(log(2.5120209874), log(29.2329569146)),
    sigma_D = c(1.5000000000, 1.2402934461),
    rho = 0.5
  )
  observed$Sigma <- list(
    keynote598_component_sigma(
      observed$sigma_R[1L], observed$sigma_D[1L], observed$rho
    ),
    keynote598_component_sigma(
      observed$sigma_R[2L], observed$sigma_D[2L], observed$rho
    )
  )

  standard <- control
  standard$mu_R <- control$mu_R + 0.5946115916
  standard$mu_D <- control$mu_D + 0.3773863423

  list(
    control = control,
    H1_calibration = standard,
    H1_observed = observed,
    targets = list(
      H1_calibration_pfs_hr = 0.69,
      H1_calibration_os_hr = 0.70,
      calibration_horizon = 31.7
    )
  )
}

keynote598_validate_arm <- function(x, label) {
  required <- c("pi", "mu_R", "sigma_R", "mu_D", "sigma_D", "rho", "Sigma")
  missing <- required[!required %in% names(x)]
  if (length(missing)) {
    stop(label, " parameters are missing: ", paste(missing, collapse = ", "))
  }
  if (length(x$pi) != 2L || any(!is.finite(x$pi)) ||
      any(x$pi <= 0) || abs(sum(x$pi) - 1) > 1e-10) {
    stop(label, " pi must contain two positive weights summing to one.")
  }
  for (nm in c("mu_R", "sigma_R", "mu_D", "sigma_D")) {
    if (length(x[[nm]]) != 2L || any(!is.finite(x[[nm]]))) {
      stop(label, " ", nm, " must contain two finite values.")
    }
  }
  if (any(x$sigma_R <= 0) || any(x$sigma_D <= 0) ||
      length(x$rho) != 1L || !is.finite(x$rho) || abs(x$rho) >= 1) {
    stop(label, " has invalid standard deviations or rho.")
  }
  if (!is.list(x$Sigma) || length(x$Sigma) != 2L) {
    stop(label, " Sigma must contain two covariance matrices.")
  }
  for (k in seq_len(2L)) {
    expected <- keynote598_component_sigma(
      x$sigma_R[k], x$sigma_D[k], x$rho
    )
    if (!isTRUE(all.equal(x$Sigma[[k]], expected, tolerance = 1e-12)) ||
        min(eigen(x$Sigma[[k]], symmetric = TRUE,
                  only.values = TRUE)$values) <= 0) {
      stop(label, " component ", k, " covariance is invalid.")
    }
  }
  invisible(TRUE)
}

keynote598_scenario_parameters <- function(
    scenario = c("H1_observed", "H1_calibration", "H0", "H1"),
    parameters = keynote598_generator_parameters()) {
  scenario <- match.arg(scenario)
  if (scenario == "H1") scenario <- "H1_calibration"
  treatment <- if (scenario == "H0") parameters$control else parameters[[scenario]]
  keynote598_validate_arm(parameters$control, "control")
  keynote598_validate_arm(treatment, "treatment")
  list(
    scenario = scenario,
    treatment = treatment,
    control = parameters$control
  )
}

keynote598_draw_times <- function(n, arm) {
  if (n == 0L) {
    return(list(
      times = matrix(numeric(0), nrow = 0L, ncol = 2L),
      component = integer(0)
    ))
  }
  component <- 1L + as.integer(stats::runif(n) > arm$pi[1L])
  z_D <- stats::rnorm(n)
  z_independent <- stats::rnorm(n)
  z_R <- arm$rho * z_D + sqrt(1 - arm$rho^2) * z_independent
  log_R <- arm$mu_R[component] + arm$sigma_R[component] * z_R
  log_D <- arm$mu_D[component] + arm$sigma_D[component] * z_D
  list(
    times = cbind(progression_t = exp(log_R), death_t = exp(log_D)),
    component = component
  )
}

data.simulation <- function(
    N.sim,
    N.max = keynote598_trial_settings()$N.max,
    mu.trt = NULL,
    Sigma.trt = NULL,
    mu.ctrl = NULL,
    Sigma.ctrl = NULL,
    cens_upper = keynote598_trial_settings()$cens_upper,
    scenario = c("H1_observed", "H1_calibration", "H0", "H1"),
    balanced_randomization =
      keynote598_trial_settings()$balanced_randomization,
    parameters = keynote598_generator_parameters()) {

  scenario <- match.arg(scenario)
  selected <- keynote598_scenario_parameters(scenario, parameters)
  N.sim <- as.integer(N.sim)
  N.max <- as.integer(N.max)
  if (N.sim < 1L || N.max < 2L) {
    stop("N.sim >= 1 and N.max >= 2 are required.")
  }
  if (!is.finite(cens_upper) || cens_upper <= 0) {
    stop("cens_upper must be a positive finite number.")
  }

  variables <- c("progression_t", "death_t", "censor_t", "group", "id",
                 "delta1", "delta2")
  out <- array(
    NA_real_, dim = c(N.max, 7L, N.sim),
    dimnames = list(NULL, Variables = variables,
                    Simulations = paste("Simulation", seq_len(N.sim)))
  )
  component_labels <- matrix(NA_integer_, nrow = N.max, ncol = N.sim)

  make_groups <- function(n) {
    if (!isTRUE(balanced_randomization)) {
      return(sample(c(0L, 1L), n, replace = TRUE))
    }
    n_trt <- n %/% 2L
    if ((n %% 2L) == 1L) {
      n_trt <- n_trt + stats::rbinom(1L, 1L, 0.5)
    }
    sample(c(rep.int(1L, n_trt), rep.int(0L, n - n_trt)),
           replace = FALSE)
  }

  for (j in seq_len(N.sim)) {
    group <- make_groups(N.max)
    idx_trt <- which(group == 1L)
    idx_ctrl <- which(group == 0L)
    event_times <- matrix(NA_real_, nrow = N.max, ncol = 2L)

    trt <- keynote598_draw_times(length(idx_trt), selected$treatment)
    ctrl <- keynote598_draw_times(length(idx_ctrl), selected$control)
    event_times[idx_trt, ] <- trt$times
    event_times[idx_ctrl, ] <- ctrl$times
    component_labels[idx_trt, j] <- trt$component
    component_labels[idx_ctrl, j] <- ctrl$component

    out[, "progression_t", j] <- event_times[, 1L]
    out[, "death_t", j] <- event_times[, 2L]
    out[, "censor_t", j] <- stats::runif(N.max, 0, cens_upper)
    out[, "group", j] <- group
    out[, "id", j] <- seq_len(N.max)
  }

  attr(out, "scenario") <- selected$scenario
  attr(out, "generator_version") <-
    keynote598_trial_settings()$generator_version
  attr(out, "component_labels") <- component_labels
  out
}

keynote598_entry_times <- function(
    N.max = keynote598_trial_settings()$N.max) {
  keynote598_trial_settings()$recruit.int *
    (seq_len(as.integer(N.max)) - 1L)
}

keynote598_admin_followup <- function(
    entry_times = keynote598_entry_times(),
    analysis_time = keynote598_trial_settings()$final_analysis_time) {
  pmax(analysis_time - entry_times, 0)
}
