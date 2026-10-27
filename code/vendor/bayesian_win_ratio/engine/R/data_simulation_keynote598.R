#' Select shared-component KEYNOTE-598 parameters for a simulation scenario
#'
#' @export
keynote598.scenario.parameters <- function(
    cfg, scenario = c("H1_calibration", "H1_observed", "H0", "H1")) {
  scenario <- match.arg(scenario)
  required <- c(
    "keynote598_control",
    "keynote598_treatment_observed",
    "keynote598_treatment_standard"
  )
  missing <- required[vapply(required, function(x) is.null(cfg[[x]]), logical(1))]
  if (length(missing)) {
    stop("KEYNOTE-598 config is missing: ", paste(missing, collapse = ", "))
  }

  if (identical(scenario, "H0")) {
    treatment <- cfg$keynote598_control
  } else if (identical(scenario, "H1_observed")) {
    treatment <- cfg$keynote598_treatment_observed
  } else {
    treatment <- cfg$keynote598_treatment_standard
  }
  list(treatment = treatment, control = cfg$keynote598_control)
}


#' Simulate KEYNOTE-598 latent progression and death times
#'
#' Each patient draws one shared mixture label K. Conditional on K, log
#' progression time and log death time follow a bivariate Normal distribution
#' with component-specific locations and standard deviations and a fixed
#' within-component correlation. PFS is derived downstream as min(R, D).
#'
#' @export
data.simulation.keynote598 <- function(
    N.sim, N.max, cens_upper, cfg,
    scenario = c("H1_calibration", "H1_observed", "H0", "H1"),
    balanced_randomization = TRUE) {

  scenario <- match.arg(scenario)
  N.sim <- as.integer(N.sim)
  N.max <- as.integer(N.max)
  if (N.sim < 1L || N.max < 2L) {
    stop("N.sim must be >= 1 and N.max must be >= 2.")
  }
  if (length(cens_upper) != 1L || !is.finite(cens_upper) ||
      cens_upper <= 0) {
    stop("cens_upper must be one positive finite value.")
  }

  parameters <- keynote598.scenario.parameters(cfg, scenario)

  validate_arm <- function(arm, label) {
    required <- c("pi", "mu_R", "sigma_R", "mu_D", "sigma_D", "rho")
    missing <- required[!required %in% names(arm)]
    if (length(missing)) {
      stop(label, " parameters are missing: ", paste(missing, collapse = ", "))
    }
    if (length(arm$pi) != 2L || any(!is.finite(arm$pi)) ||
        any(arm$pi <= 0) || abs(sum(arm$pi) - 1) > 1e-10) {
      stop(label, " pi must contain two positive weights summing to one.")
    }
    for (nm in c("mu_R", "sigma_R", "mu_D", "sigma_D")) {
      if (length(arm[[nm]]) != 2L || any(!is.finite(arm[[nm]]))) {
        stop(label, " ", nm, " must contain two finite values.")
      }
    }
    if (any(arm$sigma_R <= 0) || any(arm$sigma_D <= 0)) {
      stop(label, " component standard deviations must be positive.")
    }
    if (length(arm$rho) != 1L || !is.finite(arm$rho) ||
        abs(arm$rho) >= 1) {
      stop(label, " rho must be one finite value strictly between -1 and 1.")
    }
  }
  validate_arm(parameters$treatment, "Treatment")
  validate_arm(parameters$control, "Control")

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

  simulate_arm <- function(n, arm) {
    if (n == 0L) {
      return(list(log_times = matrix(numeric(0), ncol = 2L),
                  component = integer(0)))
    }
    component <- 1L + as.integer(stats::runif(n) > arm$pi[1L])
    z_death <- stats::rnorm(n)
    z_independent <- stats::rnorm(n)
    z_progression <- arm$rho * z_death +
      sqrt(1 - arm$rho^2) * z_independent
    log_progression <- arm$mu_R[component] +
      arm$sigma_R[component] * z_progression
    log_death <- arm$mu_D[component] + arm$sigma_D[component] * z_death
    list(
      log_times = cbind(log_progression, log_death),
      component = component
    )
  }

  myData <- array(NA_real_, dim = c(N.max, 7L, N.sim))
  component_labels <- matrix(NA_integer_, nrow = N.max, ncol = N.sim)
  for (j in seq_len(N.sim)) {
    group <- make_groups(N.max)
    idx_trt <- which(group == 1L)
    idx_ctrl <- which(group == 0L)
    log_times <- matrix(NA_real_, nrow = N.max, ncol = 2L)

    trt <- simulate_arm(length(idx_trt), parameters$treatment)
    ctrl <- simulate_arm(length(idx_ctrl), parameters$control)
    log_times[idx_trt, ] <- trt$log_times
    log_times[idx_ctrl, ] <- ctrl$log_times
    component_labels[idx_trt, j] <- trt$component
    component_labels[idx_ctrl, j] <- ctrl$component

    myData[, 1L, j] <- exp(log_times[, 1L])
    myData[, 2L, j] <- exp(log_times[, 2L])
    myData[, 3L, j] <- stats::runif(N.max, min = 0, max = cens_upper)
    myData[, 4L, j] <- group
    myData[, 5L, j] <- seq_len(N.max)
  }

  dimnames(myData) <- list(
    NULL,
    Variables = c("progression_t", "death_t", "censor_t", "group", "id",
                  "delta1", "delta2"),
    Simulations = paste("Simulation", seq_len(N.sim))
  )
  attr(myData, "component_labels") <- component_labels
  attr(myData, "generator_scenario") <- scenario
  myData
}
