#' Simulate OAK data using a two-component lognormal PFS model
#'
#' The treatment recurrence time follows a two-component lognormal mixture
#' under H1_observed or H1_calibration. Death has the same lognormal marginal
#' distribution in both treatment mixture components. The control arm also
#' follows a two-component recurrence mixture with a common locked death
#' marginal. Under H0 both randomized groups use that control-arm mixture.
#'
#' The returned array has the same layout as data.simulation(). censor_t is
#' the raw independent random-censoring duration. Administrative censoring at
#' the calendar analysis time is applied later by the analysis entries.
#'
#' @export
data.simulation.oak.mixture <- function(
    N.sim,
    N.max,
    parameters,
    scenario = c("H1_calibration", "H1_observed", "H0", "H1"),
    cens_upper = parameters$cens_upper,
    balanced_randomization = TRUE) {

  scenario <- match.arg(scenario)
  if (scenario == "H1") scenario <- "H1_observed"
  N.sim <- as.integer(N.sim)
  N.max <- as.integer(N.max)

  if (N.sim < 1L || N.max < 2L) {
    stop("N.sim must be >= 1 and N.max must be >= 2.")
  }
  if (!is.finite(cens_upper) || cens_upper <= 0) {
    stop("cens_upper must be a positive finite number.")
  }

  required <- c(
    "rho",
    "pi", "mu_R1", "sigma_R1", "mu_R2", "sigma_R2",
    "mu_D_trt", "sigma_D_trt",
    "pi_ctrl", "mu_R1_ctrl", "sigma_R1_ctrl", "mu_R2_ctrl",
    "sigma_R2_ctrl", "mu_D_ctrl", "sigma_D_ctrl"
  )
  missing_parameters <- required[!vapply(required, function(x) {
    !is.null(parameters[[x]]) && length(parameters[[x]]) == 1L &&
      is.finite(parameters[[x]])
  }, logical(1))]
  if (length(missing_parameters)) {
    stop("Missing/invalid OAK parameters: ", paste(missing_parameters, collapse = ", "))
  }
  if (parameters$rho <= -1 || parameters$rho >= 1) {
    stop("rho must lie strictly between -1 and 1.")
  }
  validate_mixture <- function(x, label) {
    if (x$pi < 0.05 || x$pi > 0.95) {
      stop(label, " pi must lie in [0.05, 0.95].")
    }
    if (x$mu_R1 >= x$mu_R2) {
      stop(label, " must satisfy mu_R1 < mu_R2.")
    }
    if (any(c(x$sigma_R1, x$sigma_R2, x$sigma_D_trt) <= 0)) {
      stop("All ", label, " log-scale standard deviations must be positive.")
    }
  }

  observed_parameters <- parameters[
    c("pi", "mu_R1", "sigma_R1", "mu_R2", "sigma_R2",
      "mu_D_trt", "sigma_D_trt")
  ]
  control_parameters <- with(parameters, list(
    pi = pi_ctrl,
    mu_R1 = mu_R1_ctrl,
    sigma_R1 = sigma_R1_ctrl,
    mu_R2 = mu_R2_ctrl,
    sigma_R2 = sigma_R2_ctrl,
    mu_D_trt = mu_D_ctrl,
    sigma_D_trt = sigma_D_ctrl
  ))
  validate_mixture(observed_parameters, "H1_observed")
  validate_mixture(control_parameters, "control")
  if (scenario == "H1_calibration") {
    calibration_required <- c("pi", "mu_R1", "sigma_R1", "mu_R2",
                              "sigma_R2", "mu_D_trt", "sigma_D_trt")
    calibration_missing <- calibration_required[!vapply(
      calibration_required,
      function(x) !is.null(parameters$calibration[[x]]) &&
        length(parameters$calibration[[x]]) == 1L &&
        is.finite(parameters$calibration[[x]]),
      logical(1)
    )]
    if (length(calibration_missing)) {
      stop("Missing/invalid calibration parameters: ",
           paste(calibration_missing, collapse = ", "))
    }
    validate_mixture(parameters$calibration, "H1_calibration")
  }

  draw_mixture <- function(n, mixture_parameters) {
    if (n == 0L) return(matrix(numeric(0), nrow = 0L, ncol = 2L))
    component_1 <- stats::runif(n) < mixture_parameters$pi
    sigma_R <- ifelse(component_1, mixture_parameters$sigma_R1,
                      mixture_parameters$sigma_R2)
    mu_R <- ifelse(component_1, mixture_parameters$mu_R1,
                   mixture_parameters$mu_R2)
    z_R <- stats::rnorm(n)
    z_D_independent <- stats::rnorm(n)
    z_D <- parameters$rho * z_R +
      sqrt(1 - parameters$rho^2) * z_D_independent
    cbind(
      recurrence_t = exp(mu_R + sigma_R * z_R),
      death_t = exp(mixture_parameters$mu_D_trt +
                      mixture_parameters$sigma_D_trt * z_D)
    )
  }

  make_groups <- function(n) {
    if (isTRUE(balanced_randomization)) {
      n_trt <- floor(n / 2)
      if ((n %% 2L) == 1L) n_trt <- n_trt + stats::rbinom(1L, 1L, 0.5)
      return(sample(c(rep.int(1, n_trt), rep.int(0, n - n_trt)), n, replace = FALSE))
    }
    stats::rbinom(n, size = 1L, prob = 0.5)
  }

  myData <- array(NA_real_, dim = c(N.max, 7L, N.sim))
  for (j in seq_len(N.sim)) {
    group <- make_groups(N.max)
    idx_trt <- which(group == 1)
    idx_ctrl <- which(group == 0)
    event_times <- matrix(NA_real_, nrow = N.max, ncol = 2L)

    if (scenario == "H1_calibration") {
      event_times[idx_trt, ] <- draw_mixture(
        length(idx_trt), parameters$calibration
      )
    } else if (scenario == "H1_observed") {
      event_times[idx_trt, ] <- draw_mixture(
        length(idx_trt), observed_parameters
      )
    } else {
      event_times[idx_trt, ] <- draw_mixture(
        length(idx_trt), control_parameters
      )
    }
    event_times[idx_ctrl, ] <- draw_mixture(
      length(idx_ctrl), control_parameters
    )

    myData[, 1L, j] <- event_times[, 1L]
    myData[, 2L, j] <- event_times[, 2L]
    myData[, 3L, j] <- stats::runif(N.max, min = 0, max = cens_upper)
    myData[, 4L, j] <- group
    myData[, 5L, j] <- seq_len(N.max)
    # delta1 and delta2 deliberately remain NA; winratio() constructs them.
  }

  dimnames(myData) <- list(
    NULL,
    Variables = c("recurrence_t", "death_t", "censor_t", "group", "id",
                  "delta1", "delta2"),
    Simulations = paste("Simulation", seq_len(N.sim))
  )
  myData
}
