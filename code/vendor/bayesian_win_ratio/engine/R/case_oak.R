#' Dispatch data generation by case configuration
#'
#' Configurations without generator_type use data.simulation().
#' OAK_A_ITT uses data.simulation.oak.mixture().
#'
.wr_oak_simulation <- function(
    N.sim, N.max, mu.trt, Sigma.trt, mu.ctrl, Sigma.ctrl, cens_upper,
    cfg = NULL,
    scenario = c("H1_calibration", "H1_observed", "H0", "H1")) {

  scenario <- match.arg(scenario)
  generator_type <- if (is.null(cfg) || is.null(cfg$generator_type)) {
    "bivariate_lognormal_v1"
  } else {
    as.character(cfg$generator_type)
  }

  if (identical(generator_type, "oak_lognormal_mixture_v1")) {
    if (is.null(cfg$oak_generator_parameters)) {
      stop("OAK config is missing oak_generator_parameters.")
    }
    oak_scenario <- if (scenario == "H1") "H1_observed" else scenario
    return(data.simulation.oak.mixture(
      N.sim = N.sim,
      N.max = N.max,
      parameters = cfg$oak_generator_parameters,
      scenario = oak_scenario,
      cens_upper = cens_upper,
      balanced_randomization = if (is.null(cfg$balanced_randomization)) TRUE else
        isTRUE(cfg$balanced_randomization)
    ))
  }

  if (!identical(generator_type, "bivariate_lognormal_v1")) {
    stop("Unsupported generator_type: ", generator_type)
  }
  data.simulation(
    N.sim = N.sim, N.max = N.max,
    mu.trt = mu.trt, Sigma.trt = Sigma.trt,
    mu.ctrl = mu.ctrl, Sigma.ctrl = Sigma.ctrl,
    cens_upper = cens_upper
  )
}

#' Fingerprint the data-generating settings used by a cache
#'
.wr_oak_fingerprint <- function(
    cfg = NULL, scenario, mu.trt, Sigma.trt, mu.ctrl, Sigma.ctrl,
    cens_upper) {

  generator_type <- if (is.null(cfg) || is.null(cfg$generator_type)) {
    "bivariate_lognormal_v1"
  } else {
    as.character(cfg$generator_type)
  }
  list(
    generator_type = generator_type,
    generator_version = if (is.null(cfg) || is.null(cfg$generator_version)) 1L else
      cfg$generator_version,
    scenario = as.character(scenario),
    mu.trt = as.numeric(mu.trt),
    Sigma.trt = as.numeric(Sigma.trt),
    mu.ctrl = as.numeric(mu.ctrl),
    Sigma.ctrl = as.numeric(Sigma.ctrl),
    cens_upper = as.numeric(cens_upper),
    balanced_randomization = if (is.null(cfg) || is.null(cfg$balanced_randomization))
      FALSE else isTRUE(cfg$balanced_randomization),
    oak_generator_parameters = if (identical(generator_type, "oak_lognormal_mixture_v1"))
      cfg$oak_generator_parameters else NULL
  )
}
