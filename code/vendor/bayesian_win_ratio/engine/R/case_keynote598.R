#' Dispatch data generation by case configuration
#'
#' The survival_colon case uses a bivariate-lognormal generator.  The
#' keynote598 case uses a shared-component two-class bivariate-lognormal
#' generator. Configurations without generator_type use data.simulation().
#'
.wr_keynote598_simulation <- function(
    N.sim, N.max, mu.trt, Sigma.trt, mu.ctrl, Sigma.ctrl, cens_upper,
    cfg = NULL,
    scenario = c("H1_calibration", "H1_observed", "H0", "H1")) {

  scenario <- match.arg(scenario)
  generator_type <- if (is.null(cfg) || is.null(cfg$generator_type)) {
    "bivariate_lognormal_v1"
  } else {
    as.character(cfg$generator_type)
  }

  if (identical(generator_type,
                "survival_colon_bivariate_lognormal_v1")) {
    return(data.simulation.survival.colon(
      N.sim = N.sim,
      N.max = N.max,
      mu.trt = mu.trt,
      Sigma.trt = Sigma.trt,
      mu.ctrl = mu.ctrl,
      Sigma.ctrl = Sigma.ctrl,
      cens_upper = cens_upper,
      balanced_randomization = if (is.null(cfg$balanced_randomization)) {
        FALSE
      } else {
        isTRUE(cfg$balanced_randomization)
      }
    ))
  }

  if (identical(generator_type,
                "keynote598_shared_component_lognormal_v1")) {
    return(data.simulation.keynote598(
      N.sim = N.sim,
      N.max = N.max,
      cens_upper = cens_upper,
      cfg = cfg,
      scenario = scenario,
      balanced_randomization = if (is.null(cfg$balanced_randomization)) {
        TRUE
      } else {
        isTRUE(cfg$balanced_randomization)
      }
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
.wr_keynote598_fingerprint <- function(
    cfg = NULL, scenario, mu.trt, Sigma.trt, mu.ctrl, Sigma.ctrl,
    cens_upper) {

  generator_type <- if (is.null(cfg) || is.null(cfg$generator_type)) {
    "bivariate_lognormal_v1"
  } else {
    as.character(cfg$generator_type)
  }

  fingerprint <- list(
    generator_type = generator_type,
    generator_version = if (is.null(cfg) || is.null(cfg$generator_version)) {
      1L
    } else {
      cfg$generator_version
    },
    scenario = as.character(scenario),
    mu.trt = as.numeric(mu.trt),
    Sigma.trt = as.numeric(Sigma.trt),
    mu.ctrl = as.numeric(mu.ctrl),
    Sigma.ctrl = as.numeric(Sigma.ctrl),
    cens_upper = as.numeric(cens_upper),
    balanced_randomization = if (is.null(cfg) ||
                                  is.null(cfg$balanced_randomization)) {
      FALSE
    } else {
      isTRUE(cfg$balanced_randomization)
    }
  )

  if (identical(generator_type,
                "keynote598_shared_component_lognormal_v1")) {
    fingerprint$keynote598_parameters <-
      keynote598.scenario.parameters(cfg, scenario)
  }
  fingerprint
}
