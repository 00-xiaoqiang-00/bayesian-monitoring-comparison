# Case interfaces preserve each generator's simulation and cache definitions.
case.data.simulation <- function(
    N.sim, N.max, mu.trt, Sigma.trt, mu.ctrl, Sigma.ctrl, cens_upper,
    cfg = NULL, scenario = c("H1_calibration", "H1_observed", "H0", "H1")) {
  scenario <- match.arg(scenario)
  type <- if (is.null(cfg) || is.null(cfg$generator_type)) "bivariate_lognormal_v1" else as.character(cfg$generator_type)
  generator <- switch(type,
    keynote598_shared_component_lognormal_v1 = .wr_keynote598_simulation,
    oak_lognormal_mixture_v1 = .wr_oak_simulation,
    survival_colon_bivariate_lognormal_v1 = .wr_survival_colon_simulation,
    bivariate_lognormal_v1 = NULL,
    stop("Unsupported generator_type: ", type))
  if (is.null(generator)) {
    return(data.simulation(N.sim=N.sim, N.max=N.max, mu.trt=mu.trt, Sigma.trt=Sigma.trt,
                          mu.ctrl=mu.ctrl, Sigma.ctrl=Sigma.ctrl, cens_upper=cens_upper))
  }
  generator(N.sim=N.sim, N.max=N.max, mu.trt=mu.trt, Sigma.trt=Sigma.trt,
            mu.ctrl=mu.ctrl, Sigma.ctrl=Sigma.ctrl, cens_upper=cens_upper,
            cfg=cfg, scenario=scenario)
}

case.generator.fingerprint <- function(
    cfg = NULL, scenario, mu.trt, Sigma.trt, mu.ctrl, Sigma.ctrl, cens_upper) {
  type <- if (is.null(cfg) || is.null(cfg$generator_type)) "bivariate_lognormal_v1" else as.character(cfg$generator_type)
  fingerprint <- switch(type,
    keynote598_shared_component_lognormal_v1 = .wr_keynote598_fingerprint,
    oak_lognormal_mixture_v1 = .wr_oak_fingerprint,
    survival_colon_bivariate_lognormal_v1 = .wr_survival_colon_fingerprint,
    bivariate_lognormal_v1 = NULL,
    stop("Unsupported generator_type: ", type))
  if (is.null(fingerprint)) {
    return(list(generator_type="bivariate_lognormal_v1", generator_version=1L,
                scenario=scenario, mu.trt=as.numeric(mu.trt), Sigma.trt=as.numeric(Sigma.trt),
                mu.ctrl=as.numeric(mu.ctrl), Sigma.ctrl=as.numeric(Sigma.ctrl),
                cens_upper=cens_upper))
  }
  fingerprint(cfg=cfg, scenario=scenario, mu.trt=mu.trt, Sigma.trt=Sigma.trt,
              mu.ctrl=mu.ctrl, Sigma.ctrl=Sigma.ctrl, cens_upper=cens_upper)
}
