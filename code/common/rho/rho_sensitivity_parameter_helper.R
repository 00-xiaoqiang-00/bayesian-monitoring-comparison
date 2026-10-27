# KM-target rho-sensitivity parameter reader and deterministic seed policy.
# This helper is shared by the packaged KM-target sensitivity entries.

rho_sensitivity_version_dir <- Sys.getenv("RHO_SENSITIVITY_VERSION_DIR", unset = "")
if (!nzchar(rho_sensitivity_version_dir)) {
  stop("RHO_SENSITIVITY_VERSION_DIR must point to the KM-target parameter directory.")
}
rho_sensitivity_version_dir <- normalizePath(rho_sensitivity_version_dir, mustWork = TRUE)
rho_parameter_file <- file.path(
  rho_sensitivity_version_dir, "rho_sensitivity_km_target_parameters.csv"
)
if (!file.exists(rho_parameter_file)) {
  stop("KM-target parameter table not found: ", rho_parameter_file)
}
rho_parameter_table <- utils::read.csv(
  rho_parameter_file, stringsAsFactors = FALSE, check.names = FALSE
)
rho_required_columns <- c(
  "case", "scenario", "arm", "rho_label", "rho", "component", "weight",
  "mu_log_R", "sigma_log_R", "mu_log_D", "sigma_log_D", "cov_log_R_D"
)
rho_missing_columns <- setdiff(rho_required_columns, names(rho_parameter_table))
if (length(rho_missing_columns)) {
  stop("KM-target parameter table is missing: ", paste(rho_missing_columns, collapse = ", "))
}

rho_case_name <- function(x) {
  x <- as.character(x)
  if (identical(x, "KEYNOTE598")) return("KEYNOTE-598")
  if (identical(x, "IMpower")) return("IMpower131")
  x
}

rho_internal_scenario <- function(x) {
  x <- as.character(x)
  if (x %in% c("H1_standard", "H1-standard", "H1_sens2", "H1_calibration")) return("H1-standard")
  if (x %in% c("H1_observed", "H1-observed", "Realcase")) return("H1-observed")
  if (identical(x, "H0")) return("H0")
  stop("Unknown internal scenario: ", x)
}

rho_get_arm <- function(case, scenario, arm, rho_level) {
  case <- rho_case_name(case)
  scenario <- rho_internal_scenario(scenario)
  if (!rho_level %in% c("low", "high")) stop("rho_level must be low or high.")
  z <- rho_parameter_table[
    rho_parameter_table$case == case &
      rho_parameter_table$scenario == scenario &
      rho_parameter_table$arm == arm &
      rho_parameter_table$rho_label == rho_level,
    , drop = FALSE
  ]
  if (!nrow(z)) stop("No KM-target parameters for ", case, "/", scenario, "/", arm, "/", rho_level)
  z <- z[order(as.integer(z$component)), , drop = FALSE]
  if (anyDuplicated(z$component) || any(!is.finite(as.numeric(z$component)))) {
    stop("Invalid or duplicated component labels for ", case, "/", scenario, "/", arm)
  }
  if (abs(sum(as.numeric(z$weight)) - 1) > 1e-8) {
    stop("Mixture weights do not sum to one for ", case, "/", scenario, "/", arm)
  }
  if (any(as.numeric(z$sigma_log_R) <= 0) || any(as.numeric(z$sigma_log_D) <= 0)) {
    stop("Non-positive log-time standard deviation in parameter table.")
  }
  rho_value <- unique(as.numeric(z$rho))
  if (length(rho_value) != 1L || !is.finite(rho_value) || abs(rho_value) >= 1) {
    stop("Invalid rho in parameter table.")
  }
  expected_cov <- rho_value * as.numeric(z$sigma_log_R) * as.numeric(z$sigma_log_D)
  if (max(abs(as.numeric(z$cov_log_R_D) - expected_cov)) > 1e-8) {
    stop("Covariance column is inconsistent with rho*sigma_R*sigma_D.")
  }
  list(
    weight = as.numeric(z$weight),
    mu_log_R = as.numeric(z$mu_log_R),
    sigma_log_R = as.numeric(z$sigma_log_R),
    mu_log_D = as.numeric(z$mu_log_D),
    sigma_log_D = as.numeric(z$sigma_log_D),
    rho = rho_value,
    covariance = as.numeric(z$cov_log_R_D),
    component = as.integer(z$component)
  )
}

rho_make_generator_parameters <- function(case, rho_level) {
  case <- rho_case_name(case)
  rho_low <- rho_get_arm(case, "H0", "control", rho_level)$rho
  if (case == "IMpower131") return(invisible(NULL))
  if (case == "survival_colon") {
    make_arm <- function(scenario, arm = "control") {
      a <- rho_get_arm(case, scenario, arm, rho_level)
      if (length(a$component) != 1L) stop("survival_colon must have one component.")
      list(
        mu = c(a$mu_log_R, a$mu_log_D),
        Sigma = matrix(
          c(a$sigma_log_R^2, a$covariance, a$covariance, a$sigma_log_D^2),
          nrow = 2L, byrow = TRUE,
          dimnames = list(c("log_R", "log_D"), c("log_R", "log_D"))
        )
      )
    }
    return(list(
      control = make_arm("H0", "control"),
      H1_calibration = make_arm("H1-standard", "treatment"),
      H1_observed = make_arm("H1-observed", "treatment"),
      targets = list(), rho = rho_low
    ))
  }
  if (case == "OAK") {
    make_arm <- function(scenario, arm) {
      a <- rho_get_arm(case, scenario, arm, rho_level)
      if (length(a$component) != 2L) stop("OAK must have two mixture components.")
      list(
        pi = a$weight[1L],
        mu_R1 = a$mu_log_R[1L], sigma_R1 = a$sigma_log_R[1L],
        mu_R2 = a$mu_log_R[2L], sigma_R2 = a$sigma_log_R[2L],
        mu_D = a$mu_log_D[1L], sigma_D = a$sigma_log_D[1L]
      )
    }
    return(list(
      rho = rho_low,
      control = make_arm("H0", "control"),
      H1_calibration = make_arm("H1-standard", "treatment"),
      H1_observed = make_arm("H1-observed", "treatment"),
      targets = list()
    ))
  }
  if (case == "KEYNOTE-598") {
    make_arm <- function(scenario, arm) {
      a <- rho_get_arm(case, scenario, arm, rho_level)
      if (length(a$component) != 2L) stop("KEYNOTE-598 must have two mixture components.")
      sigma <- lapply(seq_along(a$component), function(k) {
        matrix(c(
          a$sigma_log_R[k]^2, a$covariance[k],
          a$covariance[k], a$sigma_log_D[k]^2
        ), nrow = 2L, byrow = TRUE,
        dimnames = list(c("log_R", "log_D"), c("log_R", "log_D")))
      })
      list(
        pi = a$weight, mu_R = a$mu_log_R, sigma_R = a$sigma_log_R,
        mu_D = a$mu_log_D, sigma_D = a$sigma_log_D,
        rho = rho_low, Sigma = sigma
      )
    }
    return(list(
      control = make_arm("H0", "control"),
      H1_calibration = make_arm("H1-standard", "treatment"),
      H1_observed = make_arm("H1-observed", "treatment"),
      targets = list()
    ))
  }
  stop("Unsupported case: ", case)
}

# Each scenario uses its own arm-specific KM-target refits.


rho_make_generator_parameters_by_scenario <- function(case, rho_level) {
  case <- rho_case_name(case)
  if (case == "IMpower131") {
    make_arm <- function(scenario, arm) {
      a <- rho_get_arm(case, scenario, arm, rho_level)
      if (length(a$component) != 1L) stop("IMpower131 must have one component.")
      list(mu = c(a$mu_log_R, a$mu_log_D),
           Sigma = matrix(c(a$sigma_log_R^2, a$covariance,
                            a$covariance, a$sigma_log_D^2), nrow = 2L, byrow = TRUE))
    }
    make_bundle <- function(scenario) {
      control <- make_arm(scenario, "control")
      treatment <- if (scenario == "H0") control else make_arm(scenario, "treatment")
      list(control = control, treatment = treatment,
           H1_calibration = treatment, H1_observed = treatment)
    }
    return(list(H0 = make_bundle("H0"),
                H1_calibration = make_bundle("H1-standard"),
                H1_observed = make_bundle("H1-observed")))
  }
  if (case == "OAK") {
    make_arm <- function(scenario, arm) {
      a <- rho_get_arm(case, scenario, arm, rho_level)
      if (length(a$component) != 2L) stop("OAK must have two mixture components.")
      list(
        pi = a$weight[1L],
        mu_R1 = a$mu_log_R[1L], sigma_R1 = a$sigma_log_R[1L],
        mu_R2 = a$mu_log_R[2L], sigma_R2 = a$sigma_log_R[2L],
        mu_D = a$mu_log_D[1L], sigma_D = a$sigma_log_D[1L]
      )
    }
    make_bundle <- function(scenario) {
      control <- make_arm(scenario, "control")
      treatment <- if (scenario == "H0") control else make_arm(scenario, "treatment")
      list(rho = unique(rho_get_arm(case, scenario, "control", rho_level)$rho),
           control = control, H1_calibration = treatment,
           H1_observed = treatment, targets = list())
    }
    return(list(H0 = make_bundle("H0"),
                H1_calibration = make_bundle("H1-standard"),
                H1_observed = make_bundle("H1-observed")))
  }
  if (case == "KEYNOTE-598") {
    make_arm <- function(scenario, arm) {
      a <- rho_get_arm(case, scenario, arm, rho_level)
      if (length(a$component) != 2L) stop("KEYNOTE-598 must have two mixture components.")
      sigma <- lapply(seq_along(a$component), function(k) {
        matrix(c(a$sigma_log_R[k]^2, a$covariance[k],
                 a$covariance[k], a$sigma_log_D[k]^2), nrow = 2L, byrow = TRUE,
               dimnames = list(c("log_R", "log_D"), c("log_R", "log_D")))
      })
      list(pi = a$weight, mu_R = a$mu_log_R, sigma_R = a$sigma_log_R,
           mu_D = a$mu_log_D, sigma_D = a$sigma_log_D,
           rho = unique(a$rho), Sigma = sigma)
    }
    make_bundle <- function(scenario) {
      control <- make_arm(scenario, "control")
      treatment <- if (scenario == "H0") control else make_arm(scenario, "treatment")
      list(control = control, H1_calibration = treatment,
           H1_observed = treatment, targets = list())
    }
    return(list(H0 = make_bundle("H0"),
                H1_calibration = make_bundle("H1-standard"),
                H1_observed = make_bundle("H1-observed")))
  }
  if (case == "survival_colon") {
    make_arm <- function(scenario, arm) {
      a <- rho_get_arm(case, scenario, arm, rho_level)
      if (length(a$component) != 1L) stop("survival_colon must have one component.")
      list(mu = c(a$mu_log_R, a$mu_log_D),
           Sigma = matrix(c(a$sigma_log_R^2, a$covariance,
                            a$covariance, a$sigma_log_D^2), nrow = 2L, byrow = TRUE,
                          dimnames = list(c("log_R", "log_D"), c("log_R", "log_D"))))
    }
    make_bundle <- function(scenario) {
      control <- make_arm(scenario, "control")
      treatment <- if (scenario == "H0") control else make_arm(scenario, "treatment")
      list(control = control, H1_calibration = treatment,
           H1_observed = treatment, targets = list())
    }
    return(list(H0 = make_bundle("H0"),
                H1_calibration = make_bundle("H1-standard"),
                H1_observed = make_bundle("H1-observed")))
  }
  stop("Unsupported generator case: ", case)
}

rho_make_impower_parameters <- function(scenario, arm, rho_level) {
  a <- rho_get_arm("IMpower131", scenario, arm, rho_level)
  if (length(a$component) != 1L) stop("IMpower131 must have one component.")
  list(
    mu = c(a$mu_log_R, a$mu_log_D),
    Sigma = matrix(
      c(a$sigma_log_R^2, a$covariance, a$covariance, a$sigma_log_D^2),
      nrow = 2L, byrow = TRUE
    )
  )
}

rho_seed_component <- function(x) {
  y <- utf8ToInt(enc2utf8(as.character(x)))
  if (!length(y)) return(0L)
  as.numeric(sum(seq_along(y) * y))
}

rho_data_seed_base <- function(master_seed, case, rho_level) {
  as.integer((as.numeric(master_seed) + 1009 * rho_seed_component(case) +
                7919 * rho_seed_component(rho_level)) %% (.Machine$integer.max - 1))
}

rho_data_seed <- function(master_seed, case, scenario, rho_level) {
  as.integer((rho_data_seed_base(master_seed, case, rho_level) +
                104729 * rho_seed_component(rho_internal_scenario(scenario))) %%
               (.Machine$integer.max - 1))
}

rho_model_seed_base <- function(master_seed, case, method, rho_level) {
  as.integer((as.numeric(master_seed) + 1009 * rho_seed_component(case) +
                7919 * rho_seed_component(rho_level) +
                104729 * rho_seed_component(method)) %% (.Machine$integer.max - 1))
}

rho_parameter_summary <- function(case, rho_level) {
  out <- rho_parameter_table[
    rho_parameter_table$case == rho_case_name(case) &
      rho_parameter_table$rho_label == rho_level, , drop = FALSE
  ]
  if (!nrow(out)) stop("No parameter summary for case/rho level.")
  out
}

# Cache signatures cover the supplied source files and analysis settings.
rho_cache_context <- function(files, ...) {
  if (is.null(names(files)) || any(!nzchar(names(files))) || any(!file.exists(files))) {
    stop("Named cache source files must exist.")
  }
  paths <- normalizePath(unname(files), winslash = "/", mustWork = TRUE)
  list(source_md5 = setNames(unname(tools::md5sum(paths)), names(files)),
       settings = list(...))
}

rho_cache_is_current <- function(stats, context, n) {
  is.list(stats) && length(stats) == n &&
    identical(attr(stats, "rho_cache_context"), context)
}

rho_cache_stamp <- function(stats, context) {
  attr(stats, "rho_cache_context") <- context
  stats
}
