# Case config: OAK_A_ITT
# Arm-specific two-component lognormal recurrence models.
# H1-standard component log-scale sigmas and rho equal control.

case_id   <- "OAK_A_ITT"
case_note <- paste(
  "OAK arm A ITT; recurrence uses arm-specific two-component lognormal mixtures;",
  "H1-standard component dispersions equal control; rho=0.7"
)

# ---- trial timing and sample size ----
N.max              <- 850L
interim_cohort <- c(340L, 510L, 680L)   # 40%, 60%, 80%; N.max is final
recruit.int        <- 8.6 / 850
cens_upper         <- 109.3
final_analysis_time <- 27.9

# ---- generator selection ----
generator_type          <- "oak_lognormal_mixture_v1"
generator_version       <- 5L
balanced_randomization  <- TRUE
rho_selected            <- 0.7

oak_generator_parameters <- list(
  rho = rho_selected,
  # H1_observed treatment recurrence mixture fitted to the OAK treatment PFS curve.
  pi       = 0.2623242204350347,
  mu_R1    = 0.3645597370078664,
  sigma_R1 = 0.0806101279283579,
  mu_R2    = 1.677883933719392,
  sigma_R2 = 1.2693769377027608,
  # Locked treatment death marginal
  mu_D_trt    = 2.636754680552583,
  sigma_D_trt = 1.1709,
  # Control recurrence mixture fitted to the OAK control PFS curve.
  pi_ctrl       = 0.13123216948698738,
  mu_R1_ctrl    = 0.2666273699841142,
  sigma_R1_ctrl = 0.3000000000000000,
  mu_R2_ctrl    = 1.6274656470707884,
  sigma_R2_ctrl = 0.8489991311979127,
  # Locked control death marginal
  mu_D_ctrl    = 2.3088355183235585,
  sigma_D_ctrl = 1.0177,
  cens_upper   = cens_upper
)

# H1_calibration is a hypothetical clinically meaningful PFS-effect scenario.
# It targets PFS HR=0.70 and treatment median PFS=5.7 months while retaining
# the locked OAK treatment death marginal (protocol/design OS target HR=0.73).
# PFS=min(R,D), so these are targets for the derived PFS endpoint, not for the
# latent recurrence time R alone. Published OAK PFS HR=0.95 belongs to
# H1_observed only.
observed_os_hr                    <- 0.727
calibration_target_os_hr          <- 0.73
calibration_target_pfs_hr         <- 0.70
calibration_target_pfs_median_trt <- 5.7
calibration_target_recurrence_hr  <- NA_real_
observed_pfs_hr                   <- 0.95

# H1-standard uses the prespecified treatment locations and mixture weight.
# Its component log-scale dispersions equal control; realized PFS/OS and
# method-specific estimands follow from the complete mixture distribution.
oak_generator_parameters$calibration_target <- list(
  pfs_hr = calibration_target_pfs_hr,
  pfs_median_trt = calibration_target_pfs_median_trt,
  curve_time_scale = 1.247936722762762,
  curve_survival_power = 0.8686138491694806,
  hr_evaluation_horizon = 24.0
)
oak_generator_parameters$calibration <- with(oak_generator_parameters, list(
  pi       = 0.10183833564436579,
  mu_R1    = 0.510611455725049,
  sigma_R1 = sigma_R1_ctrl,
  mu_R2    = 1.9615174106704887,
  sigma_R2 = sigma_R2_ctrl,
  mu_D_trt = mu_D_trt,
  sigma_D_trt = sigma_D_ctrl
))

mixture_moments <- function(p, rho) {
  mean_R <- p$pi * p$mu_R1 + (1 - p$pi) * p$mu_R2
  var_R <- p$pi * (p$sigma_R1^2 + (p$mu_R1 - mean_R)^2) +
    (1 - p$pi) * (p$sigma_R2^2 + (p$mu_R2 - mean_R)^2)
  cov_RD <- rho * p$sigma_D_trt *
    (p$pi * p$sigma_R1 + (1 - p$pi) * p$sigma_R2)
  list(
    mean_R = mean_R,
    Sigma = matrix(c(var_R, cov_RD, cov_RD, p$sigma_D_trt^2),
                   nrow = 2, byrow = TRUE)
  )
}

# ---- analysis covariance ----
# winratio() uses a single bivariate-lognormal working covariance. For each
# mixture DGM, the mu/Sigma interface receives exact first and second
# moments of (log R, log D).
control_parameters <- with(oak_generator_parameters, list(
  pi = pi_ctrl,
  mu_R1 = mu_R1_ctrl,
  sigma_R1 = sigma_R1_ctrl,
  mu_R2 = mu_R2_ctrl,
  sigma_R2 = sigma_R2_ctrl,
  mu_D_trt = mu_D_ctrl,
  sigma_D_trt = sigma_D_ctrl
))
control_moments <- mixture_moments(control_parameters, rho_selected)
Sigma.ctrl <- control_moments$Sigma

# Under H0 both groups use the same control-arm recurrence mixture and death
# marginal.
Sigma_h0_trt <- Sigma.ctrl

observed_parameters <- oak_generator_parameters[
  c("pi", "mu_R1", "sigma_R1", "mu_R2", "sigma_R2",
    "mu_D_trt", "sigma_D_trt")
]
observed_moments <- mixture_moments(observed_parameters, rho_selected)
calibration_parameters <- oak_generator_parameters$calibration
calibration_moments <- mixture_moments(calibration_parameters, rho_selected)

# The mu/Sigma fields are moment summaries used by the analysis interface
# and cache metadata. Data are generated from the component parameters above.
Sigma_h1_trt <- calibration_moments$Sigma
Sigma_h1_observed_trt <- observed_moments$Sigma

# The scenario interface uses the log-time mixture means below;
# data generation uses oak_generator_parameters.
mu_h0_trt  <- c(control_moments$mean_R,
                oak_generator_parameters$mu_D_ctrl)
mu_h0_ctrl <- mu_h0_trt
mu_h1_trt  <- c(calibration_moments$mean_R,
                oak_generator_parameters$calibration$mu_D_trt)
mu_h1_ctrl <- mu_h0_ctrl
mu_h1_observed_trt <- c(observed_moments$mean_R,
                        oak_generator_parameters$mu_D_trt)
mu_h1_observed_ctrl <- mu_h0_ctrl

# ---- MCMC / simulation ----
design       <- 1L
m0           <- c(0, 0)
L0           <- diag(10^6, 2)
v0           <- 4
S0           <- diag(10^(-6), 2)
time_max     <- 120
thin_MCMC    <- 5
Niter        <- 10000
L_sim <- 10000
mc.cores <- if (.Platform$OS.type == "windows") 1L else 70L
seed <- 26073003L

# ---- decision / calibration targets ----
eta_interim        <- 1
eta_final          <- 1
enable_superiority <- FALSE
alpha <- 0.05
power_target       <- 0.80
# Calibration candidate grids.
c_F_grid <- seq(0.80, 0.99, by = 0.005)
gamma_E_grid <- seq(0.00, 4.00, by = 0.05)
rebuild_pp_cache   <- FALSE

# H1_observed is evaluated after calibration but is not used to choose the
# decision boundary.
scenarios_extra <- list(
  list(
    name = "H1_observed",
    role = "alternative",
    generator_scenario = "H1_observed",
    note = "Observed OAK treatment mixture versus locked control",
    cache_file = "calib_stats_H1_observed.rds",
    mu.trt = mu_h1_observed_trt,
    mu.ctrl = mu_h1_observed_ctrl,
    Sigma.trt = Sigma_h1_observed_trt,
    Sigma.ctrl = Sigma.ctrl
  )
)


# Analysis schedule and cache identity.
design_version <- "calibration_final_all_n_OAK_equal_H1_variance_v3"
analysis_row_rule <- "interim_n_minus_1_final_all_n_v2"
interim_fractions <- c(0.40, 0.60, 0.80)

# Formal calibration provenance. Decision parameters are selected only after
# all H0/H1 PP caches are complete.
calibration_role <- "formal_parameter_calibration"
calibration_outer_replicates <- 10000L
calibration_final_rule <- "PP_final > c_F"
calibration_futility_rule <- "PP_k <= b_k"
calibration_selection_constraints <- paste(
  "Candidate tables report H0 go probability, H1-standard power, PET and ESS;",
  "this driver generates candidates without automatic parameter recommendation;",
  "formal validation uses LOCKED_PARAMETER_MANIFEST.tsv"
)
