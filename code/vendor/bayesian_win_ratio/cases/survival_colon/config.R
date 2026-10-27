# Case config: survival_colon
# Obs versus Lev+5FU, with latent recurrence and death times generated from
# arm-specific bivariate lognormal distributions.

case_id <- "survival_colon"
case_note <- paste(
  "R survival::colon Stage III trial case; Obs is control and Lev+5FU is",
  "treatment; H1_calibration uses a common covariance and literature-anchored",
  "recurrence/death hazard ratios; H1_observed retains the fitted treatment",
  "covariance"
)

# ---- trial timing and sample size ----
N.max               <- 620L
interim_cohort <- c(248L, 372L, 496L)   # 40%, 60%, 80%; N.max is final
enrollment_duration <- 43
recruit.int         <- enrollment_duration / N.max
cens_upper          <- 1900
final_analysis_time <- 109

# ---- generator selection ----
generator_type         <- "survival_colon_bivariate_lognormal_v1"
generator_version      <- 1L
balanced_randomization <- FALSE

# ---- calibrated case parameters: c(log latent recurrence, log death) ----
mu_ctrl <- c(log(38.1268), log(64.5907))
Sigma_ctrl <- matrix(
  c(3.29767175, 2.19701400,
    2.19701400, 1.77436755),
  nrow = 2, byrow = TRUE
)

# H1_calibration changes only the treatment location parameters. This target
# reproduces cause-specific recurrence HR about 0.59 and death HR about 0.67.
mu_standard_trt <- c(log(92.6930), log(116.2198))
Sigma_standard_trt <- Sigma_ctrl

# H1_observed reproduces the curve-calibrated Lev+5FU arm, including its
# observed variance/covariance structure.
mu_observed_trt <- c(log(87.1562), log(111.4653))
Sigma_observed_trt <- matrix(
  c(4.46376802, 3.41770993,
    3.41770993, 2.92813972),
  nrow = 2, byrow = TRUE
)

calibration_target_recurrence_hr <- 0.59
calibration_target_death_hr <- 0.67
observed_population_recurrence_hr <- 0.5990
observed_population_death_hr <- 0.6240

# ---- scenario interface used by calibration/validation scripts ----
mu_h0_ctrl <- mu_ctrl
mu_h0_trt <- mu_ctrl
Sigma.ctrl <- Sigma_ctrl
Sigma_h0_trt <- Sigma_ctrl

mu_h1_ctrl <- mu_ctrl
mu_h1_trt <- mu_standard_trt
Sigma_h1_trt <- Sigma_standard_trt

mu_h1_observed_ctrl <- mu_ctrl
mu_h1_observed_trt <- mu_observed_trt
Sigma_h1_observed_trt <- Sigma_observed_trt

# ---- MCMC / simulation ----
design    <- 1L
m0        <- c(0, 0)
L0        <- diag(10^6, 2)
v0        <- 4
S0        <- diag(10^(-6), 2)
time_max  <- 600
thin_MCMC <- 5
Niter     <- 10000
L_sim <- 10000
mc.cores <- if (.Platform$OS.type == "windows") 1L else 70L
seed <- 26073004L

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

# H1_observed is evaluated after design calibration and never selects the
# decision boundary.
scenarios_extra <- list(
  list(
    name = "H1_observed",
    role = "alternative",
    generator_scenario = "H1_observed",
    note = "Observed Lev+5FU curve-calibrated distribution versus Obs control",
    cache_file = "calib_stats_H1_observed.rds",
    mu.trt = mu_h1_observed_trt,
    mu.ctrl = mu_h1_observed_ctrl,
    Sigma.trt = Sigma_h1_observed_trt,
    Sigma.ctrl = Sigma.ctrl
  )
)


# Analysis schedule and cache identity.
design_version <- "calibration_final_all_n_v2"
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
