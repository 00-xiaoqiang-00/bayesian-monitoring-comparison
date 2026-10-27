# Case config: KEYNOTE-598
# Pembrolizumab + ipilimumab (treatment) versus pembrolizumab + placebo
# (control). Latent progression R and death D use one shared component label.

case_id <- "keynote598"
case_note <- paste(
  "KEYNOTE-598 numerical curve-reproduction DGM; R means progression time;",
  "PFS is min(R,D); a single latent class K selects both R and D component",
  "parameters; mixture classes are numerical approximation components and",
  "must not be interpreted as clinical subgroups"
)

# ---- trial timing and sample size ----
N.max               <- 568L
interim_cohort <- c(228L, 340L, 454L)   # 40%, 60%, 80%; N.max is final
enrollment_duration <- 19.3
recruit.int         <- enrollment_duration / N.max
final_analysis_time <- 31.7
cens_upper          <- 141.8559446075

# With the project convention recruit.int=enrollment_duration/N.max, actual
# deterministic entry times run from 0 to 19.266021 months and final
# administrative follow-up runs from 12.433979 to 31.7 months.  cens_upper is
# a separately calibrated independent Uniform random-censoring parameter.
# The three interim cohorts are working 40/60/80% design looks, not reported
# KEYNOTE-598 interim sample sizes.

# ---- generator selection ----
generator_type         <- "keynote598_shared_component_lognormal_v1"
generator_version      <- 1L
balanced_randomization <- TRUE

component_sigma <- function(sigma_R, sigma_D, rho = 0.5) {
  matrix(
    c(
      sigma_R^2, rho * sigma_R * sigma_D,
      rho * sigma_R * sigma_D, sigma_D^2
    ),
    nrow = 2L, byrow = TRUE
  )
}

# ---- observed curve-reproduction parameters ----
# K ~ Categorical(pi); conditional on K=k,
# (log R, log D)' ~ N2(c(mu_R[k], mu_D[k]), Sigma[[k]]).
# rho is the within-component correlation of log R and log D.

keynote598_control <- list(
  pi = c(0.2070545742, 0.7929454258),
  mu_R = c(log(2.4837219024), log(16.4600817241)),
  sigma_R = c(0.4340057961, 1.5317036442),
  mu_D = c(log(2.7290953670), log(28.6890202294)),
  sigma_D = c(1.5000000000, 0.9670082477),
  rho = 0.5
)
keynote598_control$Sigma <- list(
  component_sigma(0.4340057961, 1.5000000000),
  component_sigma(1.5317036442, 0.9670082477)
)

keynote598_treatment_observed <- list(
  pi = c(0.1936579096, 0.8063420904),
  mu_R = c(log(2.0212342835), log(15.2003022720)),
  sigma_R = c(0.1500000000, 1.3810016053),
  mu_D = c(log(2.5120209874), log(29.2329569146)),
  sigma_D = c(1.5000000000, 1.2402934461),
  rho = 0.5
)
keynote598_treatment_observed$Sigma <- list(
  component_sigma(0.1500000000, 1.5000000000),
  component_sigma(1.3810016053, 1.2402934461)
)

# ---- controlled beneficial alternative used only for design calibration ----
# It retains the control weights, component SDs and rho, and shifts both
# component locations by common amounts.  Through 31.7 months this gives
# deterministic curve-level Cox HRs of 0.69 (PFS) and 0.70 (OS).
delta_R_standard <- 0.5946115916
delta_D_standard <- 0.3773863423
keynote598_treatment_standard <- keynote598_control
keynote598_treatment_standard$mu_R <-
  keynote598_control$mu_R + delta_R_standard
keynote598_treatment_standard$mu_D <-
  keynote598_control$mu_D + delta_D_standard

calibration_target_pfs_hr <- 0.69
calibration_target_os_hr  <- 0.70
calibration_horizon       <- 31.7

# ---- analysis interface ----
# The mixture itself is generated from the structured lists above.  These
# aggregate log-time moments only provide coherent initialization values to the
# single-bivariate-Normal Bayesian analysis and cache fingerprint.
mu_ctrl <- c(2.4093606993, 2.8694094224)
Sigma_ctrl <- matrix(
  c(
    2.4865535444, 1.3851046341,
    1.3851046341, 2.1160248702
  ),
  nrow = 2L, byrow = TRUE
)
mu_standard_trt <- c(3.0039722909, 3.2467957647)
Sigma_standard_trt <- Sigma_ctrl
mu_observed_trt <- c(2.3305897681, 2.9000197243)
Sigma_observed_trt <- matrix(
  c(
    2.1778491933, 1.4855762600,
    1.4855762600, 2.6166897181
  ),
  nrow = 2L, byrow = TRUE
)

mu_h0_ctrl  <- mu_ctrl
mu_h0_trt   <- mu_ctrl
Sigma.ctrl  <- Sigma_ctrl
Sigma_h0_trt <- Sigma_ctrl

mu_h1_ctrl  <- mu_ctrl
mu_h1_trt   <- mu_standard_trt
Sigma_h1_trt <- Sigma_standard_trt

mu_h1_observed_ctrl <- mu_ctrl
mu_h1_observed_trt  <- mu_observed_trt
Sigma_h1_observed_trt <- Sigma_observed_trt

# ---- Bayesian analysis / simulation settings ----
design    <- 1L
m0        <- c(0, 0)
L0        <- diag(10^6, 2)
v0        <- 4
S0        <- diag(10^(-6), 2)
time_max  <- 120
thin_MCMC <- 5
Niter     <- 10000
L_sim <- 10000
mc.cores <- if (.Platform$OS.type == "windows") 1L else 70L
seed <- 26073002L

# ---- decision / calibration targets ----
eta_interim        <- 1
eta_final          <- 1
enable_superiority <- FALSE
alpha <- 0.05
power_target       <- 0.80
c_F_grid <- seq(0.80, 0.99, by = 0.005)
gamma_E_grid <- seq(0.00, 4.00, by = 0.05)
rebuild_pp_cache   <- FALSE

# H1_observed is a real-case curve-reproduction check; it does not select the
# calibrated decision boundary.
scenarios_extra <- list(
  list(
    name = "H1_observed",
    role = "alternative",
    generator_scenario = "H1_observed",
    note = paste(
      "Observed pembrolizumab + ipilimumab curve fit versus pembrolizumab +",
      "placebo; not a beneficial calibration alternative"
    ),
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
