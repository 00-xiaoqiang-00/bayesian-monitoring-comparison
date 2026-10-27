#' Simulate the survival_colon recurrence/death case
#'
#' Generates latent recurrence and death times from arm-specific bivariate
#' lognormal distributions. Randomization is independent 1:1 unless
#' balanced_randomization is TRUE. The returned array is compatible with
#' winratio() and the posterior-probability analysis entries.
#'
#' @export
data.simulation.survival.colon <- function(
    N.sim, N.max, mu.trt, Sigma.trt, mu.ctrl, Sigma.ctrl, cens_upper,
    balanced_randomization = FALSE) {

  N.sim <- as.integer(N.sim)
  N.max <- as.integer(N.max)
  if (N.sim < 1L || N.max < 2L) {
    stop("N.sim must be >= 1 and N.max must be >= 2.")
  }
  if (!is.finite(cens_upper) || cens_upper <= 0) {
    stop("cens_upper must be a positive finite value.")
  }

  validate_arm <- function(mu, Sigma, label) {
    if (length(mu) != 2L || any(!is.finite(mu))) {
      stop(label, " mu must contain two finite values.")
    }
    if (!is.matrix(Sigma) || !identical(dim(Sigma), c(2L, 2L)) ||
        any(!is.finite(Sigma)) ||
        max(abs(Sigma - t(Sigma))) > 1e-10 ||
        min(eigen(Sigma, symmetric = TRUE, only.values = TRUE)$values) <= 0) {
      stop(label, " Sigma must be a finite symmetric positive-definite 2x2 matrix.")
    }
  }
  validate_arm(mu.trt, Sigma.trt, "Treatment")
  validate_arm(mu.ctrl, Sigma.ctrl, "Control")

  make_groups <- function(n) {
    if (!isTRUE(balanced_randomization)) {
      return(sample(c(0, 1), n, replace = TRUE))
    }
    n_trt <- n %/% 2L
    if ((n %% 2L) == 1L) n_trt <- n_trt + stats::rbinom(1L, 1L, 0.5)
    sample(c(rep.int(1, n_trt), rep.int(0, n - n_trt)), replace = FALSE)
  }

  myData <- array(NA_real_, dim = c(N.max, 7L, N.sim))
  for (j in seq_len(N.sim)) {
    group <- make_groups(N.max)
    idx_trt <- which(group == 1L)
    idx_ctrl <- which(group == 0L)
    log_times <- matrix(NA_real_, nrow = N.max, ncol = 2L)

    log_times[idx_trt, ] <- MASS::mvrnorm(
      length(idx_trt), mu = mu.trt, Sigma = Sigma.trt
    )
    log_times[idx_ctrl, ] <- MASS::mvrnorm(
      length(idx_ctrl), mu = mu.ctrl, Sigma = Sigma.ctrl
    )

    myData[, 1L, j] <- exp(log_times[, 1L])
    myData[, 2L, j] <- exp(log_times[, 2L])
    myData[, 3L, j] <- stats::runif(N.max, min = 0, max = cens_upper)
    myData[, 4L, j] <- group
    myData[, 5L, j] <- seq_len(N.max)
  }

  dimnames(myData) <- list(
    NULL,
    Variables = c("recurrence_t", "death_t", "censor_t", "group", "id",
                  "delta1", "delta2"),
    Simulations = paste("Simulation", seq_len(N.sim))
  )
  myData
}
