

#' est_cens_time
#' 
#' This function is an intermediate function used inside the winratio function to generate censoring time from estimated censoring distribution.
#' 
#' @param fit   Estimated censoring distribution.
#' @param n   A numeric value indicating number of censoring time to generate. Either n.current.trt or n.current.ctrl.
#' @returns Censoring time for each subject.
#' @examples
#' \dontrun{est_cens_time(trt_cens.est,n.current.trt)}


est_cens_time = function(fit,n){
  Cens.generated = as.numeric(stats::quantile(fit, probs = stats::runif(n) )$quantile)
  Cens.generated[is.na(Cens.generated)] = Inf # replace NAs with Inf
  return (Cens.generated)
}

#' make_pos_def
#' 
#'  This function is an intermediate function used inside the winratio function to ensure the sample covariance matrix is positive definite
#' @param mat   Sample covariance matrix
#' @param eps a small value to ensure the matrix is positive-definite
#' @returns A positive definite matrix.
#' @examples
#' \dontrun{make_pos_def(temp)}  
make_pos_def <- function(mat, eps = 1e-6) {
  if (!isSymmetric(mat)) mat <- (mat + t(mat)) / 2
  eig <- eigen(mat)
  eig$values[eig$values < eps] <- eps
  return(eig$vectors %*% diag(eig$values) %*% t(eig$vectors))
}


#' winratio
#' 
#' Computes the posterior probability of (WR > eta) at an analysis look.

#' @useDynLib BayesianWinRatio, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @exportPattern "^[[:alpha:]]+"
#' @param currentdd   The data at a interim analysis point.
#' @param n_current   Sample size at a interim analysis point.
#' @param Time_current Time at a interim analysis point.
#' @param Time_entry  Time at entry for each patient at a interim analysis point.
#' @param m0  Prior mean for mu.
#' @param L0  Prior covariance for mu.
#' @param v0  Prior degrees of freedom for Sigma.
#' @param S0  Prior scale matrix for Sigma.
#' @param time_max The upper limit for the recurrence and death time sampled from truncated normal. This will set the upper limit to to time_max rather than Inf.
#' @param eta A pre-specified lower bound of acceptable performance based on prespecified information.
#' @param lambda  Cutoff parameter.
#' @param N.max   Maximum number of patients to enroll.
#' @param thin_MCMC  Thinning degree.
#' @param Niter Number of iterations for gibbs sampler.
#' @param Sigma.ctrl A matrix representing the variance-covariance matrix of the time to each event (logarithm) for the control arm. 
#' @returns A list with the following components:\tabular{ll}{
#'    \code{probs} \tab A numeric value, which is the posterior probability of \eqn{(\widehat{WR} > eta)}.  \cr
#'    \tab \cr
#'    \code{cutoff} \tab A numeric value. cutoff = (n_current/N.max)**lambda.  \cr
#'    \tab \cr
#'    \code{WR} \tab A numeric value. WR = the average posterior win ratio.  \cr
#' }
#' @examples
#' \dontrun{
#' winratio(currentdd = currentData,n_current = n.current, Time_current = Time.current,
#' Time_entry = Time.entry,m0 = c(0,0),L0 = diag(10^6, 2),v0 = 4,S0 = diag(10^(-6), 2),
#' time_max = 20,eta = 1,lambda = 0.25,N.max = 100, thin_MCMC = 5,
#' Niter = 100000,Sigma.ctrl = matrix(c(1,0.5,0.5,1), nrow = 2, byrow = T))
#' } 

winratio = function(currentdd,n_current,Time_current,
                    Time_entry, m0,L0,v0,S0,
                    time_max,eta,lambda,N.max, thin_MCMC, Niter,Sigma.ctrl){

  # Use the supplied analysis rows: n_current - 1 at interim looks and
  # all n_current participants at the final calendar-time analysis.
  n_obs <- nrow(currentdd)
  if (!is.finite(n_obs) || n_obs < 2L) {
    stop("currentdd must contain at least two analyzed participants.")
  }
  if (length(Time_entry) < n_obs) {
    stop("Time_entry is shorter than the number of rows in currentdd.")
  }
  Time_entry_obs <- Time_entry[seq_len(n_obs)]
  
  # update actual censoring time at this interim look
  currentdd[,3] = pmin(Time_current - Time_entry_obs, currentdd[,3])
  currentdd[,6] <- as.integer(currentdd[,1] < currentdd[,2] & currentdd[,1] <= currentdd[,3])
  currentdd[,7] <-  as.integer(currentdd[,2] <= currentdd[,3])
 
 
  # update recurrence and death  and the censoring time to reflect the observed data at this interim look    
  currentdd[,1] = ifelse(currentdd[,6]==1, currentdd[,1], NA)
  currentdd[,3]  = pmin(currentdd[,2] , currentdd[,3] )
  currentdd[,2] = ifelse(currentdd[,7]==1, currentdd[,2], NA)
  
  # number enrolled to each group
  n.current.trt = sum(currentdd[,4])
  n.current.ctrl = n_obs - n.current.trt
  
  # separate data to trt and control
  currentData.trt = currentdd[currentdd[,4]==1,]
  currentData.ctrl = currentdd[currentdd[,4]==0,]
  
  # Step 1:Update posterior theta 
  ## obtain sample covariance as an input for starting value of Sigma, if NA use Sigma.ctrl
  
  temp = stats::cov( log(currentData.trt[,1:2]) , use = "pairwise.complete.obs")
  if (any(is.na(temp))){
    init_trt_Sigma = make_pos_def(Sigma.ctrl)
  }else {
    init_trt_Sigma = make_pos_def(temp)
  }   
  
  trt.post = update_theta(N_iter = Niter, dd = currentData.trt, 
                          n = n.current.trt, m0 = m0, L0 = L0, S0 = S0, v0 = v0, 
                          time_max = time_max, init_Sigma = init_trt_Sigma)

  temp = stats::cov( log(currentData.ctrl[,1:2]) , use = "pairwise.complete.obs")
  if (any(is.na(temp))){
    init_ctrl_Sigma =  make_pos_def(Sigma.ctrl)
  }else {
    init_ctrl_Sigma = make_pos_def(temp)
  }   
  ctrl.post = update_theta(N_iter = Niter, dd = currentData.ctrl, 
                           n = n.current.ctrl,m0 = m0, L0 = L0, S0 = S0, v0 = v0, 
                           time_max = time_max, init_Sigma = init_ctrl_Sigma)
  burn_MCMC = as.integer(0.3*Niter)
  idxs <- seq(burn_MCMC, Niter, by = thin_MCMC)
  M_iter = length(idxs)
  
  # Step 2: estimating censoring distribution
  trt_cens.est = survival::survfit(survival::Surv(currentData.trt[,"censor_t"], 1 - currentData.trt[,"delta2"]) ~ 1 )
  ctrl_cens.est =survival::survfit(survival::Surv(currentData.ctrl[,"censor_t"], 1 - currentData.ctrl[,"delta2"]) ~ 1 )
  
  # Step 3: sample R and D for each group, and generate censoring time for each patient
  ### initialize posterior data
  postData <-array(dim=c(n_obs,6, M_iter));  
  # Create dimension names
  dim_names <- list(
    NULL,
    Variables = c("R", "D", "C", "group",  "delta1", "delta2"),
    Iterations = paste("Iteration", 1:M_iter)
  )
  # Assign the dimension names to the array
  dimnames(postData) <- dim_names
  
  postData[1:n.current.trt,1:2,] = trt.post$post_data[,,idxs]
  postData[(n.current.trt+1):n_obs,1:2,] = ctrl.post$post_data[,,idxs]
  
  res1 = replicate(M_iter, est_cens_time(trt_cens.est,n.current.trt), simplify = F)
  res2 = replicate(M_iter, est_cens_time(ctrl_cens.est, n.current.ctrl), simplify = F)
  postData[1:n.current.trt,3,]  = matrix(unlist(res1), nrow = n.current.trt, ncol = M_iter, byrow = F)
  postData[(n.current.trt+1):n_obs,3,] =  matrix(unlist(res2), nrow = n.current.ctrl, ncol = M_iter, byrow = F)
  
  postData[1:n.current.trt,4,] = 1
  postData[(n.current.trt+1):n_obs,4,] = 0
  postData[,5,] <- as.integer(postData[,1,] < postData[,2,] & postData[,1,] <= postData[,3,])
  postData[, 6,] <- as.integer(postData[, 2,] <= postData[, 3,])
  postData[, 1,] <- ifelse( postData[, 5,]==1,  postData[, 1,], NA)
  postData[, 3,] <- pmin(postData[, 2,], postData[, 3,])
  postData[, 2,] <- ifelse( postData[, 6,]==1,  postData[, 2,], NA)
  

  # Step 4: calculate the IPCW-adjusted WR
time.trt <- time.ctrl <- surv.trt <- surv.ctrl <- vector("list", M_iter)

for (j.iter in 1:M_iter) {
  C      <- postData[, 3, j.iter] 
  delta2 <- postData[, 6, j.iter]
  
  # Estimate the treatment-arm censoring survival curve.
  # 提取 TRT 组的删失状态 (1-delta2)
  status_trt <- 1 - delta2[1:n.current.trt]
  
  # 只有当存在至少一个删失事件，且样本量大于0时才拟合
  if (length(status_trt) > 0 && sum(status_trt) > 0) {
    fit.trt <- tryCatch({
      survival::survfit(survival::Surv(C[1:n.current.trt], status_trt) ~ 1)
    }, error = function(e) NULL)
    
    if (!is.null(fit.trt)) {
      s_trt <- summary(fit.trt)
      time.trt[[j.iter]] <- s_trt$time
      surv.trt[[j.iter]] <- s_trt$surv
    } else {
      time.trt[[j.iter]] <- numeric(0)
      surv.trt[[j.iter]] <- numeric(0)
    }
  } else {
    # Empty censoring curves contribute survival probability one.
    time.trt[[j.iter]] <- numeric(0)
    surv.trt[[j.iter]] <- numeric(0)
  }

  # Estimate the control-arm censoring survival curve.
  status_ctrl <- 1 - delta2[(n.current.trt + 1):n_obs]
  
  if (length(status_ctrl) > 0 && sum(status_ctrl) > 0) {
    fit.ctrl <- tryCatch({
      survival::survfit(survival::Surv(C[(n.current.trt + 1):n_obs], status_ctrl) ~ 1)
    }, error = function(e) NULL)
    
    if (!is.null(fit.ctrl)) {
      s_ctrl <- summary(fit.ctrl)
      time.ctrl[[j.iter]] <- s_ctrl$time
      surv.ctrl[[j.iter]] <- s_ctrl$surv
    } else {
      time.ctrl[[j.iter]] <- numeric(0)
      surv.ctrl[[j.iter]] <- numeric(0)
    }
  } else {
    time.ctrl[[j.iter]] <- numeric(0)
    surv.ctrl[[j.iter]] <- numeric(0)
  }
}
  
  
  ## getting the survival probability for event times, from the censoring distribution from both groups
  ## calculate WR
  result = compare(M_iter = M_iter, n_current_ctrl = n.current.ctrl, 
                                n_current_trt = n.current.trt, postData = postData,
                                time_trt = time.trt, time_ctrl = time.ctrl, surv_trt = surv.trt, surv_ctrl = surv.ctrl)
  # Convert Inf values to NA
  result$WR[result$WR == Inf] <- NA
  WR = mean(result$WR, na.rm = T)
  probs = mean(result$WR > eta, na.rm  = T)  
  cutoff = (n_current/N.max)**lambda
  
  rm(result,postData, res1, res2)
  gc()
  
  return(list(probs = probs, cutoff = cutoff, WR = WR))
}

