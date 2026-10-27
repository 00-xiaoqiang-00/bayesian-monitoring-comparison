#include <RcppArmadillo.h>
#include <truncnorm.h>
// [[Rcpp::depends(RcppArmadillo, RcppDist)]]
using namespace Rcpp;
using namespace arma;


  
  
//' update_theta_univariate
//'
//' Performs Gibbs sampling and returns posterior mean, variance, and event times (log) for the conventional designs, following Supplementary file S1.
//' 
//' @param N_iter Number of iterations for gibbs sampler.
//' @param dd The data at a interim analysis point for either treatment or control arm. Either use currentData.trt or currentData.ctrl.
//' @param n  Sample size at a interim analysis point for either treatment or control arm. Either use n.current.trt or n.current.ctrl.
//' @param m0  Prior mean for mu.trt/mu.ctrl.
//' @param L0  Prior variance for mu.trt/mu.ctrl.
//' @param v0  v0/2 is the prior shape for Sigma.trt/Sigma.ctrl.
//' @param S0  v0*S0/2 is the prior scale for Sigma.trt/Sigma.ctrl.
//' @param time_max The upper limit for the recurrence and death time sampled from truncated normal. This will set the upper limit to to time_max rather than Inf.
//' @param init_Sigma A numeric value represents the time to event (logarithm) for the control arm.
//' @return A list of theta. theta is A matrix with N_iter rows and two columns. Row v, v = 1, ..., N_iter, contains posterior mean vector generated from the v-th iteration of the gibbs sampler.  
// [[Rcpp::export]]
List update_theta_univariate(int N_iter, NumericMatrix dd, int n, double L0,
                  double m0,  double v0,double S0, double time_max, double init_Sigma) {
  int N = dd.nrow();
  arma::mat theta(N_iter, 2, fill::none);

  // Preprocess data to remove NAs
  NumericVector temp_0 = dd( _ , 0);
  arma::vec col_data_0 = as<arma::vec>(temp_0);

  arma::vec col_data_0_lg = log(col_data_0);
  
  arma::vec col_data_0_rm = col_data_0_lg.elem(find_finite(col_data_0_lg));
 
  // Initialize theta.row(0)
  theta.row(0)  = {log(2), init_Sigma}; 

  arma::mat y_log(N_iter, n, arma::fill::none);
  
  y_log.row(0) = col_data_0_lg.t(); 
  
  // Initialize missing log-event times at the log-censoring time.

  for (int j = 0; j < n; ++j) {
  if (NumericVector::is_na(dd(j, 0))) {
    
    y_log(0, j) = log(dd(j, 1));
  }
}
  

  double Tbar = mean(y_log.row(0));
 

  double mu;double sigma;
  double mun; double t2n; double sd_for_rnorm;
  double vn; double Sn;
  
  auto y_update = [&](int j) -> NumericVector {
  return rtruncnorm(1, mu, sigma, log(dd(j, 1)), log(time_max));
};


 for (int t = 1; t < N_iter; ++t) {
   try {

      mun = (m0/L0 + n*Tbar/theta(t-1, 1))/(1/L0 + n/theta(t-1, 1));
      t2n = 1/(1/L0 + n/theta(t-1, 1));
      sd_for_rnorm = sqrt(t2n);  
      mu = R::rnorm(mun, sd_for_rnorm);
      
      vn = v0 + n;
      Sn = v0*S0 + (y_log.n_cols-1)*var(y_log.row(t-1)) +  n * (Tbar - mu) * (Tbar - mu);
      sigma = R::rgamma(vn/2, 2/Sn); // shape, scale;
      sigma = 1/sigma;
      theta.row(t) = {mu, sigma};
    

      y_log.row(t) = y_log.row(t-1);
      
      for (int j = 0; j < N; ++j) {
        if (NumericVector::is_na(dd(j, 0))) {
          y_log(t, j) = y_update(j)[0];
        }
      }
      
      Tbar = mean(y_log.row(t));
 } catch (std::exception &ex) {
      Rcpp::Rcout << "Error in iteration: " << t << std::endl;
      Rcpp::Rcout << ex.what() << std::endl;
    }
 }
  return List::create(Named("theta") = theta);
  
}

