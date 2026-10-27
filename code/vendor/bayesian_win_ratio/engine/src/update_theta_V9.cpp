#include <RcppArmadillo.h>
#include <mvnorm.h>
#include <truncnorm.h>
#include <wishart.h>
// [[Rcpp::depends(RcppArmadillo, RcppDist)]]
using namespace Rcpp;
using namespace arma;



//' update_theta
//'
//' Performs Gibbs sampling for the bivariate lognormal model, following Supplementary file S1.
//'
//' @param N_iter Number of iterations for gibbs sampler.
//' @param dd The data at a interim analysis point for either treatment or control arm. Either use currentData.trt or currentData.ctrl.
//' @param n  Sample size at a interim analysis point for either treatment or control arm. Either use n.current.trt or n.current.ctrl.
//' @param m0  Prior mean for mu.trt/mu.ctrl.
//' @param L0  Prior covariance for mu.trt/mu.ctrl.
//' @param S0  Prior scale matrix for Sigma.trt/Sigma.ctrl.
//' @param v0  Prior degrees of freedom for Sigma.trt/Sigma.ctrl.
//' @param time_max The upper limit for the recurrence and death time sampled from truncated normal. This will set the upper limit to to time_max rather than Inf.
//' @param init_Sigma An initial Sigma covariance matrix.
//' @return A list with the following components:\tabular{ll}{
//'    \code{MU} \tab A matrix with N_iter rows and two columns. Row v contains a draw of the log-time mean vector from iteration v of the Gibbs sampler. \cr
//'    \tab \cr
//'    \code{post_data} \tab A cube with n rows, two columns and N_iter slices. Slice v, v = 1, ..., N_iter, contains n pairs of event times sampled from bivariate lognormal distribution with mean MU and variance Sigma, where
//'                     Mu and Sigma are generated from the posterior normal and inverse-wishart distribution, respectively, from the v-th iteration of the gibbs sampler.  \cr
//'    \tab \cr
//'    \code{SIGMA} \tab A cube with two rows, two columns and N_iter slices. Slice v contains a draw of the covariance matrix from iteration v of the Gibbs sampler. \cr
//' }
// [[Rcpp::export]]
List update_theta(int N_iter, NumericMatrix dd, int n,
                  arma::vec m0,arma::mat L0, arma::mat S0, double v0,
                  double time_max, arma::mat init_Sigma) {
  int N = dd.nrow();

  arma::mat MU(N_iter, 2); // to store mu;
  arma::cube SIGMA(2, 2, N_iter); // to store Sigma;

  // Log-transform recurrence and death times.
  NumericVector temp_0 = dd( _ , 0);// time to recurrence
  NumericVector temp_1 = dd( _ , 1); // time to death
  arma::vec col_data_0 = as<arma::vec>(temp_0);
  arma::vec col_data_1 = as<arma::vec>(temp_1);

  arma::vec col_data_0_lg = log(col_data_0);// time to recurrence (log)
  arma::vec col_data_1_lg = log(col_data_1);// time to death (log)

  
  
  

  // Initialize theta
  MU.row(0)  = {log(2), log(4)};
  SIGMA.slice(0) = init_Sigma;

  arma::mat R_log(N_iter, n, arma::fill::none);
  arma::mat D_log(N_iter, n, arma::fill::none);


  // Assign the log of the first two columns of dd to the first row of R_log and D_log
  R_log.row(0) = col_data_0_lg.t();
  D_log.row(0) = col_data_1_lg.t();


  // Initialize missing log-event times at the log-censoring time.

  for (int j = 0; j < n; ++j) {
  if (NumericVector::is_na(dd(j, 0))) {
    
     R_log(0, j) = log(dd(j, 2));
  }

  if (NumericVector::is_na(dd(j, 1))) {
    
     D_log(0, j) = log(dd(j, 2));
  }
}


  arma::vec Tbar(2);
  Tbar(0) = mean(R_log.row(0));
  Tbar(1) = mean(D_log.row(0));



  arma::rowvec mu(2);
  arma::mat Sigma(2, 2);
  auto R_update1 = [&](int j) -> NumericVector {
  return rtruncnorm(1, mu(0), sqrt(Sigma(0, 0)), log(dd(j, 2)), log(time_max));
};

  auto R_update2 = [&](int t, int j) -> NumericVector {
  double sd = sqrt(Sigma(0, 0) - Sigma(1, 0) * Sigma(1, 0)/Sigma(1, 1) );
  double mean = mu(0) + (D_log(t, j) - mu(1)) * Sigma(0, 1) / Sigma(1, 1);
  return rtruncnorm(1, mean, sd, D_log(t, j), log(time_max));
};

auto D_update = [&](int t, int j) -> NumericVector {
  double sd = sqrt(Sigma(1, 1) - Sigma(1, 0) * Sigma(1, 0)/Sigma(0, 0) );
  double mean = mu(1) + (R_log(t, j) - mu(0)) * Sigma(0, 1) /  Sigma(1, 1);
  return rtruncnorm(1, mean, sd, log(dd(j, 2)), log(time_max));
};


arma::mat Ln(2, 2); arma::vec mun;arma::mat Sn;
arma::cube post_data(n, 2, N_iter, arma::fill::none);

 for (int t = 1; t < N_iter; ++t) {
   try {


      Ln = arma::inv(arma::inv(L0) + n*arma::inv(SIGMA.slice(t-1)));
      mun = Ln *( arma::inv(L0)*m0 + n*arma::inv(SIGMA.slice(t-1))* Tbar);
      mu = rmvnorm(1, mun, Ln);
      MU.row(t) = mu;

      arma::mat A(2, R_log.n_cols);
      A.row(0) = R_log.row(t-1) - mu(0);
      A.row(1) = D_log.row(t-1) - mu(1);

      Sn = S0 + A * A.t();
      Sn = Sn + 1e-6 * arma::eye(2, 2); // Diagonal regularization of the scale matrix.
      Sigma = riwish(v0 + n, Sn); //Sn is the scale matrix;
      SIGMA.slice(t) = Sigma;


      post_data.slice(t) = exp(rmvnorm(n, mu.t(), Sigma));

      R_log.row(t) = R_log.row(t-1);
      D_log.row(t) = D_log.row(t-1);

      for (int j = 0; j < N; ++j) {
        if (NumericVector::is_na(dd(j, 0)) && !NumericVector::is_na(dd(j, 1))) {
          R_log(t, j) = R_update2(t, j)[0];
        } else if (!NumericVector::is_na(dd(j, 0)) && NumericVector::is_na(dd(j, 1))) {
          D_log(t, j) = D_update(t, j)[0];
        } else if (NumericVector::is_na(dd(j, 0)) && NumericVector::is_na(dd(j, 1))) {
          R_log(t, j) = R_update1(j)[0];
          D_log(t, j) = D_update(t, j)[0];
        }
      }

      Tbar(0) = mean(R_log.row(t));
      Tbar(1) = mean(D_log.row(t));


 } catch (std::exception &ex) {
      Rcpp::Rcout << "Error in iteration: " << t << std::endl;
      Rcpp::Rcout << ex.what() << std::endl;
    }
 }
  return List::create(Named("MU") = MU,
                      Named("post_data") = post_data,
                      Named("SIGMA") = SIGMA);

}
