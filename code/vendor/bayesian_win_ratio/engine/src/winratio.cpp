#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]


using namespace Rcpp;
using namespace arma;


// calculateSurvival
// Evaluate the censoring survival curve at a comparison time.

double calculateSurvival(const arma::vec& timeData, const arma::vec& survData, double mytime, const arma::vec& c_time) {
    double surv = std::numeric_limits<double>::quiet_NaN();
    arma::uword n = timeData.n_elem;  
    
    if (!std::isnan(mytime) && n > 0  && arma::max(c_time) >= mytime) { // Censoring survival is undefined beyond the observed support.
         if (arma::min(timeData) > mytime) {
            surv = 1.0;
        }  else {
            arma::uword index = arma::max(arma::find(timeData <= mytime));
            surv = survData(index);  
        }
        
       
    }

    return surv;
}
    

//' compare
//'
//' This is an interim function used within the winratio function that calculates the IPCW-adjusted win ratio following Equations (1) - (3).
//' 
//' @param M_iter Number of iterations after burning and thinning.
//' @param n_current_ctrl  Sample size at a interim analysis point for the control arm. 
//' @param n_current_trt  Sample size at a interim analysis point for the treatment arm. 
//' @param postData  A cube with n rows, two columns and N_iter slices. Slice v, v = 1, ..., N_iter, contains n pairs of event times sampled from bivariate lognormal distribution with mean MU and variance Sigma, where
//'  Mu and Sigma are generated from the posterior normal and inverse-wishart distribution, respectively, from the v-th iteration of the gibbs sampler.
//' @param time_trt A field containing M_iter elements, each of which is a vector of length n_current_trt. These vectors store the censoring time generated from the estimated censoring distribution from the treatment arm.
//' @param time_ctrl A field containing M_iter elements, each of which is a vector of length n_current_ctrl. These vectors store the censoring time generated from the estimated censoring distribution from the control arm.
//' @param surv_trt A field containing M_iter elements, each of which is a vector of length n_current_trt. These vectors store the survival probabilities generated from the estimated censoring distribution from the treatment arm.
//' @param surv_ctrl A field containing M_iter elements, each of which is a vector of length n_current_ctrl. These vectors store the survival probabilities generated from the estimated censoring distribution from the control arm.
//' @return A list of WR. WR is a matrix with one row and M_iter columns. Each row contains the estimated win ratio for each iteration.
// [[Rcpp::export]]
List compare(int M_iter, int n_current_ctrl, int n_current_trt, const arma::cube& postData,
arma::field<arma::vec> time_trt, arma::field<arma::vec> time_ctrl, arma::field<arma::vec> surv_trt, arma::field<arma::vec> surv_ctrl) {

    arma::mat WR(1,M_iter, fill::none);
    
    
    for (int j_iter = 0; j_iter < M_iter; j_iter++) {
        arma::vec C = postData.slice(j_iter).col(2);
        arma::vec C_trt_vec = C.head(n_current_trt);
        arma::vec C_ctrl_vec = C.tail(n_current_ctrl);
        arma::vec delta1 = postData.slice(j_iter).col(4);
        arma::vec delta2 = postData.slice(j_iter).col(5);
        arma::vec R = postData.slice(j_iter).col(0);
        arma::vec D = postData.slice(j_iter).col(1);

        arma::cube compare_results(n_current_ctrl, 2, n_current_trt, fill::zeros);

        for (int i_trt = 0; i_trt < n_current_trt; i_trt++) {
            double delta2_trt = delta2(i_trt);
            double delta1_trt = delta1(i_trt);
            double D_trt = D(i_trt);
            double R_trt = R(i_trt);
            double C_trt = C(i_trt);
            
            for (int j_ctrl = 0; j_ctrl < n_current_ctrl; j_ctrl++) {
                double delta2_ctrl = delta2(n_current_trt + j_ctrl);
                double delta1_ctrl = delta1(n_current_trt + j_ctrl);
                double D_ctrl = D(n_current_trt + j_ctrl);
                double R_ctrl = R(n_current_trt + j_ctrl);
                double C_ctrl = C(n_current_trt + j_ctrl);
                
                if ((delta2_trt == 1 && delta2_ctrl == 1 && D_ctrl < D_trt) || 
                    (delta2_trt == 0 && delta2_ctrl == 1 && D_ctrl < C_trt)) {
                    
                    compare_results(j_ctrl, 0, i_trt) = 1.0 /  (calculateSurvival(time_ctrl(j_iter), surv_ctrl(j_iter), D_ctrl, C_ctrl_vec) * calculateSurvival(time_trt(j_iter), surv_trt(j_iter), D_ctrl, C_trt_vec));
                } else if ((delta2_trt == 1 && delta2_ctrl == 1 && D_trt < D_ctrl) || 
                           (delta2_trt == 1 && delta2_ctrl == 0 && D_trt < C_ctrl)) {
                    
                    compare_results(j_ctrl, 1, i_trt) = 1.0 / (calculateSurvival(time_ctrl(j_iter), surv_ctrl(j_iter), D_trt, C_ctrl_vec) * calculateSurvival(time_trt(j_iter), surv_trt(j_iter), D_trt, C_trt_vec));
                } else if ((delta1_trt == 1 && delta1_ctrl == 1 && R_ctrl < R_trt) || 
                           (delta1_trt == 0 && delta1_ctrl == 1 && R_ctrl < C_trt)) {
                    
                    compare_results(j_ctrl, 0, i_trt) = 1.0 / ( calculateSurvival(time_ctrl(j_iter), surv_ctrl(j_iter), R_ctrl, C_ctrl_vec) * calculateSurvival(time_trt(j_iter), surv_trt(j_iter), R_ctrl, C_trt_vec));
                } else if ((delta1_trt == 1 && delta1_ctrl == 1 && R_trt < R_ctrl) || 
                           (delta1_trt == 1 && delta1_ctrl == 0 && R_trt < C_ctrl)) {
                    
                    compare_results(j_ctrl, 1, i_trt) = 1.0 / (calculateSurvival(time_ctrl(j_iter), surv_ctrl(j_iter), R_trt, C_ctrl_vec) * calculateSurvival(time_trt(j_iter), surv_trt(j_iter), R_trt, C_trt_vec));
                }
            }
        }

        arma::cube compare_results_rm = compare_results.replace(datum::nan, 0);
        WR.col(j_iter) = arma::accu(compare_results_rm.col(0))/arma::accu(compare_results_rm.col(1));
    }
    
    return List::create(Named("WR") = WR);
    
}
            
     
     