//#include <Rcpp.h>
# include <RcppArmadillo.h>
#include <RcppEigen.h>
#include <vector>

//------------------------------------------------
// [[Rcpp::export]]
Rcpp::List predict_map_cpp(Eigen::MatrixXd data, Rcpp::List mcmc_sample,
                           Eigen::MatrixXd dist_11, Eigen::MatrixXd dist_12,
                           Eigen::MatrixXd dist_22, Rcpp::List params,
                           int inner_reps, Rcpp::List args_progress,
                           Rcpp::List args_functions, Rcpp::List args_misc);

//------------------------------------------------
// [[Rcpp::export]]
Rcpp::List wrangle_pairwise_Gst_cpp(Rcpp::List freq_array);

//------------------------------------------------
std::vector<double> calculate_pairwise_Gst(const std::vector<std::vector<std::vector<double>>> &p);
