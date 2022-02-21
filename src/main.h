#include <Rcpp.h>
#include <vector>

//------------------------------------------------
// [[Rcpp::export]]
Rcpp::List predict_map_cpp(Rcpp::List data_list, Rcpp::List mcmc_sample,
                           Rcpp::NumericMatrix dist_11, Rcpp::NumericMatrix dist_12,
                           Rcpp::NumericMatrix dist_22, Rcpp::List params,
                           int inner_reps, Rcpp::List args_progress,
                           Rcpp::List args_functions, Rcpp::List args_misc);
