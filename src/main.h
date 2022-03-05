//#include <Rcpp.h>
# include <RcppArmadillo.h>
#include <vector>

//------------------------------------------------
// [[Rcpp::export]]
Rcpp::List predict_map_cpp(arma::mat data, Rcpp::List mcmc_sample,
                           arma::mat dist_11, arma::mat dist_12,
                           arma::mat dist_22, Rcpp::List params,
                           int inner_reps, Rcpp::List args_progress,
                           Rcpp::List args_functions, Rcpp::List args_misc);

//------------------------------------------------
// [[Rcpp::export]]
arma::mat get_mean_pairwise_Gst_cpp(arma::field<arma::mat> freq_list);

//------------------------------------------------
arma::mat get_pairwise_Gst_cpp(arma::mat p);

//------------------------------------------------
// [[Rcpp::export]]
Rcpp::List predict_pairwise_Gst_cpp(arma::field<arma::mat> data_list,
                                    Rcpp::List mcmc_sample, arma::mat dist_11,
                                    Rcpp::List params, int inner_reps,
                                    Rcpp::List args_progress,
                                    Rcpp::List args_functions, Rcpp::List args_misc);

//------------------------------------------------
// [[Rcpp::export]]
Rcpp::List GeoMAPI_assign_edges_cpp(Rcpp::List args, Rcpp::List args_functions, Rcpp::List args_progress,
                                    arma::mat dist_11);

