
#include <Rcpp.h>

//------------------------------------------------
// Simulate from simple Wright-Fisher model
// [[Rcpp::export]]
Rcpp::List sim_wrightfisher_cpp(Rcpp::List args, Rcpp::List args_functions,
                                Rcpp::List args_progress);
