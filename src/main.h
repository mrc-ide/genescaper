//#include <Rcpp.h>
# include <RcppArmadillo.h>
#include <vector>

//------------------------------------------------
// [[Rcpp::export]]
arma::field<arma::cube> predict_map_cpp(arma::mat data, Rcpp::List mcmc_sample,
                                        arma::mat dist_11, arma::mat dist_12,
                                        arma::mat dist_22, Rcpp::List params,
                                        int inner_reps, Rcpp::List args_progress,
                                        Rcpp::List args_functions, Rcpp::List args_misc);

//------------------------------------------------
void draw_sigsq_mu(double &sigsq, double &mu, double gamma_0, double phi_0,
                   double alpha_0, double beta_0, int n_site, arma::mat &K_11_inv,
                   arma::vec z, arma::vec ones_site);

//------------------------------------------------
void draw_y_post(arma::vec &ret, double sigsq, double mu, arma::mat &R_11_inv,
                 arma::mat &W_inv, arma::vec &ones_site, arma::vec &z_obs, double nu);

//------------------------------------------------
void draw_z_pred(arma::mat &ret, arma::mat &R_11_inv, arma::mat &R_12,
                 arma::mat &K_22, double mu, double sigsq, arma::vec &ones_pred,
                 arma::vec &ones_site, arma::vec &y, int inner_reps);

//------------------------------------------------
void z_to_freq(arma::cube &z_cube);

//------------------------------------------------
void draw_z_null(arma::mat &ret, arma::mat &K_22, double mu, double sigsq,
                 arma::vec &ones_pred, int inner_reps);

//------------------------------------------------
// [[Rcpp::export]]
arma::field<arma::cube> null_map_cpp(arma::mat data, Rcpp::List mcmc_sample,
                                     arma::mat dist_11, arma::mat dist_22,
                                     Rcpp::List params, int inner_reps,
                                     Rcpp::List args_progress,
                                     Rcpp::List args_functions, Rcpp::List args_misc);

//------------------------------------------------
// [[Rcpp::export]]
arma::field<arma::cube> null_site_cpp(arma::mat data, Rcpp::List mcmc_sample,
                                      arma::mat dist_11, Rcpp::List params,
                                      int inner_reps, Rcpp::List args_progress,
                                      Rcpp::List args_functions, Rcpp::List args_misc);

//------------------------------------------------
// [[Rcpp::export]]
Rcpp::List GeoMAPI_assign_edges_cpp(Rcpp::List args, Rcpp::List args_functions, Rcpp::List args_progress,
                                    arma::mat dist_11);

//------------------------------------------------
// [[Rcpp::export]]
Rcpp::List post_sigsq_mu(arma::mat data, Rcpp::List mcmc_sample,
                         arma::mat dist_11, std::vector<double> true_sigsq);

