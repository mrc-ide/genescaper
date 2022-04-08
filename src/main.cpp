
#include "main.h"
#include "misc_v14.h"
#include "probability_v17.h"
#include "utils.h"

using namespace std;

//------------------------------------------------
// draw allele frequencies at new locations from predictive distribution
arma::field<arma::cube> predict_map_cpp(arma::mat data, Rcpp::List mcmc_sample,
                                        arma::mat dist_11, arma::mat dist_12,
                                        arma::mat dist_22, Rcpp::List params,
                                        int inner_reps, Rcpp::List args_progress,
                                        Rcpp::List args_functions, Rcpp::List args_misc) {
  
  // get functions
  Rcpp::Function update_progress = args_functions["update_progress"];
  
  // get misc parameters
  bool silent = rcpp_to_bool(args_misc["silent"]);
  bool pb_markdown = rcpp_to_bool(args_misc["pb_markdown"]);
  
  // get basic dimensions
  int alleles = data.n_cols;
  int n_site = dist_11.n_rows;
  int n_pred = dist_22.n_rows;
  
  // get posterior MCMC draws
  vector<double> nu_draws = rcpp_to_vector_double(mcmc_sample["nu"]);
  vector<double> lambda_draws = rcpp_to_vector_double(mcmc_sample["lambda"]);
  vector<double> omega_draws = rcpp_to_vector_double(mcmc_sample["omega"]);
  vector<double> gamma_draws = rcpp_to_vector_double(mcmc_sample["gamma"]);
  int reps = nu_draws.size();
  
  // define fixed model parameters
  double phi_0 = 0.0;
  double alpha_0 = 0.01;
  double beta_0 = 0.01;
  
  // initialise intermediate objects
  arma::vec ones_pred = arma::ones(n_pred, 1);
  arma::vec ones_site = arma::ones(n_site, 1);
  arma::vec y(n_pred);
  arma::cube z_cube(n_pred, inner_reps, alleles + 1);
  
  // initialise return object
  arma::field<arma::cube> ret(reps);
  
  // initialise progress bar
  if (!silent) {
    update_progress_cpp(args_progress, update_progress, "pb", 0, reps, !pb_markdown);
  }
  
  // loop through MCMC draws
  for (int rep_i = 0; rep_i < reps; ++rep_i) {
    double nu = nu_draws[rep_i];
    double inv_lambda = 1.0 / lambda_draws[rep_i];
    double omega = omega_draws[rep_i];
    double gamma = gamma_draws[rep_i];
    
    // allow user to exit on escape
    Rcpp::checkUserInterrupt();
    
    // initialise kernel matrices
    arma::mat R_11 = nu * exp(-pow(dist_11 * inv_lambda, omega));
    arma::mat R_11_inv = arma::inv_sympd(R_11);
    arma::mat R_12 = nu * exp(-pow(dist_12 * inv_lambda, omega));
    arma::mat R_22 = nu * exp(-pow(dist_22 * inv_lambda, omega));
    arma::mat W_11(n_site, n_site, arma::fill::eye);
    W_11 = W_11 * (1 - nu);
    arma::mat W_11_inv = arma::inv(W_11);
    arma::mat W_22(n_pred, n_pred, arma::fill::eye);
    W_22 = W_22 * (1 - nu);
    arma::mat K_11 = R_11 + W_11;
    arma::mat K_11_inv = arma::inv_sympd(K_11);
    arma::mat K_22 = R_22 + W_22;
    
    // loop through alleles
    for (int allele_i = 0; allele_i < alleles; ++allele_i) {
      
      // get data for this allele
      arma::vec z_obs = data.col(allele_i);
      
      // draw mu and sigsq from posterior
      double sigsq, mu;
      draw_sigsq_mu(sigsq, mu, gamma, phi_0, alpha_0, beta_0, n_site, K_11_inv,
                    z_obs, ones_site);
      
      // draw y from posterior
      draw_y_post(y, sigsq, mu, R_11_inv, W_11_inv, ones_site, z_obs, nu);
      
      // draw z from posterior predictive given y
      draw_z_pred(z_cube.slice(allele_i), R_11_inv, R_12, K_22, mu, sigsq,
                  ones_pred, ones_site, y, inner_reps);
      
    }  // end of allele_i loop
    
    // convert transformed z vales back to allele frequencies via stick-breaking
    // approach
    z_to_freq(z_cube);
    
    // store results of this rep
    ret(rep_i) = z_cube;
    
    // update progress bar
    if (!silent) {
      update_progress_cpp(args_progress, update_progress, "pb", rep_i + 1, reps, !pb_markdown);
    }
    
  }  // end of rep_i loop
  
  // return
  return ret;
}

//------------------------------------------------
// draw mu and sigsq (passed in by reference) from posterior
void draw_sigsq_mu(double &sigsq, double &mu, double gamma_0, double phi_0,
                   double alpha_0, double beta_0, int n_site, arma::mat &K_11_inv,
                   arma::vec z, arma::vec ones_site) {
  
  // calculate posterior parameters
  double alpha_1 = alpha_0 + 0.5 * (double)n_site;
  double gamma_1 = gamma_0 + arma::accu(K_11_inv);
  double phi_1 = (gamma_0 * phi_0 + arma::as_scalar(z.t() * K_11_inv * ones_site)) / gamma_1;
  double beta_1 = beta_0 + 0.5*(gamma_0 * sq(phi_0) - gamma_1 * sq(phi_1) + arma::as_scalar(z.t() * K_11_inv * z));
  
  // draw mu and sigsq
  sigsq = 1.0 / rgamma1(alpha_1, beta_1);
  mu = rnorm1(phi_1, pow(sigsq / gamma_1, 0.5));
}

//------------------------------------------------
// draw y from posterior distribution. Results vector passed in by reference
void draw_y_post(arma::vec &ret, double sigsq, double mu, arma::mat &R_11_inv, 
                 arma::mat &W_inv, arma::vec &ones_site, arma::vec &z_obs, double nu) {
  
  // calculate Sigma_1 and mu_1
  arma::mat Sigma_1 = sigsq* arma::inv_sympd(R_11_inv + W_inv);
  arma::vec mu_1 = Sigma_1 * (mu / sigsq * R_11_inv * ones_site + z_obs / (sigsq * (1 - nu)) );
  
  // draw y from posterior
  ret = arma::mvnrnd(mu_1, Sigma_1);
}

//------------------------------------------------
// draw z from predictive distribution given y. Results matrix passed in by reference
void draw_z_pred(arma::mat &ret, arma::mat &R_11_inv, arma::mat &R_12,
                 arma::mat &K_22, double mu, double sigsq, arma::vec &ones_pred,
                 arma::vec &ones_site, arma::vec &y, int inner_reps) {
  
  // get predictive mean and covariance. Going via the Cholesky
  // decomposition like this appears to make the result more stable in terms
  // of the final covariance mtarix being symmetric positive-definite
  arma::mat A = arma::chol(R_11_inv);
  arma::mat B = R_12.t() * A.t();
  arma::mat pred_covar = sigsq * (K_22 - B * B.t());
  arma::vec pred_mean = mu * ones_pred + R_12.t() * R_11_inv * (y - mu * ones_site);
  
  // draw multiple times from predictive distribution
  ret = arma::mvnrnd(pred_mean, pred_covar, inner_reps);
}

//------------------------------------------------
// back-transform z-values to allele frequencies. Values are overwritten in
// z_cube object, passed by reference
void z_to_freq(arma::cube &z_cube) {
  
  // get basic dimensions
  int n_pred = z_cube.n_rows;
  int inner_reps = z_cube.n_cols;
  int alleles = z_cube.n_slices - 1;
  
  // logistic transform z values back to [0,1] range
  z_cube = 1.0 / (1.0 + exp(-z_cube));
  
  // convert transformed z vales back to allele frequencies via stick-breaking
  // approach
  arma::mat stick_used(n_pred, inner_reps);
  for (int j = 0; j < alleles; ++j) {
    z_cube.slice(j) %= 1.0 - stick_used;
    stick_used += z_cube.slice(j);
  }
  z_cube.slice(alleles) = 1.0 - stick_used;
  
}

//------------------------------------------------
// draw allele frequencies at new locations from null distribution
arma::field<arma::cube> null_map_cpp(arma::mat data, Rcpp::List mcmc_sample,
                                     arma::mat dist_11, arma::mat dist_22,
                                     Rcpp::List params, int inner_reps,
                                     Rcpp::List args_progress,
                                     Rcpp::List args_functions, Rcpp::List args_misc) {
  
  // get functions
  Rcpp::Function update_progress = args_functions["update_progress"];
  
  // get misc parameters
  bool silent = rcpp_to_bool(args_misc["silent"]);
  bool pb_markdown = rcpp_to_bool(args_misc["pb_markdown"]);
  
  // get basic dimensions
  int alleles = data.n_cols;
  int n_site = dist_11.n_rows;
  int n_pred = dist_22.n_rows;
  
  // get posterior MCMC draws
  vector<double> omega_draws = rcpp_to_vector_double(mcmc_sample["omega"]);
  vector<double> lambda_draws = rcpp_to_vector_double(mcmc_sample["lambda"]);
  vector<double> nu_draws = rcpp_to_vector_double(mcmc_sample["nu"]);
  vector<double> gamma_draws = rcpp_to_vector_double(mcmc_sample["gamma"]);
  int reps = nu_draws.size();
  
  // define fixed model parameters
  double phi_0 = 0;
  double alpha_0 = 0.01;
  double beta_0 = 0.01;
  
  // initialise intermediate objects
  arma::vec ones_pred = arma::ones(n_pred, 1);
  arma::vec ones_site = arma::ones(n_site, 1);
  arma::cube z_cube(n_pred, inner_reps, alleles + 1);
  
  // initialise return object
  arma::field<arma::cube> ret(reps);
  
  // initialise progress bar
  if (!silent) {
    update_progress_cpp(args_progress, update_progress, "pb", 0, reps, !pb_markdown);
  }
  
  // loop through MCMC draws
  for (int rep_i = 0; rep_i < reps; ++rep_i) {
    double nu = nu_draws[rep_i];
    double inv_lambda = 1.0 / lambda_draws[rep_i];
    double omega = omega_draws[rep_i];
    double gamma = gamma_draws[rep_i];
    
    // allow user to exit on escape
    Rcpp::checkUserInterrupt();
    
    // initialise kernel matrices
    arma::mat R_11 = nu * exp(-pow(dist_11 * inv_lambda, omega));
    arma::mat R_22 = nu * exp(-pow(dist_22 * inv_lambda, omega));
    arma::mat W_11(n_site, n_site, arma::fill::eye);
    W_11 = W_11 * (1 - nu);
    arma::mat W_22(n_pred, n_pred, arma::fill::eye);
    W_22 = W_22 * (1 - nu);
    arma::mat K_11 = R_11 + W_11;
    arma::mat K_11_inv = arma::inv_sympd(K_11);
    arma::mat K_22 = R_22 + W_22;
    
    // loop through alleles
    for (int allele_i = 0; allele_i < alleles; ++allele_i) {
      
      // get data for this allele
      arma::vec z_obs = data.col(allele_i);
      
      // draw mu and sigsq from posterior
      double sigsq, mu;
      draw_sigsq_mu(sigsq, mu, gamma, phi_0, alpha_0, beta_0, n_site, K_11_inv,
                    z_obs, ones_site);
      
      // draw z from null distribution
      draw_z_null(z_cube.slice(allele_i), K_22, mu, sigsq, ones_pred, inner_reps);
      
    }  // end of allele_i loop
    
    // convert transformed z vales back to allele frequencies via stick-breaking
    // approach
    z_to_freq(z_cube);
    
    // store results of this rep
    ret(rep_i) = z_cube;
    
    // update progress bar
    if (!silent) {
      update_progress_cpp(args_progress, update_progress, "pb", rep_i + 1, reps, !pb_markdown);
    }
    
  }  // end of rep_i loop
  
  // return 
  return ret;
}

//------------------------------------------------
// draw z from null distribution. Results matrix passed in by reference
void draw_z_null(arma::mat &ret, arma::mat &K_22, double mu, double sigsq,
                 arma::vec &ones_pred, int inner_reps) {
  
  // draw multiple times from null distribution
  ret = arma::mvnrnd(mu * ones_pred, sigsq * K_22, inner_reps);
  
}

//------------------------------------------------
// draw allele frequencies at site locations from null distribution
arma::field<arma::cube> null_site_cpp(arma::mat data, Rcpp::List mcmc_sample,
                                      arma::mat dist_11, Rcpp::List params,
                                      int inner_reps, Rcpp::List args_progress,
                                      Rcpp::List args_functions, Rcpp::List args_misc) {
  
  // get functions
  Rcpp::Function update_progress = args_functions["update_progress"];
  
  // get misc parameters
  bool silent = rcpp_to_bool(args_misc["silent"]);
  bool pb_markdown = rcpp_to_bool(args_misc["pb_markdown"]);
  
  // get basic dimensions
  int alleles = data.n_cols;
  int n_site = dist_11.n_rows;
  
  // get posterior MCMC draws
  vector<double> omega_draws = rcpp_to_vector_double(mcmc_sample["omega"]);
  vector<double> lambda_draws = rcpp_to_vector_double(mcmc_sample["lambda"]);
  vector<double> nu_draws = rcpp_to_vector_double(mcmc_sample["nu"]);
  vector<double> gamma_draws = rcpp_to_vector_double(mcmc_sample["gamma"]);
  int reps = nu_draws.size();
  
  // define fixed model parameters
  double phi_0 = 0;
  double alpha_0 = 0.01;
  double beta_0 = 0.01;
  
  // initialise intermediate objects
  arma::vec ones_site = arma::ones(n_site, 1);
  arma::cube z_cube(n_site, inner_reps, alleles + 1);
  
  // initialise return object
  arma::field<arma::cube> ret(reps);
  
  // initialise progress bar
  if (!silent) {
    update_progress_cpp(args_progress, update_progress, "pb", 0, reps, !pb_markdown);
  }
  
  // loop through MCMC draws
  for (int rep_i = 0; rep_i < reps; ++rep_i) {
    double nu = nu_draws[rep_i];
    double inv_lambda = 1.0 / lambda_draws[rep_i];
    double omega = omega_draws[rep_i];
    double gamma = gamma_draws[rep_i];
    
    // allow user to exit on escape
    Rcpp::checkUserInterrupt();
    
    // initialise kernel matrices
    arma::mat R_11 = nu * exp(-pow(dist_11 * inv_lambda, omega));
    arma::mat W_11(n_site, n_site, arma::fill::eye);
    W_11 = W_11 * (1 - nu);
    arma::mat K_11 = R_11 + W_11;
    arma::mat K_11_inv = arma::inv_sympd(K_11);
    
    // loop through alleles
    for (int allele_i = 0; allele_i < alleles; ++allele_i) {
      
      // get data for this allele
      arma::vec z_obs = data.col(allele_i);
      
      // draw mu and sigsq from posterior
      double sigsq, mu;
      draw_sigsq_mu(sigsq, mu, gamma, phi_0, alpha_0, beta_0, n_site, K_11_inv,
                    z_obs, ones_site);
      
      // draw z from null distribution
      draw_z_null(z_cube.slice(allele_i), K_11, mu, sigsq, ones_site, inner_reps);
      
    }  // end of allele_i loop
    
    // convert transformed z vales back to allele frequencies via stick-breaking
    // approach
    z_to_freq(z_cube);
    
    // store results of this rep
    ret(rep_i) = z_cube;
    
    // update progress bar
    if (!silent) {
      update_progress_cpp(args_progress, update_progress, "pb", rep_i + 1, reps, !pb_markdown);
    }
    
  }  // end of rep_i loop
  
  // return 
  return ret;
}

//------------------------------------------------
// assign edges to grid cells based on intersection
Rcpp::List GeoMAPI_assign_edges_cpp(Rcpp::List args, Rcpp::List args_functions, Rcpp::List args_progress,
                                    arma::mat dist_11) {
  
  // load data and parameters
  vector<double> node_lon = rcpp_to_vector_double(args["node_lon"]);
  vector<double> node_lat = rcpp_to_vector_double(args["node_lat"]);
  vector<double> hex_lon = rcpp_to_vector_double(args["centroid_lon"]);
  vector<double> hex_lat = rcpp_to_vector_double(args["centroid_lat"]);
  double hex_width = rcpp_to_double(args["width"]);
  double eccentricity = rcpp_to_double(args["eccentricity"]);
  double max_dist = rcpp_to_double(args["max_dist"]);
  bool silent = rcpp_to_bool(args["silent"]);
  bool pb_markdown = rcpp_to_bool(args["pb_markdown"]);
  Rcpp::Function update_progress = args_functions["update_progress"];
  
  if (!silent) {
    print("Assigning edges to grid");
  }
  
  // get basic properties
  int n_node = node_lon.size();
  int n_hex = hex_lon.size();
  
  // store list of which edges intersect each hex
  vector<vector<int>> hex_edges(n_hex);
  
  // initialise progress bar
  if (!silent) {
    update_progress(args_progress, "pb", 0, n_hex, pb_markdown);
  }
  
  // loop through hexes
  for (int hex = 0; hex < n_hex; ++hex) {
    
    // allow user to exit on escape
    Rcpp::checkUserInterrupt();
    
    // loop through pairwise nodes
    int i = 0;
    for (int node1 = 0; node1 < (n_node - 1); ++node1) {
      for (int node2 = (node1 + 1); node2 < n_node; ++node2) {
        i++;
        
        // skip if longer than maximum allowed distance
        if (dist_11(node1, node2) > max_dist) {
          continue;
        }
        
        // determine whether ellipse intersects this hex
        bool intersects = collision_test_hex_ellipse(hex_lon[hex], hex_lat[hex], hex_width,
                                                     node_lon[node1], node_lat[node1], node_lon[node2],
                                                     node_lat[node2], eccentricity);
        
        // push back edge index if intersects
        if (intersects) {
          hex_edges[hex].push_back(i);
        }
      }
    }
    
    // update progress
    if (!silent) {
      update_progress(args_progress, "pb", hex, n_hex, pb_markdown);
    }
  }
  
  // return list
  return Rcpp::List::create(Rcpp::Named("edge_assignment") = hex_edges);
}

//------------------------------------------------
// draw sigma squared and mu from conditional posterior
Rcpp::List post_sigsq_mu(arma::mat data, Rcpp::List mcmc_sample,
                         arma::mat dist_11, std::vector<double> true_sigsq) {
  
  // get basic dimensions
  int alleles = data.n_cols;
  int n_site = dist_11.n_rows;
  
  // get posterior MCMC draws
  vector<double> nu_draws = rcpp_to_vector_double(mcmc_sample["nu"]);
  vector<double> lambda_draws = rcpp_to_vector_double(mcmc_sample["lambda"]);
  vector<double> omega_draws = rcpp_to_vector_double(mcmc_sample["omega"]);
  vector<double> gamma_draws = rcpp_to_vector_double(mcmc_sample["gamma"]);
  int reps = nu_draws.size();
  
  // define fixed model parameters
  double alpha_0 = 0.0;
  double beta_0 = 0.0;
  double phi_0 = 0.0;
  
  // initialise intermediate objects
  arma::vec ones_site = arma::ones(n_site, 1);
  
  // initialise return objects
  std::vector<std::vector<double>> sigsq(alleles, std::vector<double>(reps));
  std::vector<std::vector<double>> mu(alleles, std::vector<double>(reps));
  
  // loop through MCMC draws
  for (int rep_i = 0; rep_i < reps; ++rep_i) {
    double nu = nu_draws[rep_i];
    double inv_lambda = 1.0 / lambda_draws[rep_i];
    double omega = omega_draws[rep_i];
    double gamma = gamma_draws[rep_i];
    
    // TODO - remove
    gamma = 0.0;
    
    // allow user to exit on escape
    Rcpp::checkUserInterrupt();
    
    // initialise kernel matrices
    arma::mat R_11 = nu * exp(-pow(dist_11 * inv_lambda, omega));
    arma::mat W_11(n_site, n_site, arma::fill::eye);
    W_11 = W_11 * (1 - nu);
    arma::mat K_11 = R_11 + W_11;
    arma::mat K_11_inv = arma::inv_sympd(K_11);
    
    // loop through alleles
    for (int allele_i = 0; allele_i < alleles; ++allele_i) {
      
      // get data for this allele
      arma::vec z_obs = data.col(allele_i);
      
      // TODO - remove
      //double v = 0.01;
      //alpha_0 = true_sigsq[allele_i]*true_sigsq[allele_i] / v + 2.0;
      //beta_0 = true_sigsq[allele_i] * (alpha_0 - 1.0);
      //print(alpha_0, beta_0);
      
      // draw mu and sigsq from posterior
      draw_sigsq_mu(sigsq[allele_i][rep_i], mu[allele_i][rep_i],
                    gamma, phi_0, alpha_0, beta_0, n_site, K_11_inv, z_obs, ones_site);
      
    }  // end of allele_i loop
    
  }  // end of rep_i loop
  
  // return
  return Rcpp::List::create(Rcpp::Named("mu") = mu,
                            Rcpp::Named("sigsq") = sigsq);
}
