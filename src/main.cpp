
#include "main.h"
#include "misc_v13.h"
#include "probability_v17.h"
#include "utils.h"

using namespace std;

//------------------------------------------------
// draw allele frequencies at new locations from predictive distribution
Rcpp::List predict_map_cpp(arma::mat data, Rcpp::List mcmc_sample,
                           arma::mat dist_11, arma::mat dist_12,
                           arma::mat dist_22, Rcpp::List params,
                           int inner_reps, Rcpp::List args_progress,
                           Rcpp::List args_functions, Rcpp::List args_misc) {
  
  // get functions
  Rcpp::Function update_progress = args_functions["update_progress"];
  
  // get misc parameters
  bool pb_markdown = rcpp_to_bool(args_misc["pb_markdown"]);
  
  // get basic dimensions
  int alleles = data.n_cols;
  int n_site = dist_11.n_rows;
  int n_pred = dist_22.n_rows;
  
  // get posterior MCMC draws
  vector<double> lambda_draws = rcpp_to_vector_double(mcmc_sample["lambda"]);
  vector<double> nu_draws = rcpp_to_vector_double(mcmc_sample["nu"]);
  int reps = nu_draws.size();
  
  // get fixed model parameters
  double mu_mean = rcpp_to_double(params["mu_mean"]);
  double mu_scale = rcpp_to_double(params["mu_scale"]);
  double sigsq_mean = rcpp_to_double(params["sigsq_mean"]);
  double sigsq_var = rcpp_to_double(params["sigsq_var"]);
  double dist_power = 1.0;
  
  // reparameterise for convenience
  double phi_0 = mu_mean;
  double gamma_0 = 1.0 / mu_scale;
  double alpha_0 = sigsq_mean * sigsq_mean / sigsq_var + 2.0;
  double beta_0 = sigsq_mean * (alpha_0 - 1.0);
  
  // initialise intermediate objects
  arma::vec ones_pred = arma::ones(n_pred, 1);
  arma::vec ones_site = arma::ones(n_site, 1);
  arma::cube z_cube(n_pred, inner_reps, alleles + 1);
  
  // initialise return object
  arma::field<arma::cube> ret(reps);
  
  // loop through MCMC draws
  for (int rep_i = 0; rep_i < reps; ++rep_i) {
    double nu = nu_draws[rep_i];
    double inv_lambda = 1.0 / lambda_draws[rep_i];
    
    // allow user to exit on escape
    Rcpp::checkUserInterrupt();
    
    // update progress bar
    update_progress_cpp(args_progress, update_progress, "pb", rep_i, reps, !pb_markdown);
    
    // initialise kernel matrices
    arma::mat K_11 = nu * exp(-pow(dist_11 * inv_lambda, dist_power));
    arma::mat K_11_inv = arma::inv_sympd(K_11);
    arma::mat K_12 = nu * exp(-pow(dist_12 * inv_lambda, dist_power));
    arma::mat K_22 = nu * exp(-pow(dist_22 * inv_lambda, dist_power));
    arma::mat W(n_site, n_site, arma::fill::eye);
    W = W * (1 - nu);
    arma::mat R_11 = K_11 + W;
    arma::mat R_11_inv = arma::inv_sympd(R_11);
    
    // loop through alleles
    for (int allele_i = 0; allele_i < alleles; ++allele_i) {
      
      // get data for this allele
      arma::vec z = data.col(allele_i);
      
      // calculate posterior parameters for drawing mu and sigsq
      double gamma_1 = gamma_0 + arma::accu(R_11_inv);
      double phi_1 = (gamma_0 * phi_0 + arma::as_scalar(z.t() * R_11_inv * ones_site)) / gamma_1;
      double alpha_1 = alpha_0 + 0.5 * n_site;
      double beta_1 = beta_0 + 0.5*(gamma_0 * sq(phi_0) - gamma_1 * sq(phi_1) +  arma::as_scalar(z.t() * R_11_inv * z));
      
      // draw mu and sigsq
      double sigsq = 1.0 / rgamma1(alpha_1, beta_1);
      double mu = rnorm1(phi_1, sigsq / gamma_1);
      /*
      print("\nParams:");
      print("gamma", gamma_0, gamma_1);
      print("phi", phi_0, phi_1);
      print("alpha", alpha_0, alpha_1);
      print("beta", beta_0, beta_1);
      print("mu", mu);
      print("sigsq", sigsq);
      */
      // calculate Sigma_1 and mu_1
      arma::mat Sigma_1 = sigsq* arma::inv_sympd(K_11_inv + arma::inv(W));
      arma::vec mu_1 = Sigma_1 * (mu / sigsq * K_11_inv * ones_site + z / (sigsq*(1 - nu)) );
      
      // draw y from posterior
      arma::vec y = arma::mvnrnd(mu_1, Sigma_1);
      
      // get predictive mean and covariance. Going via the Cholesky
      // decomposition like this appears to make the result more stable in terms
      // of the final covariance mtarix being symmetric positive-definite
      arma::mat A = arma::chol(K_11_inv);
      arma::mat B = K_12.t() * A.t();
      arma::mat pred_covar = K_22 - B * B.t();
      arma::vec pred_mean = mu * ones_pred + K_12.t() * K_11_inv * (y - mu * ones_site);
      
      // draw multiple times from predictive distribution
      z_cube.slice(allele_i) = arma::mvnrnd(pred_mean, pred_covar, inner_reps);
      
    }  // end of allele_i loop
    
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
    
    // store results of this rep
    ret(rep_i) = z_cube;
    
  }  // end of rep_i loop
  
  // return list
  return Rcpp::List::create(Rcpp::Named("ret") = ret);
}

//------------------------------------------------
// calculate average Gst over a list of loci
arma::mat get_mean_pairwise_Gst_cpp(arma::field<arma::mat> freq_list) {
  
  // get basic dimensions
  int loci = freq_list.n_elem;
  int n_sites = freq_list(0).n_rows;
  
  // get mean pairwise Gst over all loci
  arma::mat mean_Gst(n_sites, n_sites);
  for (int i = 0; i < loci; ++i) {
    mean_Gst += get_pairwise_Gst_cpp(freq_list(i));
  }
  mean_Gst /= loci;
  
  return mean_Gst;
}

//------------------------------------------------
// calculate pairwise Gst given a matrix over demes and alleles
arma::mat get_pairwise_Gst_cpp(arma::mat p) {
  
  int n_demes = p.n_rows;
  arma::mat Gst(n_demes, n_demes);
  arma::vec J = arma::sum(p % p, 1);
  for (int i = 0; i < (n_demes - 1); ++i) {
    for (int j = (i + 1); j < n_demes; ++j) {
      arma::rowvec p_mean = 0.5 * (p.row(i) + p.row(j));
      double J_t = arma::sum(p_mean % p_mean);
      double J_s = 0.5 * (J(i) + J(j));
      Gst(i, j) = (J_s - J_t) / (1.0 - J_t);
    }
  }
  return Gst;
}

//------------------------------------------------
// draw pairwise Gst from GRF model
Rcpp::List predict_pairwise_Gst_cpp(arma::field<arma::mat> data_list,
                                    Rcpp::List mcmc_sample, arma::mat dist_11,
                                    Rcpp::List params, int inner_reps,
                                    Rcpp::List args_progress,
                                    Rcpp::List args_functions, Rcpp::List args_misc) {
  
  
  
  // get functions
  Rcpp::Function update_progress = args_functions["update_progress"];
  
  // get misc parameters
  bool silent = rcpp_to_bool(args_misc["silent"]);
  bool pb_markdown = rcpp_to_bool(args_misc["pb_markdown"]);
  
  // get basic dimensions
  int loci = data_list.n_elem;
  int n_site = data_list(0).n_rows;
  
  // get posterior MCMC draws
  vector<double> lambda_draws = rcpp_to_vector_double(mcmc_sample["lambda"]);
  vector<double> nu_draws = rcpp_to_vector_double(mcmc_sample["nu"]);
  int reps = nu_draws.size();
  
  // get fixed model parameters
  double mu_mean = rcpp_to_double(params["mu_mean"]);
  double mu_scale = rcpp_to_double(params["mu_scale"]);
  double sigsq_mean = rcpp_to_double(params["sigsq_mean"]);
  double sigsq_var = rcpp_to_double(params["sigsq_var"]);
  double dist_power = 1.0;
  
  // reparameterise for convenience
  double phi_0 = mu_mean;
  double gamma_0 = 1.0 / mu_scale;
  double alpha_0 = sigsq_mean * sigsq_mean / sigsq_var + 2.0;
  double beta_0 = sigsq_mean * (alpha_0 - 1.0);
  
  // initialise intermediate objects
  arma::vec ones_site = arma::ones(n_site, 1);
  
  // initialise return object
  arma::field<arma::cube> ret(reps);
  
  if (!silent) {
    print("Simulating from null model");
  }
  
  // loop through MCMC draws
  for (int rep_i = 0; rep_i < reps; ++rep_i) {
    double nu = nu_draws[rep_i];
    double inv_lambda = 1.0 / lambda_draws[rep_i];
    
    // allow user to exit on escape
    Rcpp::checkUserInterrupt();
    
    // update progress bar
    update_progress_cpp(args_progress, update_progress, "pb", rep_i, reps, !pb_markdown);
    
    // initialise kernel matrices
    arma::mat K_11 = nu * exp(-pow(dist_11 * inv_lambda, dist_power));
    arma::mat K_11_inv = arma::inv_sympd(K_11);
    arma::mat W(n_site, n_site, arma::fill::eye);
    W = W * (1 - nu);
    arma::mat R_11 = K_11 + W;
    arma::mat R_11_inv = arma::inv_sympd(R_11);
    
    // initialise return object for this rep
    arma::cube Gst_cube(n_site, n_site, inner_reps);
    
    // loop through loci
    for (int locus_i = 0; locus_i < loci; ++locus_i) {
      int alleles = data_list(locus_i).n_cols;
      
      // initialise intermediate objects
      arma::cube z_cube(n_site, inner_reps, alleles + 1);
      
      // loop through alleles
      for (int allele_i = 0; allele_i < alleles; ++allele_i) {
        
        // get data for this allele
        arma::vec z = data_list(locus_i).col(allele_i);
        
        // calculate posterior parameters for drawing mu and sigsq
        double gamma_1 = gamma_0 + arma::accu(R_11_inv);
        double phi_1 = (gamma_0 * phi_0 + arma::as_scalar(z.t() * R_11_inv * ones_site)) / gamma_1;
        double alpha_1 = alpha_0 + 0.5 * n_site;
        double beta_1 = beta_0 + 0.5*(gamma_0 * sq(phi_0) - gamma_1 * sq(phi_1) +  arma::as_scalar(z.t() * R_11_inv * z));
        /*
        print("\nParams:");
        print("gamma", gamma_0, gamma_1);
        print("phi", phi_0, phi_1);
        print("alpha", alpha_0, alpha_1);
        print("beta", beta_0, beta_1);
        */
        // draw mu and sigsq
        double sigsq = 1.0 / rgamma1(alpha_1, beta_1);
        double mu = rnorm1(phi_1, sigsq / gamma_1);
        
        // draw new data from GRF unconditional on data
        z_cube.slice(allele_i) = arma::mvnrnd(mu * ones_site, sigsq * R_11, inner_reps);
        
        // alternatively, draw from GRF conditional on data
        //arma::mat Sigma_1 = sigsq* arma::inv(K_11_inv + arma::inv(W));
        //arma::vec mu_1 = Sigma_1 * (mu / sigsq * K_11_inv * ones_site + z / (sigsq*(1 - nu)) );
        //z_cube.slice(allele_i) = arma::mvnrnd(mu_1, Sigma_1 + sigsq * W, inner_reps);
        
      }  // end of allele_i loop
      
      // logistic transform z values back to [0,1] range
      z_cube = 1.0 / (1.0 + exp(-z_cube));
      
      // convert transformed z vales back to allele frequencies via stick-breaking
      // approach
      arma::mat stick_used(n_site, inner_reps);
      for (int j = 0; j < alleles; ++j) {
        z_cube.slice(j) %= 1.0 - stick_used;
        stick_used += z_cube.slice(j);
      }
      z_cube.slice(alleles) = 1.0 - stick_used;
      
      // get Gst for all inner reps
      for (int i = 0; i < inner_reps; ++i) {
        Gst_cube.slice(i) += get_pairwise_Gst_cpp(z_cube.col(i));
      }
      
    }  // end of locus_i loop
    
    // store results of this rep
    ret(rep_i) = Gst_cube / loci;
    
  }  // end of rep_i loop
  
  // return list
  return Rcpp::List::create(Rcpp::Named("ret") = ret);
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
  
  // loop through hexes
  for (int hex = 0; hex < n_hex; ++hex) {
    
    // allow user to exit on escape
    Rcpp::checkUserInterrupt();
    
    // report progress
    if (!silent) {
      update_progress(args_progress, "pb", hex, n_hex);
    }
    
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
  }
  
  // return list
  return Rcpp::List::create(Rcpp::Named("edge_assignment") = hex_edges);
}

