
#include "main.h"
#include "misc_v13.h"
#include "probability_v17.h"
#include "utils.h"

using namespace std;

//------------------------------------------------
Rcpp::List predict_map_cpp(Eigen::MatrixXd data, Rcpp::List mcmc_sample,
                           Eigen::MatrixXd dist_11, Eigen::MatrixXd dist_12,
                           Eigen::MatrixXd dist_22, Rcpp::List params,
                           int inner_reps, Rcpp::List args_progress,
                           Rcpp::List args_functions, Rcpp::List args_misc) {
  
  // get functions
  Rcpp::Function update_progress = args_functions["update_progress"];
  
  // get misc parameters
  bool pb_markdown = rcpp_to_bool(args_misc["pb_markdown"]);
  
  // get basic dimensions
  int alleles = data.rows();
  int n_site = dist_11.rows();
  int n_pred = dist_22.rows();
  
  // get posterior MCMC draws
  vector<double> lambda_draws = rcpp_to_vector_double(mcmc_sample["lambda"]);
  vector<double> nu_draws = rcpp_to_vector_double(mcmc_sample["nu"]);
  int reps = nu_draws.size();
  
  // get model parameters
  double mu_mean = rcpp_to_double(params["mu_mean"]);
  double mu_scale = rcpp_to_double(params["mu_scale"]);
  double sigsq_mean = rcpp_to_double(params["sigsq_mean"]);
  double sigsq_var = rcpp_to_double(params["sigsq_var"]);
  
  // reparameterise for convenience
  double phi_0 = mu_mean;
  double gamma_0 = 1.0 / mu_scale;
  double alpha_0 = sigsq_mean * sigsq_mean / sigsq_var + 2.0;
  double beta_0 = sigsq_mean * (alpha_0 - 1.0);
  
  // initialise output objects
  vector<vector<vector<double>>> z_pred(inner_reps, vector<vector<double>>(alleles));
  vector<vector<vector<double>>> ret(reps * inner_reps, vector<vector<double>>(alleles + 1, vector<double>(n_pred)));
  
  // dummy objects
  vector<double> mu_1_std(n_site);
  vector<vector<double>> Sigma_1_std(n_site, vector<double>(n_site));
  vector<vector<double>> Sigma_1_chol(n_site, vector<double>(n_site));
  vector<double> y_std(n_site);
  Eigen::VectorXd y(n_site);
  vector<double> pred_mean_std(n_pred);
  vector<vector<double>> pred_covar_std(n_pred, vector<double>(n_pred));
  vector<vector<double>> pred_covar_chol(n_pred, vector<double>(n_pred));
  
  Eigen::VectorXd ones_site(n_site);
  for (int i = 0; i < ones_site.rows(); ++i) {
    ones_site(i) = 1.0;
  }
  Eigen::VectorXd ones_pred(n_pred);
  for (int i = 0; i < ones_pred.rows(); ++i) {
    ones_pred(i) = 1.0;
  }
  
  arma::vec ones_pred_arma = arma::ones(n_pred, 1);
  arma::vec ones_site_arma = arma::ones(n_site, 1);
  
  // loop through MCMC draws
  for (int rep_i = 0; rep_i < reps; ++rep_i) {
    double nu = nu_draws[rep_i];
    double inv_lambda = 1.0 / lambda_draws[rep_i];
    
    // allow user to exit on escape
    Rcpp::checkUserInterrupt();
    
    // update progress bar
    update_progress_cpp(args_progress, update_progress, "pb", rep_i, reps, !pb_markdown);
    
    // initialise kernel matrices
    Eigen::MatrixXd K_11 = nu * exp(-dist_11.array() * inv_lambda);
    Eigen::MatrixXd K_12 = nu * exp(-dist_12.array() * inv_lambda);
    Eigen::MatrixXd K_22 = nu * exp(-dist_22.array() * inv_lambda);
    Eigen::MatrixXd W = (1 - nu) * Eigen::MatrixXd::Identity(n_site, n_site);
    
    arma::mat K_11_arma(n_site, n_site);
    for (int i = 0; i < n_site; ++i) {
      for (int j = 0; j < n_site; ++j) {
        K_11_arma(i, j) = K_11(i, j);
      }
    }
    arma::mat K_12_arma(n_site, n_pred);
    for (int i = 0; i < n_site; ++i) {
      for (int j = 0; j < n_pred; ++j) {
        K_12_arma(i, j) = K_12(i, j);
      }
    }
    arma::mat K_22_arma(n_pred, n_pred);
    for (int i = 0; i < n_pred; ++i) {
      for (int j = 0; j < n_pred; ++j) {
        K_22_arma(i, j) = K_22(i, j);
      }
    }
    
    // misc
    Eigen::MatrixXd R_11 = K_11 + W;
    
    // loop through alleles
    for (int allele_i = 0; allele_i < alleles; ++allele_i) {
      
      // get data for this allele
      Eigen::VectorXd z = data.row(allele_i);
      
      // calculate posterior parameters for drawing mu and sigsq
      double gamma_1 = gamma_0 + R_11.sum();
      double phi_1 = (gamma_0 * phi_0 + z.transpose() * R_11.inverse() * ones_site) / gamma_1;
      double alpha_1 = alpha_0 + 0.5 * n_site;
      double beta_1 = beta_0 + 0.5*(gamma_0 * sq(phi_0) - gamma_1 * sq(phi_1) +  z.transpose() * R_11.inverse() * z);
      
      // draw mu and sigsq
      double sigsq = rgamma1(alpha_1, beta_1);
      double mu = rnorm1(phi_1, sigsq / gamma_1);
      
      // calculate Sigma_1 and mu_1
      Eigen::MatrixXd Sigma_1 = sigsq * (K_11.inverse() + W.inverse()).inverse();
      Eigen::VectorXd mu_1 = Sigma_1 * (mu / sigsq * K_11.inverse() * ones_site + z / (sigsq*(1 - nu)) );
      
      // convert to std
      for (int i = 0; i < n_site; ++i) {
        mu_1_std[i] = mu_1(i);
        for (int j = 0; j < n_site; ++j) {
          Sigma_1_std[i][j] = Sigma_1(i, j);
        }
      }
      cholesky(Sigma_1_chol, Sigma_1_std);
      
      // draw y from posterior
      rmnorm1(y_std, mu_1_std, Sigma_1_chol);
      for (int i = 0; i < n_site; ++i) {
        y(i) = y_std[i];
      }
      
      arma::vec y_arma(n_site);
      for (int i = 0; i < n_site; ++i) {
        y_arma(i) = y[i];
      }
      
      // calculate intermediate quantities
      //Eigen::MatrixXd S = K_12.transpose() * K_11.inverse();
      //Eigen::MatrixXd M = sigsq * K_22 - sigsq * S * K_12;
      //Eigen::MatrixXd Delta = S.transpose() * M.inverse();
      
      // calculate mean and covariance of predictive distribution
      //Eigen::MatrixXd pred_covar = (M.inverse() - Delta.transpose() * (Delta * S + Sigma_1.inverse()).inverse() * Delta).inverse();
      //Eigen::VectorXd pred_mean = pred_covar * ( Delta.transpose() * (Delta * S + Sigma_1.inverse()).inverse() * (mu * Delta * S * ones_site - mu * Delta * ones_pred + Sigma_1.inverse() * mu_1) - M.inverse()*(mu * S * ones_site - mu * ones_pred) );
      
      // force inverse of K_11 to be symmetric
      Eigen::MatrixXd K_11_inv = K_11.inverse();
      //for (int i = 0; i < (n_site - 1); ++i) {
      //  for (int j = (i + 1); j < n_site; ++j) {
      //    K_11_inv(i,j) = K_11_inv(j,i);
      //  }
      //}
      
      arma::mat K_11_inv_arma = arma::inv(K_11_arma);
      arma::mat pred_covar_arma = K_22_arma - K_12_arma.t() * K_11_inv_arma * K_12_arma;
      pred_covar_arma = (pred_covar_arma.t() + pred_covar_arma) / 2;
      
      arma::vec pred_mean_arma = mu * ones_pred_arma + K_12_arma.t() * K_11_inv_arma * (y_arma - mu * ones_site_arma);
      
      //pred_mean_arma.print();
      //Rcpp::stop("foo");
      
      // get Cholesky decomposition of K_11_inv
      Eigen::LLT<Eigen::MatrixXd> K_11_llt(K_11_inv);
      Eigen::MatrixXd K_11_chol = K_11_llt.matrixL();
      
      //Eigen::MatrixXd A = K_12.transpose() * K_11_inv * K_12;
      Eigen::MatrixXd B = K_12.transpose() * K_11_chol;
      //Eigen::MatrixXd C = B * B.transpose() + 0.00001 * Eigen::MatrixXd::Identity(n_pred, n_pred);
      Eigen::MatrixXd C = B * B.transpose();
      
      // calculate mean and covariance of predictive distribution
      Eigen::MatrixXd pred_covar = sigsq * K_22 - sigsq * K_12.transpose() * K_11_inv * K_12;
      //Eigen::MatrixXd pred_covar = sigsq * K_22 - sigsq * B * B.transpose();
      //Eigen::MatrixXd pred_covar = sigsq * K_22 - sigsq * C;// + 1e-9 * Eigen::MatrixXd::Identity(n_pred, n_pred);
      Eigen::VectorXd pred_mean = mu * ones_pred + K_12.transpose() * K_11_inv * (y - mu * ones_site);
      
      /*
      return Rcpp::List::create(Rcpp::Named("K_11") = Rcpp::wrap(K_11),
                                Rcpp::Named("K_12") = Rcpp::wrap(K_12),
                                Rcpp::Named("K_22") = Rcpp::wrap(K_22),
                                Rcpp::Named("K_11_inv") = Rcpp::wrap(K_11_inv),
                                Rcpp::Named("pred_covar") = pred_covar,
                                Rcpp::Named("K_11_inv_arma") = K_11_inv_arma,
                                Rcpp::Named("C_arma") = C_arma,
                                Rcpp::Named("pred_covar_arma") = pred_covar_arma,
                                Rcpp::Named("B") = B,
                                Rcpp::Named("C") = C,
                                Rcpp::Named("pred_mean") = pred_mean);
      */
      
      // convert to std
      /*
      for (int i = 0; i < n_pred; ++i) {
        pred_mean_std[i] = pred_mean(i);
        for (int j = 0; j < n_pred; ++j) {
          pred_covar_std[i][j] = pred_covar_arma(i, j);
        }
      }
      //print(ret_min, ret_max);
      cholesky(pred_covar_chol, pred_covar_std);
      */
      /*
      print_vector(pred_mean_std);
      print("");
      print_vector(pred_covar_std[10]);
      print("");
      print_vector(pred_covar_chol[10]);
      Rcpp::stop("");
      */
      // draw z_new repeatedly
      arma::mat norm_draws = arma::mvnrnd(pred_mean_arma, pred_covar_arma, inner_reps);
      print(norm_draws.n_rows, norm_draws.n_cols);
      Rcpp::stop("");
      
      //for (int rep_inner_i = 0; rep_inner_i < inner_reps; ++rep_inner_i) {
        //z_pred[rep_inner_i][allele_i] = pred_mean_std;
        //rmnorm1(z_pred[rep_inner_i][allele_i], pred_mean_std, pred_covar_chol);
        //print_vector(z_pred[rep_inner_i][allele_i]);
      //}
      
    }  // end of allele_i loop
    
    // transform z vales back to allele frequencies
    for (int rep_inner_i = 0; rep_inner_i < inner_reps; ++rep_inner_i) {
      int rep_counter = rep_i * inner_reps + rep_inner_i;
      for (int i = 0; i < n_pred; ++i) {
        double stick_remaining = 1.0;
        for (int j = 0; j < alleles; ++j) {
          double tmp = 1.0 / (1.0 + exp(-z_pred[rep_inner_i][j][i]));
          double p = stick_remaining * tmp;
          stick_remaining -= p;
          ret[rep_counter][j][i] = p;
        }
        ret[rep_counter][alleles][i] = stick_remaining;
      }
    }
    
  }  // end of rep_i loop
  
  // return list
  return Rcpp::List::create(Rcpp::Named("ret") = ret);
}

//------------------------------------------------
// wrangle data into shape before calculating pairwise Gst
Rcpp::List wrangle_pairwise_Gst_cpp(Rcpp::List freq_array) {
  
  // get frequency data into array
  vector<vector<vector<double>>> p = rcpp_to_array_double(freq_array);
  
  // calcualte pairwise Gst
  vector<double> Gst = calculate_pairwise_Gst(p);
  
  // return list
  return Rcpp::List::create(Rcpp::Named("Gst") = Gst);
}

//------------------------------------------------
// calculate pairwise Gst averaged over loci
vector<double> calculate_pairwise_Gst(const vector<vector<vector<double>>> &p) {
  
  // get basic data dimensions
  int loci = p.size();
  int n_demes = p[0][0].size();
  
  // initialise return object
  vector<double> Gst(0.5 * n_demes * (n_demes - 1));
  
  // initialise object for storing gene identity within demes
  vector<double> JS(n_demes);
  
  // loop over loci
  for (int l = 0; l < loci; ++l) {
    int alleles = p[l].size();
    
    // get gene identity within demes
    fill(JS.begin(), JS.end(), 0.0);
    for (int i = 0; i < n_demes; ++i) {
      for (int k = 0; k < alleles; ++k) {
        JS[i] += sq(p[l][k][i]);
      }
    }
    
    // get gene identity between demes and calculate Gst
    int ix = 0;
    for (int i = 0; i < (n_demes - 1); ++i) {
      for (int j = (i + 1); j < n_demes; ++j) {
        double JT = 0.0;
        for (int k = 0; k < alleles; ++k) {
          JT += sq(0.5*(p[l][k][i] + p[l][k][j]));
        }
        Gst[ix++] += (0.5*(JS[i] + JS[j]) - JT) / (1.0 - JT);
      }
    }
    
  }  // end loop over loci
  
  // divide Gst through by number of loci
  for (int i = 0; i < Gst.size(); ++i) {
    Gst[i] /= (double)loci;
  }
  
  return Gst;
}
