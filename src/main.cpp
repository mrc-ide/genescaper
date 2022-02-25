
#include "main.h"
#include "misc_v13.h"
#include "probability_v17.h"
#include "utils.h"

using namespace std;

//------------------------------------------------
Rcpp::List predict_map_cpp(Rcpp::List data_list, Rcpp::List mcmc_sample,
                           Rcpp::NumericMatrix dist_11, Rcpp::NumericMatrix dist_12,
                           Rcpp::NumericMatrix dist_22, Rcpp::List params,
                           int inner_reps, Rcpp::List args_progress,
                           Rcpp::List args_functions, Rcpp::List args_misc) {
  
  // get functions
  Rcpp::Function update_progress = args_functions["update_progress"];
  
  // get misc parameters
  bool pb_markdown = rcpp_to_bool(args_misc["pb_markdown"]);
  
  // get data into array
  vector<vector<double>> data_array = rcpp_to_matrix_double(data_list);
  int alleles = data_array.size();
  
  // get input dimensions
  int n_site = dist_11.ncol();
  int n_pred = dist_22.ncol();
  
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
  
  // initialise kernel matrices
  vector<vector<double>> R_11(n_site, vector<double>(n_site, 1.0));
  vector<vector<double>> R_12(n_site, vector<double>(n_pred, 1.0));
  vector<vector<double>> R_22(n_pred, vector<double>(n_pred, 1.0));
  
  // initialise output object
  vector<vector<vector<double>>> ret(reps * inner_reps, vector<vector<double>>(alleles + 1, vector<double>(n_pred)));
  
  // initialise dummy objects
  vector<vector<vector<double>>> z_pred(inner_reps, vector<vector<double>>(alleles));
  vector<vector<double>> cov_unscaled_Cholesky(n_pred, vector<double>(n_pred, 1.0));
  
  // loop through MCMC draws
  for (int rep_i = 0; rep_i < reps; ++rep_i) {
    
    // allow user to exit on escape
    Rcpp::checkUserInterrupt();
    
    // update progress bar
    update_progress_cpp(args_progress, update_progress, "pb", rep_i, reps, !pb_markdown);
    
    // populate correlation matrices based on this set of MCMC parameter draws
    for (int i = 0; i < (n_site - 1); ++i) {
      for (int j = (i + 1); j < n_site; ++j) {
        R_11[i][j] = (1 - nu_draws[rep_i]) * exp( - dist_11(i, j) / lambda_draws[rep_i]);
        R_11[j][i] = R_11[i][j];
      }
    }
    for (int i = 0; i < n_site; ++i) {
      for (int j = 0; j < n_pred; ++j) {
        R_12[i][j] = (1 - nu_draws[rep_i]) * exp( - dist_12(i, j) / lambda_draws[rep_i]);
      }
    }
    for (int i = 0; i < (n_pred - 1); ++i) {
      for (int j = (i + 1); j < n_pred; ++j) {
        R_22[i][j] = (1 - nu_draws[rep_i]) * exp( - dist_22(i, j) / lambda_draws[rep_i]);
        R_22[j][i] = R_22[i][j];
      }
    }
    
    // calculate properties of matrices that are common to all allele
    // elements
    vector<vector<double>> R_11_inv = inverse(R_11);
    double R_11_eRe = matrix_eMe(R_11_inv);
    
    // solution to t(R_12) %*% solve(R_11). Unfortunately this matrix does not
    // have an intuitive name
    vector<vector<double>> tmp_mat_1(n_pred, vector<double>(n_site));
    for (int i = 0; i < n_pred; ++i) {
      for (int j = 0; j < n_site; ++j) {
        for (int k = 0; k < n_site; ++k) {
          tmp_mat_1[i][j] += R_12[k][i] * R_11_inv[k][j];
        }
      }
    }
    
    // solution to t(R_12) %*% solve(R_11) %*% rep(1, n_site). Unfortunately
    // this vector does not have an intuitive name
    vector<double> tmp_vec_1(n_pred);
    for (int i = 0; i < n_pred; ++i) {
      tmp_vec_1[i] = sum(tmp_mat_1[i]);
    }
    
    // covariance of predictive distribution before scaling by sigsq, and
    // Cholesky decomposition of this matrix
    vector<vector<double>> cov_unscaled = R_22;
    for (int i = 0; i < n_pred; ++i) {
      for (int j = 0; j < n_pred; ++j) {
        for (int k = 0; k < n_site; ++k) {
          cov_unscaled[i][j] -= tmp_mat_1[i][k] * R_12[k][j];
        }
      }
    }
    cholesky(cov_unscaled_Cholesky, cov_unscaled);
    
    // loop through alleles
    for (int allele_i = 0; allele_i < alleles; ++allele_i) {
      
      // get data for this allele
      vector<double> z = data_array[allele_i];
      
      // perform matrix operations specific to this allele 
      double R_11_zRe = matrix_xMe(z, R_11_inv);
      double R_11_zRz = matrix_xMx(z, R_11_inv);
      
      // calculate posterior parameters for drawing mu and sigsq
      double gamma_1 = gamma_0 + R_11_eRe;
      double phi_1 = (gamma_0 * phi_0 + R_11_zRe) / gamma_1;
      double alpha_1 = alpha_0 + 0.5 * n_site;
      double beta_1 = beta_0 + 0.5*(gamma_0 * sq(phi_0) - gamma_1 * sq(phi_1) + R_11_zRz);
      
      // draw mu and sigsq
      double sigsq = rgamma1(alpha_1, beta_1);
      double mu = rnorm1(phi_1, sigsq / gamma_1);
      
      // calculate mean of predictive distribution
      vector<double> pred_mean(n_pred);
      for (int i = 0; i < n_pred; ++i) {
        pred_mean[i] = mu - mu*tmp_vec_1[i];
        for (int j = 0; j < n_site; ++j) {
          pred_mean[i] += tmp_mat_1[i][j] * z[j];
        }
      }
      
      // draw z_new repeatedly
      double sigma = pow(sigsq, 0.5);
      for (int rep_inner_i = 0; rep_inner_i < inner_reps; ++rep_inner_i) {
        rmnorm1(z_pred[rep_inner_i][allele_i], pred_mean, cov_unscaled_Cholesky, sigma);
      }
      
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
