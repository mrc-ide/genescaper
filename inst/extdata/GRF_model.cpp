#include <Rcpp.h>
using namespace std;

//------------------------------------------------
// helper function for printing a single value or series of values
template<typename TYPE>
void print(TYPE x) {
  Rcpp::Rcout << x << "\n";
}

template<typename TYPE, typename... Args>
void print(TYPE first, Args... args) {
  Rcpp::Rcout << first << " ";
  print(args...);
}

//------------------------------------------------
// density of gamma(shape,rate) distribution
double dgamma1(double x, double shape, double rate, bool return_log) {
  return(R::dgamma(x, shape, 1.0 / rate, return_log));
}

//------------------------------------------------
// density of beta(shape1,shape2) distribution
double dbeta1(double x, double shape1, double shape2, bool return_log) {
  return R::dbeta(x, shape1, shape2, return_log);
}

//------------------------------------------------
// calculate Cholesky decomposition of positive definite matrix sigma
void cholesky(vector<vector<double>> &chol, const vector<vector<double>> &sigma) {
  
  for (int i = 0; i < int(sigma.size()); ++i) {
    for (int j = 0; j < (i+1); ++j) {
      chol[i][j] = sigma[i][j];
      if (i == j) {
        if (i > 0) {
          for (int k = 0; k < i; ++k) {
            chol[i][i] -= chol[i][k]*chol[i][k];
          }
        }
        chol[i][i] = sqrt(chol[i][i]);
      } else {
        if (j > 0) {
          for (int k = 0; k < j; ++k) {
            chol[i][j] -= chol[i][k]*chol[j][k];
          }
        }
        chol[i][j] /= chol[j][j];
      }
    }
  }
  
}

//------------------------------------------------
// return (in log space) determinant of positive definite matrix. chol is the
// Cholesky decomposition of the target matrix.
double log_determinant(const vector<vector<double>> &chol) {
  double ret = 0;
  for (int i = 0; i < int(chol.size()); i++) {
    ret += 2*log(chol[i][i]);
  }
  return ret;
}

//------------------------------------------------
// return diagonal matrix of dimension d with elements x
template<class TYPE>
std::vector<std::vector<TYPE>> diag(int d, TYPE x) {
  std::vector<std::vector<TYPE>> ret(d, std::vector<TYPE>(d));
  for (int i = 0; i < d; i++) {
    ret[i][i] = x;
  }
  return ret;
}

//------------------------------------------------
// find the inverse of the matrix M by Gauss-Jordan elimination. Note that M is
// modified inside the function, and therefore cannot be passed in by reference.
vector<vector<double>> inverse(vector<vector<double>> M) {
  
  int d = int(M.size());
  vector<vector<double>> IM = diag(d, 1.0);
  
  double temp;
  for (int j = 0; j < d; j++) {
    for (int i = j; i < d; i++) {
      if (i == j) {
        temp = M[i][i];
        for (int k = 0; k < d; k++) {
          M[i][k] /= temp;
          IM[i][k] /= temp;
        }
      } else {
        if (M[i][j]!=0) {
          temp = M[i][j];
          for (int k = 0; k < d; k++) {
            M[i][k] -= temp * M[j][k];
            IM[i][k] -= temp * IM[j][k];
          }
        }
      }
    }
  }
  for (int j = (d - 1); j > 0; j--) {
    for (int i = (j - 1); i >= 0; i--) {
      temp = M[i][j];
      for (int k = 0; k < d; k++) {
        M[i][k] -= temp * M[j][k];
        IM[i][k] -= temp * IM[j][k];
      }
    }
  }
  
  return IM;
}

//------------------------------------------------
// sum a matrix
template<class TYPE>
TYPE matrix_sum(vector<vector<TYPE>> x) {
  TYPE ret = 0;
  for (int i = 0; i < x.size(); ++i) {
    for (int j = 0; j < x[i].size(); ++j) {
      ret += x[i][j];
    }
  }
  return ret;
}

// ------------------------------------------------------------------
// [[Rcpp::export]]
SEXP loglike(Rcpp::NumericVector params, Rcpp::List data, Rcpp::List misc) {
  
  // get model parameters
  double nu = params["nu"];
  double log_lambda = params["log_lambda"];
  double lambda = std::exp(log_lambda);
  double inv_lambda = 1.0 / lambda;
  double omega = params["omega"];
  double gamma = params["gamma"];
  
  // define fixed parameters
  double phi_0 = 0.0;
  double alpha_0 = 0.01;
  double beta_0 = 0.01;
  
  // TODO - remove
  alpha_0 = 0.0;
  beta_0 = 0.0;
  gamma = 0.0;
  phi_0 = 0.0;
  
  // TODO - remove
  std::vector<double> true_sigsq = Rcpp::as<vector<double>>(misc["true_sigsq"]);
  
  // get distance matrix between all sites
  Rcpp::NumericMatrix site_dist = misc["site_dist"];
  int n_site = misc["n_site"];
  
  // create spatial kernel. Note, this has 1s on the diagonal
  vector<vector<double>> K(n_site, vector<double>(n_site, 1.0));
  for (int i = 0; i < (n_site - 1); ++i) {
    for (int j = (i + 1); j < n_site; ++j) {
      double d = site_dist(i, j);
      K[i][j] = nu * exp(-pow(d * inv_lambda, omega));
      K[j][i] = K[i][j];
    }
  }
  
  // calculate intermediate quantities that apply to all locus-allele
  // combinations
  vector<vector<double>> K_inv = inverse(K);
  double K_inv_sum = matrix_sum(K_inv);
  vector<vector<double>> K_chol(n_site, vector<double>(n_site));
  cholesky(K_chol, K);
  double K_logdet = log_determinant(K_chol);
  
  // sum loglikelihood over all locus-allele combinations
  double ret = 0.0;
  for (int combo_i = 0; combo_i < data.size(); ++combo_i) {
    
    // get data
    vector<double> z = Rcpp::as< vector<double> >(data[combo_i]);
    
    // calculate intermediate quantities that apply to this combination
    double K_inv_zsum = 0.0;
    double K_inv_zsq = 0.0;
    for (int i = 0; i < n_site; ++i) {
      for (int j = 0; j < n_site; ++j) {
        K_inv_zsum += z[i] * K_inv[i][j];
        K_inv_zsq += z[i] * z[j] * K_inv[i][j];
      }
    }
    
    // TODO - remove
    //double v = 0.01;
    //alpha_0 = true_sigsq[combo_i]*true_sigsq[combo_i] / v + 2.0;
    //beta_0 = true_sigsq[combo_i] * (alpha_0 - 1.0);
    
    double gamma_1 = gamma + K_inv_sum;
    double alpha_1 = alpha_0 + (double)n_site / 2.0;
    double phi_1 = (gamma * phi_0 + K_inv_zsum) / gamma_1;
    double beta_1 = beta_0 + 0.5*(gamma*phi_0*phi_0 - gamma_1*phi_1*phi_1 + K_inv_zsq);
    
    //ret += alpha_0*log(beta_0) - alpha_1*log(beta_1) + lgamma(alpha_1) - lgamma(alpha_0) + 0.5*log(gamma) - 0.5*log(gamma_1) - 0.5*K_logdet;
    ret += - alpha_1*log(beta_1) - 0.5*log(gamma_1) - 0.5*K_logdet;
    
  }
  
  // catch non-finite return value and replace with arbitrary small value
  if (isnan(ret)) {
    ret = -DBL_MAX / 100;
  }
  
  // return as SEXP
  return Rcpp::wrap(ret);
}

// ------------------------------------------------------------------
// [[Rcpp::export]]
SEXP logprior(Rcpp::NumericVector params, Rcpp::List misc) {
  
  // get parameters
  double nu = params["nu"];
  double log_lambda = params["log_lambda"];
  double lambda = exp(log_lambda);
  //double omega = params["omega"]; // (uniform prior does not depend on parameter value)
  //double gamma = params["gamma"]; // (uniform prior does not depend on parameter value)
  
  double nu_shape1 = misc["nu_shape1"];
  double nu_shape2 = misc["nu_shape2"];
  double lambda_shape = misc["lambda_shape"];
  double lambda_rate = misc["lambda_rate"];
  
  // calculate logprior
  double ret = dgamma1(lambda, lambda_shape, lambda_rate, true) +
    dbeta1(nu, nu_shape1, nu_shape2, true) - log(3.0 - 1.0) - log(10.0 - 0.1);
  
  // TODO - remove
  //ret = 0;
  
  // add adjustment factors for transformations
  ret += log_lambda;
  
  // return as SEXP
  return Rcpp::wrap(ret);
}

// ------------------------------------------------------------------
// NOTE: Do not edit this function name
// [[Rcpp::export]]  
SEXP create_xptr(string function_name) {  
  typedef SEXP (*funcPtr_likelihood)(Rcpp::NumericVector params, Rcpp::List data, Rcpp::List misc);  
  typedef SEXP (*funcPtr_prior)(Rcpp::NumericVector params, Rcpp::List misc);  
  
  // NOTE: If your loglikelihood function is not called "loglike" please edit:
  if (function_name == "loglike"){
    return(Rcpp::XPtr<funcPtr_likelihood>(new funcPtr_likelihood(&loglike)));
  } 
  // NOTE: If your logprior function is not called "logprior" please edit:
  if (function_name == "logprior"){
    return(Rcpp::XPtr<funcPtr_prior>(new funcPtr_prior(&logprior)));
  } 
  
  Rcpp::stop("cpp function %i not found", function_name);
}
