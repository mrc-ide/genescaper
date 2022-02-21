#include <Rcpp.h>
using namespace std;

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
  double log_lambda = params["log_lambda"];
  double lambda = std::exp(log_lambda);
  double inv_lambda = 1.0 / lambda;
  double nu = params["nu"];
  
  // get fixed parameters
  double mu_mean = misc["mu_mean"];
  double mu_scale = misc["mu_scale"];
  double sigsq_mean = misc["sigsq_mean"];
  double sigsq_var = misc["sigsq_var"];
  
  // reparameterise for convenience
  double phi_0 = mu_mean;
  double gamma_0 = 1.0 / mu_scale;
  double alpha_0 = sigsq_mean * sigsq_mean / sigsq_var + 2.0;
  double beta_0 = sigsq_mean * (alpha_0 - 1.0);
  
  // get distance matrix between all sites
  Rcpp::NumericMatrix site_dist = misc["site_dist"];
  int n_site = misc["n_site"];
  
  // create correlation matrix
  vector<vector<double>> R(n_site, vector<double>(n_site, 1.0));
  for (int i = 0; i < (n_site - 1); ++i) {
    for (int j = (i + 1); j < n_site; ++j) {
      double d = site_dist(i, j);
      R[i][j] = (1 - nu) * exp(-d * inv_lambda);
      R[j][i] = R[i][j];
    }
  }
  
  // calculate intermediate quantities that apply to all combos
  vector<vector<double>> R_inv = inverse(R);
  double R_inv_sum = matrix_sum(R_inv);
  vector<vector<double>> R_chol(n_site, vector<double>(n_site));
  cholesky(R_chol, R);
  double R_logdet = log_determinant(R_chol);
  double gamma_1 = gamma_0 + R_inv_sum;
  double alpha_1 = alpha_0 + (double)n_site / 2.0;
  
  // sum loglikelihood over all locus&allele combos
  double ret = 0.0;
  for (int combo_i = 0; combo_i < data.size(); ++combo_i) {
    
    // get data
    vector<double> z = Rcpp::as< vector<double> >(data[combo_i]);
    
    // calculate intermediate quantities that apply to this combo
    double R_inv_xsum = 0.0;
    double R_inv_xsq = 0.0;
    for (int i = 0; i < n_site; ++i) {
      for (int j = 0; j < n_site; ++j) {
        R_inv_xsum += z[i] * R_inv[i][j];
        R_inv_xsq += z[i] * z[j] * R_inv[i][j];
      }
    }
    double phi_1 = (gamma_0 * phi_0 + R_inv_xsum) / gamma_1;
    double beta_1 = beta_0 + 0.5*(gamma_0*phi_0*phi_0 - gamma_1*phi_1*phi_1 + R_inv_xsq);
    
    // calculate multivariate normal likelihood
    ret += alpha_0*log(beta_0) - alpha_1*log(beta_1) + lgamma(alpha_1) - lgamma(alpha_0) + 0.5*log(gamma_0) -0.5*log(gamma_1) - 0.5*R_logdet;
    
  }
  
  // return as SEXP
  return Rcpp::wrap(ret);
}

// ------------------------------------------------------------------
// [[Rcpp::export]]
SEXP logprior(Rcpp::NumericVector params, Rcpp::List misc) {
  
  // get parameters
  double log_lambda = params["log_lambda"];
  double lambda = exp(log_lambda);
  double nu = params["nu"];
  
  double lambda_mean = misc["lambda_mean"];
  double lambda_var = misc["lambda_var"];
  double nu_shape1 = misc["nu_shape1"];
  double nu_shape2 = misc["nu_shape2"];
  
  // reparameterise for convenience
  double lambda_shape = lambda_mean * lambda_mean / lambda_var;
  double lambda_rate = lambda_mean / lambda_var;
  
  // calculate logprior
  double ret = dgamma1(lambda, lambda_shape, lambda_rate, true) +
    dbeta1(nu, nu_shape1, nu_shape2, true);
  
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
