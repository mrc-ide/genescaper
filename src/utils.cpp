
#include "utils.h"

using namespace std;

//------------------------------------------------
// evaluate the matrix operation t(e) %*% M %*% e, where e is a vector of 1s of
// the appropriate length. Not that this is equivalent to the sum of all
// elements in the matrix.
double matrix_eMe(const vector<vector<double>> &M) {
  double ret = 0;
  for (int i = 0; i < M.size(); ++i) {
    for (int j = 0; j < M[i].size(); ++j) {
      ret += M[i][j];
    }
  }
  return ret;
}

//------------------------------------------------
// evaluate the matrix operation t(x) %*% M %*% e, where x is a vector of
// the appropriate length and e is a vector of 1s.
double matrix_xMe(const vector<double> &x, const vector<vector<double>> &M) {
  double ret = 0;
  for (int i = 0; i < M.size(); ++i) {
    for (int j = 0; j < M[i].size(); ++j) {
      ret += x[i] * M[i][j];
    }
  }
  return ret;
}

//------------------------------------------------
// evaluate the matrix operation t(x) %*% M %*% x, where x is a vector of
// the appropriate length.
double matrix_xMx(const vector<double> &x, const vector<vector<double>> &M) {
  double ret = 0;
  for (int i = 0; i < M.size(); ++i) {
    for (int j = 0; j < M[i].size(); ++j) {
      ret += x[i] * x[j] * M[i][j];
    }
  }
  return ret;
}
