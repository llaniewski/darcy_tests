#define ARMA_USE_SUPERLU 1
#define ARMA_64BIT_WORD 1

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;

// [[Rcpp::export]]
vec arma_spsolve(umat locations, vec values, vec rhs) {
  sp_mat A(locations, values);
  vec x=rhs;
  bool status = spsolve(x, A, rhs);
  if (status == false)  { Rcpp::Rcout << "no solution" << endl; }
  return x;
}
