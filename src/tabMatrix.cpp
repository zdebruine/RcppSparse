#include "../inst/include/RcppSparse.h"

//[[Rcpp::depends(RcppClock)]]
#include "RcppClock.h"

//[[Rcpp::export]]
Rcpp::S4 as_tabMatrix(RcppSparse::Matrix A) {
  RcppSparse::tabMatrix A_(A);
  return A_.wrap();
}

//[[Rcpp::export]]
Rcpp::S4 as_dgCMatrix(RcppSparse::tabMatrix A) {
  return A.as_Matrix().wrap();
}

//[[Rcpp::export]]
Rcpp::NumericMatrix print_tabMatrix(RcppSparse::tabMatrix A, unsigned int ncol = 10, unsigned int nrow = 10) {
  Rcpp::NumericMatrix res(nrow, ncol);
  if (ncol > A.cols()) ncol = A.cols();
  for (unsigned int i = 0; i < ncol; ++i) {
    for (RcppSparse::tabMatrix::InnerIterator it(A, i); it; ++it) {
      if ((unsigned int)it.row() < nrow) {
        res(it.row(), it.col()) = it.value();
      }
    }
  }
  return res;
}

//[[Rcpp::export]]
std::vector<int> tabMatrix_colSums(RcppSparse::tabMatrix A) {
  std::vector<int> sums(A.cols());
  #ifdef _OPENMP
  #pragma omp parallel for
  #endif
  for (unsigned int col = 0; col < A.cols(); ++col)
    for (RcppSparse::tabMatrix::InnerIterator it(A, col); it; ++it)
      sums[col] += it.value();
  return sums;
}

//[[Rcpp::export]]
std::vector<int> dgCMatrix_colSums(RcppSparse::Matrix A) {
  std::vector<int> sums(A.cols());
  #ifdef _OPENMP
  #pragma omp parallel for
  #endif
  for (unsigned int col = 0; col < A.cols(); ++col)
    for (RcppSparse::Matrix::InnerIterator it(A, col); it; ++it)
      sums[col] += it.value();
  return sums;
}

//[[Rcpp::export]]
Rcpp::NumericVector tabVector_as_numeric(Rcpp::IntegerVector x, int length) {
  Rcpp::NumericVector result(length);
  if (length > 0) {
    int cur_value = std::abs(x[0]);
    for (int i = 1; i < x.size(); ++i) {
      if(x[i] < 0) {
        cur_value = std::abs(x[i]);
        ++i;
        result[x[i]] = cur_value;
      } else {
        result[x[i]] = cur_value;
      }
    }
  }
  return result;
}

//[[Rcpp::export]]
Rcpp::S4 tabMatrix_column(RcppSparse::tabMatrix A, int col){
  RcppSparse::tabVector v = A.col(col);
  return v.wrap();
}