#include "../inst/include/RcppSparse.h"

//' Simple RcppSparse::Matrix example
//' 
//' Compute column sums of a sparse matrix
//' 
//' @param A an object of class \code{dgCMatrix}
//' @examples
//' library(Matrix)
//' A <- rsparsematrix(nrow = 10, ncol = 5, density = 0.5)
//' columnSums(A)
//' 
//' # Relevant C++ code:
//' \dontrun{
//' ```{Rcpp}
//' Rcpp::NumericVector columnSums(RcppSparse::Matrix& A) {
//'      Rcpp::NumericVector sums(A.cols());
//'      for (size_t col = 0; col < A.cols(); ++col)
//'           for (RcppSparse::Matrix::InnerIterator it(A, col); it; ++it)
//'                sums(col) += it.value();
//'      return sums;
//' }
//' ```
//' }
//[[Rcpp::export]]
Rcpp::NumericVector columnSums(RcppSparse::Matrix& A) {
    Rcpp::NumericVector sums(A.cols());
    for (size_t col = 0; col < A.cols(); ++col)
        for (RcppSparse::Matrix::InnerIterator it(A, col); it; ++it)
            sums(col) += it.value();
    return sums;
}