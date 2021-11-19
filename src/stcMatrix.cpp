#include "../inst/include/RcppSparse.h"

//[[Rcpp::depends(RcppClock)]]
#include "RcppClock.h"

//[[Rcpp::export]]
Rcpp::S4 as_stcMatrix(RcppSparse::Matrix A) {
    RcppSparse::stcMatrix A_(A);
    return A_.wrap();
}

//[[Rcpp::export]]
Rcpp::S4 as_dgCMatrix(RcppSparse::stcMatrix A){
  return A.as_Matrix().wrap();
}

//[[Rcpp::export]]
void rowit(RcppSparse::stcMatrix A){
  Rprintf("%5s %5s %5s\n", "row", "col", "val");
  for(int i = 0; i < A.rows(); ++i)
    for(RcppSparse::stcMatrix::RowIterator it(A, i); it; ++it)
      Rprintf("%5i %5i %5i\n", it.row(), it.col(), it.value());
}