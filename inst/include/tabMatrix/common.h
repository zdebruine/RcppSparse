#ifndef TABMATRIX_COMMON_H
#define TABMATRIX_COMMON_H

#include <RcppCommon.h>

// forward declare classes
namespace tabMatrix {
class cscMatrix;
template <int RTYPE>
class tabVector;
template <int RTYPE>
class tabMatrix;
}  // namespace tabMatrix

// forward declare Rcpp::as<> Exporter
namespace Rcpp {

namespace traits {

template <>
class Exporter<tabMatrix::cscMatrix>;

template <int RTYPE>
class Exporter<tabMatrix::tabVector<RTYPE>>;

template <int RTYPE>
class Exporter<tabMatrix::tabMatrix<RTYPE>>;
}  // namespace traits
}  // namespace Rcpp

#include <Rcpp.h>

//[[Rcpp::plugins(openmp)]]
#ifdef _OPENMP
#include <omp.h>
#endif

// typedefs for convenience
namespace tabMatrix {
using itabVector = tabVector<INTSXP>;
using dtabVector = tabVector<REALSXP>;
using itabMatrix = tabMatrix<INTSXP>;
using dtabMatrix = tabMatrix<REALSXP>;
}  // namespace tabMatrix

#endif