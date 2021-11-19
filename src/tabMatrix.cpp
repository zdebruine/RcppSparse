#include "../inst/include/tabMatrix.h"

// -- TABVECTOR METHODS

// COERCION METHODS

//[[Rcpp::export]]
Rcpp::IntegerVector Rcpp_itabVector_to_Vector(tabMatrix::itabVector x) {
    return x.as_Vector();
}

//[[Rcpp::export]]
Rcpp::NumericVector Rcpp_dtabVector_to_Vector(tabMatrix::dtabVector x) {
  return x.as_Vector();
}

//[[Rcpp::export]]
Rcpp::S4 Rcpp_tabVector_to_isparseVector(tabMatrix::tabVector<INTSXP> x) {
    Rcpp::S4 s("isparseVector");
    Rcpp::IntegerVector x_, i_;
    for (tabMatrix::tabVector<INTSXP>::Iterator it(x); it; ++it) {
        x_.push_back(it.value());
        i_.push_back(it.row() + 1);
    }
    s.slot("length") = x.length;
    s.slot("x") = x_;
    s.slot("i") = i_;
    return s;
}

//[[Rcpp::export]]
Rcpp::S4 Rcpp_isparseVector_to_itabVector(const Rcpp::S4 x) {
    Rcpp::IntegerVector i_ = x.slot("i");
    for (int& i__ : i_)
        i__ -= 1;
    tabMatrix::itabVector v = tabMatrix::itabVector(x.slot("x"), i_, x.slot("length"));
    return v.as_S4();
}

// SUBSETTING

//[[Rcpp::export]]
Rcpp::S4 Rcpp_subset_itabvector_vector(tabMatrix::itabVector x, Rcpp::IntegerVector ind) {
    return x(ind).as_S4();
}

//[[Rcpp::export]]
int Rcpp_sum_itabvector(tabMatrix::itabVector x) {
    return x.sum();
}

//[[Rcpp::export]]
Rcpp::LogicalVector Rcpp_negate(tabMatrix::itabVector x) {
    return x.is_zero();
}

//[[Rcpp::export]]
int Rcpp_crossprod_itabVector_itabVector(tabMatrix::itabVector v1, tabMatrix::itabVector v2) {
    return v1.crossprod(v2);
}

//[[Rcpp::export]]
int Rcpp_crossprod_itabVector_IntegerVector(tabMatrix::itabVector v1, Rcpp::IntegerVector v2) {
    return v1.crossprod(v2);
}

//[[Rcpp::export]]
int sum1(tabMatrix::itabVector v1) {
    int v = 0;
    for (tabMatrix::itabVector::Iterator it(v1); it; ++it)
        v += it.value();
    return v;
}

//[[Rcpp::export]]
int sum2(tabMatrix::itabVector v1) {
    int v = 0;
    int curr_value, prev_index = 0;
    tabMatrix::itabVector::Iterator it(v1);
    while (it) {
        curr_value = it.value();
        it.nextValue();
        v += curr_value * (it.index() - 1 - prev_index);
        prev_index = it.index();
    }
    return v;
}

// -- TABMATRIX METHODS

//[[Rcpp::export]]
Rcpp::S4 as_itabMatrix_from_cscMatrix(tabMatrix::cscMatrix A) {
    tabMatrix::itabMatrix A_(A);
    return A_.as_S4();
}

//[[Rcpp::export]]
Rcpp::S4 as_dgCMatrix_from_itabMatrix(tabMatrix::itabMatrix A) {
    return A.as_cscMatrix().wrap();
}