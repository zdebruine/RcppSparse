#ifndef TABMATRIX_RCPP_EXPORTER_H
#define TABMATRIX_RCPP_EXPORTER

namespace Rcpp {
namespace traits {

template <>
class Exporter<tabMatrix::cscMatrix> {
    Rcpp::NumericVector x_;
    Rcpp::IntegerVector i, p, Dim;

   public:
    Exporter(SEXP x) {
        Rcpp::S4 s(x);
        if (!s.hasSlot("x") || !s.hasSlot("p") || !s.hasSlot("i") || !s.hasSlot("Dim"))
            throw std::invalid_argument("Cannot construct tabMatrix::cscMatrix from this S4 object");
        x_ = s.slot("x");
        i = s.slot("i");
        p = s.slot("p");
        Dim = s.slot("Dim");
    }

    tabMatrix::cscMatrix get() {
        return tabMatrix::cscMatrix(x_, i, p, Dim);
    }
};

template <int RTYPE>
class Exporter<tabMatrix::tabVector<RTYPE>> {
    Rcpp::Vector<RTYPE> x_;
    int length;

   public:
    Exporter(SEXP x) {
        Rcpp::S4 s(x);
        if (!s.hasSlot("x") || !s.hasSlot("length"))
            throw std::invalid_argument("Cannot construct tabMatrix::tabVector from this S4 object");
        x_ = s.slot("x");
        length = s.slot("length");
    }

    tabMatrix::tabVector<RTYPE> get() {
        return tabMatrix::tabVector<RTYPE>(x_, length);
    }
};

template <int RTYPE>
class Exporter<tabMatrix::tabMatrix<RTYPE>> {
    Rcpp::IntegerVector Dim, p;
    Rcpp::Vector<RTYPE> x_;

   public:
    Exporter(SEXP x) {
        Rcpp::S4 s(x);
        if (!s.hasSlot("x") || !s.hasSlot("p") || !s.hasSlot("Dim"))
            throw std::invalid_argument("Cannot construct tabMatrix::tabMatrix from this S4 object");
        x_ = s.slot("x");
        p = s.slot("p");
        Dim = s.slot("Dim");
    }

    tabMatrix::tabMatrix<RTYPE> get() {
        return tabMatrix::tabMatrix<RTYPE>(x_, p, Dim);
    }
};
}  // namespace traits
}  // namespace Rcpp

#endif