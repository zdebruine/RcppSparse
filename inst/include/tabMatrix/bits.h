#ifndef TABMATRIX_BITS_H
#define TABMATRIX_BITS_H

#ifndef TABMATRIX_COMMON_H
#include "common.h"
#endif

template <int RTYPE>
void reorder(std::vector<int> idx, Rcpp::Vector<RTYPE>& v) {
    Rcpp::Vector<RTYPE> result(idx.size());
    for (size_t i = 0; i < idx.size(); ++i)
        result[i] = v[idx[i]];
    v = result;
}

void reorder(std::vector<int> idx, std::vector<int>& v) {
    std::vector<int> result(idx.size());
    for (size_t i = 0; i < idx.size(); ++i)
        result[i] = v[idx[i]];
    v = result;
}

template <typename T>
void reorder(std::vector<int> idx, std::vector<int>& v) {
    std::vector<T> result(idx.size());
    for (size_t i = 0; i < idx.size(); ++i)
        result[i] = v[idx[i]];
    v = result;
}

template <int RTYPE>
Rcpp::Vector<RTYPE> join(Rcpp::Vector<RTYPE> a, Rcpp::Vector<RTYPE> b) {
    Rcpp::Vector<RTYPE> c(a.size() + b.size());
    for (int i = 0; i < a.size(); ++i)
        c(i) = a(i);
    for (int j = a.size(), k = 0; k < b.size(); ++j, ++k)
        c(j) = b(k);
    return c;
}

Rcpp::IntegerVector find_true(Rcpp::LogicalVector v) {
    // Rcpp doesn't have handy features like .resize, .reserve, etc.
    Rcpp::IntegerVector res(v.size());
    int j = 0;
    for (int i = 0; i < v.size(); ++i) {
        if (v[i]) {
            res[j] = i;
            ++j;
        }
    }
    Rcpp::IntegerVector result(j + 1);
    for (int i = 0; i < (j + 1); ++i)
        result[i] = res[i];
    return result;
}

#endif