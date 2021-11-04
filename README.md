# RcppSparse

- [Documentation](https://zdebruine.github.io/RcppSparse/articles/Documentation.html)
- [Example](https://github.com/zdebruine/RcppSparse/blob/main/src/example.cpp)
- Install: `install.packages("RcppSparse")`

RcppSparse provides a seamless Rcpp object class for R sparse matrix objects. The `RcppSparse::Matrix` class can directly import `Matrix::dgCMatrix-class` objects from R without any copying by simply using base Rcpp types (`IntegerVector` and `NumericVector`).

The result is a constant, by-reference, zero-copy Rcpp sparse matrix class. In contrast, RcppArmadillo and RcppEigen are deep copies.

Install RcppSparse from CRAN with `install.packages("RcppSparse")` and then load the header library in your C++ file:

```{Rcpp}
//[[Rcpp::depends(RcppSparse)]]
#include RcppSparse.h // this line before ever including Rcpp.h
```

Now you're set to roll!  Here's a simple program to compute column sums of a sparse matrix:

```{Rcpp}
//[[Rcpp::export]]
Rcpp::NumericVector columnSums(RcppSparse::Matrix& A) {
    Rcpp::NumericVector sums(A.cols());
    for (size_t col = 0; col < A.cols(); ++col)
        for (RcppSparse::Matrix::InnerIterator it(A, col); it; ++it)
            sums(col) += it.value();
    return sums;
}
```

To use this function from the R end, just make sure you pass a `dgCMatrix` into the function:

```{R}
library(Matrix)
A <- rsparsematrix(nrow = 10, ncol = 10, density = 0.1)
class(A) # this is a "dgCMatrix"
columnSums(A)
```

The idea of a zero-copy sparse matrix class is discussed on the [Rcpp gallery](https://gallery.rcpp.org/articles/sparse-matrix-class/).
