# Rcpp::dgCMatrix

See https://gallery.rcpp.org/articles/sparse-matrix-class/ for a vignette about this matrix. This repository and class is no longer under active development. Consider using the RcppML sparse matrix class instead:

R code:
```{R}
install.packages(RcppML)
library(RcppML)
```

Then, in your C++ file:

```Rcpp
#include <RcppML.hpp>
```

## Overview
_Rcpp::dgCMatrix_ is Rcpp class for compressed sparse column (CSC) matrices with seamless zero-copy reference passing between R and C++. Usage is very similar to the Armadillo _SpMat_ class, including row and column iterators.

* [Class structure](#class-structure)
* [Rcpp::as and Rcpp::wrap](#Rcppas-and-Rcppwrap)
* [Constructors](#constructors)
* [Attributes](#attributes)
* [Copy](#copy)
* [Sparse iterators](#sparse-iterators)
* [Element access](#element-access)
* [Marginal totals](#marginal-totals)
* [Coefficient-wise operations](#coefficient-wise-operations)
* [Scalar functions](#scalar-functions)

**Speed**: _Rcpp::dgCMatrix_ sparse iterators are at least as fast as Armadillo, and only slightly slower than Eigen ([**benchmarks**](https://github.com/zdebruine/RcppSparse/wiki/Microbenchmarks:-RcppArmadillo-and-RcppEigen)).  _Rcpp::dgCMatrix_ is almost universally faster than equivalent R "Matrix" package operations ([**benchmarks**](https://github.com/zdebruine/RcppSparse/wiki/Microbenchmarks:--Matrix-R-package)).

## Why to use
RcppArmadillo and RcppEigen make deep copies of a sparse matrix when converting from R double-type to C++ type. _Rcpp::dgcMatrix_ passes by reference, not copy. This means no time or memory is wasted doing a deep copy.

## When to use/not to use
* Use when you need to access an R dgCMatrix in C++ in _double_ type and do not need to modify non-zero values
* Use when the overhead of a deep copy of your sparse matrix is unacceptable or a copy of your very large matrix will not fit in memory
* Do not use when you need to modify zero-valued indices or perform linear algebra operations on the matrix itself
* Do not use when you need a type other than _double_ (i.e. _float_), but note that the cost of computing in _double_ rather than _float_ may be balanced out by the overhead of copying the R _double_ type to _float_ (and possibly back again to R _double_).

## How to use
Include `<RcppSparse.h>` in place of `<Rcpp.h>` in your .cpp file, or simply copy the contents of `<RcppSparse.h>` into your .cpp file.

## Class structure
The _Rcpp::dgCMatrix_ exactly corresponds to the dgCMatrix S4 class structure in the "Matrix" R package. There are three primary vectors:
* `i`: row pointers for nonzeros (_Rcpp::IntegerVector_)
* `p`: column pointers for row pointers (_Rcpp::IntegerVector_)
* `x`: nonzeros (_Rcpp::NumericVector_)

**Note:** These vectors are public, but do not manipulate them directly unless you know what you are doing, as you may invalidate the S4 object.

## Rcpp::as and Rcpp::wrap
`Rcpp::as` and `Rcpp::wrap` methods for `Rcpp::dgCMatrix` ensure seamless conversion of any _dgCMatrix_ passed to/from R to Rcpp.

**C++:**
```{Cpp}
//[[Rcpp::export]]
Rcpp::dgCMatrix square_mat(Rcpp::dgCMatrix& mat){ // Seamless conversion to/from R
     return mat.square();                         // changes the R object too!
}
```
**R:**
```{R}
library(Matrix)
mat <- rsparsematrix(1000, 1000, 0.1)
sum1 <- sum(mat)
square_mat(mat)
sum(mat) / sum1
# [1] 2
```

## Constructors
In addition to passing a dgCMatrix from R to C++ via Rcpp, _Rcpp::dgCMatrix_ may be constructed in C++. You must specify `i`, `p`, and `x`, and the number of rows, values which correspond to the respective slots in an R `Matrix::dgCMatrix`.

* `dgCMatrix(i, p, x, nrow)` as _(IntegerVector, IntegerVector, NumericVector, int)_
* `dgCMatrix(i, p, x, nrow, ncol)` as _(IntegerVector, IntegerVector, NumericVector, int, int)_
* `dgCMatrix(i, p, x, nrow, Dimnames)` as _(IntegerVector, IntegerVector, NumericVector, int, List)_
* `dgCMatrix(i, p, x, nrow, ncol, Dimnames)` as _(IntegerVector, IntegerVector, NumericVector, int, int, List)_
* `dgCMatrix(S4)` containing slot attributes for _i_, _p_, _x_, _Dim_, and _Dimnames_ as _(IntegerVector, IntegerVector, NumericVector, IntegerVector, List)_

**Note**: Dimnames must be an _Rcpp::List_ of two _Rcpp::CharacterVector_

**C++**
```{Cpp}
//[[Rcpp::export]]
Rcpp::dgCMatrix constructor_example(){
  Rcpp::IntegerVector i = {0, 2, 0, 1, 1};
  Rcpp::IntegerVector p = {0, 0, 1, 2, 4, 5};
  Rcpp::NumericVector x = {0.41, 0.35, 0.84, 0.37, 0.26};
  int nrow = 5;
  // Rcpp::dgCMatrix mat(i, p, x, nrow);
  
  // dimnames may also be specified as a list of two character vectors
  Rcpp::CharacterVector rownames = {"r1", "r2", "r3", "r4", "r5"};
  Rcpp::CharacterVector colnames = {"c1", "c2", "c3", "c4", "c5"};
  Rcpp::List dimnames = Rcpp::List::create(rownames, colnames);
  Rcpp::dgCMatrix mat(i, p, x, nrow, dimnames);

  return mat;  
}
```

**R**
```{R}
> constructor_example()
5 x 5 sparse Matrix of class "dgCMatrix"
   c1   c2   c3   c4   c5
r1  . 0.41 .    0.84 .   
r2  . .    .    0.37 0.26
r3  . .    0.35 .    .   
r4  . .    .    .    .   
r5  . .    .    .    .   
```

## Attributes
* `.rows()`, `.nrow()`, `.n_rows` number of rows
* `.cols()`, `.ncol()`, `.n_cols` number of columns
* `.n_nonzero()` number of non-zero values

## .copy()
Because Rcpp objects are references to R objects existing in memory, any manipulation of an `Rcpp::dgCMatrix` that has been passed to C++ via `Rcpp::as` will also affect the R object it points to. `.copy()` creates a deep copy of the object (as is done in `Rcpp::as` for RcppArmadillo and RcppEigen sparse matrix classes) so that operations on the returned object will no longer alter the original R object.

```{Cpp}
//[[Rcpp::export]]
Rcpp::dgCMatrix square_mat(Rcpp::dgCMatrix& mat){
     Rcpp::dgCMatrix mat2 = mat.copy();           // consumes memory and time
     return mat2.square();                        // does not change the original R object!
}
```

## Sparse iterators
Forward iterators traverse over non-zero elements within the specified range with either read-only (copy) or read/write (reference) access.

### Iterator classes:
* `::iterator` forward iterator traversing all non-zero elements, passes reference to each element for read/write access
* `::const_iterator` passes copy of value of each element for read-only access.
* `::col_iterator` read/write access to all non-zero elements of a given column
* `::const_col_iterator` read-only access to all non-zero elements of a given column
* `::row_iterator` read/write access to all non-zero elements of a given row (not as efficient as col_iterator)
* `::const_row_iterator` read-only to all non-zero elements of a given row

### Iterator constructors:
* `.begin()`, `.end()`: returns read/write iterators at the first and last non-zero indices
* `.const_begin()`, `.const_end()` returns const iterators at the first and last non-zero indices
* `.begin_col(int)`, `.end_col(int)` returns read/write iterators at the first and last non-zero indices in a given column
* `.const_begin_col(int)`, `.const_end_col(int)` returns const iterators at the first and last non-zero indices in a given column
* `.begin_row(int)`, `.end_row(int)` returns read/write iterators at the first and last non-zero indices in a given row
* `.const_begin_row(int)`, `.const_end_row(int)` returns const iterators at the first and last non-zero indices in a given row

### Member functions:
* `.value()` or `*`: non-zero value at index
* `.row()`: row at index (fast)
* `.col()`: col at index (slow if not using a column iterator)
(applies to all iterator classes)

### Operators:
* `!=`, `<` logical comparison of two iterators of the same class
* `++` post-increment

**Note**:
* Setting values to zero will not alter indexing. The matrix pointers are not adjusted.
* Indices with zero values cannot be accessed as in RcppArmadillo or RcppEigen sparse matrix classes.
* Generally avoid row iterators if possible, although this implementation is not inferior to RcppArmadillo or RcppEigen.
* There is no support for mixing-and-matching constant and non-constant iterators

```{Cpp}
//[[Rcpp::export]]
void print_mat(Rcpp::dgCMatrix& mat){     
     for(Rcpp::dgCMatrix::const_iterator it = mat.const_begin(); it != X.const_end(); ++it){
          Rprintf("row: %5d, col: %5d, val: %5d \n", it.row(), it.col(), (*it));
     }
}

// example of read-only access: dot product between a sparse column and dense vector
//[[Rcpp::export]]
double sparse_dot(Rcpp::dgCMatrix& mat, Rcpp::NumericVector& vec, int col){
     double dot_product = 0;
     for(Rcpp::dgCMatrix::const_col_iterator it = mat.const_col_begin(col); it < mat.const_col_end(col); it++){
          dot_product += (*it) * vec[it.row()];
     }
     return dot_product;
}

// example of write access: multiply each sparse column by a dense vector
//[[Rcpp::export]]
void multiply_by_vec(Rcpp::dgCMatrix& mat, Rcpp::NumericVector& vec){
     for(Rcpp::dgCMatrix::iterator it = mat.begin(); it < mat.end(); it++)
          (*it) *= vec[it.row()];
}
```

## Element Access
### Matrix Subviews
`(i, j)` where `i` and `j` may be any combination of _int_ or _IntegerVector_ giving valid indices within the sparse matrix.

The above returns value copies of non-contiguous subviews copied to a _double_, _NumericVector_, or _NumericMatrix_, including zeros.

### Dense Marginal Views
`.col(int)`, `.row(int)`, `.cols(IntegerVector)`, `.rows(IntegerVector)`
Returns a _NumericVector_ or _NumericMatrix_ corresponding to the marginal values, including zeros.

If you need write access to non-zero values in subviews, use sparse iterators. If this does not meet your needs, consider RcppArmadillo and RcppEigen.

### Direct Access
Elements may be directly accessed (and modified) by manipulating the underlying vectors:
*`i` or `.innerIndexPtr()` is a referenced _IntegerVector_ giving pointers to row indices
*`p` or `.outerIndexPtr()` is a referenced _IntegerVector_ giving pointers in `i` corresponding to the first non-zero element in each column.
*`x` or `.nonzeros()` is a referenced _NumericVector_ giving non-zero values corresponding to rows in `i` and columns delineated by `p`.

Note that the above organization scheme is identical to the `Matrix::dgCMatrix` class in R.

```{Cpp}
Rcpp::IntegerVector indices = {0, 2, 0, 1, 1};
A(1, 5);
A(1, indices);
A(indices, indices);
A.rows(indices);

// direct access (example, square all non-zero values in row 0)
for(int ptr = 0; ptr < mat.n_nonzero(); ++ptr)
     if(i[ptr] == 0) x[ptr] *= x[ptr];
```

## Marginal totals
`.colSums()`, `.colMeans()`, `.rowSums()`, `.rowMeans()`

These member functions use iterators to sum marginal values and return _NumericVector_.

```{Cpp}
//[[Rcpp::export]]
Rcpp::NumericVector Rcpp_rowSums(Rcpp::dgCMatrix& mat){
     return mat.rowSums();
}

// which is equivalent to:
Rcpp::NumericVector Rcpp_rowSums(Rcpp::dgCMatrix& mat){
     Rcpp::NumericVector sums(mat.rows());
     for(Rcpp::dgCMatrix::const_iterator it = mat.const_begin(); it < mat.const_end(); it++)
          sums[it.row()] += it.value();
     return sums;
}
```

## Coefficient-wise operations
`.abs()`, `.ceil()`, `.floor()`, `.round()`, `.sqrt()`, `.square()`, `.trunc()`

## Scalar functions
`.sum()`, `.mean()`

## .crossprod()
Returns the cross-product of the matrix as _NumericMatrix_ (matrix of column-wise dot products). If OpenMP is defined, column-wise parallelization is employed.
