# RcppSparse

This is a lightweight read-only `Rcpp::dgCMatrix` class which references a corresponding R dgCMatrix. This achieves zero-copy/zero-overhead access and conversion between R and C++ (and vice versa).

Compare to RcppEigen `SparseMatrix` and RcppArmadillo `SpMat` which are deep copies, not references to R objects.


**WHAT IT IS**

The dgCMatrix class is simply three vectors:
- `i`: row pointers for nonzeros (`Rcpp::IntegerVector`)
-  `p`: column pointers for row pointers (`Rcpp::IntegerVector`)
-  `x`: nonzeros (`Rcpp::NumericVector`)
And a handy `const_iterator` class complete with methods for element access or subviews.


**WHEN TO USE**

- Use when you only need **read-only** access to elements, columns, rows, subviews, or iterators.
- Use when a copy of your very large matrix will not fit in memory


**WHEN NOT TO USE**

- You need **write access**, linear algebra operations, or element-wise operations.
- You need a type other than "double", such as "float".*
*Note that the cost of computing in "double" may be balanced out by the overhead of copying R doubles to C++ floats (and possibly back again to R doubles).


**HOW TO USE**

Include this header in your R package or in your C++ file.

**DOCUMENTATION**

See the wiki!
