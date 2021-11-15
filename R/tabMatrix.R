#' Tabulated sparse matrix
#' 
#' @details 
#' The Tabulated Matrix (\code{tabMatrix}) and tabulated vector (\code{tabVector}) classes hold discrete positive values in a compressed form.
#' 
#' \code{tabMatrix} is intended for:
#' 
#' 1. Compressed storage of CSC sparse matrices with discrete positive values.
#' 2. Column-wise operations in C++ against scalars, dense vectors, or matrices (not sparse-sparse operations).
#'
#' **Important:** Tabulated vectors are sorted firstly by value and secondarily by row index, thus column-major random access is not be ordered by row. This 
#' does not affect the performance of dense-sparse operations because row access is stochastic regardless of how sparse values are sorted. However, 
#' this does affect the performance of sparse-sparse operations because iterator coupling (i.e. forward traversal iterators) is no longer efficient.
#'
#' Prefer \code{tabMatrix} for dense-sparse operations and for efficient storage, prefer \code{dgCMatrix} if any sparse-sparse operations are necessary.
#'
#' Coercion methods to/from \code{dgCMatrix} or \code{matrix} are fully parallelized.
#'
#' Random access iterators to values in tabulated vectors are equally as efficient as \code{dgCMatrix} and \code{sparseVector} equivalents.
#'
#' Operations and subsetting \code{tabMatrix} is not currently supported in R. The corresponding Rcpp class in the tabMatrix header library provides a 
#' convenient column iterator class for traversal over non-zero values.
#'
#' @methods
#' Coerce to/from \code{dgCMatrix} using \code{as(matrix, "tabMatrix")}. Any other coercion attempts to coerce to \code{dgCMatrix} either directly or 
#' via \code{matrix}.
#'
#' * \code{t()}: Transposition is supported, but not recommended over transposition of \code{dgCMatrix}.
#' 
#' @slot p pointers to indices in \code{x} corresponding to a new column
#' @slot Dim integer vector of two giving number of rows and columns
#' @slot Dimnames list of two giving rownames and colnames
#' @slot x A vector containing values followed by row indices, ordered firstly by column, secondly by increasing value, and thirdly by increasing row index.
#' @slot factors list for value mappings, etc.
#' @md
setClass("tabMatrix", 
         representation(p = "integer", Dim = "integer", Dimnames = "list", x = "integer", factors = "list"),
         prototype(Dim = c(0L, 0L), Dimnames = list(NULL, NULL), p = NA_integer_, x = NA_integer_, factors = list()))

setAs("dgCMatrix", "tabMatrix", function(from, to) {
  # check if matrix contains any "NA" values
  if(any(is.na(from@x))) stop("matrix contains NA values")
  
  # check if values in from@x are integers
  len <- length(from@x)
  if(len > 100) len <- 100
  if(sum(as.integer(from@x[1:len])) != sum(from@x[1:len])) stop("matrix values are not integer counts")

  # check that from@x is a numeric vector, for compatibility with Rcpp::NumericVector
  if(class(from@x) != "numeric") from@x <- as.numeric(from@x)

  # check that the minimum value in from@x is not less than 0
  if(min(from@x) < 0) stop("negative values are not permitted")

  # run C++ conversion routine
  to <- as_tabMatrix(from)
  to@Dimnames <- from@Dimnames
  to@factors <- from@factors
  to
})

setAs("matrix", "tabMatrix", function(from, to) as(as(from, "dgCMatrix"), "tabMatrix"))
setAs("tabMatrix", "matrix", function(from, to) as(as(from, "dgCMatrix"), "matrix"))

setAs("tabMatrix", "dgCMatrix", function(from, to) {
  to <- as_dgCMatrix(from)
  to@Dimnames <- from@Dimnames
  to@factors <- from@factors
  ind <- which.min(to@x) - 1
  to@x <- to@x[1:ind]
  to@i <- to@i[1:ind]
  to
})

setMethod("print", "tabMatrix", function(x, nrow = 10, ncol = 10){
  p <- print_tabMatrix(x, ncol, nrow)
  cat(x@Dim[1],"x",x@Dim[2],"sparse Matrix of class \"tabMatrix\"\n\n")
  rownames(p) <- x@Dimnames[[1]][1:nrow]
  print.table(p, zero.print = ".")
  if(nrow < x@Dim[1] && ncol < x@Dim[2]) {
    cat("\nSuppressing", x@Dim[1] - nrow, "rows and", x@Dim[2] - ncol,"columns\n")
  } else if(nrow < x@Dim[1]){
    cat("\nSuppressing", x@Dim[1] - nrow, "rows\n")
  } else if(ncol < x@Dim[2]){
    cat("\nSuppressing", x@Dim[2] - ncol, "columns\n")
  }
  invisible(x)
})

setMethod("t", "tabMatrix", function(x){
  x <- as(t(as(x, "dgCMatrix")), "tabMatrix")
})

setMethod("dim", "tabMatrix", function(x) x@Dim)

setClass("tabVector",
  representation(length = "integer", x = "integer"),
  prototype(length = 0L, x = NA_integer_))

setAs("tabVector", "numeric", function(from, to) {
  tabVector_as_numeric(from@x, from@length)
})

setMethod("[", signature("tabMatrix"), function(x, i = NULL, j = NULL) {
  if(!is.null(i)) cat(i)
  if(!is.null(j)) cat(j)

# [i, j] where i and j are scalar -> integer
# [i, j] where i or j is IntegerVector -> tabVector
# [i, j] where i and j are IntegerVector -> tabMatrix

})

setMethod("[", c("tabMatrix", "integer", "integer", "ANY"), function(x, i, j, ...) {
  if(length(i) == 1 && length(j) == 1){
    # return scalar
  } else if(length(i) == 1){
    # return tabVector, look up a row
  } else if(length(j) == 1){
    # return tabVector, look up a column
  } else if(length(i) > 1 && length(j) > 1){
    # return tabMatrix
  } else stop("subsetting method not implemented")
})

setMethod("[", c("tabMatrix", "missing", "integer", "ANY"), function(x, i, j, ...) {
  if(length(j) == 1){
    # return scalar
  } else {
    # return tabVector
  }
})

setMethod("[", c("tabMatrix", "integer", "missing", "ANY"), function(x, i, j, ...) {
  if(length(i) == 1){
    # return scalar
  } else {
    # return tabVector
  }
})
