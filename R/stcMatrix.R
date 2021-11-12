#' The "stcMatrix" class
#' 
#' Sparse Tabulated Counts Matrix
#' 
#' @details 
#' The \code{stcMatrix} class is a lightweight column-major storage format for CSC matrices with highly redundant discrete values.
#' 
#' As redundancy among values increases, the size of an \code{stcMatrix} approaches half that of the corresponding \code{dgCMatrix}.
#' 
#' \code{stcMatrix} is ideal when memory is a limitation, when uncompressed file size must be minimized, and for applications requiring
#' only random-access column-major constant iterator access. Random access column iterators are conforable to \code{dgCMatrix}, but
#' operations that benefit from ordered access to rows (i.e. cross-products, transposition) are inefficient and not supported.
#' 
#' Some operations are much faster in an \code{stcMatrix} compared to \code{dgCMatrix}, such as column sums and certain column-wise operations, 
#' because of improved vectorization.
#' 
#' The only supported method on the R side is \code{as()} to/from a \code{dgCMatrix}. Note that the value of \code{slot("x")} must be a non-negative integer vector.
#' 
#' @slot Dim vector giving rows and columns in matrix
#' @slot p vector of pointers to indices in \code{x} corresponding to column starts
#' @slot x vector of values and row indices, values are multiplied by -1 and followed by all row indices in that column corresponding to each value. Values and row indices are non-decreasing.
#' @md
setClass("stcMatrix", representation(Dim = "integer", p = "integer", x = "integer"))

setAs("dgCMatrix", "stcMatrix", function(from, to) {
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
  as_stcMatrix(from)
})

setAs("stcMatrix", "dgCMatrix", function(from, to) {
  to <- as_dgCMatrix(from)
  ind <- which.min(to@x) - 1
  to@x <- to@x[1:ind]
  to@i <- to@i[1:ind]
  to
})