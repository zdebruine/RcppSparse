#' Tabulated sparse vector format
#' 
#' @details
#' 
#' A tabVector efficiently represents sparse and redundant data, for example, counts or ratings. A tabVector cannot contain values less than zero. 
#' The length of the vector (in dense form) is given in the \code{length} slot, while the vector values and row indices are given in the \code{x} slot.
#' 
#' The \code{x} slot contains values (indicated by \code{-} signs) followed by row indices corresponding to that value. 
#' Row indices are always sorted in increasing order between values, and must be unique within the vector.
#' Values are usually unique, but uniqueness is not required. Each value instance must be followed by at least one row index.
#' 
#' A tabMatrix is simply a collection of concatenated tabVectors, one for each column, with a pointer vector giving the start index for each column.
#' 
#' An extensive Rcpp "tabMatrix" class is implemented, and supports most operations available in R as well iterators for random access within a vector and traversal through the intersection between vectors. 
#' The Rcpp header library and class is documented in the package vignette.
#' 
#' @section tabVector Methods:
#' 
#' Coercion:
#'      "as.numeric", "as.vector", "as.integer", "as.matrix"
#'      "as(from, to)" between "itabVector" or "dtabVector" and any of "sparseVector", "integer", "numeric", "matrix". Coercing "numeric" to "integer" will "floor" values.
#'
#' Subset:
#'      "[" Subsetting is implemented for indices given as a scalar, dense numeric/integer/logical vector, or "nsparseVector".
#'      "[<-" Writing is supported, but coerces to/from "sparseVector"
#' 
#' Arithmetic:
#' 
#'      "+" "-" "/" "*" between tabVector and scalar/dense vector returns dense
#'      "*" between two tabVectors returns tabVector ("+" "/" "-" not implemented)
#'      "%*%" "crossprod(x, y = NULL)" between tabVectors and/or dense/tabVector returns appropriate scalar
#' 
#' Compare:
#' 
#'      "==" ">" "<" "!=" "<=" ">=" between a tabVector and a scalar/dense vector or between two tabVectors returns \code{nsparseVector}
#'      "all.equal" between two tabVectors or tabVector/dense vector
#' 
#' Math:
#' 
#'      For \code{itabVector} returning \code{dtabVector}: "sqrt", "log", "log10", "log2", "exp"
#'      For \code{dtabVector} returning \code{itabVector}: "ceiling", "floor", "trunc", "round" (digits = 0)
#'      For \code{dtabVector} returning \code{dtabVector}: "round" (digits > 0), "sqrt", "log", "log10", "log2", "acos", "acosh", "asin", "asinh", "atan", "atanh", "exp", "cos", "cosh", "cospi", "sin", "sinh", "sinpi", "tan", "tanh", "tanpi"
#'
#' Summary:
#' 
#'      "mean", "max", "min", "range", "prod", "sum", "any", "all"
#' 
#' Distance:
#' 
#'      \code{dist(x, y, method = "euclidean")}: Compute distance between two tabVectors, or a tabVector and a dense vector. Implemented methods are "euclidean", "manhattan", "intersection", "cosine", "jaccard", and "kldivergence".
#'      \code{cor(x, y)}, \code{cosine(x, y)}: Pearson correlation or cosine similarity between two tabVectors or a dense vector/tabVector pair.
#'
#' Other:
#'      "!" return dense logical vector giving zero-valied indices
#'      "length" "head" "tail", "show", "print"
#' 
setClass("tabVector",
         representation(length = "integer", x = "integer"),
         prototype(length = 0L, x = NA_integer_))

# coerce to/from isparseVector/integer
setAs("tabVector", "integer", function(from, to){ Rcpp_tabVector_to_integer(from) })
setAs("tabVector", "isparseVector", function(from, to) { Rcpp_tabVector_to_isparseVector(from) })
setAs("isparseVector", "tabVector", function(from, to) { Rcpp_isparseVector_to_tabVector(from) })
setAs("integer", "tabVector", function(from, to) { as(as(as(from, "sparseVector"), "isparseVector"), "tabVector") })

# coerce to/from other common types using the above routines
setAs("tabVector", "sparseVector", function(from, to) { as(as(from, "isparseVector"), "sparseVector") })
setAs("tabVector", "numeric", function(from, to) { as(as(from, "integer"), "numeric") })
setAs("tabVector", "logical", function(from, to) { as(from, "integer") != 0 })
setAs("tabVector", "matrix", function(from, to) { as(as(from, "integer"), "matrix") })
setAs("tabVector", "vector", function(from, to) { as(as(from, "numeric"), "vector") })
setAs("tabVector", "array", function(from, to) { as(as(from, "numeric"), "array") })
setAs("sparseVector", "tabVector", function(from, to) { as(as(from, "isparseVector"), "tabVector") })
setAs("numeric", "tabVector", function(from, to) { as(as(from, "integer"), "tabVector") })
setAs("matrix", "tabVector", function(from, to) { as(as(as(from, "vector"), "integer"), "tabVector") })
setAs("vector", "tabVector", function(from, to) { as(as(from, "integer"), "tabVector") })
setAs("array", "tabVector", function(from, to) { as(as(as(from, "numeric"), "integer"), "tabVector") })
setMethod("as.integer", "tabVector", function(x){ as(x, "integer") })
setMethod("as.vector", "tabVector", function(x){ as(x, "numeric") })
setMethod("as.numeric", "tabVector", function(x) { as(x, "numeric") })
setMethod("as.matrix", "tabVector", function(x) { as(x, "matrix") })
setMethod("as.logical", "tabVector", function(x) { as(x, "logical") })
setMethod("as.array", "tabVector", function(x) { as(x, "array") })

setMethod("[", "tabVector", function(x, i) {
  if(length(i) == 1) {
    if(i < 0) {
      # remove a single value
    } else {
      # get a single value
      Rcpp_subset_tabvector_int(x, i)
    }
  } else {
    if(all(i < 0)) {
      # remove several values
    } else if(all(i > 0)){
      # get several values
      Rcpp_subset_tabvector_vector(x, i)
    } else stop("mixed sign subsetting is not implemented for tabVector")
  }
})

setMethod("length", "tabVector", function(x) x@length)

setMethod("print", "tabVector", function(x){
  cat("sparse tabulated vector (length =", length(x), ") of class \"tabVector\"\n")
  len <- getOption("width")
  if(len > length(x)) len <- length(x)
  # subset v to len
  x <- as(x, "integer")
  x <- x[1:len]
  print.table(x, zero.print = ".")
})

setMethod("show", "tabVector", function(object) { print(object) })

rtabvector <- function(length = 100, density = 0.1, rand.x = function(n) sample(1:5, n, replace = TRUE)){
  v <- rand.x(length)
  v <- v * sample(0:1, length, prob = c(1 - density, density), replace = T)
  v <- as.integer(v)
  as(v, "tabVector")
}

setMethod("!", "tabVector", function(x){
  Rcpp_negate(x)
})

# NEXT
setMethod("+", c(e1="tabVector", e2="integer"), function(e1, e2){
    if(length(e2 == 1)){
      new(Class = "tabVector", x = sapply(e1, function(x) { ifelse(x > 0, x + e2, x) }), length = length(e1))
    } else if(length(e2) == length(e1)){
        # new(Class = "tabVector", x = ?, length = length(e1))
    } else {
      stop("length of LHS and RHS of + operation were not equal")
    }
})

setMethod("all.equal", c("tabVector", "tabVector"), function(target, current){
  if(target@length == current@length){
    if(length(target@x) == length(current@x)){
      if(all(target@x == current@x)){
        return(TRUE)
      }
    }
  }
  return(FALSE)
})

setMethod("all.equal", c("tabVector", "integer"), function(target, current) all.equal(target, as(current, "tabVector")))
setMethod("all.equal", c("tabVector", "numeric"), function(target, current) all.equal(target, as(current, "tabVector")))
setMethod("all.equal", c("tabVector", "matrix"), function(target, current) all.equal(target, as(current, "tabVector")))
setMethod("all.equal", c("tabVector", "sparseVector"), function(target, current) all.equal(target, as(current, "tabVector")))
setMethod("all.equal", c("tabVector", "array"), function(target, current) all.equal(target, as(current, "tabVector")))

setMethod("head", "tabVector", function(x, n = 6L) { x[1:n] })
setMethod("tail", "tabVector", function(x, n = 6L) { x[(length(x) - n): length(x)] })
setMethod("sum", "tabVector", function(x) { Rcpp_sum_tabvector(x) })
setMethod("mean", "tabVector", function(x) { sum(x) / length(x) })

setMethod("Compare", c("integer", "tabVector"), function(e1, e2){
  Compare(e1, as(e2, "sparseVector"))
})

setMethod("crossprod", c("tabVector", "tabVector"), function(x, y) {
  Rcpp_crossprod(x, y)
})

setMethod("crossprod", c("tabVector", "integer"), function(x, y) {
  Rcpp_crossprod2(x, y)
})

# methods - [ [<- * / &

# methods:
#  [<- assign to subset
#  %*% crossprod (s.t)
#  all.equal test all values for equality
#  anyNA test for any NA values
# Arith + - / *
# c
# Logic
# Math
# Math2
# Ops
# rep
# tcrossprod
# dot