#' RcppSparse
#'
#' @description
#' An Rcpp sparse matrix class header file for iterator-based read-only access to \code{Matrix::`dgCMatrix-class`} objects in R.
#'
#' @details 
#' The RcppSparse class is a zero-copy by-reference object constructed from \code{Rcpp::IntegerVector} and \code{Rcpp::NumericVector} 
#' that is to objects of class \code{dgCMatrix} as \code{Rcpp::NumericMatrix} is to objects of class \code{"matrix"}.
#' 
#' The class is equipped with iterators that allow for efficient marginal traversal of values in columns, rows, or ranges of a column.
#' Performance of the iterators exceeds that of \code{RcppArmadillo::SpMat} and approaches that of \code{Eigen::SparseMatrix}.
#' 
#' The class iterators and member functions in \code{RcppSparse::Matrix} align with method counterparts in the \code{Eigen::SparseMatrix<double>}, such
#' that templated functions accepting either class can be used without further specialization.
#' 
#' Note: This class is read-only. While access to the \code{i}, \code{p}, \code{Dim}, and \code{x} vectors is public, there are no methods 
#' for inserting or erasing values.
#' 
#' See code for example \code{\link{columnSums}} to get started, and the package vignette for documentation.
#' 
#' @docType package
#' @name RcppSparse
#' @author Zach DeBruine
#' @md
#' @seealso \code{\link{columnSums}}
#'
NULL