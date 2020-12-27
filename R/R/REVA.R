#' REVA: testing whther a sclalar value differs less for coler points  
#'
#' We hava a scalar value for a set of points and we have a distance between the points.
#' Null hypothesis: the scalar does not depend ont the location, so the values differ between distant points as strons as between close points
#' Main function: 
#' 
#' @section REVA functions:
#'	REVA.test calculates the p-value corresponding to the null
#'	REVA.user.distance calculate distance as Kendall tau and calculate the p-value	
#'	REVA.user.distance.Geneset does REVA.user.distance for a gene set

#' @docType package
#' @name REVA 
#' @useDynLib REVA
#' @importFrom Rcpp evalCpp
# these two are Rcpp - scecific ivocations 
NULL
