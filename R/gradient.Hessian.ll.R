#' Hessian and Gradient for a Probit model
#'
#' The function computes the Hessian and Gradient of the log-likelihood under a Probit model
#'
#' @param theta parameter vector of a logit model
#' @param Events vector of observed response events (sucesses). First row of the matrix created with \code{get.vec.from.mat} (see vignette for details)
#' @param Risk vector of observed events (sucesses and failure). Second row of the matrix created with \code{get.vec.from.mat} (see vignette for details)
#'
#' @return A list of two elements (H, G) where H anf G are the Hessian and Gradient of the log-likelihood under a Probit model
#'
#' @keywords internal
#' @author Steffen Ventz  \email{ventzer@@yahoo.de}
gradient.Hessian.ll = function(theta, Events, Risk){

  CDF   = pnorm(theta)
  cCDF  = 1-CDF
  p     = dnorm(theta)
  p_CDF = p/CDF
  p_cCDF= p/cCDF
  N     = (Risk-Events)
  G     = Events*p_CDF -N*p_cCDF
  H     = -Events*((p_CDF)^2 +theta*(p_CDF)) +N*(theta*p_cCDF -(p_cCDF)^2)

  return(list(G=G, H=diag(H))) }
