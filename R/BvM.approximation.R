
#' Posterior for the Bayesian Probit Model
#'
#' The function compute the posterior mode of a Bayesian probit model for Gaussion prior using a Neton-Rapson Algorithm
#'
#' @param theta.mu starting value of the algorithm
#' @param m1H.prior the negative Hessian of the log-prior density, which is the negative of the inverse prior covariance matrix
#' @param Events A vector of length \code{nrow(m1H.prior)} of sucesses
#' @param Risk A vector of length \code{nrow(m1H.prior)} of sucesses plus failure
#' @param maxir Maximum number of iterations of the algorithm
#' @param eps Stoping criteria, stop if the L1 distance between the parameter vector for two succesive iterations is below \code{eps}.
#'
#' @return
#' The function returns a lsit of two object \code{theta.mode} and \code{Asymp.Var}, which are the posterior mode and the approximate Covatiance matrix of the estimated posterior mode.
#'
#' @author Steffen Ventz  \email{ventzer@@yahoo.de}
#' @keywords internal



BvM.approximation = function(theta.mu=rep(0,length(Events)), m1H.prior, Events, Risk, maxir=50, eps=10^-5){
  i         = 0
  theta     = theta.mu
  continue  = T

  ## Newton-Rapson algorithm
  while(continue){
    i          = i+1
    theta.old  = theta
    G_and_H    = gradient.Hessian.ll(as.vector(theta), Events, Risk)
    G.log.post = G_and_H$G + m1H.prior %*% theta
    H.log.post = G_and_H$H + m1H.prior
    theta      = theta - solve(H.log.post) %*% G.log.post
    continue   = (sum(abs(theta-theta.old)) > eps) && (i <= maxir)  }

  Asymp.Var  = -solve(m1H.prior+gradient.Hessian.ll(as.vector(theta), Events, Risk)$H)

  return(list(theta.mode=as.vector(theta), Asymp.Var=Asymp.Var))
}
