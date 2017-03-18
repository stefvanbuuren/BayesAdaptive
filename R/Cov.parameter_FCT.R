
#' Covariance Matrix for a the Gaussian Prior
#'
#' Covariance Matrix for the Gaussian Prior for a Bayesian Hierarchical Probit Regression Model
#'
#' @param eligibility.array
#'  Eligibility array, a two or three dimensional array of zeros and ones.
#' The entry \code{eligibility.array[d,a,m]} is one if patients with disease \code{d} in subpopulation \code{m} are eligible for
#' treatment with agent \code{a}.
#'
#' @param Var.vec Vector of five length five of variances of the parameter
#'  \eqn{\alpha_{d}}, \eqn{\alpha_{d,m}}, \eqn{\beta_{a}},  \eqn{\beta_{a,m}} and \eqn{\beta_{d,a,m}}.
#'
#' @return
#' Returns the covariance matrix of a Gausian prior for one module with variance components \code{var.alpha}, \code{var.beta} and \code{var.gamma}.
#' @examples
#' eligibility.array      = array( dim=c(4,4,2))
#' eligibility.array[,,1] = cbind( c(0,1,0,1), rep(0,4), c(0,1,0,1), c(0,0,0,1))
#' eligibility.array[,,2] = cbind( c(1,1,1,1), c(1,0,1,1), c(0,1,1,0), c(1,0,0,1))
#' dimnames(eligibility.array) = list( paste0("disease_", 1:4),
#' c("control", paste0("agent", 1:3)),
#' paste0("marker", 1:2))
#'
#' ## Exa_1: one module only
#' Cov.parameter_FCT(eligibility.array[,,2], c(0, .8, 0, .1, .05))
#'
#' ## Exa_2: two subpopulations
#' Cov.parameter_FCT(eligibility.array, c(.1, .8, .1, .1, .05))
#'
#' @export
#'
#' @author Steffen Ventz  \email{ventzer@@yahoo.de}



Cov.parameter_FCT = function(eligibility.array, Var.vec=c(0, .8, 0, .1, .05) ){

 Dim    = dim(eligibility.array)

 if(length(Dim)==2){
  Dim        = c(Dim,1)
  e          = eligibility.array
  eligibility.array      = array(dim=Dim)
  eligibility.array[,,1] = e

 }else{ if(length(Dim)<2) stop("specify proper two or three dimensional eligibility array") }

 C.mat1 = array(0, dim = c(dim(eligibility.array), dim(eligibility.array)) )
 C.mat  = matrix(0, prod(Dim), prod(Dim))
 Name   = rep(NA, nrow(C.mat))

for(m in 1:Dim[3]) for(mp in 1:Dim[3]) for(a in 1:Dim[2]) for(ap in 1:Dim[2]) for(d in 1:Dim[1]) for(dp in 1:Dim[1]){

      ## control and experimental agents (alpha part)
      if(d==dp & eligibility.array[d,a,m] & eligibility.array[dp,ap,mp]==1 )
      {           C.mat1[d,a,m, dp,ap,mp] = Var.vec[1]
       if(m==mp)  C.mat1[d,a,m, dp,ap,mp] = Var.vec[2] + C.mat1[d,a,m, dp,ap,mp]
      }

      ## experimental agents (beta part)
      if(a==ap & a>1 & eligibility.array[d,a,m] & eligibility.array[dp,ap,mp]==1 )
      {           C.mat1[d,a,m, dp,ap,mp] = Var.vec[3] + C.mat1[d,a,m, dp,ap,mp]
       if(m==mp)
       {          C.mat1[d,a,m, dp,ap,mp] = Var.vec[4] + C.mat1[d,a,m, dp,ap,mp]
        if(d==dp) C.mat1[d,a,m, dp,ap,mp] = Var.vec[5] + C.mat1[d,a,m, dp,ap,mp]
       }
      }
}

 for(a in 1:Dim[2]) for(d in 1:Dim[1]) for(m in 1:Dim[3]) {
                                          id         = id.fct(a=a,m=m,d=d, Dim)
                                          C.vec      = NULL
  for(ap in 1:Dim[2]) for(mp in 1:Dim[3]) C.vec      = c(C.vec, C.mat1[d,a,m, ,ap,mp])
                                          C.mat[,id] = C.vec
                                 if(a==1) Name[id]   = paste0("C_dm", d, m)
                                 if(a >1) Name[id]   = paste0("A_adm", a-1, d, m)
 }

rownames(C.mat) = colnames(C.mat) = Name
C.rowsum        = rowSums(C.mat) == 0
C.colsum        = colSums(C.mat) == 0
C.mat           = C.mat[!C.rowsum,]
C.mat           = C.mat[,!C.colsum]

return(C.mat)}

