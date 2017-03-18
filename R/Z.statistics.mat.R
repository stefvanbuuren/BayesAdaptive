#' Array-valued Z-statistics
#'
#' Computes the coordinate-wise Z-statistics
#'
#' The function takes the array of observed events and number of observed outcomes,
#' and computes an array of Z-statistics for Binomial data.
#'
#' @param N.0
#' Matrix of observed outcome for the standard of care for all disease by subpopulations.
#'
#' @param N.a
#' Array of observed outcome for the experimental agents for all disease by subpopulations.
#'
#' @param R.0
#' Matrix of observed response/succeses for the standard of care for all disease by subpopulations.
#'
#' @param R.a
#' Array of observed response/succeses for the experimental agents for all disease by subpopulations.
#'
#' @param One.sample
#' If \code{TRUE}, use z-statistics for  the onse-sided  null hypothesis \eqn{ p_{d,a,m} \le } \code{p.historical[d,m]}.
#'
#' @param p.historical
#' An matrix  of probabilityies corresponding to the null respons probability \eqn{p_{d,a,m} \le } \code{p.historical[d,m]}.
#'
#' @author Steffen Ventz  \email{ventzer@@yahoo.de}
#' @keywords internal
#'

# # @examples
# ### example: 2 cancer, 2 experimental agents and 2 modules, active control
# eligibility.array = array(1, c(2, 5, 2))
# Dim               = dim(eligibility.array)
# Resp_Events       = Resp_RiskPop = array(0, Dim)
# eligibility.array[2,2,2] = 0
# eligibility.array[1,,1] = 0
#
#
# ### response rates
# rate                     = array(NA, Dim)
# rate[,,1]                = matrix(c(.3, .3, .5, .3, .3, .3, .3, .2, .4, .3), 2, 5, byrow=TRUE)
# rate[,,2]                = matrix(c(.3, .2, .4, .3, .3, .3, .3, .2, .4, .3), 2, 5, byrow=TRUE)
#
# ### generate data
# Resp_Events = Resp_RiskPop = array(dim=Dim)
#
# for(d in 1:Dim[1]) for(a in 1:Dim[2])
# Resp_RiskPop[d,a,] = rpois(Dim[3], 150)
# for(d in 1:Dim[1]) for(a in 1:Dim[2])
# Resp_Events[d,a,] = rbinom(Dim[3], Resp_RiskPop[d,a,], rate[d,a,])
#
# Resp_Events[eligibility.array==0]  = NA
# Resp_RiskPop[eligibility.array==0] = NA
#
# Z.statistics.arr(N.0=Resp_RiskPop[,1,], N.a=Resp_RiskPop[,-1,],
#                  R.0=Resp_Events[,1,], R.a=Resp_Events[,-1,])
#




Z.statistics.arr = function(N.0=NULL, N.a, R.0=NULL, R.a, One.sample=FALSE, p.historical=NULL){

if(!One.sample){      Dim           = dim(R.a)

 if(length(Dim)>2){   N.0.mat       = R.0.mat = array(0, dim=Dim)
  for(m in 1:Dim[3]){ N.0.mat[,,m]  = N.0[,m]
                      R.0.mat[,,m]  = R.0[,m] }

 }else{               N.0.mat       = matrix(N.0, Dim[1], Dim[2])
                      R.0.mat       = matrix(R.0, Dim[1], Dim[2])
 }




  No.event.0 = N.0.mat ==0
  No.event.a = N.a     ==0

  P.a        = R.a/N.a
  P.0        = R.0.mat/N.0.mat
  P          = (R.a + R.0.mat)/(N.0.mat + N.a)
  Z.DN       = sqrt(P*(1-P) * (1/N.a +1/N.0.mat) )
  Z          = (P.a-P.0) / Z.DN

  Z[Z.DN == 0]  = 0
  Z[No.event.0] = 0
  Z[No.event.a] = 0

}else{
 if(is.null(p.historical))
  stop("response probability array p.0.array missing")

  p.0.array  = array(dim=dim(R.a));  for(m in 1:dim(R.a)[3]) p.0.array[,,m] = matrix(p.historical[,m], dim(R.a)[1], dim(R.a)[2])

  Z          = (R.a-N.a*p.0.array)/ sqrt(N.a*p.0.array*(1-p.0.array))
  Z[N.a==0]  = 0
}


return(Z)}
