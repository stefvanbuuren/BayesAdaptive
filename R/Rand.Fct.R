#' Response Adaptive Randomization
#'
#' Response adaptive randomization
#'
#' @param Rand.Pr
#' Posterior pandomization statistics, created as output with \code{update_Posterior}.
#'
#' @param Resp_RiskPop
#' Array of number of observed outcome for the control and experimental agents for all disease by subpopulations.
#'
#' @param Accrual
#' Array of number of patients randimized to the control and experimental agents for all disease by subpopulations.
#'
#' @param Active
#' Array of "logical" variables which indicate the activity status
#' for the control and experimental agents for all disease by subpopulations.
#'
#' @param rand.vec
#' A list of randomization parameter to  compute the power function \eqn{h_{d,m}(N) = a_{d,m} N^b }. The elements are:
#' (1) \code{a} array of identical slopes for each experimental agents by disease and subpolulation.
#' (2) \code{b} The exponent \eqn{b}.
#' (3) \code{c} The exponent for the control arms.
#' (4) \code{N.star} The minimum number of patients to be randomized to each combination of arm/disease/subpopulation.
#'
#' @param Design
#' A vector of length 2 with specifies the experimental design of the trial.
#' For \code{Design[1] = 1, 2} the function returns
#' functions of the posterior for a Phase II subpopulation-finding (\code{Design[1] = 1})
#' and subpopulation-stratified design (\code{Design[1] = 2}).
#' For \code{Design[2] = 1, 2, 3} specifise the presence of an (1) active,
#' (2) historical control arm, and (3) the absence of a control arm.
#'
#'
#' @examples
#'
#' ### example1 : 2 cancer, 4 experimental agents and 2 modules, active control
#'  eligibility.array        = array(1, c(2, 5, 2))
#'  Resp_Events              = Resp_RiskPop = array(0, c(2, 5, 2))
#'  eligibility.array[2,2,2] = 0
#'  eligibility.array[1,,1]  = 0
#'  Dim                      = dim(eligibility.array)
#'  rate                     = array(dim=Dim)
#'  rate[,,1] = rate[,,2]    = matrix(c(.3, .3, .5, .3, .3, .3, .3, .2, .4, .3), 2, 5, byrow=TRUE)
#'  Resp_Events              = Resp_RiskPop = array(0, dim=Dim)
#'
#'  for(d in 1:Dim[1]) for(a in 1:Dim[2])
#'  Resp_RiskPop[d,a,] = rpois(Dim[3], 150)
#'  for(d in 1:Dim[1]) for(a in 1:Dim[2])
#'  Resp_Events[d,a,] = rbinom(Dim[3], Resp_RiskPop[d,a,], rate[d,a,])
#'
#'  Resp_Events[eligibility.array==0]  = 0
#'  Resp_RiskPop[eligibility.array==0] = 0
#'
#' Prior      = Prior.list.Fct(eligibility.array,  Var.vec=c(1, 2, .2, .1, .05))
#' Posterior1 = updatePosterior(Resp_Events, Resp_RiskPop, Prior, c(2,1))
#' E.sample.s = cbind(c(10, 20), c(30, 40))
#' A          = array(dim=c(2,4,2))
#' for(m in 1:Dim[3]) for(d in 1:Dim[1]) A[d,,m] = E.sample.s[d,m]
#' rand.vec   = list(a = 4/A^(-log(4)/log(.5)), b= -log(4)/log(.5) , c=0.25, N.star=5)
#'
#' Rand.Fct(Posterior1$Rand.Pr, Resp_RiskPop, Accrual=Resp_RiskPop,
#'         Active=eligibility.array, rand.vec, Design=c(2,1))
#'
#' ### example1 : 2 cancer, 5 experimental agents and 2 modules, no active or historical control
#' eligibility.array[,1,]=0
#' Prior                 = Prior.list.Fct(eligibility.array,  Var.vec=c(0, 0, 0, .1, .05))
#' p.historical          = matrix(.3, 2, 2)
#' alpha.historical      = qnorm(p.historical)
#' Posterior2 = updatePosterior(Resp_Events, Resp_RiskPop, Prior, c(1,3),
#'                              alpha.historical, p.historical)
#' Rand.Fct(Posterior2$Rand.Pr, Resp_RiskPop, Accrual=Resp_RiskPop,
#'         Active=eligibility.array, rand.vec, Design=c(1,3))
#'
#'
#' @export
#' @author Steffen Ventz  \email{ventzer@@yahoo.de}
#'


## randomization function
Rand.Fct = function(Rand.Pr, Resp_RiskPop, Accrual, Active, rand.vec, Design){

### add power-function to randomization probabilities "rand_prob"
Resp_RiskPop[!Active] = 0
not.active         = !Active
Dim                = dim(Resp_RiskPop)
p                  = array(0, dim=Dim)
n.total            = array(dim=Dim-c(0,1,0)); for(m in 1:Dim[3]) n.total[,,m] = rowSums(Resp_RiskPop[,,m], na.rm = TRUE)
n.Active           = sapply(1:Dim[3], function(m) rowSums(Active[,,m]))                     ## nr of active+eligible experimental arms by disease

if(Dim[2]>2){   n.max = sapply(1:Dim[3], function(m) apply(Accrual[,-1,m], 1, max, na.rm=TRUE))
}else       {   n.max = Accrual[,2,]        }

# ## remove the maximum and normalize
# for(m in 1:Dim[3]) p[,-1,m]     = p[,-1,m] - matrix(apply( matrix(p[,-1,m], Dim[1], Dim[2]-1), 1,  max, na.rm = TRUE), Dim[1], Dim[2]-1, FALSE)
#                   p[,-1,]      = exp(p[,-1,])

                   p[,-1, ]     = Rand.Pr ^ (rand.vec$a * n.total ^ rand.vec$b)
                   p[!Active]   = 0



for(m in 1:Dim[3]){row.sum      =  rowSums(p[,,m])
                   NULL.rows    =  row.sum == 0
                   p[,,m]       =  p[,,m] / row.sum
                p[NULL.rows,,m] = 0                       }  ## if no agent active for disease d, use only control



## put control & experimental arm together
if(Design[2] ==1){   p[,1,] = exp(rand.vec$c*(n.max-Resp_RiskPop[,1,])) /(n.Active-1)
              p[not.active] = 0

  for(m in 1:Dim[3]){  R.sum  = rowSums(p[,,m])
                       p[,,m] = p[,,m] / matrix(R.sum, nrow=Dim[1], ncol=Dim[2])
                    p[R.sum==0,,m] = 0 }
}



for(m in 1:Dim[3]){
 if(is.null(dim(rand.vec$N.star)))
 {     p[,,m]        = (p[,,m] + 10^-2)* exp(apply(rand.vec$N.star     -Accrual[,,m], 2, function(x) pmax.int(x, 0)))
 }else p[,,m]        = (p[,,m] + 10^-2)* exp(apply(rand.vec$N.star[,,m]-Accrual[,,m], 2, function(x) pmax.int(x, 0)))

 p[not.active] = 0
 p[,,m]        = p[,,m] / matrix(rowSums(p[,,m]), nrow=Dim[1], ncol=Dim[2])
 p[not.active] = 0
}


return(p)}

