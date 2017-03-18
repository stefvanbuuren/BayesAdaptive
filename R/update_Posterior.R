#' Compute posterior for a hierarchical trial design
#' 
#' Compoute functions of the posterior require for the adaptive randomization and for applying early the decision rules.
#' 
#' @details
#' For \code{Design[1] = 1} or \code{Design[1] = 2} the function returns 
#' functions of the posterior for a Phase II Biomarker-finding or  Biomarker-stratified design
#' For \code{Design[2] = 1, 2, 3} specifise the presence of an (1) active, (2) historical control arm, and (3) the absence of a control arm.  
#' 
#' 
#' 
#' 
#' @param Resp_Events 
#' Array of observed response/succeses for the control and experimental agents for all disease by subpopulations.
#' 
#' @param Resp_RiskPop 
#' Array of observed outcome for the control and experimental agents for all disease by subpopulations.
#' 
#' @param Prior
#' A list of prior parameter generated with \code{Prior.list.Fct}
#' 
#' @param Design 
#' A vector of length 2 with specifies the experimental design of the trial.
#' For \code{Design[1] = 1} or \code{Design[1] = 2} the function returns 
#' functions of the posterior for a Phase II Biomarker-finding or Biomarker-stratified design
#' For \code{Design[2] = 1, 2, 3} specifise the presence of an (1) active, (2) historical control arm, and (3) the absence of a control arm.  
#' 
#' @param alpha.historical
#' A matrix of the inverse normal of the historical/null response probabilities by disease (rows) and subpopulations (columns), i.e.
#' \code{pnorm(alpha.historical[d,m])} should be equal to \code{p.historical[d,m]}.
#' 
#' @param p.historical
#'  A matrix of  historical/null response probabilities by disease (rows) and subpopulations (columns).
#'  Note that it is required that \code{pnorm(alpha.historical[d,m])} should be equal to \code{p.historical[d,m]}.
#' 
#' @import mnormt
#' 
#' @return
#' A list of three elements;  
#' (i) \code{Rand.Pr} array of non-amplified randomization statistics,
#' (ii) \code{Fut.Pr} array of  futility statistics,
#' (iii) \code{Eff.stat} array of  efficacy statistics.
#' 
#' @author Steffen Ventz  \email{ventzer@@yahoo.de}
#' 
#' @examples
#' ### example1 : 2 cancer, 5 experimental agents and 2 modules, active control
#'  eligibility.array        = array(1, c(2, 5, 2))
#'  Resp_Events              = Resp_RiskPop = array(0, c(2, 5, 2))
#'  eligibility.array[2,2,2] = 0
#'  eligibility.array[1,,1]  = 0
#'  Dim                      = dim(eligibility.array)
#'  
#'  ### response rates
#'  rate                      = array(NA, Dim)
#'  rate[,,1]                 = matrix(c(.3, .3, .5, .3, .3, .3, .3, .2, .4, .3), 2, 5, byrow=TRUE)
#'  rate[,,2]                 = matrix(c(.3, .2, .4, .3, .3, .3, .3, .2, .4, .3), 2, 5, byrow=TRUE)
#'  
#'  ### generate data
#'  Resp_Events = Resp_RiskPop = array(0, dim=Dim)
#'  
#'  for(d in 1:Dim[1]) for(a in 1:Dim[2])  
#'  Resp_RiskPop[d,a,] = rpois(Dim[3], 150)
#'  for(d in 1:Dim[1]) for(a in 1:Dim[2])   
#'  Resp_Events[d,a,] = rbinom(Dim[3], Resp_RiskPop[d,a,], rate[d,a,])
#'  
#'  Resp_Events [eligibility.array==0] = 0
#'  Resp_RiskPop[eligibility.array==0] = 0
#'      
#'  Prior = Prior.list.Fct(eligibility.array=eligibility.array, 
#'               Var.vec = c(alpha1=1, alpha2=2, beta1=.2, beta2=.1, beta3=.05))
#'  updatePosterior(Resp_Events, Resp_RiskPop, Prior, Design = c(2,1)) ## SSD
#'  updatePosterior(Resp_Events, Resp_RiskPop, Prior, Design = c(1,1)) ## SFD 
#'  
#'  
#' ### example2 : 2 cancer, 5 experimental agents and 2 modules, no/historical control
#'  eligibility.array[,1,]             = 0
#'  
#'  Prior = Prior.list.Fct(eligibility.array=eligibility.array, 
#'                  Var.vec = c(alpha1=0, alpha2=0, beta1=.2, beta2=.1, beta3=.05))
#'  Resp_Events[eligibility.array==0]  = Resp_RiskPop[eligibility.array==0] =0
#'                    
#' 
#'  updatePosterior(Resp_Events, Resp_RiskPop, Prior, Design = c(1,3),
#'                   alpha.historical   = matrix(qnorm(.3), nrow=Dim[1], ncol=Dim[3]),
#'                   p.historical       = matrix(.3, nrow=Dim[1], ncol=Dim[3]) ) ## SFD no control
#' 
#'  updatePosterior(Resp_Events, Resp_RiskPop,Prior, Design = c(2,2),
#'                   alpha.historical   = matrix(qnorm(.3), nrow=Dim[1], ncol=Dim[3]),
#'                   p.historical       = matrix(.3, nrow=Dim[1], ncol=Dim[3]) )
#' @export



updatePosterior = function( Resp_Events, 
                            Resp_RiskPop, 
                            Prior, 
                            Design = c(2,3),
                            alpha.historical=NULL,
                            p.historical=NULL){
 
 if(Design[2]!=1) if(is.null(alpha.historical) || is.null(p.historical)) 
  stop("specify alpha.historical and p.historical without active control agent") 
 
 if(is.null(Prior)) stop("speficy prior parameter with ")
 
 ###### Newton-Rapson Algorithm                              Prior$eligibility.vec
 E_and_R       = get.vec.from.mat(Resp_Events, Resp_RiskPop, Prior$eligibility.vec$e.vec.indicator)
 post.mode.Cov = BvM.approximation(Prior$Prior_list$theta, Prior$Prior_list$m1H.prior, E_and_R[,1], E_and_R[,2])
 
 
 ###### put posterior mode and variance of theta[d,a,m] into the array mode.theta.array[d,a,m] and var.theta.array[d,a,m]
 Dim                    = dim(Resp_Events)
 nr.theta               = length(post.mode.Cov$theta.mode)
 theta.var              = diag(post.mode.Cov$Asymp.Var)
 mode.theta.array       = array(0, dim=Dim)
 Cov.theta.with.control = Cov.theta.with.control = var.theta.array = array(0, dim=Dim)    
 
 
 ## theta parameter ####################################################################################
 for(j in 1:nr.theta){mode.theta.array[Prior$eligibility.vec$e.vec.d.a.m[1,j],
                                       Prior$eligibility.vec$e.vec.d.a.m[2,j],
                                       Prior$eligibility.vec$e.vec.d.a.m[3,j]] =  post.mode.Cov$theta.mode[j]
 
                      var.theta.array[Prior$eligibility.vec$e.vec.d.a.m[1,j],
                                      Prior$eligibility.vec$e.vec.d.a.m[2,j],
                                      Prior$eligibility.vec$e.vec.d.a.m[3,j]]  = theta.var[j] }
 

  ## historical response
  if(Design[2]!=1) mode.theta.array[,1,]  = alpha.historical
    
  ## for an active control: covariance betw control and experimental
  if(Design[2]==1){  for(i in 1: Prior$eligibility.vec$nr.agent ){
    
    control.ind = which( Prior$eligibility.vec$control.var.index[1,] == Prior$eligibility.vec$agent.var.index[1,i] &
                         Prior$eligibility.vec$control.var.index[3,] == Prior$eligibility.vec$agent.var.index[3,i] )
                              
    Cov.theta.with.control[ Prior$eligibility.vec$agent.var.index[1,i],
                            Prior$eligibility.vec$agent.var.index[2,i],
                            Prior$eligibility.vec$agent.var.index[3,i]] =  post.mode.Cov$Asymp.Var[control.ind, Prior$eligibility.vec$nr.control + i]    

   }}
 
 
  ## treatment effect ####################################################################################
                    beta.mode  = lapply(1:Dim[3], function(m) mode.theta.array[,-1,m] - mode.theta.array[,1,m])
  if(Design[2]==1){ beta.var   = lapply(1:Dim[3], function(m)  var.theta.array[,-1,m] + var.theta.array[,1,m] - 2*Cov.theta.with.control[,-1,m]) 
  }else             beta.var   = lapply(1:Dim[3], function(m)  var.theta.array[,-1,m])

 
 
  ## randomization statistics ####################################################################################
  if(Design[2]!=3){  Rand.Pr.l   = lapply(1:Dim[3], function(m) pnorm(0, mean=beta.mode[[m]], sd= sqrt(beta.var[[m]]), lower.tail=FALSE))                 
  }else{             Rand.Pr.l   = lapply(1:Dim[3], function(m) apply(matrix(1:(Dim[2]-1)),1, function(a) sapply(1:Dim[1], function(d){                     
                                    if(Prior$eligibility.vec$eligibility.TRUE[d,a+1,m]){
                                           ac      = which(Prior$eligibility.vec$eligibility.TRUE[d,-1,m])
                                           ac      = ac[ac!=a]
                                           Cov.mat = diag(beta.var[[m]][d,ac]) + beta.var[[m]][d,a]
                                           Fut.p   = mnormt::sadmvn(lower  = rep(-Inf, length(ac)), 
                                                                    upper  = rep(0, length(ac)), 
                                                                    mean   = beta.mode[[m]][d,ac]-beta.mode[[m]][d,a], varcov = Cov.mat)[1]
                                    }else{ Fut.p = 0 }
                                           return(Fut.p) }))) }
 
                    Rand.Pr      = array(dim=Dim-c(0,1,0))
 for(m in 1:Dim[3]) Rand.Pr[,,m] = Rand.Pr.l[[m]]
 
 
 
 
 ## futility statistics ####################################################################################
 if(Design[1]==2){    Fut.Pr      = Rand.Pr                   ## Subpopulation-Stratified design  
 }else{               Fut.Pr      = array(dim=Dim-c(0,1,0))
   for(m in 1:Dim[3]) Fut.Pr[,,m] = matrix(apply(Rand.Pr[,,m], 2, max), Dim[1], Dim[2]-1, byrow=TRUE)
 }

 
 
 
 ## efficacy statistics ####################################################################################
  if(Design[2]==1){ Eff.stat.single = Z.statistics.arr(N.0= Resp_RiskPop[,1,], N.a=Resp_RiskPop[,-1,], R.0=Resp_Events[,1,], R.a=Resp_Events[,-1,])
  }else             Eff.stat.single = Z.statistics.arr(N.a=Resp_RiskPop[,-1,], R.a=Resp_Events[,-1,], One.sample = TRUE, p.historical = p.historical)
  
if(Dim[3]==1){      Eff.stat        = array(dim=Dim-c(0,1,0))
                    Eff.stat[,,1]   = Eff.stat.single
}else               Eff.stat        = Eff.stat.single

return(list(Rand.Pr  = Rand.Pr, 
            Fut.Pr   = Fut.Pr,
            Eff.stat = Eff.stat))}


